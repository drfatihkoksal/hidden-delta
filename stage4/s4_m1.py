"""HIDDEN-DELTA Aşama 4 — M1: medyan vuru üzerinde 1D ResNet (ANALYSIS_PLAN §4). Yalnızca MGH.

Girdi: baskın küme medyan vurusu, 12 × 500 (500 Hz, fiducial −400…+600 ms), mV.
Arama alanı (genişletilmez): derinlik {4, 8} blok × mod {ince ayar, linear probe} = 4 yapılandırma.
Her dış kat için: (1) SSL — maskelenmiş yeniden yapılandırma; havuz = MGH'nin T2 dışı çıkarılmış EKG'leri, o katın
test hastaları hariç; (2) ince ayar (iç doğrulama: eğitim kümelerinin %10'u, erken durdurma) ya da linear probe
(SSL gömmesi üzerinde L2 lojistik, C iç CV ile). Seçim: kat dışı (OOF) birleşik AUROC.
Final: seçilen yapılandırma tüm MGH ile (SSL: tüm havuz; ince ayar: katların en iyi epoch medyanı kadar).
Kalibrasyon: Platt, seçilen yapılandırmanın OOF tahminleriyle. Eşik: OOF kalibre tahminlerde kontrollerin 95. persentili.
Hız: bf16 karma hassasiyet ve cuDNN algoritma seçimi (eğitim bit düzeyinde tekrarlanabilir değil; final ağırlıklar özetlenir).
Augmentasyon: rastgele filtre zinciri (nedensel ya da sıfır fazlı Butterworth alçak geçiren 40/100/150 Hz; 50/60 Hz notch;
ya da hiçbiri), genlik ölçekleme 0.9–1.1, ±4 örnek kaydırma, 5 µV gürültü.

Çıktı: stage4/out/m1_mgh.md (agregat); ~/data/hidden-delta/stage4/m1/ (ağırlıklar, OOF tahminler, kalibrasyon).
Çalıştırma: python3 stage4/s4_m1.py
"""
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy import signal
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_auc_score

P4 = Path("/home/fatih/data/hidden-delta/stage4")
OUTD = P4 / "m1"
OUT = Path(__file__).resolve().parent / "out"
SEED = 20261119
FS, L, NPAD = 500, 500, 1000
DEV = torch.device("cuda")
SSL_EPOCHS, FT_MAX_EPOCHS, PATIENCE = 15, 40, 6
CONFIGS = [(4, "finetune"), (4, "probe"), (8, "finetune"), (8, "probe")]
AMP = dict(device_type="cuda", dtype=torch.bfloat16)   # hız: bf16 karma hassasiyet (RTX 5090 tensor çekirdekleri)
warnings.filterwarnings("ignore", category=FutureWarning)


def seed_all(s):
    np.random.seed(s)
    torch.manual_seed(s)
    torch.cuda.manual_seed_all(s)
    # hız için cuDNN algoritma seçimi; eğitim bit düzeyinde tekrarlanabilir değildir (final ağırlıklar dondurulup özetlenir)
    torch.backends.cudnn.deterministic = False
    torch.backends.cudnn.benchmark = True


# ------------------------------------------------------------------ filtre augmentasyonu (frekans alanında, GPU)
def filter_bank():
    w = np.fft.rfftfreq(NPAD, 1 / FS)
    bank = [np.ones_like(w, complex)]
    for fc in (40, 100, 150):
        sos = signal.butter(2, fc, btype="low", fs=FS, output="sos")
        _, h = signal.sosfreqz(sos, worN=w, fs=FS)
        bank += [h, np.abs(h) ** 2]                       # nedensel, sıfır fazlı (ileri-geri)
    for f0 in (50, 60):
        b, a = signal.iirnotch(f0, 30, fs=FS)
        _, h = signal.freqz(b, a, worN=w, fs=FS)
        bank += [h, np.abs(h) ** 2]
    return torch.tensor(np.stack(bank), dtype=torch.complex64, device=DEV)


BANK = None


def augment(x, gen):
    n = x.shape[0]
    idx = torch.randint(0, BANK.shape[0], (n,), device=DEV, generator=gen)
    xp = F.pad(x, (0, NPAD - L), mode="replicate")
    x = torch.fft.irfft(torch.fft.rfft(xp, dim=-1) * BANK[idx][:, None, :], n=NPAD, dim=-1)[..., :L]
    x = x * (0.9 + 0.2 * torch.rand(n, 1, 1, device=DEV, generator=gen))
    sh = int(torch.randint(-4, 5, (1,), generator=gen, device=DEV))
    x = torch.roll(x, sh, dims=-1)
    return x + 0.005 * torch.randn(x.shape, device=DEV, generator=gen)


# ------------------------------------------------------------------ model
class Block(nn.Module):
    def __init__(self, cin, cout, stride):
        super().__init__()
        self.c1 = nn.Conv1d(cin, cout, 7, stride, 3, bias=False)
        self.b1 = nn.BatchNorm1d(cout)
        self.c2 = nn.Conv1d(cout, cout, 7, 1, 3, bias=False)
        self.b2 = nn.BatchNorm1d(cout)
        self.sk = nn.Sequential() if (cin == cout and stride == 1) else nn.Sequential(
            nn.Conv1d(cin, cout, 1, stride, bias=False), nn.BatchNorm1d(cout))

    def forward(self, x):
        return F.relu(self.b2(self.c2(F.relu(self.b1(self.c1(x))))) + self.sk(x))


class Encoder(nn.Module):
    def __init__(self, depth):
        super().__init__()
        self.stem = nn.Sequential(nn.Conv1d(12, 32, 15, 1, 7, bias=False), nn.BatchNorm1d(32), nn.ReLU())
        chans = [32, 64, 128, 128]
        per = depth // 4
        layers, cin = [], 32
        for i, c in enumerate(chans):
            for j in range(per):
                layers.append(Block(cin, c, 2 if (j == 0 and i > 0) else 1))
                cin = c
        self.blocks = nn.Sequential(*layers)
        self.dim = cin

    def forward(self, x):
        return self.blocks(self.stem(x))                 # B × 128 × T'

    def embed(self, x):
        return self.forward(x).mean(-1)


class SSL(nn.Module):
    def __init__(self, enc):
        super().__init__()
        self.enc = enc
        self.dec = nn.Sequential(nn.Conv1d(enc.dim, 64, 5, 1, 2), nn.ReLU(), nn.Conv1d(64, 12, 5, 1, 2))

    def forward(self, x):
        h = self.enc(x)
        return self.dec(F.interpolate(h, size=L, mode="linear", align_corners=False))


class Clf(nn.Module):
    def __init__(self, enc):
        super().__init__()
        self.enc, self.head = enc, nn.Linear(enc.dim, 1)

    def forward(self, x):
        return self.head(self.enc.embed(x)).squeeze(-1)


def mask(x, gen):
    """Her örnekte 4 rastgele 40 ms'lik zaman bölümü ve rastgele 2 derivasyon sıfırlanır."""
    m = torch.ones_like(x)
    n = x.shape[0]
    for _ in range(4):
        st = torch.randint(0, L - 20, (n,), device=DEV, generator=gen)
        ar = torch.arange(L, device=DEV)[None, :]
        m = m * (~((ar >= st[:, None]) & (ar < st[:, None] + 20)))[:, None, :].float()
    leads = torch.randint(0, 12, (n, 2), device=DEV, generator=gen)
    for k in range(2):
        m[torch.arange(n, device=DEV), leads[:, k], :] = 0
    return x * m, m


# ------------------------------------------------------------------ eğitim
def pretrain(X, depth, seed):
    seed_all(seed)
    gen = torch.Generator(device=DEV).manual_seed(seed)
    model = SSL(Encoder(depth)).to(DEV)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-2)
    steps = SSL_EPOCHS * int(np.ceil(len(X) / 256))
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=1e-3, total_steps=steps)
    Xt = torch.tensor(X, device=DEV)
    model.train()
    for ep in range(SSL_EPOCHS):
        perm = torch.randperm(len(Xt), device=DEV, generator=gen)
        for i in range(0, len(Xt), 256):
            xb = augment(Xt[perm[i:i + 256]], gen)
            xm, m = mask(xb, gen)
            with torch.autocast(**AMP):
                rec = model(xm)
            loss = (((rec.float() - xb) ** 2) * (1 - m)).sum() / (1 - m).sum().clamp(min=1)
            opt.zero_grad()
            loss.backward()
            opt.step()
            sched.step()
    return model.enc


@torch.no_grad()
def predict(model_or_enc, X, embed=False):
    model_or_enc.eval()
    out = []
    for i in range(0, len(X), 512):
        xb = torch.tensor(X[i:i + 512], device=DEV)
        with torch.autocast(**AMP):
            o = model_or_enc.embed(xb) if embed else model_or_enc(xb)
        out.append(o.float().cpu().numpy())
    return np.concatenate(out)


def finetune(enc, Xtr, ytr, Xva, yva, seed, epochs=None):
    """epochs verilirse erken durdurma yok (final model); verilmezse iç doğrulamayla en iyi epoch."""
    seed_all(seed)
    gen = torch.Generator(device=DEV).manual_seed(seed)
    model = Clf(enc).to(DEV)
    opt = torch.optim.AdamW([{"params": model.enc.parameters(), "lr": 3e-4},
                             {"params": model.head.parameters(), "lr": 1e-3}], weight_decay=1e-2)
    pos_w = torch.tensor((1 - ytr.mean()) / ytr.mean(), device=DEV)
    Xt, yt = torch.tensor(Xtr, device=DEV), torch.tensor(ytr, dtype=torch.float32, device=DEV)
    best, best_ep, best_state, wait = -1, 0, None, 0
    n_ep = epochs or FT_MAX_EPOCHS
    for ep in range(1, n_ep + 1):
        model.train()
        perm = torch.randperm(len(Xt), device=DEV, generator=gen)
        for i in range(0, len(Xt), 128):
            j = perm[i:i + 128]
            with torch.autocast(**AMP):
                logit = model(augment(Xt[j], gen))
            loss = F.binary_cross_entropy_with_logits(logit.float(), yt[j], pos_weight=pos_w)
            opt.zero_grad()
            loss.backward()
            opt.step()
        if epochs is None:
            auc = roc_auc_score(yva, predict(model, Xva))
            if auc > best:
                best, best_ep, wait = auc, ep, 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                wait += 1
                if wait >= PATIENCE:
                    break
    if epochs is None:
        model.load_state_dict(best_state)
    return model, best_ep


def probe_fit(Ztr, ytr):
    return LogisticRegressionCV(Cs=10, cv=5, scoring="roc_auc", max_iter=5000, random_state=SEED).fit(Ztr, ytr)


def main():
    global BANK
    t0 = time.time()
    OUTD.mkdir(parents=True, exist_ok=True)
    BANK = filter_bank()
    dev = pd.read_parquet(P4 / "I0001_dev.parquet")
    Xd = np.load(P4 / "I0001_dev_mb.npy").astype(np.float32)
    ssl = pd.read_parquet(P4 / "I0001_ssl.parquet")
    Xs = np.load(P4 / "I0001_ssl_mb.npy", mmap_mode="r")
    y = dev.y.to_numpy()
    oof = {c: np.full(len(dev), np.nan) for c in CONFIGS}
    best_eps = {c: [] for c in CONFIGS}
    rng = np.random.default_rng(SEED)
    for f in range(5):
        te = dev.fold.to_numpy() == f
        tr_idx = np.where(~te)[0]
        tr_sets = dev.set_id.to_numpy()[tr_idx]
        va_sets = set(rng.choice(np.unique(tr_sets), int(0.1 * len(np.unique(tr_sets))), replace=False))
        va = np.array([s in va_sets for s in tr_sets])
        test_pids = set(dev.pid[te])
        Xpool = np.asarray(Xs[(~ssl.pid.isin(test_pids)).to_numpy()], dtype=np.float32)
        for depth in (4, 8):
            enc = pretrain(Xpool, depth, SEED + 10 * f + depth)
            torch.save(enc.state_dict(), OUTD / f"ssl_fold{f}_d{depth}.pt")
            # linear probe
            Ztr, Zte = predict(enc, Xd[tr_idx], embed=True), predict(enc, Xd[te], embed=True)
            pr = probe_fit(Ztr, y[tr_idx])
            oof[(depth, "probe")][te] = pr.decision_function(Zte)
            # ince ayar (SSL ağırlıklarından başlar)
            enc2 = Encoder(depth).to(DEV)
            enc2.load_state_dict(torch.load(OUTD / f"ssl_fold{f}_d{depth}.pt"))
            model, ep = finetune(enc2, Xd[tr_idx[~va]], y[tr_idx[~va]], Xd[tr_idx[va]], y[tr_idx[va]], SEED + f)
            oof[(depth, "finetune")][te] = predict(model, Xd[te])
            best_eps[(depth, "finetune")].append(ep)
            print(f"kat {f} derinlik {depth}: probe AUROC {roc_auc_score(y[te], oof[(depth, 'probe')][te]):.3f}, "
                  f"ince ayar {roc_auc_score(y[te], oof[(depth, 'finetune')][te]):.3f} (epoch {ep})  "
                  f"{time.time() - t0:.0f} s", flush=True)
    rows = []
    for c in CONFIGS:
        rows.append({"derinlik": c[0], "mod": c[1], "OOF AUROC": round(roc_auc_score(y, oof[c]), 4),
                     "kat AUROC'ları": ", ".join(f"{roc_auc_score(y[dev.fold == f], oof[c][dev.fold == f]):.3f}"
                                                 for f in range(5)),
                     "en iyi epoch'lar": ", ".join(map(str, best_eps[c])) or "—"})
    tab = pd.DataFrame(rows)
    best = max(CONFIGS, key=lambda c: roc_auc_score(y, oof[c]))
    # Platt kalibrasyonu ve eşik: seçilen yapılandırmanın OOF skorları
    platt = LogisticRegression(C=1e6, max_iter=1000).fit(oof[best].reshape(-1, 1), y)
    p_oof = platt.predict_proba(oof[best].reshape(-1, 1))[:, 1]
    thr = float(np.quantile(p_oof[y == 0], 0.95))
    sens = float((p_oof[y == 1] > thr).mean())
    pd.DataFrame({"set_id": dev.set_id, "pid": dev.pid, "k": dev.k, "y": y, "fold": dev.fold,
                  **{f"oof_d{c[0]}_{c[1]}": oof[c] for c in CONFIGS}, "p_M1": p_oof}).to_parquet(OUTD / "oof_m1.parquet")
    # final model: seçilen yapılandırma, tüm MGH
    depth, mode = best
    enc = pretrain(np.asarray(Xs, dtype=np.float32), depth, SEED + 100 + depth)
    if mode == "probe":
        Z = predict(enc, Xd, embed=True)
        pr = probe_fit(Z, y)
        torch.save(enc.state_dict(), OUTD / "final_encoder.pt")
        np.savez(OUTD / "final_probe.npz", coef=pr.coef_, intercept=pr.intercept_)
    else:
        n_ep = int(np.median(best_eps[best]))
        model, _ = finetune(enc, Xd, y, None, None, SEED + 200, epochs=n_ep)
        torch.save(model.state_dict(), OUTD / "final_model.pt")
    meta = {"config": {"depth": depth, "mode": mode}, "platt": [float(platt.coef_[0, 0]), float(platt.intercept_[0])],
            "threshold_spec95": thr, "oof_sens_at_thr": sens, "best_epochs": best_eps.get(best),
            "seed": SEED, "ssl_epochs": SSL_EPOCHS, "runtime_s": round(time.time() - t0)}
    (OUTD / "final_meta.json").write_text(json.dumps(meta, indent=1))
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# M1 — MGH 5 katlı çapraz doğrulama (geliştirme; hipotez testi değil)\n",
          tab.to_markdown(index=False),
          f"\nSeçilen yapılandırma: derinlik {depth}, {mode} (OOF AUROC en yüksek).",
          f"Platt kalibrasyonu OOF skorlarla; %95 özgüllük eşiği (OOF kontrollerin 95. persentili) {thr:.4f}; "
          f"bu eşikte OOF duyarlılık {sens:.3f}.",
          f"Final model tüm MGH ile eğitildi. Süre {meta['runtime_s']} s.\n",
          "Not: geliştirme sayılarıdır; birincil test MIMIC'te dondurulmuş modelle (ANALYSIS_PLAN §1)."]
    (OUT / "m1_mgh.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
