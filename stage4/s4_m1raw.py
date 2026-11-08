"""HIDDEN-DELTA Aşama 4 — duyarlılık analizi §8.7: ham 10 s girişli M1 (yalnızca MGH'de geliştirme).

M1'in seçilen yapılandırması korunur (derinlik 8, linear probe); yalnızca girdi değişir: taban çizgisi düzeltilmiş
(medyan filtre) 10 s, 12 derivasyon, 500 Hz (12 × 5000). Başlangıç katmanı 2 adımla örnekleme düşürür. SSL: maskelenmiş
yeniden yapılandırma (8 × 100 ms bölüm + 2 derivasyon), aynı filtre augmentasyonu (8192 noktalı FFT). 5 katlı CV (kat
başına SSL, katın test hastaları hariç) ve final model (tüm MGH). Platt kalibrasyonu OOF skorlarla.

Veri: `build` modu ham dizileri üretir (dev ve SSL havuzu, float16). Çalıştırma:
  python3 stage4/s4_m1raw.py build      # ~birkaç dk, ~6 GB
  python3 stage4/s4_m1raw.py train      # GPU
"""
import json
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy import signal
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage2"))
import hdsig  # noqa: E402
import s4_m1 as M1  # noqa: E402
import statlib as S  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
P2 = Path("/home/fatih/data/hidden-delta/stage2/I0001")
P4 = Path("/home/fatih/data/hidden-delta/stage4")
OUTD = P4 / "m1raw"
OUT = Path(__file__).resolve().parent / "out"
SEED = 20261121
FS, LR, NPADR = 500, 5000, 8192
DEPTH = 8
DEV = M1.DEV
AMP = M1.AMP


def load_raw(rec):
    """Taban çizgisi düzeltilmiş 12 × 5000 (float16); 500 Hz ve 10 s değilse None."""
    try:
        fs, X = hdsig.load_record(rec)
    except Exception:  # noqa: BLE001
        return None
    if fs != FS or X.shape[1] < LR:
        return None
    return hdsig.remove_baseline(X[:, :LR], fs).astype(np.float16)


def raw_array(recs, workers=28):
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(load_raw, recs, chunksize=16))
    out = np.zeros((len(recs), 12, LR), np.float16)
    ok = np.zeros(len(recs), bool)
    for i, r in enumerate(res):
        if r is not None:
            out[i], ok[i] = r, True
    return out, ok


def build():
    jobs = pd.read_parquet(P2 / "jobs.parquet")
    rec = dict(zip(jobs.k, jobs.rec))
    for name in ("dev", "ssl"):
        t = pd.read_parquet(P4 / f"I0001_{name}.parquet")
        arr, ok = raw_array([rec[k] for k in t.k])
        np.save(P4 / f"I0001_{name}_raw.npy", arr)
        np.save(P4 / f"I0001_{name}_raw_ok.npy", ok)
        print(f"{name}: {ok.sum():,}/{len(ok):,} ham kayıt", flush=True)


# ------------------------------------------------------------------ model (M1 bloklarıyla, girdi uzunluğundan bağımsız)
class EncoderRaw(M1.Encoder):
    def __init__(self, depth):
        super().__init__(depth)
        self.stem = nn.Sequential(nn.Conv1d(12, 32, 15, 2, 7, bias=False), nn.BatchNorm1d(32), nn.ReLU())


class SSLRaw(nn.Module):
    def __init__(self, enc):
        super().__init__()
        self.enc = enc
        self.dec = nn.Sequential(nn.Conv1d(enc.dim, 64, 5, 1, 2), nn.ReLU(), nn.Conv1d(64, 12, 5, 1, 2))

    def forward(self, x):
        return self.dec(F.interpolate(self.enc(x), size=x.shape[-1], mode="linear", align_corners=False))


def bank_raw():
    w = np.fft.rfftfreq(NPADR, 1 / FS)
    bank = [np.ones_like(w, complex)]
    for fc in (40, 100, 150):
        _, h = signal.sosfreqz(signal.butter(2, fc, btype="low", fs=FS, output="sos"), worN=w, fs=FS)
        bank += [h, np.abs(h) ** 2]
    for f0 in (50, 60):
        b, a = signal.iirnotch(f0, 30, fs=FS)
        _, h = signal.freqz(b, a, worN=w, fs=FS)
        bank += [h, np.abs(h) ** 2]
    return torch.tensor(np.stack(bank), dtype=torch.complex64, device=DEV)


BANK = None


def augment(x, gen):
    n = x.shape[0]
    idx = torch.randint(0, BANK.shape[0], (n,), device=DEV, generator=gen)
    xp = F.pad(x, (0, NPADR - LR), mode="replicate")
    x = torch.fft.irfft(torch.fft.rfft(xp, dim=-1) * BANK[idx][:, None, :], n=NPADR, dim=-1)[..., :LR]
    x = x * (0.9 + 0.2 * torch.rand(n, 1, 1, device=DEV, generator=gen))
    x = torch.roll(x, int(torch.randint(-4, 5, (1,), generator=gen, device=DEV)), dims=-1)
    return x + 0.005 * torch.randn(x.shape, device=DEV, generator=gen)


def mask(x, gen):
    n = x.shape[0]
    m = torch.ones_like(x)
    ar = torch.arange(LR, device=DEV)[None, :]
    for _ in range(8):
        st = torch.randint(0, LR - 50, (n,), device=DEV, generator=gen)
        m = m * (~((ar >= st[:, None]) & (ar < st[:, None] + 50)))[:, None, :].float()
    leads = torch.randint(0, 12, (n, 2), device=DEV, generator=gen)
    for k in range(2):
        m[torch.arange(n, device=DEV), leads[:, k], :] = 0
    return x * m, m


def pretrain(X, seed, epochs=15, bs=128):
    M1.seed_all(seed)
    gen = torch.Generator(device=DEV).manual_seed(seed)
    model = SSLRaw(EncoderRaw(DEPTH)).to(DEV)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-2)
    steps = epochs * int(np.ceil(len(X) / bs))
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=1e-3, total_steps=steps)
    Xt = torch.tensor(X, device=DEV)                    # float16 GPU'da; adımda float32'ye
    model.train()
    for _ in range(epochs):
        perm = torch.randperm(len(Xt), device=DEV, generator=gen)
        for i in range(0, len(Xt), bs):
            xb = augment(Xt[perm[i:i + bs]].float(), gen)
            xm, m = mask(xb, gen)
            with torch.autocast(**AMP):
                rec = model(xm)
            loss = (((rec.float() - xb) ** 2) * (1 - m)).sum() / (1 - m).sum().clamp(min=1)
            opt.zero_grad()
            loss.backward()
            opt.step()
            sched.step()
    del Xt
    torch.cuda.empty_cache()
    return model.enc


@torch.no_grad()
def embed(enc, X, bs=256):
    enc.eval()
    out = []
    for i in range(0, len(X), bs):
        with torch.autocast(**AMP):
            out.append(enc.embed(torch.tensor(X[i:i + bs], device=DEV).float()).float().cpu().numpy())
    return np.concatenate(out)


def train():
    global BANK
    t0 = time.time()
    OUTD.mkdir(parents=True, exist_ok=True)
    BANK = bank_raw()
    dev = pd.read_parquet(P4 / "I0001_dev.parquet")
    ssl = pd.read_parquet(P4 / "I0001_ssl.parquet")
    Xd, okd = np.load(P4 / "I0001_dev_raw.npy"), np.load(P4 / "I0001_dev_raw_ok.npy")
    Xs = np.load(P4 / "I0001_ssl_raw.npy", mmap_mode="r")
    oks = np.load(P4 / "I0001_ssl_raw_ok.npy")
    y = dev.y.to_numpy()
    oof = np.full(len(dev), np.nan)
    for f in range(5):
        te = (dev.fold.to_numpy() == f) & okd
        tr = (dev.fold.to_numpy() != f) & okd
        pool = oks & ~ssl.pid.isin(set(dev.pid[dev.fold == f])).to_numpy()
        enc = pretrain(np.asarray(Xs[pool]), SEED + f)
        torch.save(enc.state_dict(), OUTD / f"ssl_fold{f}.pt")
        pr = LogisticRegressionCV(Cs=10, cv=5, scoring="roc_auc", max_iter=5000, random_state=SEED, n_jobs=-1).fit(
            embed(enc, Xd[tr]), y[tr])
        oof[te] = pr.decision_function(embed(enc, Xd[te]))
        print(f"kat {f}: AUROC {S.auroc(y[te], oof[te]):.3f}  {time.time() - t0:.0f} s", flush=True)
    m = ~np.isnan(oof)
    auc, ci = S.auc_boot(y[m], oof[m], dev.set_id.to_numpy()[m], 2000, SEED)
    platt = LogisticRegression(C=1e6, max_iter=1000).fit(oof[m].reshape(-1, 1), y[m])
    enc = pretrain(np.asarray(Xs[oks]), SEED + 100)
    torch.save(enc.state_dict(), OUTD / "final_encoder.pt")
    pr = LogisticRegressionCV(Cs=10, cv=5, scoring="roc_auc", max_iter=5000, random_state=SEED, n_jobs=-1).fit(
        embed(enc, Xd[okd]), y[okd])
    np.savez(OUTD / "final_probe.npz", coef=pr.coef_, intercept=pr.intercept_)
    meta = {"config": {"depth": DEPTH, "mode": "probe", "input": "raw 10 s, 12 × 5000, 500 Hz"},
            "platt": [float(platt.coef_[0, 0]), float(platt.intercept_[0])], "oof_auroc": auc, "oof_ci": ci.tolist(),
            "n_dev": int(m.sum()), "seed": SEED, "runtime_s": round(time.time() - t0)}
    (OUTD / "final_meta.json").write_text(json.dumps(meta, indent=1))
    pd.DataFrame({"k": dev.k, "y": y, "set_id": dev.set_id, "oof_raw": oof}).to_parquet(OUTD / "oof_m1raw.parquet")
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# §8.7 duyarlılık: ham 10 s girişli M1 — MGH 5 katlı CV (geliştirme)\n",
          f"Yapılandırma M1 ile aynı (derinlik 8, linear probe); girdi 12 × 5000. OOF AUROC {auc:.3f} (%95 GA {ci[0]:.3f}–{ci[1]:.3f}); "
          f"{int(m.sum()):,} hasta. Final model tüm MGH ile. Süre {meta['runtime_s']} s."]
    (OUT / "m1raw_mgh.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    {"build": build, "train": train}[sys.argv[1]]()
