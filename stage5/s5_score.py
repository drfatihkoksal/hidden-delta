"""HIDDEN-DELTA Aşama 5 — dondurulmuş modellerle skorlama, oklüzyon ve integrated gradients (plan §1, §5.1–5.2, §7.5).

Resmî çalıştırma (MIMIC, I0006): önce `FROZEN_MANIFEST.sha256` doğrulanır; çıktı zaten varsa çalışmaz (bir kez).
Prova (`--rehearsal`): çıktı ~/data/hidden-delta/stage5/rehearsal/<site>/. Test setlerinde (MIMIC, I0006) prova,
gerçek model yerine rastgele ağırlıklı bir modelle yapılır; böylece dondurmadan önce hiçbir gerçek test skoru üretilmez.
MGH (I0001) provası gerçek modelle yapılır (geliştirme verisi).

Çıktı: scores.parquet (k; M1 logit/olasılık; M0, M-demo, M2 logit; oklüzyon logitleri; atıf payları), MIMIC'te h3_scores.parquet.
Çalıştırma: python3 stage5/s5_score.py MIMIC [--rehearsal]
"""
import argparse
import subprocess

import numpy as np
import pandas as pd
import torch

import common5 as C

TEST = {"MIMIC", "I0006"}
MANIFEST = C.ROOT / "FROZEN_MANIFEST.sha256"


def verify_manifest():
    if not MANIFEST.exists():
        raise SystemExit("FROZEN_MANIFEST.sha256 yok: önce dondurma (s5_freeze.py).")
    r = subprocess.run(["sha256sum", "-c", "--quiet", str(MANIFEST)], cwd=C.ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("Manifest doğrulanamadı:\n" + r.stdout + r.stderr)


def scramble(fz, seed):
    """Prova: gerçek model yerine rastgele ağırlıklar (test setlerinde gerçek skor üretilmez)."""
    torch.manual_seed(seed)
    for m in fz.enc.modules():
        if hasattr(m, "reset_parameters"):
            m.reset_parameters()
    rng = np.random.default_rng(seed)
    fz.w = (rng.normal(size=fz.w.shape) / np.sqrt(len(fz.w))).astype(np.float32)
    for p in fz.jm.values():
        p["lr"]["coef"] = rng.normal(size=len(p["lr"]["coef"])).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site", choices=["MIMIC", "I0006", "I0001"])
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    src = C.P5 / a.site
    out = (C.P5 / "rehearsal" / a.site) if a.rehearsal else src
    if not a.rehearsal:
        if a.site not in TEST:
            raise SystemExit("Resmî çalıştırma yalnızca MIMIC ve I0006 için.")
        verify_manifest()
        if (out / "scores.parquet").exists():
            raise SystemExit("scores.parquet zaten var: skorlama bir kez yapılır.")
    out.mkdir(parents=True, exist_ok=True)
    fz = C.Frozen()
    if a.rehearsal and a.site in TEST:
        scramble(fz, C.SEED + 7)
    idx = pd.read_parquet(src / "index.parquet")
    mb = np.load(src / "index_mb.npy")
    Z = fz.embed(mb)
    logit = Z @ fz.w + fz.b
    sc = pd.DataFrame({"k": idx.k, "logit_M1": logit, "p_M1": fz.m1_prob(logit),
                       "logit_M0": fz.lin("M0", idx), "logit_demo": fz.lin("M-demo", idx), "logit_M2": fz.lin("M2", idx, Z)})
    # oklüzyon (§5.1): başlangıç ±40 ms, ST-T (+40…+400 ms), aynı genişlikte rastgele pencereler
    w_on, w_st = zip(*[C.windows(o, f) for o, f in zip(idx.on_low_ms, idx.off_ms)])
    r_on = C.random_windows([b - a_ for a_, b in w_on], C.SEED + 11)
    r_st = C.random_windows([b - a_ for a_, b in w_st], C.SEED + 12)
    for name, wins in (("occ_onset", w_on), ("occ_stt", w_st), ("occ_rand_onset", r_on), ("occ_rand_stt", r_st)):
        sc[name] = fz.embed(C.masked(mb, wins)) @ fz.w + fz.b
    # integrated gradients (§5.2): mutlak atfın pencerelere düşen payı
    attr = np.abs(fz.integrated_gradients(mb))
    tot = attr.sum(axis=(1, 2))
    share = lambda wins: np.array([attr[i, :, a_:b].sum() for i, (a_, b) in enumerate(wins)]) / tot
    sc["share_onset"], sc["share_stt"] = share(w_on), share(w_st)
    p_win = [(C.ms2s(p0), C.ms2s(p1)) if np.isfinite(p0) and np.isfinite(p1) else (0, 0)
             for p0, p1 in zip(idx.p_on_ms, idx.p_end_ms)]
    sc["share_p"] = share(p_win)
    if (src / "index_raw.npy").exists():                         # §8.7: ham 10 s girişli M1
        import s4_m1raw as R
        enc = R.EncoderRaw(R.DEPTH).to(C.DEV)
        enc.load_state_dict(torch.load(C.P4 / "m1raw" / "final_encoder.pt"))
        pr = np.load(C.P4 / "m1raw" / "final_probe.npz")
        w, b0 = pr["coef"].ravel(), float(pr["intercept"][0])
        if a.rehearsal and a.site in TEST:
            torch.manual_seed(C.SEED + 8)
            for m in enc.modules():
                if hasattr(m, "reset_parameters"):
                    m.reset_parameters()
            w = np.random.default_rng(C.SEED + 8).normal(size=w.shape) / np.sqrt(len(w))
        raw, rok = np.load(src / "index_raw.npy"), np.load(src / "index_raw_ok.npy")
        lr = R.embed(enc, raw) @ w + b0
        sc["logit_M1raw"] = np.where(rok, lr, np.nan)
    if (src / "agg.parquet").exists():                           # §8.4: hasta düzeyinde ortalama/maksimum
        agg = pd.read_parquet(src / "agg.parquet")
        agg.assign(logit_M1=fz.m1_logit(np.load(src / "agg_mb.npy"))).to_parquet(out / "agg_scores.parquet", index=False)
    if a.site == "MIMIC" and (src / "h3.parquet").exists():
        h3 = pd.read_parquet(src / "h3.parquet")
        hl = fz.m1_logit(np.load(src / "h3_mb.npy"))
        h3.assign(logit_M1=hl, p_M1=fz.m1_prob(hl)).to_parquet(out / "h3_scores.parquet", index=False)
    sc.to_parquet(out / "scores.parquet", index=False)
    print(f"{a.site}{' (prova)' if a.rehearsal else ''}: {len(sc):,} EKG skorlandı → {out}")


if __name__ == "__main__":
    main()
