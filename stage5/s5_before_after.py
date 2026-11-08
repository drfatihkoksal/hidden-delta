"""HIDDEN-DELTA — kayıt sonrası keşifsel analiz (2026-11-08): indeks EKG preeksite EKG'den önce mi, sonra mı?

Birincil K1 kümelerinde vakalar indeks EKG'nin zamanına göre ayrılır:
- önce: indeks EKG hastanın ilk preeksite EKG'sinden önce;
- sonra: ilk preeksite EKG'den sonra (aralarda preeksite EKG'ler olabilir); ayrıca en yakın preeksite EKG'ye ≤30 gün ve
  >30 gün (kardiyak bellek en çok ilk haftalarda beklenir).
Her tabaka kendi eşleşmiş kontrolleriyle analiz edilir. AUROC farkı (sonra − önce) için tabakalar içinde ayrı ayrı küme
bootstrap'ı (2,000). Preeksite EKG zamanı bilinmeyen vakalar (işlenmiş veride T2 yok) tabakalara girmez. Siteler: MGH (M1 OOF, geliştirme), Emory ve MIMIC (dondurulmuş M1, resmî skorlar). Yeni skorlama yok.

Çalıştırma: python3 stage5/s5_before_after.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage4"))
import statlib as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
D = Path("/home/fatih/data/hidden-delta")
NB, SEED = 2000, 20261108


def mgh():
    dev = pd.read_parquet(D / "stage4/I0001_dev.parquet", columns=["pid", "k", "t", "y", "set_id"])
    oof = pd.read_parquet(D / "stage4/m1/oof_m1.parquet", columns=["k", "oof_d8_probe"])
    e = pd.read_parquet(D / "stage1/heedb_I0001_ecg_labels.parquet", columns=["pid", "t", "t2_cons"])
    t2 = e[e.t2_cons.fillna(False).astype(bool) & e.pid.isin(set(dev.pid[dev.y == 1]))]
    first = t2.groupby("pid").t.min()
    times = t2.groupby("pid").t.apply(lambda s: s.to_numpy())
    d = dev.merge(oof, on="k").rename(columns={"oof_d8_probe": "score"})
    c = d[d.y == 1].copy()
    c["before"] = c.t < c.pid.map(first)
    c["days"] = [np.min(np.abs((tt - t) / np.timedelta64(1, "D"))) for p, t in zip(c.pid, c.t.to_numpy()) for tt in [times[p]]]
    return d, c[["set_id", "before", "days"]]


def test_site(s):
    sets = pd.read_parquet(D / "stage5" / s / "sets.parquet")
    sets = sets[sets.analysis == "primary_k1"]
    idx = pd.read_parquet(D / "stage5" / s / "index.parquet", columns=["pid", "k", "before_first_t2", "days_to_t2"])
    sc = pd.read_parquet(D / "stage5" / s / "scores.parquet", columns=["k", "logit_M1"])
    d = sets.merge(sc, on="k").assign(y=lambda z: (z.grp == "case").astype(int)).rename(columns={"logit_M1": "score"})
    c = d[d.y == 1].merge(idx, on=["pid", "k"])
    c = c.assign(before=c.before_first_t2.astype(bool), days=c.days_to_t2)
    c = c[c.days.notna()]                     # preeksite EKG zamanı işlenmiş veride yoksa önce/sonra belirsiz: dışarıda
    return d, c[["set_id", "before", "days"]]


def auc(x, seed):
    if x.y.sum() < 10:
        return {"n_case": int(x.y.sum()), "n_ctrl": int((1 - x.y).sum()), "auroc": None, "ci": None}
    a, ci = S.auc_boot(x.y, x.score, x.set_id, NB, seed)
    return {"n_case": int(x.y.sum()), "n_ctrl": int((1 - x.y).sum()), "auroc": float(a), "ci": [float(ci[0]), float(ci[1])]}


def diff(xa, xb, seed):
    """AUROC(xb) − AUROC(xa); tabakalar içinde ayrı küme bootstrap'ı."""
    rng = np.random.default_rng(seed)
    parts = []
    for x in (xa, xb):
        x = x.reset_index(drop=True)
        parts.append((x.y.to_numpy(), x.score.to_numpy(), list(x.groupby("set_id").indices.values())))
    est = S.auroc(parts[1][0], parts[1][1]) - S.auroc(parts[0][0], parts[0][1])
    vals = []
    for _ in range(NB):
        a = []
        for y, s, ind in parts:
            ix = np.concatenate([ind[j] for j in rng.integers(0, len(ind), len(ind))])
            a.append(S.auroc(y[ix], s[ix]) if 0 < y[ix].sum() < len(ix) else np.nan)
        vals.append(a[1] - a[0])
    return {"diff": float(est), "ci": np.nanpercentile(vals, [2.5, 97.5]).tolist()}


def main():
    R, md = {}, ["# Keşifsel (kayıt sonrası): indeks EKG ilk preeksite EKG'den önce mi, sonra mı?\n",
                 "| Site | Tabaka | Vaka | Kontrol | AUROC (%95 GA) |", "|---|---|---|---|---|"]
    for name, (d, c) in (("MGH", mgh()), ("Emory", test_site("I0006")), ("MIMIC", test_site("MIMIC"))):
        sub = lambda sets: d[d.set_id.isin(sets)]
        before, after = set(c.set_id[c.before]), set(c.set_id[~c.before])
        near = set(c.set_id[~c.before & (c.days <= 30)])
        far = set(c.set_id[~c.before & (c.days > 30)])
        r = {"önce": auc(sub(before), 1), "sonra": auc(sub(after), 2), "sonra, ≤30 gün": auc(sub(near), 3),
             "sonra, >30 gün": auc(sub(far), 4), "fark sonra − önce": diff(sub(before), sub(after), 5)}
        R[name] = r
        for k in ("önce", "sonra", "sonra, ≤30 gün", "sonra, >30 gün"):
            v = r[k]
            a = f"{v['auroc']:.3f} ({v['ci'][0]:.3f}–{v['ci'][1]:.3f})" if v["auroc"] is not None else "— (az vaka)"
            md.append(f"| {name} | {k} | {v['n_case']} | {v['n_ctrl']} | {a} |")
        f = r["fark sonra − önce"]
        md.append(f"| {name} | **fark, sonra − önce** | | | {f['diff']:+.3f} ({f['ci'][0]:+.3f} – {f['ci'][1]:+.3f}) |")
    md.append("\nMGH: M1 OOF tahminleri (geliştirme). Emory ve MIMIC: dondurulmuş M1. Her tabaka kendi eşleşmiş kontrolleriyle.")
    out = ROOT / "stage5" / "out"
    (out / "before_after.md").write_text("\n".join(md) + "\n")
    (out / "before_after.json").write_text(json.dumps(R, indent=1, ensure_ascii=False))
    print("\n".join(md))


if __name__ == "__main__":
    main()
