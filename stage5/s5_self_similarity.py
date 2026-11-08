"""HIDDEN-DELTA — kayıt sonrası referans analizi (2026-11-08): 5.3'ün "aynı kişi" düzeyi.

5.3, indeks EKG'nin ilk 20 ms profilini hastanın kendi preeksite EKG'lerinin ilk 20 ms eğim profiliyle karşılaştırır;
permütasyon başka hastalarla karşılaştırır. Aynı kişinin herhangi iki EKG'si kalp konumu ve elektrot yerleşimi nedeniyle
birbirine benzediğinden, bu referans aynı ölçüyü aynı kişinin iki preeksitasyonsuz EKG'si arasında hesaplar:
- A: indeks EKG'nin ilk 20 ms profili (5.3 ile aynı: `C.initial_profile`, 8 derivasyon);
- B: aynı kişinin indeks dışındaki preeksitasyonsuz EKG'lerinin (§8.4 kümesi) ilk 20 ms eğim profili, EKG başına
  normalize edilip ortalanır (5.3'teki T2 profiliyle aynı);
- ortalama Pearson r ve kişiler arası permütasyon (5.3 ile aynı fonksiyon, 10,000 permütasyon); kontrollerde ve vakalarda.
Yeni skorlama yok; model kullanılmaz. Eğimler `hdsig.beat_features` ile, aynı pipeline koduyla hesaplanır.

Çalıştırma: python3 stage5/s5_self_similarity.py
"""
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import common5 as C
import s5_analysis as A

SLOPE = [f"slope20_{x}" for x in C.L8]


def feats(mb):
    try:
        f = C.hdsig.beat_features(mb.astype(float), C.FS)
        return [f[c] for c in SLOPE]
    except Exception:  # noqa: BLE001
        return [np.nan] * len(SLOPE)


def site(s):
    P = C.P5 / s
    sets = pd.read_parquet(P / "sets.parquet")
    sets = sets[sets.analysis == "primary_k1"]
    idx = pd.read_parquet(P / "index.parquet", columns=["pid", "k", "on_low_ms"])
    mb_idx = np.load(P / "index_mb.npy")
    pos = {k: i for i, k in enumerate(pd.read_parquet(P / "index.parquet", columns=["k"]).k)}
    agg = pd.read_parquet(P / "agg.parquet")
    mb_agg = np.load(P / "agg_mb.npy")
    with ProcessPoolExecutor(28) as ex:
        sl = np.array(list(ex.map(feats, list(mb_agg), chunksize=32)), float)
    agg = pd.concat([agg, pd.DataFrame(sl, columns=SLOPE)], axis=1)
    ix = sets.merge(idx, on=["pid", "k"], how="left")
    other = agg.merge(ix[["pid", "k"]].rename(columns={"k": "k_index"}), on="pid")
    other = other[(other.k != other.k_index)].dropna(subset=SLOPE)
    prof = other.groupby("pid")[SLOPE].apply(
        lambda g: pd.Series(np.nanmean(g.to_numpy() / np.linalg.norm(g.to_numpy(), axis=1, keepdims=True), 0), index=SLOPE))
    R = {"site": s, "note": "kayıt sonrası referans; resmî sonuçlardan sonra yazıldı"}
    rng = np.random.default_rng(C.SEED + 700)
    for grp, lab in (("k1", "controls"), ("case", "cases")):
        g = ix[(ix.grp == grp) & ix.pid.isin(prof.index)]
        Aa = np.array([C.initial_profile(mb_idx[pos[k]], o) for k, o in zip(g.k, g.on_low_ms)])
        obs, null, p = A.perm_agreement(Aa, prof.loc[g.pid].to_numpy(), rng)
        R[lab] = {"n": len(g), "mean_r": obs, "perm_mean_r": null, "p_one_sided": p,
                  "n_other_ecg_median": float(other[other.pid.isin(set(g.pid))].groupby("pid").size().median())}
    a = json.loads((P / "analysis.json").read_text())["H2"]["5.3_delta"]
    R["registered_5_3"] = {"n": a["n"], "mean_r": a["mean_r"], "perm_mean_r": a["perm_mean_r"]}
    (P / "self_similarity.json").write_text(json.dumps(R, indent=1, ensure_ascii=False))
    return R


def main():
    res = {s: site(s) for s in ("MIMIC", "I0006")}
    md = ["# Kayıt sonrası referans: 5.3 ölçüsü aynı kişinin iki preeksitasyonsuz EKG'si arasında\n",
          "| Site | Karşılaştırma | n | Ortalama r | Permütasyon r |", "|---|---|---|---|---|"]
    for s, R in res.items():
        n = C.SITE_NAME[s]
        md.append(f"| {n} | Vaka: indeks T0 – kendi T2'si (kayıtlı 5.3) | {R['registered_5_3']['n']} | "
                  f"{R['registered_5_3']['mean_r']:.3f} | {R['registered_5_3']['perm_mean_r']:.3f} |")
        for lab, t in (("cases", "Vaka: indeks T0 – kendi diğer T0'ları"), ("controls", "Kontrol: indeks T0 – kendi diğer T0'ları")):
            md.append(f"| {n} | {t} | {R[lab]['n']} | {R[lab]['mean_r']:.3f} | {R[lab]['perm_mean_r']:.3f} |")
    (C.OUT / "self_similarity.md").write_text("\n".join(md) + "\n")
    out = {s: R for s, R in res.items()}
    (C.OUT / "self_similarity.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print("\n".join(md))


if __name__ == "__main__":
    main()
