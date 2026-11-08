"""HIDDEN-DELTA — kayıt sonrası keşifsel analiz (2026-11-08): cinsiyete göre AUROC.

Kümeler cinsiyet üzerinde kesin eşleştirildiği için her küme tek cinsiyetlidir; AUROC kadın ve erkek kümelerde ayrı
hesaplanır, fark (kadın − erkek) tabakalar içinde ayrı küme bootstrap'ıyla (2,000). Cinsiyet kaynak veritabanlarındaki
kayıtlı değerdir (HEEDB meta verisi, MIMIC-IV `patients.gender`). MGH: M1 OOF (geliştirme); Emory ve MIMIC: dondurulmuş M1.
Yeni skorlama yok.

Çalıştırma: python3 stage5/s5_by_sex.py
"""
import json
from pathlib import Path

import pandas as pd

import s5_before_after as B

ROOT = Path(__file__).resolve().parents[1]
D = B.D


def mgh():
    dev = pd.read_parquet(D / "stage4/I0001_dev.parquet", columns=["k", "y", "set_id", "male"])
    oof = pd.read_parquet(D / "stage4/m1/oof_m1.parquet", columns=["k", "oof_d8_probe"])
    return dev.merge(oof, on="k").rename(columns={"oof_d8_probe": "score"})


def test_site(s):
    sets = pd.read_parquet(D / "stage5" / s / "sets.parquet")
    sets = sets[sets.analysis == "primary_k1"]
    idx = pd.read_parquet(D / "stage5" / s / "index.parquet", columns=["pid", "k", "male"])
    sc = pd.read_parquet(D / "stage5" / s / "scores.parquet", columns=["k", "logit_M1"])
    d = sets.merge(idx, on=["pid", "k"]).merge(sc, on="k").assign(y=lambda z: (z.grp == "case").astype(int))
    return d.rename(columns={"logit_M1": "score"})


def main():
    R, md = {}, ["# Keşifsel (kayıt sonrası): cinsiyete göre AUROC\n", "| Site | Grup | Vaka | Kontrol | AUROC (%95 GA) |",
                 "|---|---|---|---|---|"]
    for name, d in (("MGH", mgh()), ("Emory", test_site("I0006")), ("MIMIC", test_site("MIMIC"))):
        sex = d.groupby("set_id").male.transform("first")
        f, m = d[sex == 0], d[sex == 1]
        r = {"kadın": B.auc(f, 1), "erkek": B.auc(m, 2), "fark kadın − erkek": B.diff(m, f, 3)}
        R[name] = r
        for k in ("kadın", "erkek"):
            v = r[k]
            md.append(f"| {name} | {k} | {v['n_case']} | {v['n_ctrl']} | {v['auroc']:.3f} ({v['ci'][0]:.3f}–{v['ci'][1]:.3f}) |")
        x = r["fark kadın − erkek"]
        md.append(f"| {name} | **fark, kadın − erkek** | | | {x['diff']:+.3f} ({x['ci'][0]:+.3f} – {x['ci'][1]:+.3f}) |")
    out = ROOT / "stage5" / "out"
    (out / "by_sex.md").write_text("\n".join(md) + "\n")
    (out / "by_sex.json").write_text(json.dumps(R, indent=1, ensure_ascii=False))
    print("\n".join(md))


if __name__ == "__main__":
    main()
