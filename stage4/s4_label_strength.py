"""HIDDEN-DELTA — keşifsel (kayıt sonrası, 2026-11-08): vaka etiketinin gücü, sıkı T0 ve AUROC.

Soru: dış testlerin AUROC'u MGH'den yüksek; vaka etiketlerinin daha kesin olması bunu açıklıyor mu, yoksa daha kesin
etiketli vakalarda sınırda (algoritmanın kaçırdığı minimal preeksitasyonlu) T0'lar mı daha sık?
- Vaka hastası başına etiket gücü: herhangi bir konsensüs T2 EKG'sinde üç kaynağın üçü de preeksitasyon diyor mu (3/3);
  ≥2 konsensüs T2 EKG var mı.
- MGH ve Emory vakalarında bu oranlar.
- Birleştirilmiş OOF AUROC ile kat içi AUROC ortalaması (MGH).
- AUROC, etiket gücü × sıkı T0 tabakalarında: MGH (M1 OOF, derinlik 8 probe; geliştirme) ve Emory (dondurulmuş M1;
  resmî skorlar, kayıt sonrası keşifsel). Tabaka: o tabakadaki vakalar ve eşleşmiş kontrolleri. Sıkı T0 (plan 7.1,
  HEEDB): QRS < 110 ms ve PR − 42 ≥ 140 ms; s5_analysis ile aynı biçimde vakası ve ≥1 kontrolü sıkı olan kümeler.

Çalıştırma: python3 stage4/s4_label_strength.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import statlib as S

D = Path("/home/fatih/data/hidden-delta")
OUT = Path(__file__).resolve().parent / "out"
NB = 2000
J = {}                                                  # JSON çıktısı (agregat sonuçlar)


def strength(site, case_pids):
    e = pd.read_parquet(D / "stage1" / f"heedb_{site}_ecg_labels.parquet", columns=["pid", "t2_cons", "pre_v24", "pre_sw", "pre_ph"])
    e = e[e.t2_cons.fillna(False).astype(bool) & e.pid.isin(case_pids)].copy()
    e["n_src"] = e[["pre_v24", "pre_sw", "pre_ph"]].fillna(False).astype(int).sum(1)
    g = e.groupby("pid").agg(n_t2=("n_src", "size"), max_src=("n_src", "max"))
    g["all3"], g["multi"] = g.max_src == 3, g.n_t2 >= 2
    return g


def auc_row(x, lab, seed, site=None):
    if x.y.sum() < 10:
        return f"| {lab} | {int(x.y.sum())} | {int((1 - x.y).sum())} | — (az vaka) |"
    a, ci = S.auc_boot(x.y, x.score, x.set_id, NB, seed)
    if site:
        J.setdefault(site, {})[lab] = {"n_case": int(x.y.sum()), "n_ctrl": int((1 - x.y).sum()), "auroc": float(a),
                                       "ci": [float(ci[0]), float(ci[1])]}
    return f"| {lab} | {int(x.y.sum())} | {int((1 - x.y).sum())} | {a:.3f} ({ci[0]:.3f}–{ci[1]:.3f}) |"


def strict_sets(x):
    s = x[x.strict]
    return s[s.groupby("set_id").y.transform("sum").eq(1) & s.groupby("set_id").y.transform("count").gt(1)]


def grid(o, g, title, site):
    """o: y, set_id, pid, score, strict. g: vaka etiket gücü."""
    o = o.merge(g[["all3"]], left_on="pid", right_index=True, how="left")
    md = [f"\n### {title}\n", "| Tabaka | Vaka | Kontrol | AUROC (%95 GA, küme bootstrap) |", "|---|---|---|---|"]
    case_all3 = o.loc[o.y == 1].set_index("set_id").all3
    for lab, sets in (("Tüm vakalar", set(o.set_id)), ("3/3 kaynak", set(case_all3.index[case_all3 == True])),  # noqa: E712
                      ("2/3 kaynak", set(case_all3.index[case_all3 == False]))):  # noqa: E712
        x = o[o.set_id.isin(sets)]
        md.append(auc_row(x, lab, 1, site))
        md.append(auc_row(strict_sets(x), f"{lab}, sıkı T0", 2, site))
    return md


def main():
    dev = pd.read_parquet(D / "stage4" / "I0001_dev.parquet")
    oof = pd.read_parquet(D / "stage4" / "m1" / "oof_m1.parquet").dropna(subset=["oof_d8_probe"])
    gm = strength("I0001", set(dev.pid[dev.y == 1]))
    P5 = D / "stage5" / "I0006"
    st = pd.read_parquet(P5 / "sets.parquet")
    st = st[st.analysis == "primary_k1"]
    ge = strength("I0006", set(st.pid[st.grp == "case"]))
    md = ["# Keşifsel: vaka etiketinin gücü, sıkı T0 ve AUROC (kayıt sonrası)\n",
          "| Site | Vaka | 3/3 kaynak | ≥2 konsensüs T2 |", "|---|---|---|---|"]
    md += [f"| {n} | {len(g)} | {g.all3.mean():.1%} | {g.multi.mean():.1%} |" for n, g in (("MGH", gm), ("Emory", ge))]
    pooled = S.auroc(oof.y, oof.oof_d8_probe)
    fold = np.mean([S.auroc(x.y, x.oof_d8_probe) for _, x in oof.groupby("fold")])
    md += ["", f"MGH OOF AUROC: birleştirilmiş {pooled:.3f}, kat içi ortalama {fold:.3f}.", "",
           "| MGH tabakası | Vaka | AUROC (%95 GA, küme bootstrap) |", "|---|---|---|"]
    o = oof.merge(gm[["all3", "multi"]], left_on="pid", right_index=True, how="left")
    for col, lab in (("all3", "3/3 kaynak"), ("multi", "≥2 konsensüs T2")):
        for v in (True, False):
            x = o[o.set_id.isin(set(o.set_id[(o.y == 1) & (o[col] == v)]))]
            a, ci = S.auc_boot(x.y, x.oof_d8_probe, x.set_id, NB, 1)
            md.append(f"| {lab}: {'evet' if v else 'hayır'} | {int(x.y.sum())} | {a:.3f} ({ci[0]:.3f}–{ci[1]:.3f}) |")

    # ---- etiket gücü × sıkı T0
    m = oof[["k", "pid", "set_id", "y", "oof_d8_probe"]].merge(dev[["k", "qrs_ms", "pr_ms"]], on="k", how="left")
    m = m.assign(score=m.oof_d8_probe, strict=(m.qrs_ms < 110) & ((m.pr_ms - 42) >= 140))
    md += grid(m, gm, "MGH — M1 OOF (geliştirme)", "MGH")
    idx = pd.read_parquet(P5 / "index.parquet", columns=["pid", "k", "strict_t0"])
    sc = pd.read_parquet(P5 / "scores.parquet", columns=["k", "logit_M1"])
    e = st.merge(idx, on=["pid", "k"], how="left").merge(sc, on="k", how="left")
    e = e.assign(y=(e.grp == "case").astype(int), score=e.logit_M1, strict=e.strict_t0.fillna(False).astype(bool))
    md += grid(e, ge, "Emory — dondurulmuş M1 (test seti; kayıt sonrası keşifsel)", "Emory")
    J["share_all3"] = {"MGH": float(gm.all3.mean()), "Emory": float(ge.all3.mean())}
    J["n_cases"] = {"MGH": len(gm), "Emory": len(ge)}
    md.append(f"\nSıkı T0 oranı (vaka / kontrol): MGH {m.strict[m.y == 1].mean():.1%} / {m.strict[m.y == 0].mean():.1%}; "
              f"Emory {e.strict[e.y == 1].mean():.1%} / {e.strict[e.y == 0].mean():.1%}.")
    for n, x, g in (("MGH", m, gm), ("Emory", e, ge)):
        c = x[x.y == 1].merge(g[["all3"]], left_on="pid", right_index=True)
        md.append(f"{n} vakalarında sıkı T0 oranı: 3/3 {c.strict[c.all3].mean():.1%}, 2/3 {c.strict[~c.all3].mean():.1%}.")
    OUT.mkdir(exist_ok=True)
    (OUT / "label_strength.md").write_text("\n".join(md) + "\n")
    (OUT / "label_strength.json").write_text(json.dumps(J, indent=1, ensure_ascii=False))
    print("\n".join(md))


if __name__ == "__main__":
    main()
