"""HIDDEN-DELTA Aşama 4 — M0 (özellik modeli) ve M-demo (yaş, cinsiyet, KH); MGH 5 katlı çapraz doğrulama + final.

M0: ANALYSIS_PLAN §4.1'deki sabit özellikler + yaş + cinsiyet (`m0lib`). Eksik değerler doldurulmaz: eksik özellikli
hastalar o modelden çıkar (sayı raporlanır). Değerlendirme: kat dışı (OOF) tahminlerle hasta düzeyinde AUROC;
%95 GA eşleştirilmiş küme bootstrap'ı (2,000). Final: tüm MGH ile; tasarım ve katsayılar JSON'a yazılır.

Çıktı: stage4/out/m0_mgh.md (agregat); ~/data/hidden-delta/stage4/I0001_oof_m0.parquet, final/{M0,M-demo}.json.
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

import m0lib as M

warnings.filterwarnings("ignore", category=FutureWarning)
P4 = Path("/home/fatih/data/hidden-delta/stage4")
OUT = Path(__file__).resolve().parent / "out"
N_BOOT = 2000


def oof(dev, cont, binv):
    pred = pd.Series(np.nan, index=dev.index)
    for f in sorted(dev.fold.unique()):
        tr, te = dev[dev.fold != f], dev[dev.fold == f]
        prm = M.fit_design(tr, cont, binv)
        m = M.fit_lr(M.design(tr, prm), tr.y)
        pred[te.index] = m.decision_function(M.design(te, prm))
    return pred


def boot(d, fn, rng):
    g = {s: x for s, x in d.groupby("set_id")}
    sets = d.set_id.unique()
    out = []
    for _ in range(N_BOOT):                       # eşleştirilmiş küme bootstrap'ı (yerine koyarak)
        b = pd.concat([g[s] for s in rng.choice(sets, len(sets))])
        if b.y.nunique() == 2:
            out.append(fn(b))
    return np.percentile(out, [2.5, 97.5])


def main():
    rng = np.random.default_rng(M.SEED)
    dev = pd.read_parquet(P4 / "I0001_dev.parquet")
    (P4 / "final").mkdir(parents=True, exist_ok=True)
    res, out = [], dev[["set_id", "pid", "k", "y", "fold"]].copy()
    for name, (cont, binv) in M.MODELS.items():
        ok = dev[cont + binv].notna().all(axis=1)
        d = dev[ok].copy()
        d["p"] = oof(d, cont, binv)
        out.loc[d.index, f"p_{name}"] = d.p
        ci = boot(d, lambda b: roc_auc_score(b.y, b.p), rng)
        res.append({"model": name, "hasta": len(d), "vaka": int(d.y.sum()), "çıkan (eksik özellik)": int((~ok).sum()),
                    "AUROC": round(roc_auc_score(d.y, d.p), 3), "%95 GA": f"{ci[0]:.3f}–{ci[1]:.3f}"})
        prm = M.fit_design(d, cont, binv)                          # final: tüm MGH
        prm["lr"] = M.lr_params(M.fit_lr(M.design(d, prm), d.y))
        (P4 / "final" / f"{name}.json").write_text(json.dumps(prm, indent=1))
    both = out.dropna(subset=["p_M0", "p_M-demo"])
    diff = roc_auc_score(both.y, both["p_M0"]) - roc_auc_score(both.y, both["p_M-demo"])
    ci = boot(both, lambda b: roc_auc_score(b.y, b["p_M0"]) - roc_auc_score(b.y, b["p_M-demo"]), rng)
    out.to_parquet(P4 / "I0001_oof_m0.parquet", index=False)
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# M0 ve M-demo — MGH 5 katlı çapraz doğrulama (geliştirme; hipotez testi değil)\n",
          pd.DataFrame(res).to_markdown(index=False),
          f"\nM0 − M-demo AUROC farkı (ortak {len(both):,} hasta): {diff:+.3f}, küme bootstrap %95 GA "
          f"{ci[0]:+.3f} – {ci[1]:+.3f}.\n",
          "Final modeller tüm MGH ile eğitildi (`~/data/hidden-delta/stage4/final/`).\n",
          "Not: bu sayılar geliştirme setindendir; birincil test MIMIC'te dondurulmuş modelle yapılır (ANALYSIS_PLAN §1)."]
    (OUT / "m0_mgh.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
