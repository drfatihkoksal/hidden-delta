"""HIDDEN-DELTA Aşama 4 — M2: M1 kodlayıcı gömmesi + M0 özellikleri, L2 lojistik (ANALYSIS_PLAN §4). Yalnızca MGH.

M1 gömmesi: M1'in seçilen derinliğindeki SSL kodlayıcısının 128 boyutlu ortalama havuzlanmış çıktısı. Kat başına o katın
SSL kodlayıcısı kullanılır (SSL havuzu katın test hastalarını içermez; etiket kullanmaz). Final: tüm MGH SSL kodlayıcısı
(M1 final ile aynı tohum) + tüm MGH ile lojistik. Eksik M0 özelliği olan hastalar çıkar (M0 ile aynı kural).

Çıktı: stage4/out/m2_mgh.md; ~/data/hidden-delta/stage4/final/M2.json, m1/final_ssl_encoder_m2.pt, I0001_oof_m2.parquet.
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

import m0lib as M
import s4_m1 as M1
import statlib as S

warnings.filterwarnings("ignore", category=FutureWarning)
P4 = Path("/home/fatih/data/hidden-delta/stage4")
OUT = Path(__file__).resolve().parent / "out"
N_BOOT = 2000


def emb_design(Z, df, prm):
    return np.hstack([Z, M.design(df, prm)])


def main():
    M1.BANK = M1.filter_bank()
    rng = np.random.default_rng(M.SEED + 2)
    meta = json.loads((P4 / "m1" / "final_meta.json").read_text())
    depth = meta["config"]["depth"]
    dev = pd.read_parquet(P4 / "I0001_dev.parquet")
    Xd = np.load(P4 / "I0001_dev_mb.npy").astype(np.float32)
    cont, binv = M.MODELS["M0"]
    ok = dev[cont + binv].notna().all(axis=1).to_numpy()
    pred = np.full(len(dev), np.nan)
    for f in range(5):
        enc = M1.Encoder(depth).to(M1.DEV)
        enc.load_state_dict(torch.load(P4 / "m1" / f"ssl_fold{f}_d{depth}.pt"))
        Z = M1.predict(enc, Xd, embed=True)
        tr, te = ok & (dev.fold.to_numpy() != f), ok & (dev.fold.to_numpy() == f)
        prm = M.fit_design(dev[tr], cont, binv)
        Ztr, Zte = Z[tr], Z[te]
        zm, zs = Ztr.mean(0), np.where(Ztr.std(0) > 0, Ztr.std(0), 1)
        m = M.fit_lr(emb_design((Ztr - zm) / zs, dev[tr], prm), dev.y[tr])
        pred[te] = m.decision_function(emb_design((Zte - zm) / zs, dev[te], prm))
        print(f"kat {f} bitti", flush=True)
    d = dev[ok].assign(p=pred[ok])
    oof_m0 = pd.read_parquet(P4 / "I0001_oof_m0.parquet")[["k", "p_M0"]]
    oof_m1 = pd.read_parquet(P4 / "m1" / "oof_m1.parquet")[["k", "p_M1"]]
    d = d.merge(oof_m0, on="k").merge(oof_m1, on="k")
    rows = []
    for i, (col, name) in enumerate((("p", "M2"), ("p_M1", "M1"), ("p_M0", "M0"))):
        auc, ci = S.auc_boot(d.y, d[col], d.set_id, N_BOOT, M.SEED + i)
        rows.append({"model": name, "hasta (ortak)": len(d), "OOF AUROC": round(auc, 3), "%95 GA": f"{ci[0]:.3f}–{ci[1]:.3f}"})
    dd, dci, _ = S.auc_diff_boot(d.y, d.p, d.p_M1, d.set_id, N_BOOT, M.SEED + 9)
    print("bootstrap bitti", flush=True)
    # final M2
    if meta["config"]["mode"] == "probe":                     # M1 final kodlayıcısı zaten SSL kodlayıcısı: aynısı kullanılır
        enc = M1.Encoder(depth).to(M1.DEV)
        enc.load_state_dict(torch.load(P4 / "m1" / "final_encoder.pt"))
    else:
        enc = M1.pretrain(np.asarray(np.load(P4 / "I0001_ssl_mb.npy", mmap_mode="r"), dtype=np.float32), depth,
                          M1.SEED + 100 + depth)
        torch.save(enc.state_dict(), P4 / "m1" / "final_ssl_encoder_m2.pt")
    Z = M1.predict(enc, Xd, embed=True)[ok]
    prm = M.fit_design(dev[ok], cont, binv)
    zm, zs = Z.mean(0), np.where(Z.std(0) > 0, Z.std(0), 1)
    prm["z_mean"], prm["z_scale"], prm["depth"] = zm.tolist(), zs.tolist(), depth
    prm["lr"] = M.lr_params(M.fit_lr(emb_design((Z - zm) / zs, dev[ok], prm), dev.y[ok]))
    (P4 / "final" / "M2.json").write_text(json.dumps(prm, indent=1))
    d[["set_id", "pid", "k", "y", "fold", "p"]].rename(columns={"p": "p_M2"}).to_parquet(P4 / "I0001_oof_m2.parquet")
    md = ["# M2 — MGH 5 katlı çapraz doğrulama (geliştirme; hipotez testi değil)\n",
          f"M1 gömmesi: derinlik {depth} SSL kodlayıcısı (128 boyut) + M0 özellikleri; ortak hastalarda karşılaştırma.\n",
          pd.DataFrame(rows).to_markdown(index=False),
          f"\nM2 − M1 AUROC farkı: {dd:+.3f} "
          f"(küme bootstrap %95 GA {dci[0]:+.3f} – {dci[1]:+.3f}).\n",
          "Final M2 tüm MGH ile eğitildi (`final/M2.json`; kodlayıcı: M1 probe ise `m1/final_encoder.pt`)."]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "m2_mgh.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
