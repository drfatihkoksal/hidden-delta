"""HIDDEN-DELTA Aşama 4 — MGH geliştirme veri setinin hazırlanması (ANALYSIS_PLAN §3–4).

Birim: eşleştirilmiş kümelerdeki hasta (vaka + K1), indeks EKG başına bir satır. Katlar: 5, eşleştirilmiş kümeye göre
(küme bölünmez), sabit tohum. Ayrık dalga bayrağı (§6) medyan vurudan hesaplanır.
SSL havuzu: MGH'de çıkarılmış tüm başarılı vaka-T0, K1 ve K2 EKG'leri (T2 EKG'leri hiçbir aşamada kullanılmaz;
concept §7). Her EKG hastasının katı taşınır; M1 kendi test katındaki hastaları SSL'den çıkarır.

Çıktı (~/data/hidden-delta/stage4/): I0001_dev.parquet, I0001_dev_mb.npy, I0001_ssl.parquet, I0001_ssl_mb.npy.
Çalıştırma: python3 stage4/s4_data.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "stage2"))
sys.path.insert(0, str(ROOT / "stage3"))
import hdsig  # noqa: E402
from s3_match import load_heedb  # noqa: E402

SITE = "I0001"
P2 = Path(f"/home/fatih/data/hidden-delta/stage2/{SITE}")
P3 = Path("/home/fatih/data/hidden-delta/stage3")
P4 = Path("/home/fatih/data/hidden-delta/stage4")
SEED = 20261117
K = 5
FS = 500


def load_beats(keys):
    """İstenen EKG'lerin baskın küme medyan vurusu (N×12×500), istenen sırayla."""
    want = {k: i for i, k in enumerate(keys)}
    out = np.full((len(keys), 12, 500), np.nan, np.float32)
    for f in sorted(P2.glob("beats_*.npz")):
        z = np.load(f)
        ids = z["ids"]
        sel = [(j, want[k]) for j, k in enumerate(ids) if k in want]
        if sel:
            src, dst = zip(*sel)
            out[list(dst)] = z["mb0"][list(src)]
    return out


def wave_row(mb, on_low_ms):
    on = int(round(hdsig.PRE_S * FS + on_low_ms * FS / 1000))
    flag, d = hdsig.discrete_wave(mb.astype(float), FS, on)
    return dict(wave=float(flag), wave_amp=d.get("amp", np.nan), wave_fwhm=d.get("fwhm", np.nan),
                wave_peak_ms=d.get("peak_ms", np.nan), wave_nleads=d.get("n_leads", 0))


def main():
    P4.mkdir(parents=True, exist_ok=True)
    df, _ = load_heedb(SITE)
    m = pd.read_parquet(P3 / f"{SITE}_matched.parquet")
    m = m[m.pool == "k1"]
    dev = m.merge(df, on=["pid", "k"], how="left")
    dev["y"] = (dev.grp == "case").astype(int)
    sets = np.array(sorted(dev.set_id.unique()))
    rng = np.random.default_rng(SEED)
    fold_of = dict(zip(rng.permutation(sets), np.arange(len(sets)) % K))
    dev["fold"] = dev.set_id.map(fold_of)
    dev = dev.sort_values(["fold", "set_id", "y", "k"], ascending=[True, True, False, True]).reset_index(drop=True)
    mb = load_beats(dev.k.tolist())
    assert not np.isnan(mb[:, :, 0]).any(), "eksik medyan vuru"
    dev = pd.concat([dev, pd.DataFrame([wave_row(mb[i], dev.on_low_ms[i]) for i in range(len(dev))])], axis=1)
    dev["septal_q_absent"] = (dev[["q_I", "q_aVL", "q_V5", "q_V6"]] == 0).all(axis=1).astype(float)
    dev["v1_rs_c"] = dev.v1_r / np.maximum(dev.v1_s, 0.01)          # S yoksa tanımlı kalsın (impütasyon değil)
    dev["male"] = (dev.sex == "Male").astype(float)
    keep = [c for c in dev.columns if not c.startswith("c2_")]
    dev[keep].to_parquet(P4 / f"{SITE}_dev.parquet", index=False)
    np.save(P4 / f"{SITE}_dev_mb.npy", mb)

    # SSL havuzu: T2 hariç tüm başarılı EKG'ler; hasta katı (dev'de yoksa −1)
    ssl = df[df.role.isin(["case_t0", "k1", "k2"]) & ~df.t2_prim][["pid", "k", "role"]]
    pid_fold = dev.drop_duplicates("pid").set_index("pid").fold
    ssl = ssl.assign(fold=ssl.pid.map(pid_fold).fillna(-1).astype(int)).sort_values("k").reset_index(drop=True)
    np.save(P4 / f"{SITE}_ssl_mb.npy", load_beats(ssl.k.tolist()))
    ssl.to_parquet(P4 / f"{SITE}_ssl.parquet", index=False)

    print(f"geliştirme: {len(dev):,} hasta ({dev.y.sum():,} vaka, {(1 - dev.y).sum():,} kontrol), "
          f"{dev.set_id.nunique():,} küme; kat başına vaka {dev.groupby('fold').y.sum().tolist()}")
    print(f"ayrık dalga: vaka %{100 * dev.wave[dev.y == 1].mean():.1f}, kontrol %{100 * dev.wave[dev.y == 0].mean():.1f}")
    print(f"SSL havuzu: {len(ssl):,} EKG; rol {ssl.role.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
