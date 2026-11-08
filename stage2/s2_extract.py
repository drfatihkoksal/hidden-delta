"""HIDDEN-DELTA Aşama 2 — toplu özellik ve medyan vuru çıkarımı.

İş listesi (HEEDB sitesi için):
  case   herhangi bir T2 tanımına göre vaka olan hastaların (stage1 cohort) 500 Hz T0 ve T2 EKG'lerinin tamamı
  k1     K1 havuzu: hiçbir kaynakta preeksitasyon/kısa PR kodu yok, I45.6/426.7 ve SVT kodu yok, ≥2 EKG;
         rastgele N_K1 hasta × 1 T0 EKG (overread'li, 500 Hz). Eşleştirme Aşama 3'te bu havuzdan yapılır.
  k2     SVT kodu var, I45.6/426.7 yok, preeksitasyon kodu yok; rastgele N_K2 hasta × 1 T0 EKG (500 Hz)

Çıktı (~/data/hidden-delta/stage2/<site>/): parça başına features_XXX.parquet ve beats_XXX.npz
(ids, mb0 = baskın küme medyan vurusu, mb1 = ikinci küme; 12×500, 500 Hz, fiducial −400…+600 ms).
Var olan parçalar atlanır (kaldığı yerden sürer). Çalıştırma: python3 stage2/s2_extract.py I0001
"""
import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage0"))
import hdsig  # noqa: E402
from common import HEEDB, MIMIC_ECG, PRIVATE  # noqa: E402

P1 = Path("/home/fatih/data/hidden-delta/stage1")
P2 = Path("/home/fatih/data/hidden-delta/stage2")
SEED = 20261112
N_K1, N_K2 = 20000, 3000
THIN = {"I0001": 40, "I0006": 4}   # K1 havuzunu deterministik inceltme (Emory küçük)
SHARD = 2000
L = int((hdsig.PRE_S + hdsig.POST_S) * 500)


def read_fs(rec):
    try:
        with open(rec + ".hea") as fh:
            return int(float(fh.readline().split()[2].split("/")[0]))
    except (OSError, IndexError, ValueError):
        return None


def sample_one_per_patient(df, n, rng, rec_of):
    """Hasta başına bir EKG, 500 Hz olanlardan; fs başlıktan okunur."""
    df = df.sample(frac=1.0, random_state=int(rng.integers(2**31))).drop_duplicates("pid")
    keep = []
    for row in df.itertuples():
        if read_fs(rec_of(row.k)) == 500:
            keep.append(row.Index)
            if len(keep) == n:
                break
    return df.loc[keep]


def build_jobs(site):
    rng = np.random.default_rng(SEED)
    con = duckdb.connect()
    rec_of = lambda k: str(HEEDB / site / "WFDB" / k)
    case = con.execute(f"""SELECT pid, k, CASE WHEN t2_cons OR t2_ph OR t2_v24 THEN 'case_t2' ELSE 'case_t0' END AS role
        FROM '{P1 / f'heedb_{site}_ecg_labels.parquet'}'
        WHERE fs = 500 AND (t0_lab OR t2_cons OR t2_ph OR t2_v24)
          AND pid IN (SELECT pid FROM '{P1 / f'heedb_{site}_cohort.parquet'}' WHERE case_cons OR case_ph OR case_v24)""").df()
    flags = f"'{PRIVATE / f'heedb_{site}_ecg_flags.parquet'}'"
    icdpre, icdsvt = f"'{PRIVATE / f'heedb_{site}_icd_pre.parquet'}'", f"'{PRIVATE / f'heedb_{site}_icd_svt.parquet'}'"
    con.execute(f"""CREATE VIEW pp AS SELECT pid, count(*) n_ecg, bool_or(pre) any_pre, bool_or(spr) any_spr
                    FROM {flags} GROUP BY pid""")
    # büyük havuzları deterministik incelt, sonra rastgele hasta başına bir EKG
    k1 = con.execute(f"""SELECT pid, k FROM {flags} WHERE t0 AND has_ph AND (hash(k || '{SEED}') % {THIN[site]}) = 0
        AND pid IN (SELECT pid FROM pp WHERE NOT any_pre AND NOT any_spr AND n_ecg >= 2)
        AND pid NOT IN (SELECT pid FROM {icdpre}) AND pid NOT IN (SELECT pid FROM {icdsvt})""").df()
    k2 = con.execute(f"""SELECT pid, k FROM {flags} WHERE t0 AND has_ph
        AND pid IN (SELECT pid FROM {icdsvt}) AND pid NOT IN (SELECT pid FROM {icdpre})
        AND pid IN (SELECT pid FROM pp WHERE NOT any_pre)""").df()
    k1 = sample_one_per_patient(k1, N_K1, rng, rec_of).assign(role="k1")
    k2 = sample_one_per_patient(k2, N_K2, rng, rec_of).assign(role="k2")
    jobs = pd.concat([case, k1, k2], ignore_index=True).drop_duplicates("k")
    jobs = jobs.sample(frac=1.0, random_state=SEED).reset_index(drop=True)   # parçalara dengeli dağılım
    jobs["rec"] = [rec_of(k) for k in jobs.k]
    return jobs


def build_jobs_mimic():
    """MIMIC: herhangi T2 tanımına göre vaka olan hastaların T0 ve T2 EKG'leri; ablasyon kodlu hastaların T0/T2'leri;
    K1 (hiç preeksitasyon/kısa PR ifadesi yok, I45.6/426.7, SVT ve ablasyon kodu yok, ≥2 EKG) rastgele N_K1 hasta × 1 T0;
    K2 hastalarının (stage1) her biri × 1 T0. MIMIC'in tamamı 500 Hz. Etiketli analiz yapılmaz; yalnızca sinyal."""
    rng = np.random.default_rng(SEED)
    con = duckdb.connect()
    con.execute(f"CREATE VIEW f AS SELECT * FROM '{PRIVATE / 'mimic_ecg_flags.parquet'}'")
    con.execute(f"""CREATE VIEW rl AS SELECT study_id::BIGINT study, path
                    FROM read_csv('{MIMIC_ECG / 'record_list.csv'}', all_varchar=true)""")
    con.execute(f"CREATE VIEW lab AS SELECT * FROM '{P1 / 'mimic_ecg_labels.parquet'}'")
    con.execute(f"CREATE VIEW coh AS SELECT * FROM '{P1 / 'mimic_cohort.parquet'}'")
    for v, fname in (("icdpre", "mimic_icd_pre"), ("icdsvt", "mimic_icd_svt"), ("abl", "mimic_abl")):
        con.execute(f"CREATE VIEW {v} AS SELECT * FROM '{PRIVATE / (fname + '.parquet')}'")
    q = lambda where, role: con.execute(f"""SELECT f.sid AS pid, f.study::VARCHAR AS k, {role} AS role, rl.path
        FROM f JOIN rl USING (study) WHERE {where}""").df()
    case = q("""(f.t0 OR f.pre) AND f.sid IN (SELECT sid FROM coh WHERE case_any)""",
             "CASE WHEN f.pre THEN 'case_t2' ELSE 'case_t0' END")
    abl = q("""(f.t0 OR f.pre) AND f.sid IN (SELECT sid FROM coh WHERE d_abl IS NOT NULL)""",
            "CASE WHEN f.pre THEN 'abl_t2' ELSE 'abl_t0' END")
    con.execute("CREATE VIEW pp AS SELECT sid, count(*) n_ecg, bool_or(pre) any_pre, bool_or(spr) any_spr FROM f GROUP BY sid")
    k1 = q("""f.t0 AND f.sid IN (SELECT sid FROM pp WHERE NOT any_pre AND NOT any_spr AND n_ecg >= 2)
              AND f.sid NOT IN (SELECT sid FROM icdpre) AND f.sid NOT IN (SELECT sid FROM icdsvt)
              AND f.sid NOT IN (SELECT sid FROM abl)""", "'k1'")
    k2 = q("f.t0 AND f.sid IN (SELECT sid FROM coh WHERE k2)", "'k2'")
    one = lambda d, n: d.sample(frac=1.0, random_state=int(rng.integers(2**31))).drop_duplicates("pid").head(n)
    jobs = pd.concat([case, abl[~abl.k.isin(case.k)], one(k1, N_K1), one(k2, None)], ignore_index=True)
    jobs = jobs.drop_duplicates("k").sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    jobs["rec"] = [str(MIMIC_ECG / p) for p in jobs.pop("path")]
    return jobs


def work(rec):
    t = time.time()
    try:
        feat, mbs = hdsig.process_record(rec)
        feat["status"], feat["error"] = "ok", ""
    except Exception as e:  # noqa: BLE001
        feat, mbs = {"status": "fail", "error": f"{type(e).__name__}: {str(e)[:60]}"}, {}
    feat["proc_s"] = round(time.time() - t, 3)
    blank = np.full((12, L), np.nan, np.float32)
    return feat, mbs.get(0, blank), mbs.get(1, blank)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site", choices=["I0001", "I0006", "MIMIC"])
    ap.add_argument("--workers", type=int, default=28)
    a = ap.parse_args()
    out = P2 / a.site
    out.mkdir(parents=True, exist_ok=True)
    jpath = out / "jobs.parquet"
    if jpath.exists():
        jobs = pd.read_parquet(jpath)
    else:
        jobs = build_jobs_mimic() if a.site == "MIMIC" else build_jobs(a.site)
        jobs.to_parquet(jpath, index=False)
    print("iş sayısı:", jobs.role.value_counts().to_dict(), flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        for s in range(0, len(jobs), SHARD):
            fpath, bpath = out / f"features_{s // SHARD:03d}.parquet", out / f"beats_{s // SHARD:03d}.npz"
            if fpath.exists() and bpath.exists():
                continue
            part = jobs.iloc[s:s + SHARD]
            res = list(ex.map(work, part.rec, chunksize=8))
            feats = pd.DataFrame([r[0] for r in res])
            feats.insert(0, "k", part.k.values)
            feats.insert(1, "pid", part.pid.values)
            feats.insert(2, "role", part.role.values)
            np.savez(bpath, ids=part.k.values.astype(str), mb0=np.stack([r[1] for r in res]),
                     mb1=np.stack([r[2] for r in res]))
            feats.to_parquet(fpath, index=False)
            done = min(s + SHARD, len(jobs))
            print(f"{done:,}/{len(jobs):,}  ok %{100 * (feats.status == 'ok').mean():.1f}  "
                  f"{time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
