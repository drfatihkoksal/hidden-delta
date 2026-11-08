"""HIDDEN-DELTA — önceki IAB çalışmalarında okunmuş EKG'lerin hastaları (test setlerinden dışlama için).

Gözlenen QRS öncesi dalga, yazarın Emory ve MIMIC EKG'lerini okuduğu IAB çalışmalarında görüldü (AUTO-IAB,
https://osf.io/s7tgr; LEAD-IAB okuma setleri; P-SPLIT Emory P işaret seti). Bu hastalar Emory ve MIMIC test
setlerinden eşleştirmeden önce çıkarılır (ANALYSIS_PLAN v1.2, §13).

Kaynaklar: okuma anahtarları (sentetik olanlar hariç). Çıktı: ~/data/hidden-delta/stage3/iab_read_patients.parquet
(dataset, pid) ve agregat sayılar. Çalıştırma: python3 stage3/s3_iab_overlap.py
"""
import re
from pathlib import Path

import duckdb
import pandas as pd

H = Path.home() / "projects"
P0 = Path("/home/fatih/data/hidden-delta/stage0")
P3 = Path("/home/fatih/data/hidden-delta/stage3")
MIMIC_KEYS = [f"LEAD-IAB/out/{k}/key.csv" for k in
              ("adjudication", "validation_x", "validation_v3", "round1_reread", "reread_v3", "review_tp_only",
               "validation_v4", "reread")]
EMORY_KEYS = ["P-SPLIT/out/reading_emory/gizli_anahtar/key.csv", "AUTO-IAB/out/validation/gizli_anahtar/key.csv",
              "AUTO-IAB/out/b5_reread/gizli_anahtar/key.csv", "AUTO-IAB/out/calibration_practice/gizli_anahtar/key.csv"]


def nk(f):
    f = re.sub(r"\s", "", str(f))
    f = re.sub(r"^\./|^/", "", f)
    f = re.sub(r"^.*?WFDB/", "", f)
    return re.sub(r"\.hea$|\.dat$", "", f)


def main():
    studies = set()
    for p in MIMIC_KEYS:
        studies |= set(pd.read_csv(H / p).study_id.astype("int64"))
    files = set()
    for p in EMORY_KEYS:
        files |= {nk(x) for x in pd.read_csv(H / p).file}
    con = duckdb.connect()
    fm = con.execute(f"SELECT sid AS pid, study FROM '{P0 / 'mimic_ecg_flags.parquet'}'").df()
    fe = con.execute(f"SELECT pid, k FROM '{P0 / 'heedb_I0006_ecg_flags.parquet'}'").df()
    pm = set(fm[fm.study.isin(studies)].pid)
    pe = set(fe[fe.k.isin(files)].pid)
    out = pd.concat([pd.DataFrame({"dataset": "MIMIC", "pid": sorted(pm)}),
                     pd.DataFrame({"dataset": "I0006", "pid": sorted(pe)})], ignore_index=True)
    P3.mkdir(parents=True, exist_ok=True)
    out.to_parquet(P3 / "iab_read_patients.parquet", index=False)
    print(f"MIMIC: okunan çalışma {len(studies):,}, hasta {len(pm):,} | Emory: okunan EKG {len(files):,} "
          f"(flags'te bulunan {int(fe.k.isin(files).sum()):,}), hasta {len(pe):,}")


if __name__ == "__main__":
    main()
