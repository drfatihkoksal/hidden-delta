"""HIDDEN-DELTA Aşama 0 — ortak yollar, ifade kalıpları ve kod kümeleri.

DUA kuralı: bu klasöre (stage0/out) yalnızca agregat sayılar yazılır.
Hasta düzeyindeki ara listeler PRIVATE altına gider ve makineden çıkmaz.
"""
from pathlib import Path

import duckdb

V = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c")

# MIMIC-IV 3.1, MIMIC-IV-ECG 1.0, MIMIC-IV-Note
MIMIC_ECG = V / "mimic-iv-ecg/1.0"
MIMIC_HOSP = V / "mimiciv/3.1/hosp"
MIMIC_NOTE = V / "mimicNote/mimic-iv-note/note/discharge.csv"

# HEEDB 5.0 (I0001 = MGH, I0006 = Emory)
HEEDB = Path("/home/fatih/datasets/heedb/ECG")
HEEDB_SITES = {"I0001": "MGH", "I0006": "Emory"}

OUT = Path(__file__).resolve().parent / "out"          # agregat
PRIVATE = Path("/home/fatih/data/hidden-delta/stage0")  # hasta düzeyi, paylaşılmaz
DUCKDB_TMP = Path("/home/fatih/data/hidden-delta/duckdb_tmp")

# --- MIMIC makine raporu kalıpları (büyük harfe çevrilmiş rapor metninde, RE2) ---
RX_PRE = r"WOLFF|WPW|PRE-?EXCIT|DELTA WAVE"          # preeksitasyon ifadesi (katı)
RX_SPR = r"SHORT PR"                                  # kısa PR (çoğu WPW değil; ayrı tutulur)
RX_SINUS = r"SINUS"
RX_AFL = r"ATRIAL FIBRILLATION|FLUTTER"
RX_PACED = r"PACE|PACING|PACED|PACER"
RX_WIDE = r"BUNDLE BRANCH BLOCK|IVCD|INTRAVENTRICULAR (CONDUCTION|BLOCK)"
RX_INCOMPLETE_BBB = r"INCOMPLETE (RIGHT|LEFT) BUNDLE BRANCH BLOCK"

# --- 12SL kodları (HEEDB, diagnoses_dictionary.csv; iki merkezde aynı sözlük) ---
C_PRE = [300, 302, 303, 304]   # WPWA, WPWB, ALTWPW, WPW
C_ALT = [303]                  # "with fusion or intermittent ventricular pre-excitation (wpw)"
C_SPR = [102]                  # "with short pr"
C_SINUS = [19, 21, 22, 23, 24]
C_AFL = [161, 162, 273, 288]
C_PACED = [183, 184, 185, 186, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 326, 1669]
C_WIDE = [440, 442, 460, 478, 479, 480, 481, 482, 487]

# --- ICD ---
# Tanı: preeksitasyon I45.6 / 426.7; SVT I47.1x (I47.11 hariç) / 427.0
# Ablasyon (yolağa özgü değil): ICD-10-PCS 02583ZZ, ICD-9-CM 37.34
ABL_ICD10 = "02583ZZ"
ABL_ICD9 = "3734"

AGE_MIN = 18
POST_ABL_EARLY_DAYS = 90   # ablasyon sonrası 0–3 ay: kardiyak bellek olası


def connect(threads: int = 16, memory: str = "40GB") -> duckdb.DuckDBPyConnection:
    DUCKDB_TMP.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"SET threads={threads}; SET memory_limit='{memory}'; SET temp_directory='{DUCKDB_TMP}'")
    con.execute("SET preserve_insertion_order=false")
    return con


def md_table(rows, headers) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for r in rows:
        out.append("| " + " | ".join(f"{x:,}" if isinstance(x, int) else str(x) for x in r) + " |")
    return "\n".join(out)
