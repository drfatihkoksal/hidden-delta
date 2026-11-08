"""HIDDEN-DELTA Aşama 0 — HEEDB sayımları (I0001 = MGH, I0006 = Emory).

Kaynaklar: 12SL_diagnoses (v24 yeniden işlenmiş + kayıt anındaki yazılım ve hekim kodları),
metadata (hasta, zaman, yaş), ICD_codes (I45.6/426.7, I47.1/427.0), WFDB başlıkları (yalnızca adaylar için fs).
HEEDB'de prosedür kodu yok → ablasyon (H3) burada sayılamaz.

Çıktı: stage0/out/heedb_stage0.{md,json} (agregat) ve PRIVATE/heedb_*.parquet (hasta düzeyi).
Çalıştırma: python3 stage0/s0_heedb.py
"""
import json
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from common import (AGE_MIN, C_AFL, C_ALT, C_PACED, C_PRE, C_SINUS, C_SPR, C_WIDE, HEEDB,
                    HEEDB_SITES, OUT, PRIVATE, connect, md_table)


def one(con, sql):
    return con.execute(sql).fetchone()


def has(col, codes):
    return f"coalesce(list_has_any({col}, {codes}), false)"


def read_fs(path):
    try:
        with open(path) as f:
            parts = f.readline().split()
        return parts[2].split("/")[0] if len(parts) > 2 else "parse_error"
    except FileNotFoundError:
        return "missing"


def run_site(con, site, name):
    base = HEEDB / site
    S = {"site": site, "name": name}
    t0 = time.time()
    rd = lambda f: f"read_csv('{base / f}', all_varchar=true, header=true)"

    con.execute(f"""CREATE OR REPLACE TABLE v24 AS SELECT nk(FileName) k,
        bool_or({has('cl(codes)', C_PRE)}) pre, bool_or({has('cl(codes)', C_ALT)}) alt,
        bool_or({has('cl(codes)', C_SPR)}) spr, bool_or({has('cl(codes)', C_SINUS)}) sinus,
        bool_or({has('cl(codes)', C_AFL)}) afl, bool_or({has('cl(codes)', C_PACED)}) paced,
        bool_or({has('cl(codes)', C_WIDE)}) wide
      FROM {rd('12SL_diagnoses/diagnoses_v24.csv')} GROUP BY 1""")
    con.execute(f"""CREATE OR REPLACE TABLE acq AS SELECT nk(FileName) k,
        bool_or({has('cl(codes_software)', C_PRE)}) pre_sw, bool_or({has('cl(codes_physician)', C_PRE)}) pre_ph,
        bool_or({has('cl(codes_software)', C_ALT)} OR {has('cl(codes_physician)', C_ALT)}) alt,
        bool_or({has('cl(codes_software)', C_SPR)} OR {has('cl(codes_physician)', C_SPR)}) spr,
        bool_or({has('cl(codes_software)', C_SINUS)} OR {has('cl(codes_physician)', C_SINUS)}) sinus,
        bool_or({has('cl(codes_software)', C_AFL)} OR {has('cl(codes_physician)', C_AFL)}) afl,
        bool_or({has('cl(codes_software)', C_PACED)} OR {has('cl(codes_physician)', C_PACED)}) paced,
        bool_or({has('cl(codes_software)', C_WIDE)} OR {has('cl(codes_physician)', C_WIDE)}) wide,
        bool_or(len(cl(codes_physician)) > 0) has_ph
      FROM {rd('12SL_diagnoses/diagnoses_acquisition.csv')} GROUP BY 1""")
    # Emory metadata'da kimlik "123456789.0" biçiminde, ICD'de "123456789" → BIGINT'e normalize.
    # AgeAtAcquisition gün cinsinden (90 yaşta kesilmiş) → yıl.
    con.execute(f"""CREATE OR REPLACE TABLE meta AS SELECT pidn(BDSPPatientID) pid, nk(FileName) k,
        try_cast(ECGAcquisitionTime AS TIMESTAMP) t, try_cast(AgeAtAcquisition AS DOUBLE) / 365.25 age
      FROM {rd('metadata/metadata.csv')}""")
    con.execute(f"""CREATE OR REPLACE TABLE ecg AS
      SELECT m.pid, m.k, m.t, m.age,
        (v.k IS NOT NULL) has_v24, (a.k IS NOT NULL) has_acq, coalesce(a.has_ph, false) has_ph,
        coalesce(v.pre, false) pre_v24, coalesce(a.pre_sw, false) pre_sw, coalesce(a.pre_ph, false) pre_ph,
        (coalesce(v.pre, false) OR coalesce(a.pre_sw, false) OR coalesce(a.pre_ph, false)) pre,
        (coalesce(v.alt, false) OR coalesce(a.alt, false)) alt,
        (coalesce(v.spr, false) OR coalesce(a.spr, false)) spr,
        (coalesce(v.sinus, false) OR coalesce(a.sinus, false)) sinus,
        (coalesce(v.afl, false) OR coalesce(a.afl, false)) afl,
        (coalesce(v.paced, false) OR coalesce(a.paced, false)) paced,
        (coalesce(v.wide, false) OR coalesce(a.wide, false)) wide
      FROM meta m LEFT JOIN v24 v USING (k) LEFT JOIN acq a USING (k)""")
    con.execute(f"""ALTER TABLE ecg ADD COLUMN t0 BOOLEAN""")
    con.execute(f"""UPDATE ecg SET t0 = (NOT pre AND NOT spr AND sinus AND NOT afl AND NOT paced AND NOT wide
                                       AND coalesce(age, 0) >= {AGE_MIN})""")
    S["load_s"] = round(time.time() - t0)

    S["n_ecg"], S["n_pat"] = one(con, "SELECT count(*), count(DISTINCT pid) FROM ecg")
    S["n_rows_v24"], S["n_rows_acq"] = one(con, "SELECT (SELECT count(*) FROM v24), (SELECT count(*) FROM acq)")
    S["coverage"] = one(con, """SELECT count(*) FILTER (WHERE has_v24), count(*) FILTER (WHERE has_acq),
        count(*) FILTER (WHERE has_ph), count(*) FILTER (WHERE age IS NULL), count(*) FILTER (WHERE t IS NULL) FROM ecg""")
    S["date_range"] = [str(x) for x in one(con, "SELECT min(t), max(t) FROM ecg")]
    S["t0"] = one(con, "SELECT count(*), count(DISTINCT pid) FROM ecg WHERE t0")

    S["codes"] = [(lab,) + one(con, f"SELECT count(*), count(DISTINCT pid) FROM ecg WHERE {cond}") for lab, cond in [
        ("Preeksitasyon kodu (300/302/303/304), herhangi kaynak", "pre"),
        ("  v24 yeniden işlenmiş 12SL", "pre_v24"),
        ("  kayıt anı yazılım", "pre_sw"),
        ("  kayıt anı hekim overread", "pre_ph"),
        ("ALTWPW (303: fusion/intermittan preeksitasyon)", "alt"),
        ("Yalnızca kısa PR (102), preeksitasyon kodu yok", "spr AND NOT pre"),
    ]]

    con.execute("""CREATE OR REPLACE TABLE pat AS SELECT pid, count(*) n_ecg,
        count(*) FILTER (WHERE pre) n_pre, count(*) FILTER (WHERE pre_ph) n_pre_ph, count(*) FILTER (WHERE alt) n_alt,
        count(*) FILTER (WHERE spr AND NOT pre) n_spr_only, count(*) FILTER (WHERE t0) n_t0,
        min(t) FILTER (WHERE pre) first_pre, min(t) FILTER (WHERE t0) first_t0, min(age) FILTER (WHERE t0) age_first_t0
      FROM ecg GROUP BY pid""")

    # ICD (tanı); HEEDB'de I47.1 alt kodları yok
    con.execute(f"""CREATE OR REPLACE TABLE icd AS
      SELECT pidn(BDSPPatientID) pid, DIAGNOSIS_ICD10_CD code FROM {rd('ICD_codes/icd10_codes.csv')}
        WHERE DIAGNOSIS_ICD10_CD LIKE 'I45.6%' OR DIAGNOSIS_ICD10_CD LIKE 'I47.1%'
      UNION ALL
      SELECT pidn(BDSPPatientID) pid, DIAGNOSIS_ICD9_CD code FROM {rd('ICD_codes/icd9_codes.csv')}
        WHERE DIAGNOSIS_ICD9_CD LIKE '426.7%' OR DIAGNOSIS_ICD9_CD LIKE '427.0%'""")
    con.execute("""CREATE OR REPLACE TABLE icd_pre AS SELECT DISTINCT pid FROM icd
                   WHERE code LIKE 'I45.6%' OR code LIKE '426.7%'""")
    con.execute("""CREATE OR REPLACE TABLE icd_svt AS SELECT DISTINCT pid FROM icd
                   WHERE (code LIKE 'I47.1%' AND code NOT LIKE 'I47.11%') OR code LIKE '427.0%'""")

    con.execute("""CREATE OR REPLACE TABLE cand AS SELECT p.*, (first_t0 < first_pre) t0_before_first_pre,
        (pid IN (SELECT pid FROM icd_pre)) icd_pre FROM pat p WHERE n_pre >= 1 AND n_t0 >= 1""")
    S["pat_pre_any"] = one(con, "SELECT count(*), count(*) FILTER (WHERE n_ecg>=2) FROM pat WHERE n_pre>=1")
    S["cand"] = one(con, """SELECT count(*), sum(n_t0), sum(n_pre),
        count(*) FILTER (WHERE t0_before_first_pre), count(*) FILTER (WHERE n_t0>=2),
        count(*) FILTER (WHERE n_pre_ph>=1), count(*) FILTER (WHERE n_alt>=1), count(*) FILTER (WHERE icd_pre),
        median(age_first_t0) FROM cand""")
    S["icd_pre"] = one(con, """SELECT count(*), count(p.pid), count(*) FILTER (WHERE p.n_pre>=1),
        count(*) FILTER (WHERE coalesce(p.n_pre,0)=0 AND p.n_t0>=1)
        FROM icd_pre i LEFT JOIN pat p USING (pid)""")
    S["pool_union"] = one(con, "SELECT count(*) FROM pat WHERE n_t0>=1 AND (n_pre>=1 OR pid IN (SELECT pid FROM icd_pre))")[0]
    S["spr_only_pool"] = one(con, "SELECT count(*) FROM pat WHERE n_spr_only>0 AND n_pre=0 AND n_t0>=1")[0]
    S["k2"] = one(con, """SELECT count(*), sum(p.n_t0), count(*) FILTER (WHERE p.n_spr_only>0)
        FROM icd_svt s JOIN pat p USING (pid)
        WHERE p.n_pre=0 AND p.n_t0>=1 AND s.pid NOT IN (SELECT pid FROM icd_pre)""")
    S["svt_all"] = one(con, "SELECT count(*) FROM icd_svt")[0]

    # Örnekleme hızı: yalnızca katı adayların preeksitasyon + T0 EKG'leri
    keys = [r[0] for r in con.execute("""SELECT k FROM ecg WHERE pid IN (SELECT pid FROM cand) AND (pre OR t0)""").fetchall()]
    with ThreadPoolExecutor(32) as ex:
        fs = Counter(ex.map(read_fs, [base / "WFDB" / f"{k}.hea" for k in keys]))
    S["fs_cand"] = sorted(fs.items(), key=lambda x: -x[1])

    con.execute(f"COPY (SELECT * FROM cand) TO '{PRIVATE / f'heedb_{site}_case_candidates.parquet'}' (FORMAT parquet)")
    # Aşama 1 kontrol örneklemi için tam EKG bayrak tablosu ve ICD listeleri
    con.execute(f"COPY (SELECT * FROM ecg) TO '{PRIVATE / f'heedb_{site}_ecg_flags.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT pid FROM icd_pre) TO '{PRIVATE / f'heedb_{site}_icd_pre.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT pid FROM icd_svt) TO '{PRIVATE / f'heedb_{site}_icd_svt.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT * FROM pat WHERE coalesce(n_pre,0)=0 AND n_t0>=1 AND pid IN (SELECT pid FROM icd_pre)) TO '{PRIVATE / f'heedb_{site}_icd_only_pool.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT * FROM ecg WHERE pid IN (SELECT pid FROM cand UNION SELECT pid FROM icd_pre)) TO '{PRIVATE / f'heedb_{site}_candidate_ecgs.parquet'}' (FORMAT parquet)")
    S["site_s"] = round(time.time() - t0)
    return S


def site_md(S):
    c, icd, k2, cov = S["cand"], S["icd_pre"], S["k2"], S["coverage"]
    return "\n".join([
        f"## {S['name']} ({S['site']})\n",
        f"- EKG: {S['n_ecg']:,}; hasta: {S['n_pat']:,}; tarih aralığı (kaydırılmış): {S['date_range'][0][:10]} – {S['date_range'][1][:10]}",
        f"- 12SL v24 kodu olan EKG: {cov[0]:,}; kayıt anı kodu olan: {cov[1]:,}; hekim overread'i olan: {cov[2]:,}; yaşı eksik: {cov[3]:,}; zamanı eksik: {cov[4]:,}",
        f"- T0 adayı EKG: {S['t0'][0]:,} ({S['t0'][1]:,} hasta)\n",
        md_table(S["codes"], ["Tanım", "EKG", "Hasta"]), "",
        md_table([
            ("≥1 preeksitasyon kodlu EKG'si olan hasta", S["pat_pre_any"][0]),
            ("  bunlardan ≥2 EKG'si olan", S["pat_pre_any"][1]),
            ("**Aday (katı): ≥1 preeksitasyon EKG + ≥1 T0 adayı EKG**", c[0]),
            ("  adayların T0 adayı EKG sayısı", int(c[1] or 0)),
            ("  adayların preeksitasyon EKG sayısı", int(c[2] or 0)),
            ("  ilk preeksitasyon EKG'sinden önce T0'ı olan (prospektif varyant)", c[3]),
            ("  ≥2 T0 adayı EKG'si olan", c[4]),
            ("  hekim overread'inde preeksitasyon olan", c[5]),
            ("  ALTWPW (303) kodlu EKG'si olan", c[6]),
            ("  I45.6/426.7 kodu da olan", c[7]),
            ("I45.6/426.7 kodlu hasta", icd[0]),
            ("  bunlardan EKG'si olan", icd[1]),
            ("  bunlardan ≥1 preeksitasyon EKG'si olan", icd[2]),
            ("  **yalnızca ICD (hiç preeksitasyon kodu yok) + ≥1 T0 adayı**", icd[3]),
            ("Uzman okuma havuzu: (katı aday ∪ yalnızca-ICD)", S["pool_union"]),
            ("Yalnızca kısa PR + ≥1 T0 adayı (geniş havuz)", S["spr_only_pool"]),
            ("SVT kodlu hasta (I47.1 / 427.0)", S["svt_all"]),
            ("K2 adayı: SVT kodu, preeksitasyon kodu ve I45.6/426.7 yok, ≥1 T0 adayı", k2[0]),
            ("  K2 adaylarının T0 adayı EKG sayısı", int(k2[1] or 0)),
            ("  kısa PR kodlu EKG'si olan", k2[2]),
        ], ["Grup", "Hasta"]),
        f"\nKatı adayların ilk T0 EKG'sindeki medyan yaş: {round(c[8], 1)}.\n",
        "Katı adayların preeksitasyon + T0 EKG'lerinde örnekleme hızı (WFDB başlığı):\n",
        md_table([(str(k), v) for k, v in S["fs_cand"]], ["fs (Hz)", "EKG"]),
        f"\nSüre: {S['site_s']} s.\n"])


def main():
    PRIVATE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    con = connect(threads=24, memory="40GB")
    con.execute(r"""CREATE MACRO nk(f) AS regexp_replace(regexp_replace(regexp_replace(regexp_replace(
        f, '\s', '', 'g'), '^\./|^/', ''), '^WFDB/', ''), '\.hea$', '')""")
    con.execute("CREATE MACRO pidn(x) AS try_cast(try_cast(x AS DOUBLE) AS BIGINT)")
    con.execute("""CREATE MACRO cl(c) AS list_filter(list_transform(string_split(coalesce(c, ''), ','),
        x -> try_cast(trim(x) AS INTEGER)), x -> x IS NOT NULL)""")
    res = [run_site(con, s, n) for s, n in HEEDB_SITES.items()]
    md = ["# Aşama 0 — HEEDB sayımları\n",
          "HEEDB 5.0. Yalnızca agregat. Preeksitasyon kodları 12SL 300/302/303/304; kaynaklar v24 yeniden işlenmiş "
          "12SL, kayıt anı yazılım ve hekim overread. T0 adayı: hiçbir kaynakta preeksitasyon/kısa PR kodu yok, "
          f"sinüs kodu var, AF/flutter/pace/dal bloğu/IVCD kodu yok, yaş ≥{AGE_MIN}. HEEDB'de prosedür kodu ve "
          "WFDB başlığında filtre bilgisi yok.\n"] + [site_md(S) for S in res]
    (OUT / "heedb_stage0.md").write_text("\n".join(md))
    (OUT / "heedb_stage0.json").write_text(json.dumps(res, indent=1, default=str))
    print("\n".join(md))


if __name__ == "__main__":
    main()
