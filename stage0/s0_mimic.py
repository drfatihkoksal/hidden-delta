"""HIDDEN-DELTA Aşama 0 — MIMIC sayımları.

Kaynaklar: MIMIC-IV-ECG 1.0 (machine_measurements), MIMIC-IV 3.1 hosp (patients, diagnoses_icd,
procedures_icd), MIMIC-IV-Note (discharge; yalnızca K2 alt tiplemesi için, metin hiçbir çıktıya girmez).

Çıktı: stage0/out/mimic_stage0.{md,json} (agregat) ve PRIVATE/mimic_*.parquet (hasta düzeyi).
Çalıştırma: python3 stage0/s0_mimic.py
"""
import json
import time

from common import (ABL_ICD9, ABL_ICD10, AGE_MIN, MIMIC_ECG, MIMIC_HOSP, MIMIC_NOTE, OUT,
                    POST_ABL_EARLY_DAYS, PRIVATE, RX_AFL, RX_INCOMPLETE_BBB, RX_PACED, RX_PRE,
                    RX_SINUS, RX_SPR, RX_WIDE, connect, md_table)

# Taburcu özeti kalıpları (RE2, büyük/küçük harf duyarsız). Olumsuzlama ayıklanmaz: kaba sayım.
NX_AVNRT = r"(?i)\bAVNRT\b|AV[ -]?NODAL RE-?ENTR|ATRIO-?VENTRICULAR NODAL RE-?ENTR|SLOW[ -]PATHWAY"
NX_AP = r"(?i)CONCEALED (ACCESSORY )?(PATHWAY|BYPASS)|\bAVRT\b|ORTHODROMIC|ACCESSORY PATHWAY|BYPASS TRACT"
NX_WPW = r"(?i)WOLFF|\bWPW\b|PRE-?EXCITATION"
NX_EPS = r"(?i)\bEP STUDY\b|ELECTROPHYSIOLOG(Y|IC) STUDY|\bEPS\b|ABLATION"


def one(con, sql, params=None):
    return con.execute(sql, params or []).fetchone()


def main():
    t_start = time.time()
    PRIVATE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    con = connect()
    R = {}

    # ---------------- EKG tablosu ----------------
    cols = ", ".join(f"report_{i}" for i in range(18))
    con.execute(f"""
    CREATE TABLE ecg0 AS
    WITH raw AS (
      SELECT subject_id::BIGINT sid, study_id::BIGINT study, cart_id,
             try_cast(ecg_time AS TIMESTAMP) t,
             upper(concat_ws(' | ', {cols})) rep, bandwidth, filtering,
             CASE WHEN try_cast(qrs_end AS INT) < 29999 AND try_cast(qrs_onset AS INT) < 29999
                  THEN try_cast(qrs_end AS INT) - try_cast(qrs_onset AS INT) END qrsd
      FROM read_csv('{MIMIC_ECG / "machine_measurements.csv"}', all_varchar=true, header=true, sample_size=-1)
    )
    SELECT r.*, p.anchor_age::INT + year(r.t) - p.anchor_year::INT AS age,
      regexp_matches(rep, '{RX_PRE}') pre,
      regexp_matches(rep, '{RX_SPR}') spr,
      regexp_matches(rep, '{RX_SINUS}') sinus,
      regexp_matches(rep, '{RX_AFL}') afl,
      regexp_matches(rep, '{RX_PACED}') paced,
      (coalesce(qrsd, 0) >= 120
       OR regexp_matches(regexp_replace(rep, '{RX_INCOMPLETE_BBB}', '', 'g'), '{RX_WIDE}')) wide
    FROM raw r
    LEFT JOIN read_csv('{MIMIC_HOSP / "patients.csv.gz"}', all_varchar=true) p ON r.sid = p.subject_id::BIGINT
    """)
    con.execute(f"""
    CREATE TABLE ecg AS SELECT *,
      (NOT pre AND NOT spr AND sinus AND NOT afl AND NOT paced AND NOT wide AND coalesce(age, 0) >= {AGE_MIN}) t0
    FROM ecg0""")
    con.execute("DROP TABLE ecg0")

    con.execute("""
    CREATE TABLE pat AS SELECT sid, count(*) n_ecg,
      count(*) FILTER (WHERE pre) n_pre,
      count(*) FILTER (WHERE spr AND NOT pre) n_spr_only,
      count(*) FILTER (WHERE t0) n_t0,
      min(t) FILTER (WHERE pre) first_pre, max(t) FILTER (WHERE pre) last_pre,
      min(t) FILTER (WHERE t0) first_t0, max(t) FILTER (WHERE t0) last_t0,
      min(t) first_ecg, max(t) last_ecg
    FROM ecg GROUP BY sid""")

    # ---------------- ICD ve prosedürler ----------------
    con.execute(f"""CREATE TABLE dx AS SELECT subject_id::BIGINT sid, icd_code, icd_version::INT v
                    FROM read_csv('{MIMIC_HOSP / "diagnoses_icd.csv.gz"}', all_varchar=true)""")
    con.execute("""CREATE TABLE icd_pre AS SELECT DISTINCT sid FROM dx
                   WHERE (v=10 AND icd_code LIKE 'I456%') OR (v=9 AND icd_code LIKE '4267%')""")
    con.execute("""CREATE TABLE icd_svt AS SELECT DISTINCT sid FROM dx
                   WHERE (v=10 AND icd_code LIKE 'I471%' AND icd_code <> 'I4711') OR (v=9 AND icd_code LIKE '4270%')""")
    con.execute(f"""CREATE TABLE abl AS SELECT subject_id::BIGINT sid, min(try_cast(chartdate AS DATE)) d_abl,
                      count(DISTINCT chartdate) n_abl_days
                    FROM read_csv('{MIMIC_HOSP / "procedures_icd.csv.gz"}', all_varchar=true)
                    WHERE (icd_version='10' AND icd_code='{ABL_ICD10}') OR (icd_version='9' AND icd_code='{ABL_ICD9}')
                    GROUP BY 1""")

    # ---------------- 0. Genel ----------------
    R["n_ecg"], R["n_pat"] = one(con, "SELECT count(*), count(DISTINCT sid) FROM ecg")
    R["n_ecg_age_missing"] = one(con, "SELECT count(*) FROM ecg WHERE age IS NULL")[0]
    R["ecg_date_range"] = [str(x) for x in one(con, "SELECT min(t), max(t) FROM ecg")]
    R["n_t0_ecg"], R["n_t0_pat"] = one(con, "SELECT count(*), count(DISTINCT sid) FROM ecg WHERE t0")

    # ---------------- 1. Regex yakalama ----------------
    r1 = []
    for lab, cond in [("Concept regex (preeksitasyon ifadesi VEYA kısa PR)", "pre OR spr"),
                      ("Preeksitasyon ifadesi (katı: WPW/Wolff/pre-excit/delta)", "pre"),
                      ("Yalnızca kısa PR (preeksitasyon ifadesi yok)", "spr AND NOT pre")]:
        n_e, n_p = one(con, f"SELECT count(*), count(DISTINCT sid) FROM ecg WHERE {cond}")
        r1.append((lab, n_e, n_p))
    R["regex"] = r1
    R["n_pat_spr_never_pre"] = one(con, "SELECT count(*) FROM pat WHERE n_spr_only>0 AND n_pre=0")[0]

    # ---------------- 2. İntermittan aday havuzları ----------------
    con.execute("""CREATE TABLE cand AS SELECT p.*,
        (first_t0 < first_pre) t0_before_first_pre,
        (p.sid IN (SELECT sid FROM icd_pre)) icd_pre
      FROM pat p WHERE n_pre >= 1 AND n_t0 >= 1""")
    R["cand_strict"] = one(con, """SELECT count(*), sum(n_t0), sum(n_pre),
        count(*) FILTER (WHERE t0_before_first_pre), count(*) FILTER (WHERE n_t0>=2), count(*) FILTER (WHERE icd_pre)
        FROM cand""")
    R["pat_pre_any"] = one(con, "SELECT count(*), count(*) FILTER (WHERE n_ecg>=2) FROM pat WHERE n_pre>=1")
    R["icd_pre"] = one(con, """SELECT count(*), count(p.sid), count(*) FILTER (WHERE p.n_pre>=1),
        count(*) FILTER (WHERE coalesce(p.n_pre,0)=0 AND p.n_t0>=1)
        FROM icd_pre i LEFT JOIN pat p USING (sid)""")
    R["spr_only_pool"] = one(con, "SELECT count(*), sum(n_t0) FROM pat WHERE n_spr_only>0 AND n_pre=0 AND n_t0>=1")
    R["pool_union"] = one(con, """SELECT count(*) FROM pat p WHERE n_t0>=1 AND
        (n_pre>=1 OR sid IN (SELECT sid FROM icd_pre))""")[0]
    # Yaşa göre: aday hastalar ilk T0 EKG'sinde
    R["cand_age"] = one(con, """SELECT median(age), min(age), max(age) FROM
        (SELECT e.sid, min(e.age) age FROM ecg e JOIN cand c USING (sid) WHERE e.t0 GROUP BY 1)""")

    # ---------------- 3. Ablasyon ----------------
    con.execute(f"""
    CREATE TABLE pabl AS
    SELECT a.sid, a.d_abl, a.n_abl_days,
      (a.sid IN (SELECT sid FROM icd_pre)) icdpre,
      count(e.t) n_ecg_all,
      count(e.t) FILTER (WHERE e.t::DATE < a.d_abl) n_ecg_before,
      count(e.t) FILTER (WHERE e.pre AND e.t::DATE < a.d_abl) n_pre_before,
      count(e.t) FILTER (WHERE e.t::DATE > a.d_abl) n_ecg_after,
      count(e.t) FILTER (WHERE e.t0 AND e.t::DATE > a.d_abl AND e.t::DATE <= a.d_abl + INTERVAL {POST_ABL_EARLY_DAYS} DAY) n_t0_early,
      count(e.t) FILTER (WHERE e.t0 AND e.t::DATE > a.d_abl + INTERVAL {POST_ABL_EARLY_DAYS} DAY) n_t0_late,
      count(e.t) FILTER (WHERE e.pre AND e.t::DATE > a.d_abl) n_pre_after
    FROM abl a LEFT JOIN ecg e ON e.sid = a.sid
    GROUP BY a.sid, a.d_abl, a.n_abl_days""")
    con.execute("CREATE TABLE pabl_ap AS SELECT * FROM pabl WHERE icdpre OR n_pre_before > 0")
    R["abl_any"] = one(con, "SELECT count(*), count(*) FILTER (WHERE n_ecg_all>0) FROM pabl")
    R["abl_ap"] = one(con, """SELECT count(*),
        count(*) FILTER (WHERE n_pre_before>0),
        count(*) FILTER (WHERE icdpre AND n_pre_before=0),
        count(*) FILTER (WHERE n_ecg_before>0 AND n_ecg_after>0),
        count(*) FILTER (WHERE n_pre_before>0 AND (n_t0_early+n_t0_late)>0),
        count(*) FILTER (WHERE n_pre_before>0 AND n_t0_early>0),
        count(*) FILTER (WHERE n_pre_before>0 AND n_t0_late>0),
        count(*) FILTER (WHERE n_pre_after>0)
        FROM pabl_ap""")
    R["abl_ap_in_cand"] = one(con, "SELECT count(*) FROM pabl_ap WHERE sid IN (SELECT sid FROM cand)")[0]

    # ---------------- 4. Filtre / cart envanteri ----------------
    R["filter_all"] = con.execute("""SELECT coalesce(filtering,'NA'), coalesce(bandwidth,'NA'), count(*) n
        FROM ecg GROUP BY 1,2 ORDER BY n DESC""").fetchall()
    R["filter_cand"] = con.execute("""SELECT coalesce(filtering,'NA'), coalesce(bandwidth,'NA'), count(*) n
        FROM ecg WHERE sid IN (SELECT sid FROM cand) AND (pre OR t0) GROUP BY 1,2 ORDER BY n DESC""").fetchall()
    R["carts"] = one(con, """SELECT count(DISTINCT cart_id),
        (SELECT count(DISTINCT cart_id) FROM ecg WHERE sid IN (SELECT sid FROM cand) AND (pre OR t0)) FROM ecg""")

    # ---------------- 5. K2 (SVT kontrolü) ----------------
    con.execute("""CREATE TABLE k2 AS SELECT p.sid, p.n_t0, (p.n_spr_only>0) any_spr,
        (p.sid IN (SELECT sid FROM abl)) any_abl
      FROM icd_svt s JOIN pat p USING (sid)
      WHERE p.n_pre=0 AND p.n_t0>=1 AND s.sid NOT IN (SELECT sid FROM icd_pre)""")
    R["k2"] = one(con, """SELECT count(*), sum(n_t0), count(*) FILTER (WHERE any_spr), count(*) FILTER (WHERE any_abl) FROM k2""")
    R["svt_all"] = one(con, "SELECT count(*) FROM icd_svt")[0]

    t_note = time.time()
    con.execute(f"""
    CREATE TABLE k2note AS
    SELECT subject_id::BIGINT sid, count(*) n_notes,
      bool_or(regexp_matches(text, $1)) avnrt, bool_or(regexp_matches(text, $2)) ap,
      bool_or(regexp_matches(text, $3)) wpw, bool_or(regexp_matches(text, $4)) eps
    FROM read_csv('{MIMIC_NOTE}', header=true, all_varchar=true, quote='"', escape='"',
                  columns={{'note_id':'VARCHAR','subject_id':'VARCHAR','hadm_id':'VARCHAR','note_type':'VARCHAR',
                            'note_seq':'VARCHAR','charttime':'VARCHAR','storetime':'VARCHAR','text':'VARCHAR'}})
    WHERE subject_id::BIGINT IN (SELECT sid FROM k2)
    GROUP BY 1""", [NX_AVNRT, NX_AP, NX_WPW, NX_EPS])
    R["note_seconds"] = round(time.time() - t_note)
    R["k2_notes"] = one(con, """SELECT count(*),
        count(*) FILTER (WHERE avnrt AND NOT ap), count(*) FILTER (WHERE ap AND NOT avnrt),
        count(*) FILTER (WHERE avnrt AND ap), count(*) FILTER (WHERE NOT avnrt AND NOT ap),
        count(*) FILTER (WHERE wpw), count(*) FILTER (WHERE eps),
        count(*) FILTER (WHERE avnrt AND NOT ap AND NOT wpw), count(*) FILTER (WHERE ap AND NOT avnrt AND NOT wpw)
        FROM k2 LEFT JOIN k2note USING (sid) WHERE k2note.sid IS NOT NULL""")

    # ---------------- Hasta düzeyi listeler (PRIVATE) ----------------
    con.execute(f"COPY (SELECT c.*, a.d_abl, a.n_pre_before, a.n_t0_early, a.n_t0_late, a.n_pre_after FROM cand c LEFT JOIN pabl_ap a USING (sid)) TO '{PRIVATE / 'mimic_case_candidates.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT p.* FROM pat p WHERE coalesce(n_pre,0)=0 AND n_t0>=1 AND sid IN (SELECT sid FROM icd_pre)) TO '{PRIVATE / 'mimic_icd_only_pool.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT * FROM pabl_ap) TO '{PRIVATE / 'mimic_pathway_ablation.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT k2.*, n.n_notes, n.avnrt, n.ap, n.wpw, n.eps FROM k2 LEFT JOIN k2note n USING (sid)) TO '{PRIVATE / 'mimic_k2.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT sid, study, t, pre, spr, t0, sinus, afl, paced, wide, age, filtering, bandwidth, cart_id FROM ecg WHERE sid IN (SELECT sid FROM cand UNION SELECT sid FROM icd_pre UNION SELECT sid FROM pabl_ap)) TO '{PRIVATE / 'mimic_candidate_ecgs.parquet'}' (FORMAT parquet)")

    # Aşama 1 kontrol örneklemi için tam EKG bayrak tablosu (rapor metni hariç) ve ICD/ablasyon listeleri
    con.execute(f"COPY (SELECT sid, study, t, age, pre, spr, t0, sinus, afl, paced, wide, qrsd, filtering, bandwidth, cart_id FROM ecg) TO '{PRIVATE / 'mimic_ecg_flags.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT sid FROM icd_pre) TO '{PRIVATE / 'mimic_icd_pre.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT sid FROM icd_svt) TO '{PRIVATE / 'mimic_icd_svt.parquet'}' (FORMAT parquet)")
    con.execute(f"COPY (SELECT * FROM abl) TO '{PRIVATE / 'mimic_abl.parquet'}' (FORMAT parquet)")

    R["runtime_s"] = round(time.time() - t_start)
    write_report(R)


def write_report(R):
    cs, icd, ab, k2, kn = R["cand_strict"], R["icd_pre"], R["abl_ap"], R["k2"], R["k2_notes"]
    md = [f"# Aşama 0 — MIMIC sayımları\n",
          f"MIMIC-IV-ECG 1.0 + MIMIC-IV 3.1 hosp + MIMIC-IV-Note. Yalnızca agregat. Çalışma süresi {R['runtime_s']} s.\n",
          f"- EKG: {R['n_ecg']:,}; hasta: {R['n_pat']:,}; EKG tarih aralığı (kaydırılmış): {R['ecg_date_range'][0][:10]} – {R['ecg_date_range'][1][:10]}",
          f"- Yaşı hesaplanamayan EKG: {R['n_ecg_age_missing']:,}",
          f"- T0 adayı EKG (preeksitasyon/kısa PR ifadesi yok, sinüs, AF/flutter/pace yok, QRS <120 ms ve dal bloğu/IVCD yok, yaş ≥{AGE_MIN}): {R['n_t0_ecg']:,} EKG, {R['n_t0_pat']:,} hasta\n",
          "## 1. Rapor regex'i\n",
          md_table(R["regex"], ["Tanım", "EKG", "Hasta"]),
          f"\nYalnızca kısa PR ifadesi olup hiç preeksitasyon ifadesi olmayan hasta: {R['n_pat_spr_never_pre']:,}. "
          "MIMIC makine raporlarında preeksitasyonla ilgili seri karşılaştırma ifadesi (\"no longer evident\" vb.) yok.\n",
          "## 2. İntermittan aday havuzları (makine raporu düzeyinde; uzman doğrulaması öncesi)\n",
          md_table([
              ("≥1 preeksitasyon ifadeli EKG olan hasta", R["pat_pre_any"][0]),
              ("  bunlardan ≥2 EKG'si olan", R["pat_pre_any"][1]),
              ("**Aday (katı): ≥1 preeksitasyon EKG + ≥1 T0 adayı EKG**", cs[0]),
              ("  adayların T0 adayı EKG sayısı", int(cs[1] or 0)),
              ("  adayların preeksitasyon EKG sayısı", int(cs[2] or 0)),
              ("  ilk preeksitasyon EKG'sinden önce T0'ı olan (prospektif varyant)", cs[3]),
              ("  ≥2 T0 adayı EKG'si olan", cs[4]),
              ("  I45.6/426.7 kodu da olan", cs[5]),
              ("I45.6/426.7 kodlu hasta", icd[0]),
              ("  bunlardan EKG'si olan", icd[1]),
              ("  bunlardan ≥1 preeksitasyon EKG'si olan", icd[2]),
              ("  **yalnızca ICD (hiç preeksitasyon ifadesi yok) + ≥1 T0 adayı**", icd[3]),
              ("Uzman okuma havuzu: (katı aday ∪ yalnızca-ICD)", R["pool_union"]),
              ("Yalnızca kısa PR + ≥1 T0 adayı (geniş havuz, çoğu WPW değil)", R["spr_only_pool"][0]),
          ], ["Grup", "Hasta"]),
          f"\nKatı adayların ilk T0 EKG'sindeki yaş: medyan {R['cand_age'][0]}, aralık {R['cand_age'][1]}–{R['cand_age'][2]}.\n",
          "## 3. Ablasyon (H3)\n",
          f"Ablasyon kodu (02583ZZ / 37.34; yolağa özgü değil) olan hasta: {R['abl_any'][0]:,}; EKG'si olan: {R['abl_any'][1]:,}.\n",
          md_table([
              ("**Yolak ablasyonu vekili** (ablasyon + [I45.6/426.7 ya da ablasyon öncesi preeksitasyon EKG'si])", ab[0]),
              ("  ablasyon öncesi preeksitasyon EKG'si olan", ab[1]),
              ("  yalnızca ICD ile tanımlanan", ab[2]),
              ("  ablasyon öncesi ve sonrası EKG'si olan", ab[3]),
              ("  öncesi preeksitasyon EKG + sonrası T0 adayı", ab[4]),
              (f"    sonrası T0 0–{POST_ABL_EARLY_DAYS} gün (bellek olası)", ab[5]),
              (f"    sonrası T0 >{POST_ABL_EARLY_DAYS} gün", ab[6]),
              ("  ablasyon sonrası preeksitasyon EKG'si olan (başarısız/nüks?)", ab[7]),
              ("  aynı zamanda katı intermittan aday olan", R["abl_ap_in_cand"]),
          ], ["Grup", "Hasta"]),
          "\n## 4. Filtre / cart envanteri\n",
          f"Farklı cart_id: tümünde {R['carts'][0]}, katı adayların preeksitasyon + T0 EKG'lerinde {R['carts'][1]}.\n",
          "Tüm EKG'ler:\n", md_table(R["filter_all"], ["filtering", "bandwidth", "EKG"]),
          "\nKatı adayların preeksitasyon + T0 EKG'leri:\n", md_table(R["filter_cand"], ["filtering", "bandwidth", "EKG"]),
          "\n## 5. K2 — SVT kontrolü\n",
          f"SVT kodlu hasta (I47.1x, I47.11 hariç / 427.0): {R['svt_all']:,}.\n",
          md_table([
              ("K2 adayı: SVT kodu, preeksitasyon ifadesi ve I45.6/426.7 yok, ≥1 T0 adayı", k2[0]),
              ("  T0 adayı EKG sayısı", int(k2[1] or 0)),
              ("  kısa PR ifadeli EKG'si olan", k2[2]),
              ("  ablasyon kodu olan", k2[3]),
              ("Taburcu özeti olan K2 hastası", kn[0]),
              ("  yalnızca AVNRT kalıbı", kn[1]),
              ("  yalnızca aksesuar yol / AVRT kalıbı", kn[2]),
              ("  ikisi birden", kn[3]),
              ("  hiçbiri", kn[4]),
              ("  notta WPW/preeksitasyon sözü geçen (kontaminasyon adayı)", kn[5]),
              ("  notta EPS/ablasyon sözü geçen", kn[6]),
              ("  **K2a adayı: yalnızca AVNRT, WPW sözü yok**", kn[7]),
              ("  **K2b adayı: yalnızca aksesuar yol/AVRT, WPW sözü yok**", kn[8]),
          ], ["Grup", "Hasta"]),
          "\nNot kalıpları olumsuzlamayı ayıklamaz (\"no evidence of AVNRT\" da sayılır); alt tip ancak okuma ile kesinleşir.\n"]
    (OUT / "mimic_stage0.md").write_text("\n".join(md))
    (OUT / "mimic_stage0.json").write_text(json.dumps(R, indent=1, default=str))
    print("\n".join(md))


if __name__ == "__main__":
    main()
