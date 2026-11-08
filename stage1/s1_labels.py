"""HIDDEN-DELTA Aşama 1 — etiket tanımının uygulanması (concept.md §5.2), kesin G1 sayımı ve
kaynaklar arası uyum. Manuel okuma yok.

HEEDB  T2 (birincil) = üç kaynaktan ≥2'sinde preeksitasyon kodu (v24 12SL, kayıt anı yazılım, hekim
       overread). T2 (duyarlılık) = yalnızca overread; yalnızca v24 12SL. T0 = overread var + hiçbir kaynakta preeksitasyon/kısa PR
       kodu yok + sinüs, AF/flutter/pace/dal bloğu yok, yaş ≥18. İndeks (T0) EKG'ler yalnızca 500 Hz.
MIMIC  T2 (birincil) = kesin preeksitasyon ifadesi (Possible/Consider/Borderline içermeyen);
       T2 (duyarlılık) = herhangi preeksitasyon ifadesi. T0 = Aşama 0 tanımı.
       Sıkı T0 (§8 etiket gürültüsü denetimi) = T0 + makine PR ≥140 ms + QRS <110 ms.

Girdi: Aşama 0 PRIVATE çıktıları (stage0/s0_*.py önce çalışmış olmalı).
Çıktı: stage1/out/labels_summary.{md,json} (agregat); ~/data/hidden-delta/stage1/*.parquet (hasta düzeyi).
Çalıştırma: python3 stage1/s1_labels.py
"""
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage0"))
from common import (HEEDB, HEEDB_SITES, MIMIC_ECG, POST_ABL_EARLY_DAYS, PRIVATE, RX_PRE,  # noqa: E402
                    connect, md_table)
from s0_heedb import read_fs  # noqa: E402

P1 = Path("/home/fatih/data/hidden-delta/stage1")
OUT = Path(__file__).resolve().parent / "out"
RX_SOFT = r"POSSIBLE|CONSIDER|BORDERLINE"
PR_MIN_STRICT, QRS_MAX_STRICT = 140, 110


def q(name):
    return f"'{PRIVATE / name}'"


def one(con, sql):
    return con.execute(sql).fetchone()


def kappa2(n11, n10, n01, n00):
    n = n11 + n10 + n01 + n00
    po = (n11 + n00) / n
    pe = ((n11 + n10) * (n11 + n01) + (n01 + n00) * (n10 + n00)) / n**2
    return (po - pe) / (1 - pe), 2 * n11 / max(1, 2 * n11 + n10 + n01)


def fleiss3(patterns):
    """patterns: {(a,b,c): count} ikili üç değerlendirici."""
    n = sum(patterns.values())
    pbar = sum(c * (k * (k - 1) + (3 - k) * (2 - k)) / 6 for (a, b, cc), c in patterns.items() for k in [a + b + cc]) / n
    p1 = sum(c * (a + b + cc) for (a, b, cc), c in patterns.items()) / (3 * n)
    pe = p1**2 + (1 - p1) ** 2
    return (pbar - pe) / (1 - pe)


# ------------------------------------------------------------------ MIMIC
def run_mimic(con):
    R = {}
    cols = ", ".join(f"report_{i}" for i in range(18))
    con.execute(f"""CREATE OR REPLACE TABLE mm AS SELECT subject_id::BIGINT sid, study_id::BIGINT study,
        list_bool_or(list_transform(string_split(upper(concat_ws(' | ', {cols})), ' | '),
             s -> regexp_matches(s, '{RX_PRE}') AND NOT regexp_matches(s, '{RX_SOFT}'))) pre_def,
        CASE WHEN try_cast(p_onset AS INT) < 29999 AND try_cast(qrs_onset AS INT) < 29999
             THEN try_cast(qrs_onset AS INT) - try_cast(p_onset AS INT) END pr
      FROM read_csv('{MIMIC_ECG / 'machine_measurements.csv'}', all_varchar=true, sample_size=-1)""")
    con.execute(f"""CREATE OR REPLACE TABLE me AS SELECT f.*, m.pr,
        (f.pre AND m.pre_def) t2_def, f.pre t2_any,
        (f.t0 AND m.pr >= {PR_MIN_STRICT} AND f.qrsd < {QRS_MAX_STRICT}) t0_strict
      FROM {q('mimic_ecg_flags.parquet')} f JOIN mm m USING (sid, study)""")
    con.execute("""CREATE OR REPLACE TABLE mp AS SELECT sid,
        count(*) FILTER (WHERE t2_def) n_t2_def, count(*) FILTER (WHERE t2_any) n_t2_any,
        count(*) FILTER (WHERE t0) n_t0, count(*) FILTER (WHERE t0_strict) n_t0s,
        min(t) FILTER (WHERE t2_def) first_t2_def, min(t) FILTER (WHERE t2_any) first_t2_any,
        min(t) FILTER (WHERE t0) first_t0, min(t) FILTER (WHERE t0_strict) first_t0s
      FROM me GROUP BY sid""")
    R["t0_pr_available"] = one(con, "SELECT count(*) FILTER (WHERE pr IS NOT NULL), count(*) FROM me WHERE t0")
    for lab, t2 in [("def", "def"), ("any", "any")]:
        R[f"case_{lab}"] = one(con, f"""SELECT count(*), sum(n_t0),
            count(*) FILTER (WHERE first_t0 < first_t2_{t2}), count(*) FILTER (WHERE n_t0s >= 1),
            count(*) FILTER (WHERE first_t0s < first_t2_{t2})
            FROM mp WHERE n_t2_{t2} >= 1 AND n_t0 >= 1""")
    # H3: yolak ablasyonu vekili, kesin T2 ile
    con.execute(f"""CREATE OR REPLACE TABLE mabl AS
        SELECT a.sid, a.d_abl, (a.sid IN (SELECT sid FROM {q('mimic_icd_pre.parquet')})) icdpre,
          count(e.t) FILTER (WHERE e.t2_def AND e.t::DATE < a.d_abl) n_t2def_before,
          count(e.t) FILTER (WHERE e.t2_any AND e.t::DATE < a.d_abl) n_t2any_before,
          count(e.t) FILTER (WHERE e.t0 AND e.t::DATE > a.d_abl AND e.t::DATE <= a.d_abl + INTERVAL {POST_ABL_EARLY_DAYS} DAY) n_t0_early,
          count(e.t) FILTER (WHERE e.t0 AND e.t::DATE > a.d_abl + INTERVAL {POST_ABL_EARLY_DAYS} DAY) n_t0_late
        FROM {q('mimic_abl.parquet')} a LEFT JOIN me e ON e.sid = a.sid GROUP BY a.sid, a.d_abl""")
    R["h3"] = one(con, """SELECT
        count(*) FILTER (WHERE n_t2def_before > 0 AND n_t0_early + n_t0_late > 0),
        count(*) FILTER (WHERE n_t2def_before > 0 AND n_t0_early > 0),
        count(*) FILTER (WHERE n_t2def_before > 0 AND n_t0_late > 0),
        count(*) FILTER (WHERE n_t2any_before > 0 AND n_t0_early + n_t0_late > 0)
        FROM mabl""")
    k2 = f"{q('mimic_k2.parquet')}"
    R["k2"] = one(con, f"""SELECT count(*), count(*) FILTER (WHERE p.n_t0s >= 1),
        count(*) FILTER (WHERE k.avnrt AND NOT k.ap AND NOT k.wpw), count(*) FILTER (WHERE k.avnrt AND NOT k.ap AND NOT k.wpw AND p.n_t0s >= 1),
        count(*) FILTER (WHERE k.ap AND NOT k.avnrt AND NOT k.wpw)
        FROM {k2} k JOIN mp p USING (sid)""")
    R["filter_case_def"] = con.execute("""SELECT coalesce(filtering,'NA'), coalesce(bandwidth,'NA'), count(*) n FROM me
        WHERE (t0 OR t2_def) AND sid IN (SELECT sid FROM mp WHERE n_t2_def >= 1 AND n_t0 >= 1)
        GROUP BY 1, 2 ORDER BY n DESC""").fetchall()
    # hasta düzeyi çıktılar
    con.execute(f"""COPY (SELECT p.*, (n_t2_def >= 1 AND n_t0 >= 1) case_def, (n_t2_any >= 1 AND n_t0 >= 1) case_any,
          (sid IN (SELECT sid FROM {k2})) k2,
          a.d_abl, a.n_t2def_before, a.n_t0_early, a.n_t0_late
        FROM mp p LEFT JOIN mabl a USING (sid)
        WHERE (n_t2_any >= 1 AND n_t0 >= 1) OR sid IN (SELECT sid FROM {k2}) OR a.sid IS NOT NULL)
        TO '{P1 / 'mimic_cohort.parquet'}' (FORMAT parquet)""")
    con.execute(f"""COPY (SELECT sid, study, t, age, t2_def, t2_any, t0, t0_strict, pr, qrsd, filtering, bandwidth, cart_id
        FROM me WHERE sid IN (SELECT sid FROM mp WHERE n_t2_any >= 1 AND n_t0 >= 1)
           OR sid IN (SELECT sid FROM {k2}) OR sid IN (SELECT sid FROM mabl))
        TO '{P1 / 'mimic_ecg_labels.parquet'}' (FORMAT parquet)""")
    return R


# ------------------------------------------------------------------ HEEDB
# T2 varyantları: birincil = üç kaynaktan ≥2'si (uzlaşı); duyarlılık = yalnızca overread, yalnızca v24 12SL
T2_VARIANTS = {"cons": "(pre_v24::INT + pre_sw::INT + pre_ph::INT) >= 2", "ph": "pre_ph", "v24": "pre_v24"}
T2_LABEL = {"cons": "≥2 kaynak uzlaşısı (birincil)", "ph": "hekim overread (duyarlılık)", "v24": "v24 12SL (duyarlılık)"}


def run_heedb(con, site):
    S = {"site": site, "name": HEEDB_SITES[site]}
    t2cols = ", ".join(f"{expr} AS t2_{v}" for v, expr in T2_VARIANTS.items())
    con.execute(f"""CREATE OR REPLACE TABLE he AS SELECT pid, k, t, age, has_v24, has_acq, has_ph,
        pre_v24, pre_sw, pre_ph, pre, alt, spr, {t2cols}, (t0 AND has_ph) AS t0_lab
      FROM {q(f'heedb_{site}_ecg_flags.parquet')}""")
    # kaynaklar arası uyum: üç kaynağın da olduğu EKG'ler
    pat = {tuple(int(x) for x in r[:3]): r[3] for r in con.execute("""SELECT pre_v24, pre_sw, pre_ph, count(*)
        FROM he WHERE has_v24 AND has_acq AND has_ph GROUP BY 1, 2, 3""").fetchall()}
    S["n_kappa"] = sum(pat.values())
    S["patterns"] = {"".join(map(str, k)): v for k, v in sorted(pat.items())}

    def pair(i, j):
        c = lambda a, b: sum(v for k, v in pat.items() if k[i] == a and k[j] == b)
        return kappa2(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
    S["kappa"] = {"v24–yazılım": pair(0, 1), "v24–overread": pair(0, 2), "yazılım–overread": pair(1, 2)}
    S["fleiss"] = fleiss3(pat)

    any_t2 = " OR ".join(f"t2_{v}" for v in T2_VARIANTS)
    keys = [r[0] for r in con.execute(f"""SELECT k FROM he WHERE (t0_lab OR {any_t2}) AND pid IN
        (SELECT pid FROM he GROUP BY pid HAVING bool_or({any_t2}) AND bool_or(t0_lab))""").fetchall()]
    t_fs = time.time()
    with ThreadPoolExecutor(32) as ex:
        raw = list(ex.map(read_fs, [HEEDB / site / "WFDB" / f"{k}.hea" for k in keys]))
    fs = [int(x) if str(x).isdigit() else None for x in raw]
    S["fs_s"] = round(time.time() - t_fs)
    S["fs_missing"] = sum(1 for x in fs if x not in (250, 500))
    con.register("fsmap", __import__("pandas").DataFrame({"k": keys, "fs": fs}))
    con.execute("CREATE OR REPLACE TABLE hl AS SELECT he.*, fsmap.fs FROM he JOIN fsmap USING (k)")
    per_v = ",\n        ".join(
        f"count(*) FILTER (WHERE t2_{v}) n_t2_{v}, count(*) FILTER (WHERE t2_{v} AND fs = 500) n_t2_{v}_500, "
        f"min(t) FILTER (WHERE t2_{v}) first_t2_{v}" for v in T2_VARIANTS)
    con.execute(f"""CREATE OR REPLACE TABLE hp AS SELECT pid, {per_v},
        count(*) FILTER (WHERE t0_lab) n_t0, count(*) FILTER (WHERE t0_lab AND fs = 500) n_t0_500,
        min(t) FILTER (WHERE t0_lab AND fs = 500) first_t0_500
      FROM hl GROUP BY pid""")
    for v in T2_VARIANTS:
        S[f"case_{v}"] = one(con, f"""SELECT
            count(*) FILTER (WHERE n_t2_{v} >= 1 AND n_t0 >= 1),
            count(*) FILTER (WHERE n_t2_{v} >= 1 AND n_t0_500 >= 1),
            sum(n_t0_500) FILTER (WHERE n_t2_{v} >= 1 AND n_t0_500 >= 1),
            count(*) FILTER (WHERE n_t2_{v} >= 1 AND n_t0_500 >= 1 AND first_t0_500 < first_t2_{v}),
            count(*) FILTER (WHERE n_t2_{v}_500 >= 1 AND n_t0_500 >= 1)
            FROM hp""")
    S["overlap"] = one(con, """SELECT count(*) FILTER (WHERE n_t2_cons >= 1 AND n_t2_ph >= 1),
        count(*) FILTER (WHERE n_t2_ph >= 1 AND n_t2_cons = 0), count(*) FILTER (WHERE n_t2_v24 >= 1 AND n_t2_cons = 0)
        FROM hp WHERE n_t0_500 >= 1""")
    # K2 (fs kontrolü yok): SVT kodu, I45.6/426.7 yok, hiçbir kaynakta preeksitasyon yok, ≥1 T0
    S["k2"] = one(con, f"""SELECT count(*) FROM (SELECT pid, bool_or(pre) any_pre, count(*) FILTER (WHERE t0_lab) n_t0
        FROM he WHERE pid IN (SELECT pid FROM {q(f'heedb_{site}_icd_svt.parquet')})
          AND pid NOT IN (SELECT pid FROM {q(f'heedb_{site}_icd_pre.parquet')}) GROUP BY pid)
        WHERE NOT any_pre AND n_t0 >= 1""")[0]
    case_cols = ", ".join(f"(n_t2_{v} >= 1 AND n_t0_500 >= 1) case_{v}" for v in T2_VARIANTS)
    con.execute(f"""COPY (SELECT *, {case_cols} FROM hp WHERE n_t0 >= 1)
        TO '{P1 / f'heedb_{site}_cohort.parquet'}' (FORMAT parquet)""")
    con.execute(f"COPY (SELECT * FROM hl) TO '{P1 / f'heedb_{site}_ecg_labels.parquet'}' (FORMAT parquet)")
    return S


# ------------------------------------------------------------------ rapor
def report(M, H):
    cd, ca = M["case_def"], M["case_any"]
    pr_av = M["t0_pr_available"]
    md = ["# Aşama 1 — etiketler, kesin G1 ve kaynaklar arası uyum\n",
          "concept.md §5.2 tanımları; manuel okuma yok. Yalnızca agregat.\n",
          "## G1 özeti\n",
          md_table([
              ("MIMIC (dış test) — kesin T2", cd[0], cd[2], cd[3]),
              ("MIMIC — herhangi T2 (duyarlılık)", ca[0], ca[2], ca[3]),
          ] + [(f"{s['name']} — {T2_LABEL[v]}, 500 Hz T0", s[f"case_{v}"][1], s[f"case_{v}"][3], "—")
               for s in H for v in T2_VARIANTS],
              ["Tanım", "Vaka hastası", "T2'den önce T0'ı olan", "≥1 sıkı T0'ı olan"]),
          f"\nG1 eşikleri: geliştirme ≥100 (MGH, uzlaşı: {H[0]['case_cons'][1]:,}), dış test ≥40 (MIMIC kesin: {cd[0]:,}).\n",
          "## MIMIC ayrıntı\n",
          md_table([
              ("Kesin T2 vaka hastası", cd[0]), ("  T0 EKG sayısı", int(cd[1] or 0)),
              ("  ilk T2'den önce T0'ı olan (prospektif)", cd[2]),
              (f"  ≥1 sıkı T0 (PR ≥{PR_MIN_STRICT} ms, QRS <{QRS_MAX_STRICT} ms)", cd[3]),
              ("  ilk T2'den önce sıkı T0'ı olan", cd[4]),
              ("Herhangi T2 vaka hastası (duyarlılık)", ca[0]), ("  ≥1 sıkı T0", ca[3]),
              ("H3: kesin T2 ablasyon öncesi + ablasyon sonrası T0", M["h3"][0]),
              (f"  sonrası T0 0–{POST_ABL_EARLY_DAYS} gün", M["h3"][1]), (f"  sonrası T0 >{POST_ABL_EARLY_DAYS} gün", M["h3"][2]),
              ("H3, herhangi T2 ile", M["h3"][3]),
              ("K2 hastası (≥1 T0)", M["k2"][0]), ("  ≥1 sıkı T0", M["k2"][1]),
              ("  K2a (yalnızca AVNRT kalıbı)", M["k2"][2]), ("    ≥1 sıkı T0", M["k2"][3]),
              ("  K2b (yalnızca aksesuar yol kalıbı)", M["k2"][4]),
          ], ["Grup", "Hasta"]),
          f"\nT0 EKG'lerinde makine PR ölçümü olan: {pr_av[0]:,} / {pr_av[1]:,}.\n",
          "Kesin T2 vakalarının T2 + T0 EKG'lerinde filtre ayarı:\n",
          md_table(M["filter_case_def"], ["filtering", "bandwidth", "EKG"]), ""]
    for s in H:
        rows = []
        for v in T2_VARIANTS:
            c = s[f"case_{v}"]
            bold = "**" if v == "cons" else ""
            rows += [(f"{bold}T2 = {T2_LABEL[v]} + ≥1 500 Hz T0{bold}", c[1]),
                     ("  (fs kısıtı öncesi)", c[0]), ("  500 Hz T0 EKG sayısı", int(c[2] or 0)),
                     ("  ilk T2'den önce 500 Hz T0'ı olan (prospektif)", c[3]),
                     ("  ≥1 500 Hz T2'si de olan (H2 delta vektörü için)", c[4])]
        o = s["overlap"]
        rows += [("Uzlaşı ve overread vakası birlikte", o[0]), ("Yalnızca overread'e göre vaka (uzlaşı değil)", o[1]),
                 ("Yalnızca v24'e göre vaka (uzlaşı değil)", o[2]),
                 ("K2 hastası (≥1 overread'li T0; fs kontrol edilmedi)", s["k2"])]
        md += [f"## {s['name']} ({s['site']})\n", md_table(rows, ["Grup", "Hasta"]),
               f"\nBaşlığı okunamayan / beklenmeyen fs: {s['fs_missing']:,} EKG.\n",
               f"**Kaynaklar arası uyum** (üç kaynağın da olduğu {s['n_kappa']:,} EKG, preeksitasyon kodu):\n",
               md_table([(k, f"{v[0]:.3f}", f"{v[1]:.3f}") for k, v in s["kappa"].items()],
                        ["Çift", "Cohen κ", "Pozitif uyum"]),
               f"\nFleiss κ (üç kaynak): {s['fleiss']:.3f}. Desen sayıları (v24, yazılım, overread): {s['patterns']}.\n"]
    (OUT / "labels_summary.md").write_text("\n".join(md))
    (OUT / "labels_summary.json").write_text(json.dumps({"mimic": M, "heedb": H}, indent=1, default=str))
    print("\n".join(md))


def main():
    t0 = time.time()
    P1.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    con = connect(threads=24, memory="40GB")
    M = run_mimic(con)
    H = [run_heedb(con, s) for s in HEEDB_SITES]
    report(M, H)
    print(f"\nToplam süre {round(time.time() - t0)} s")


if __name__ == "__main__":
    main()
