"""HIDDEN-DELTA Aşama 3 — vaka-kontrol eşleştirmesi (concept.md §5.3 "Eşleştirme uygulaması"; ANALYSIS_PLAN §2–3).

Vaka: birincil tanım (HEEDB: üç kaynaktan ≥2'si; MIMIC: kesin ifade). İndeks EKG: hasta başına bir rastgele başarılı
500 Hz T0 (yaş bilinen; MIMIC'te ablasyon öncesi). Kontrol: K1 (1:5) ve K2 (1:2) havuzlarının çıkarılmış tek T0'ı.
MIMIC'te ayrıca H4 için AVNRT (K2a) hastaları K1'e 1:5 eşleştirilir.
Kesin tabakalar: cinsiyet, KH bandı, EKG sayısı bandı; HEEDB'de dönem (WFDB klasör yılı, 5 yıllık bant);
MIMIC'te cihaz filtre ayarı. Kaliperler: yaş ±3 yıl, log(EKG sayısı) ±0.25. Açgözlü, yerine koymadan.

Test setlerinde (I0006, MIMIC) yalnızca denge raporlanır; vaka-kontrol özellik karşılaştırması yapılmaz
(ANALYSIS_PLAN §2: dondurmaya kadar etiketli analiz yok). MGH'de betimsel karşılaştırma da yazılır.

Çıktı: stage3/out/match_<site>.md (agregat); ~/data/hidden-delta/stage3/<site>_matched.parquet (hasta düzeyi).
Çalıştırma: python3 stage3/s3_match.py [I0001|I0006|MIMIC]
"""
import argparse
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

P0 = Path("/home/fatih/data/hidden-delta/stage0")
P1 = Path("/home/fatih/data/hidden-delta/stage1")
P2R = Path("/home/fatih/data/hidden-delta/stage2")
P3 = Path("/home/fatih/data/hidden-delta/stage3")
HEEDB = Path("/home/fatih/datasets/heedb/ECG")
MIMIC_HOSP = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp")
OUT = Path(__file__).resolve().parent / "out"
SEED = 20261116
CALIPER = 3.0
CALIPER_LOGN = 0.25   # sapma (2026-11-08): EKG sayısı dengesi için log(EKG sayısı) kaliperi (±%28)
RATIO = {"k1": 5, "k2": 2}
N_BOOT = 2000
NAME = {"I0001": "MGH", "I0006": "Emory", "MIMIC": "MIMIC"}
COMPARE_ALLOWED = {"I0001"}


def band(x, edges, labels):
    return pd.cut(x, bins=edges, labels=labels, right=False).astype(str)


def sex_heedb(con, site, col):
    meta = HEEDB / site / "metadata" / "metadata.csv"
    # deterministik: hasta başına en sık cinsiyet; eşitlikte alfabetik ilk
    return con.execute(f"""SELECT pid, arg_max(sex, cnt * 1000000 - rk) sex FROM (
          SELECT pid, sex, count(*) cnt, dense_rank() OVER (PARTITION BY pid ORDER BY sex) rk FROM (
            SELECT try_cast(try_cast(BDSPPatientID AS DOUBLE) AS BIGINT) pid, {col} sex
            FROM read_csv('{meta}', all_varchar=true)
            WHERE try_cast(try_cast(BDSPPatientID AS DOUBLE) AS BIGINT) IN (SELECT pid FROM pids) AND {col} IS NOT NULL)
          GROUP BY pid, sex) GROUP BY pid""").df()


def load_heedb(site):
    con = duckdb.connect()
    P2 = P2R / site
    con.execute(f"""CREATE TABLE f AS SELECT * FROM read_parquet('{P2}/features_*.parquet', union_by_name=true)
                    WHERE status = 'ok'""")
    con.execute(f"""CREATE TABLE lab AS SELECT k, t, age, t0_lab, t2_cons FROM '{P1 / f'heedb_{site}_ecg_labels.parquet'}'""")
    con.execute(f"""CREATE TABLE flags AS SELECT pid, k, t, age FROM '{P0 / f'heedb_{site}_ecg_flags.parquet'}'""")
    con.execute("CREATE TABLE necg AS SELECT pid, count(*) n_ecg FROM flags GROUP BY pid")
    cases = con.execute(f"SELECT pid FROM '{P1 / f'heedb_{site}_cohort.parquet'}' WHERE case_cons").df().pid
    ypart = 2 if site == "I0001" else 1                      # MGH: S0001/<yıl>/..., Emory: <yıl>/...
    df = con.execute(f"""SELECT f.*, fl.t, fl.age, n.n_ecg,
            try_cast(split_part(f.k, '/', {ypart}) AS INT) AS fyear, year(fl.t) AS shifted_year,
            l.t0_lab, l.t2_cons
        FROM f JOIN flags fl USING (k) JOIN necg n ON n.pid = f.pid LEFT JOIN lab l USING (k)""").df()
    con.register("pids", pd.DataFrame({"pid": df.pid.unique()}))
    sex = sex_heedb(con, site, "SexDSC" if site == "I0001" else "Sex")
    sex["sex"] = sex.sex.str.capitalize()                     # Emory: MALE/FEMALE → Male/Female
    df = df.merge(sex, on="pid", how="left")
    df["t0_idx"] = (df.role == "case_t0") & df.t0_lab.fillna(False)
    df["t2_prim"] = df.t2_cons.fillna(False)
    df["filt"] = "NA"
    return df, set(cases)


def load_mimic():
    con = duckdb.connect()
    con.execute(f"""CREATE TABLE f AS SELECT * FROM read_parquet('{P2R / 'MIMIC'}/features_*.parquet', union_by_name=true)
                    WHERE status = 'ok'""")
    con.execute(f"CREATE TABLE flags AS SELECT * FROM '{P0 / 'mimic_ecg_flags.parquet'}'")
    con.execute("CREATE TABLE necg AS SELECT sid, count(*) n_ecg FROM flags GROUP BY sid")
    con.execute(f"CREATE TABLE lab AS SELECT study, t2_def FROM '{P1 / 'mimic_ecg_labels.parquet'}'")
    cases = con.execute(f"SELECT sid FROM '{P1 / 'mimic_cohort.parquet'}' WHERE case_def").df().sid
    df = con.execute(f"""SELECT f.*, fl.t, fl.age, n.n_ecg, NULL::INT AS fyear, year(fl.t) AS shifted_year,
            fl.t0 AS t0_flag, coalesce(l.t2_def, false) AS t2_prim,
            coalesce(fl.filtering, 'NA') || ' | ' || coalesce(fl.bandwidth, 'NA') AS filt,
            a.d_abl
        FROM f JOIN flags fl ON fl.study = f.k::BIGINT JOIN necg n ON n.sid = f.pid
        LEFT JOIN lab l ON l.study = f.k::BIGINT
        LEFT JOIN '{P0 / 'mimic_abl.parquet'}' a ON a.sid = f.pid""").df()
    sex = con.execute(f"""SELECT subject_id::BIGINT pid, CASE gender WHEN 'M' THEN 'Male' WHEN 'F' THEN 'Female' END sex
        FROM read_csv('{MIMIC_HOSP / 'patients.csv.gz'}', all_varchar=true)""").df()
    df = df.merge(sex, on="pid", how="left")
    # indeks T0: ablasyon öncesi (concept §5.3)
    pre_abl = df.d_abl.isna() | (df.t < pd.to_datetime(df.d_abl))
    df["t0_idx"] = (df.role == "case_t0") & df.t0_flag.fillna(False) & pre_abl
    k2 = pd.read_parquet(P1.parent / "stage0" / "mimic_k2.parquet")
    k2a = set(k2[(k2.avnrt == True) & (k2.ap != True) & (k2.wpw != True)].sid)  # noqa: E712
    df["k2a"] = df.pid.isin(k2a)
    return df, set(cases)


def strata(d):
    d = d.copy()
    d["era"] = (d.fyear // 5 * 5).astype("Int64").astype(str)
    d["hr_band"] = band(d.hr, [0, 60, 80, 100, 1e9], ["<60", "60–79", "80–99", "≥100"])
    # sapma (2026-11-08): ilk sürümde ≥13 tek banttı; EKG sayısı dengesi hedefi tutmadı (SMD 0.30) → üst bantlar inceltildi
    d["necg_band"] = band(d.n_ecg, [0, 4, 7, 13, 25, 50, 1e9], ["2–3", "4–6", "7–12", "13–24", "25–49", "≥50"])
    d["stratum"] = d.sex.fillna("NA") + "|" + d.era + "|" + d.hr_band + "|" + d.necg_band + "|" + d.filt
    return d


def match(cases, pool, ratio, rng):
    pool = pool.copy()
    used = np.zeros(len(pool), bool)
    by = {s: np.where(pool.stratum.values == s)[0] for s in pool.stratum.unique()}
    rows = []
    for i in rng.permutation(len(cases)):
        c = cases.iloc[i]
        idx = by.get(c.stratum, np.array([], int))
        idx = idx[~used[idx]]
        if len(idx) == 0:
            continue
        d_age = np.abs(pool.age.values[idx] - c.age)
        d_n = np.abs(np.log(pool.n_ecg.values[idx]) - np.log(c.n_ecg))
        okm = (d_age <= CALIPER) & (d_n <= CALIPER_LOGN)
        ok = idx[okm]
        if len(ok) == 0:
            continue
        score = d_age[okm] / CALIPER + d_n[okm] / CALIPER_LOGN      # iki boyutlu, kaliperlere göre ölçekli
        take = ok[np.argsort(score, kind="stable")[:ratio]]
        used[take] = True
        rows.append(pd.DataFrame({"set_id": c.pid, "pid": pool.pid.values[take], "k": pool.k.values[take]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["set_id", "pid", "k"])


def smd(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return float((a.mean() - b.mean()) / s) if s > 0 else 0.0


def balance(case, ctrl, site):
    rows = []
    vars_ = [("age", "yaş"), ("hr", "KH"), ("n_ecg", "EKG sayısı")]
    if site != "MIMIC":
        vars_.append(("fyear", "dönem (klasör yılı)"))
    vars_.append(("shifted_year", "kaydırılmış yıl"))
    for v, name in vars_:
        rows.append((name, round(case[v].mean(), 2), round(ctrl[v].mean(), 2), round(smd(case[v], ctrl[v]), 3)))
    lg = lambda d: np.log(d.n_ecg.astype(float))
    rows.append(("log(EKG sayısı)", round(lg(case).mean(), 2), round(lg(ctrl).mean(), 2), round(smd(lg(case), lg(ctrl)), 3)))
    male = lambda d: (d.sex == "Male").astype(float)
    rows.append(("erkek (%)", round(100 * male(case).mean(), 1), round(100 * male(ctrl).mean(), 1),
                 round(smd(male(case), male(ctrl)), 3)))
    return rows


FEATS = [("septal_q_yok", "Septal q yok (I, aVL, V5, V6)", True), ("qv6", "V6'da q var", True),
         ("slow_init_ms", "Yavaş başlangıç (ms)", False), ("qrs_ms", "QRS süresi (ms)", False),
         ("pr_ms", "PR (ms; +42 ms sapmalı)", False), ("v1_rs", "V1 R/S", False),
         ("qrst_angle", "QRS-T açısı (°)", False), ("pqlow_mag", "QRS öncesi genlik, düşük eşik (mV)", False),
         ("pqhigh_mag", "QRS öncesi genlik, görünür başlangıç (mV)", False)]


def matched_compare(case, ctrl, rng):
    """Eşleştirilmiş kümeler düzeyinde: vaka değeri − kümedeki kontrollerin ortalaması; küme bootstrap'ı."""
    rows = []
    for col, name, binary in FEATS:
        cm = ctrl.groupby("set_id")[col].mean()
        cv = case.set_index("pid")[col].reindex(cm.index)
        d = (cv - cm).dropna().to_numpy()
        boots = [rng.choice(d, len(d)).mean() for _ in range(N_BOOT)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        f = (lambda x: f"{100 * x:.1f}%") if binary else (lambda x: f"{x:.3g}")
        fd = (lambda x: f"{100 * x:+.1f} puan") if binary else (lambda x: f"{x:+.3g}")
        rows.append((name, f(np.nanmean(cv)), f(np.nanmean(ctrl[col])), fd(d.mean()), f"[{fd(lo)}, {fd(hi)}]",
                     round(smd(cv.dropna(), ctrl[col].dropna()), 3)))
    return rows


def section(md, title, cases_all, pool, cases_m, ctrl, ratio, site):
    sizes = ctrl.groupby("set_id").size()
    md += [f"## {title} (1:{ratio})\n",
           f"Havuz {len(pool):,}; eşleşen {len(cases_m):,} / {len(cases_all):,} "
           f"(%{100 * len(cases_m) / max(1, len(cases_all)):.1f}); kontrol {len(ctrl):,}; tam (1:{ratio}) küme "
           f"{int((sizes == ratio).sum()):,}.\n",
           "**Denge** (SMD; hedef |SMD| < 0.1):\n",
           "| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |",
           "|---|---|---|---|---|---|---|"]
    for (n1, a1, b1, s1), (_, a2, b2, s2) in zip(balance(cases_all, pool, site), balance(cases_m, ctrl, site)):
        md.append(f"| {n1} | {a1} | {b1} | {s1} | {a2} | {b2} | {s2} |")
    md.append("")
    return dict(pool=len(pool), matched=len(cases_m), of=len(cases_all), controls=len(ctrl),
                full_sets=int((sizes == ratio).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("site", nargs="?", default="I0001", choices=["I0001", "I0006", "MIMIC"])
    a = ap.parse_args()
    site = a.site
    rng = np.random.default_rng(SEED)
    df, case_ids = load_mimic() if site == "MIMIC" else load_heedb(site)
    n_iab = 0
    if site in ("I0006", "MIMIC"):                            # ANALYSIS_PLAN v1.2: IAB okunmuş hastalar dışlanır
        iab = pd.read_parquet(P3 / "iab_read_patients.parquet")
        iab_p = set(iab.loc[iab.dataset == site, "pid"])
        n_iab = int(df.pid.isin(iab_p).groupby(df.pid).any().sum())
        df = df[~df.pid.isin(iab_p)]
        case_ids = case_ids - iab_p
    df["septal_q_yok"] = (df[["q_I", "q_aVL", "q_V5", "q_V6"]] == 0).all(axis=1).astype(float)
    df["qv6"] = (df.q_V6 > 0).astype(float)
    first_t2 = df[df.t2_prim].groupby("pid").t.min().rename("first_t2")
    df = df.merge(first_t2, on="pid", how="left")
    df = df.sort_values("k").reset_index(drop=True)          # DuckDB satır sırası sabit değil → örnekleme öncesi sırala
    df = strata(df)
    t0 = df[df.t0_idx & df.pid.isin(case_ids) & df.age.notna()]
    t0 = t0.assign(prospective=t0.t < t0.first_t2)
    case_idx = t0.sample(frac=1.0, random_state=int(rng.integers(2**31))).drop_duplicates("pid")
    pro_pids = set(t0.loc[t0.prospective, "pid"])
    compare = site in COMPARE_ALLOWED
    res, md = {}, [f"# Aşama 3 — eşleştirme ({NAME[site]})\n",
                   "concept.md §5.3 \"Eşleştirme uygulaması\"; ANALYSIS_PLAN §2–3. Yalnızca agregat.\n",
                   f"Birincil vaka hastası (uygun indeks T0'ı olan): {len(case_idx):,}; ilk T2'den önce T0'ı olan: {len(pro_pids):,}.\n"]
    if not compare:
        md.append("Test seti: yalnızca denge raporlanır; vaka-kontrol özellik karşılaştırması dondurmadan sonra.\n")
        md.append(f"Önceki IAB çalışmalarında EKG'si okunmuş {n_iab:,} hasta eşleştirmeden önce çıkarıldı (ANALYSIS_PLAN v1.2).\n")
    out_sets = []
    for role in ("k1", "k2"):
        pool = df[(df.role == role) & df.age.notna()]
        m = match(case_idx, pool, RATIO[role], rng)
        ctrl = pool.merge(m[["set_id", "k"]], on="k")
        matched_cases = case_idx[case_idx.pid.isin(m.set_id)]
        res[role] = section(md, role.upper(), case_idx, pool, matched_cases, ctrl, RATIO[role], site)
        if compare:
            md += ["**Eşleşme sonrası betimsel karşılaştırma** (vaka − küme kontrol ortalaması; küme bootstrap %95 GA):\n",
                   "| Özellik | Vaka | Kontrol | Fark | %95 GA | SMD |", "|---|---|---|---|---|---|"]
            md += [f"| {a} | {b} | {c} | {d} | {e} | {f} |" for a, b, c, d, e, f in matched_compare(matched_cases, ctrl, rng)]
            md.append("")
        out_sets.append(pd.concat([matched_cases[["pid", "k"]].assign(set_id=matched_cases.pid, grp="case",
                                                                      prospective=matched_cases.pid.isin(pro_pids)),
                                   ctrl[["pid", "k", "set_id"]].assign(grp=role, prospective=False)]).assign(pool=role))
    if site == "MIMIC":                                          # H4: AVNRT (K2a) → K1, 1:5
        k2a_idx = df[(df.role == "k2") & df.k2a & df.age.notna()].drop_duplicates("pid")
        pool = df[(df.role == "k1") & df.age.notna()]
        m = match(k2a_idx, pool, RATIO["k1"], rng)
        ctrl = pool.merge(m[["set_id", "k"]], on="k")
        k2a_m = k2a_idx[k2a_idx.pid.isin(m.set_id)]
        res["k2a_vs_k1"] = section(md, "H4: AVNRT (K2a) → K1", k2a_idx, pool, k2a_m, ctrl, RATIO["k1"], site)
        out_sets.append(pd.concat([k2a_m[["pid", "k"]].assign(set_id=k2a_m.pid, grp="k2a", prospective=False),
                                   ctrl[["pid", "k", "set_id"]].assign(grp="k1", prospective=False)]).assign(pool="k2a_vs_k1"))
    note = ("HEEDB'de cart/filtre bilgisi yok; dönem, klasör yılından vekil olarak alındı." if site != "MIMIC"
            else "MIMIC'te dönem tabakası yok (kaydırılmış tarih); cihaz filtre ayarı kesin tabakadır. İndeks T0 ablasyon öncesidir.")
    prefix = "Betimseldir; karar ve hipotez testi Aşama 4'te. " if compare else ""
    md.append(f"**Not.** {prefix}{note}")
    P3.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.concat(out_sets, ignore_index=True).to_parquet(P3 / f"{site}_matched.parquet", index=False)
    (OUT / f"match_{site}.md").write_text("\n".join(md) + "\n")
    (OUT / f"match_{site}.json").write_text(json.dumps(res, indent=1))
    print("\n".join(md))


if __name__ == "__main__":
    main()
