"""HIDDEN-DELTA Aşama 5 — test girdilerinin hazırlanması (skor yok; yalnızca kovaryat ve sinyal özellikleri).

Kümeler (`sets.parquet`, sütun `analysis`):
  primary_k1, primary_k2      Aşama 3 eşleşmeleri (birincil vaka tanımı; IAB okunmuş hastalar çıkarılmış)
  k2a_vs_k1                   MIMIC, H4 (AVNRT → K1)
  prospective                 birincil vakalar; indeks = ilk T2'den önceki en erken T0; K1 1:5 (plan §8.2)
  var_any (MIMIC)             herhangi preeksitasyon ifadesi (plan §8.1)
  var_ph, var_v24 (HEEDB)     yalnızca overread / yalnızca v24 12SL (plan §8.1)
Eşleştirme kuralları Aşama 3 ile aynıdır (`s3_match.match`, `strata`).

İndeks EKG tablosu (`index.parquet` + `index_mb.npy`): özellikler, türetilmiş M0 özellikleri, ayrık dalga, sıkı T0
bayrağı (plan §7.1), yapısal kalp hastalığı kodu (§7.7), son T2'ye uzaklık (§5.4), filtre grubu (MIMIC).
T2 tablosu (`t2.parquet`): birincil vakaların T2 EKG'lerinde 8 derivasyonluk ilk 20 ms eğim ve QRS integral profili (§5.3–5.4).
MIMIC: `h3.parquet` + `h3_mb.npy` (yolak ablasyonu vekili; ablasyon öncesi/sonrası T0'lar, §7.5).
Emory: `alt.parquet` (12SL 303 kodlu, iki morfoloji kümeli EKG'ler; §7.3).

Çıktı: ~/data/hidden-delta/stage5/<site>/. Çalıştırma: python3 stage5/s5_prepare.py MIMIC I0006 I0001
"""
import sys

import duckdb
import numpy as np
import pandas as pd

import common5 as C
import s3_match as S3

STRUCT10 = r"^(I42|I43|Q2[0-8]|I44\.?[4-7]|I45\.?[0-2])"
STRUCT9 = r"^(425|74[5-7]|426\.?[2-5])"
MIMIC_MM = "/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0/machine_measurements.csv"


def load_site(site):
    df, case_ids = S3.load_mimic() if site == "MIMIC" else S3.load_heedb(site)
    if site in ("I0006", "MIMIC"):
        iab = pd.read_parquet(C.P3 / "iab_read_patients.parquet")
        drop = set(iab.loc[iab.dataset == site, "pid"])
        df, case_ids = df[~df.pid.isin(drop)], case_ids - drop
    con = duckdb.connect()
    if site == "MIMIC":                                       # varyant T2: herhangi ifade
        df["t2_any"] = df.role.eq("case_t2") | df.role.eq("abl_t2")
        var_cases = {"var_any": set(con.execute(f"SELECT sid FROM '{C.P1 / 'mimic_cohort.parquet'}' WHERE case_any").df().sid)}
        var_t2 = {"var_any": "t2_any"}
    else:
        lab = con.execute(f"SELECT k, t2_ph, t2_v24 FROM '{C.P1 / f'heedb_{site}_ecg_labels.parquet'}'").df()
        df = df.merge(lab, on="k", how="left")
        for c in ("t2_ph", "t2_v24"):
            df[c] = df[c].fillna(False).astype(bool)
        coh = con.execute(f"SELECT pid, case_ph, case_v24 FROM '{C.P1 / f'heedb_{site}_cohort.parquet'}'").df()
        var_cases = {"var_ph": set(coh.pid[coh.case_ph]), "var_v24": set(coh.pid[coh.case_v24])}
        var_t2 = {"var_ph": "t2_ph", "var_v24": "t2_v24"}
    if site in ("I0006", "MIMIC"):
        for v in var_cases:
            var_cases[v] -= drop
    df = df.sort_values("k").reset_index(drop=True)
    return df, case_ids, var_cases, var_t2


def first_t2(df, col):
    return df[df[col]].groupby("pid").t.min()


def index_cases(df, case_ids, t2col, rng, prospective=False):
    ft2 = first_t2(df, t2col)
    t0 = df[df.t0_idx & df.pid.isin(case_ids) & df.age.notna()].copy()
    t0["first_t2v"] = t0.pid.map(ft2)
    if prospective:
        t0 = t0[t0.t < t0.first_t2v].sort_values(["pid", "t", "k"]).drop_duplicates("pid")
    else:
        t0 = t0.sample(frac=1.0, random_state=int(rng.integers(2**31))).drop_duplicates("pid")
    return t0


def match_sets(cases, pool, ratio, rng, label):
    m = S3.match(cases, pool, ratio, rng)
    ctrl = pool.merge(m[["set_id", "k"]], on="k")
    mc = cases[cases.pid.isin(m.set_id)]
    return pd.concat([mc[["pid", "k"]].assign(set_id=mc.pid, grp="case"),
                      ctrl[["pid", "k", "set_id"]].assign(grp="k1")]).assign(analysis=label)


def machine_intervals_mimic(studies):
    """Makine PR ve QRS süresi (ms) — sıkı T0 denetimi (plan §7.1)."""
    con = duckdb.connect()
    con.register("st", pd.DataFrame({"study": [int(x) for x in studies]}))
    return con.execute(f"""SELECT study_id::BIGINT::VARCHAR AS k,
          CASE WHEN try_cast(p_onset AS INT) < 29999 AND try_cast(qrs_onset AS INT) < 29999
               THEN try_cast(qrs_onset AS INT) - try_cast(p_onset AS INT) END AS m_pr,
          CASE WHEN try_cast(qrs_end AS INT) < 29999 AND try_cast(qrs_onset AS INT) < 29999
               THEN try_cast(qrs_end AS INT) - try_cast(qrs_onset AS INT) END AS m_qrs
        FROM read_csv('{MIMIC_MM}', all_varchar=true)
        WHERE study_id::BIGINT IN (SELECT study FROM st)""").df()


def structural_patients(site, pids):
    con = duckdb.connect()
    con.register("pp", pd.DataFrame({"pid": list(pids)}))
    if site == "MIMIC":
        q = f"""SELECT DISTINCT subject_id::BIGINT pid FROM
                read_csv('/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp/diagnoses_icd.csv.gz', all_varchar=true)
                WHERE subject_id::BIGINT IN (SELECT pid FROM pp) AND
                ((icd_version = '10' AND regexp_matches(icd_code, '{STRUCT10}')) OR
                 (icd_version = '9' AND regexp_matches(icd_code, '{STRUCT9}')))"""
        return set(con.execute(q).df().pid)
    base = f"/home/fatih/datasets/heedb/ECG/{site}/ICD_codes"
    out = set()
    for f, col, rx in (("icd10_codes.csv", "DIAGNOSIS_ICD10_CD", STRUCT10), ("icd9_codes.csv", "DIAGNOSIS_ICD9_CD", STRUCT9)):
        out |= set(con.execute(f"""SELECT DISTINCT try_cast(try_cast(BDSPPatientID AS DOUBLE) AS BIGINT) pid
            FROM read_csv('{base}/{f}', all_varchar=true)
            WHERE try_cast(try_cast(BDSPPatientID AS DOUBLE) AS BIGINT) IN (SELECT pid FROM pp)
              AND regexp_matches({col}, '{rx}')""").df().pid)
    return out


def main(site):
    out = C.P5 / site
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(C.SEED + {"MIMIC": 1, "I0006": 2, "I0001": 3}[site])
    df, case_ids, var_cases, var_t2 = load_site(site)
    df = S3.strata(df)
    pool = df[(df.role == "k1") & df.age.notna()]
    sets = [pd.read_parquet(C.P3 / f"{site}_matched.parquet").assign(analysis=lambda d: d["pool"].map(
        {"k1": "primary_k1", "k2": "primary_k2", "k2a_vs_k1": "k2a_vs_k1"}))[["pid", "k", "set_id", "grp", "analysis"]]]
    sets.append(match_sets(index_cases(df, case_ids, "t2_prim", rng, prospective=True), pool, 5, rng, "prospective"))
    for v, cset in var_cases.items():
        sets.append(match_sets(index_cases(df, cset, var_t2[v], rng), pool, 5, rng, v))
    sets = pd.concat(sets, ignore_index=True)
    sets.to_parquet(out / "sets.parquet", index=False)

    # indeks EKG'leri
    idx = df[df.k.isin(set(sets.k))].drop_duplicates("k").reset_index(drop=True)
    mb = C.load_beats(site, idx.k.tolist())
    idx = C.add_derived(idx, mb)
    ft2 = first_t2(df, "t2_prim")
    t2_times = df[df.t2_prim].groupby("pid").t.apply(list)
    idx["days_to_t2"] = [min((abs((t - x).days) for x in t2_times[p]), default=np.nan) if p in t2_times.index else np.nan
                         for p, t in zip(idx.pid, idx.t)]
    idx["before_first_t2"] = idx.t < idx.pid.map(ft2)
    if site == "MIMIC":
        mi = machine_intervals_mimic(idx.k)
        idx = idx.merge(mi, on="k", how="left")
        idx["strict_t0"] = (idx.m_pr >= 140) & (idx.m_qrs < 110)
    else:
        idx["strict_t0"] = (idx.qrs_ms < 110) & ((idx.pr_ms - 42) >= 140)
    idx["structural"] = idx.pid.isin(structural_patients(site, set(idx.pid)))
    idx.to_parquet(out / "index.parquet", index=False)
    np.save(out / "index_mb.npy", mb)

    # T2 profilleri (birincil vakalar)
    cases = set(sets.loc[(sets.grp == "case") & sets.analysis.str.startswith("primary"), "pid"])
    t2 = df[df.t2_prim & df.pid.isin(cases)].reset_index(drop=True)
    t2mb = C.load_beats(site, t2.k.tolist())
    ok = ~np.isnan(t2mb[:, 0, 0])
    t2, t2mb = t2[ok].reset_index(drop=True), t2mb[ok]
    prof = pd.DataFrame([C.integral_profile(t2mb[i], t2.on_low_ms[i], t2.off_ms[i]) for i in range(len(t2))],
                        columns=[f"qint_{x}" for x in C.L8])
    t2 = pd.concat([t2[["pid", "k", "t"] + [f"slope20_{x}" for x in C.L8]], prof], axis=1)
    t2.to_parquet(out / "t2.parquet", index=False)

    msg = [f"{site}: kümeler {sets.groupby('analysis').set_id.nunique().to_dict()}",
           f"indeks EKG {len(idx):,}; T2 profili {len(t2):,} EKG ({t2.pid.nunique():,} vaka)"]
    if site == "MIMIC":                                          # H3
        coh = duckdb.connect().execute(f"""SELECT sid AS pid, d_abl FROM '{C.P1 / 'mimic_cohort.parquet'}'
                                           WHERE d_abl IS NOT NULL AND n_t2def_before > 0""").df()
        h3 = df[df.pid.isin(set(coh.pid)) & df.role.isin(["case_t0", "abl_t0"]) & df.t0_flag.fillna(False)].drop(columns=["d_abl"]).merge(coh, on="pid")
        h3["days_from_abl"] = (h3.t - pd.to_datetime(h3.d_abl)).dt.days
        h3 = h3.reset_index(drop=True)
        np.save(out / "h3_mb.npy", C.load_beats(site, h3.k.tolist()))
        h3[["pid", "k", "t", "days_from_abl"]].to_parquet(out / "h3.parquet", index=False)
        msg.append(f"H3: {h3.pid.nunique()} hasta, {len(h3)} EKG")
    if site == "I0006":                                          # §7.3
        alt = duckdb.connect().execute(f"""SELECT k FROM '{C.P0 / 'heedb_I0006_ecg_flags.parquet'}' WHERE alt""").df()
        a = df[df.k.isin(set(alt.k)) & df.c2_qrs_ms.notna()][["pid", "k", "qrs_ms", "c2_qrs_ms", "pqlow_mag", "c2_pqlow_mag",
                                                               "on_low_ms", "c2_on_low_ms"]]
        a.to_parquet(out / "alt.parquet", index=False)
        msg.append(f"ALTWPW iki kümeli EKG: {len(a)}")
    print(" | ".join(msg))





# ------------------------------------------------------------------ §8.4: hasta düzeyinde ortalama/maksimum için girdiler
def prepare_agg(site):
    """Birincil K1 kümelerinde: vakanın tüm T0'ları; her kontrolden vakanın T0 sayısı kadar rastgele T0 (yerine koymadan;
    kontrolün daha az T0'ı varsa hepsi). Çıkarılmamış kontrol EKG'leri Aşama 2 koduyla (etiketsiz) işlenir.
    Çıktı: agg.parquet (set_id, pid, grp, k) + agg_mb.npy."""
    from concurrent.futures import ProcessPoolExecutor
    import s2_extract as E
    out = C.P5 / site
    rng = np.random.default_rng(C.SEED + 40 + {"MIMIC": 1, "I0006": 2, "I0001": 3}[site])
    df, _, _, _ = load_site(site)
    sets = pd.read_parquet(out / "sets.parquet")
    sets = sets[sets.analysis == "primary_k1"]
    cases = sets[sets.grp == "case"]
    case_t0 = df[df.t0_idx & df.pid.isin(set(cases.pid))].sort_values("k")
    n_t0 = case_t0.groupby("pid").size()
    con = duckdb.connect()
    if site == "MIMIC":
        fl = con.execute(f"SELECT sid AS pid, study::VARCHAR AS k FROM '{C.P0 / 'mimic_ecg_flags.parquet'}' WHERE t0").df()
        rl = con.execute("""SELECT study_id::BIGINT::VARCHAR AS k, path FROM
            read_csv('/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0/record_list.csv', all_varchar=true)""").df()
        rec_of = dict(zip(rl.k, "/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0/" + rl.path))
        ok500 = lambda k: True
    else:
        fl = con.execute(f"SELECT pid, k FROM '{C.P0 / f'heedb_{site}_ecg_flags.parquet'}' WHERE t0 AND has_ph").df()
        rec_of = None
        ok500 = lambda k: E.read_fs(str(E.HEEDB / site / "WFDB" / k)) == 500
    ctrl = sets[sets.grp == "k1"]
    fl = fl[fl.pid.isin(set(ctrl.pid))].sort_values(["pid", "k"])
    by_pid = {p: g.k.tolist() for p, g in fl.groupby("pid")}
    rows = [case_t0[["pid", "k"]].assign(set_id=case_t0.pid, grp="case")]
    for r in ctrl.itertuples():
        cand = [k for k in rng.permutation(by_pid.get(r.pid, [r.k])).tolist()]
        need = int(n_t0.get(r.set_id, 1))
        pick = []
        for k in cand:
            if ok500(k):
                pick.append(k)
            if len(pick) == need:
                break
        rows.append(pd.DataFrame({"pid": r.pid, "k": pick or [r.k], "set_id": r.set_id, "grp": "k1"}))
    agg = pd.concat(rows, ignore_index=True)
    have = set(C.load_features(site, agg.k.unique()).k)
    extra = sorted(set(agg.k) - have)
    recs = [rec_of[k] if rec_of else str(E.HEEDB / site / "WFDB" / k) for k in extra]
    with ProcessPoolExecutor(28) as ex:
        res = list(ex.map(E.work, recs, chunksize=8))
    ok_extra = {k: r[1] for k, r in zip(extra, res) if r[0].get("status") == "ok"}
    keys = agg.k.tolist()
    mb = C.load_beats(site, keys)
    for i, k in enumerate(keys):
        if k in ok_extra:
            mb[i] = ok_extra[k]
    good = ~np.isnan(mb[:, 0, 0])
    agg, mb = agg[good].reset_index(drop=True), mb[good]
    agg.to_parquet(out / "agg.parquet", index=False)
    np.save(out / "agg_mb.npy", mb)
    print(f"{site} §8.4: {len(agg):,} EKG (vaka {int((agg.grp == 'case').sum()):,}, kontrol {int((agg.grp == 'k1').sum()):,}); "
          f"yeni çıkarılan {len(ok_extra):,}/{len(extra):,}")


# ------------------------------------------------------------------ §8.7: ham 10 s girişli M1 için indeks EKG dizileri
def prepare_raw(site):
    """İndeks EKG'lerinin taban çizgisi düzeltilmiş ham 10 s dizileri (etiketsiz). Çıktı: index_raw.npy, index_raw_ok.npy."""
    import s4_m1raw as R
    import s2_extract as E
    out = C.P5 / site
    idx = pd.read_parquet(out / "index.parquet")
    if site == "MIMIC":
        rl = duckdb.connect().execute("""SELECT study_id::BIGINT::VARCHAR AS k, path FROM
            read_csv('/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0/record_list.csv', all_varchar=true)""").df()
        rec = dict(zip(rl.k, "/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0/" + rl.path))
        recs = [rec[k] for k in idx.k]
    else:
        recs = [str(E.HEEDB / site / "WFDB" / k) for k in idx.k]
    arr, ok = R.raw_array(recs)
    np.save(out / "index_raw.npy", arr)
    np.save(out / "index_raw_ok.npy", ok)
    print(f"{site} §8.7: {ok.sum():,}/{len(ok):,} ham indeks kaydı")


if __name__ == "__main__":
    args = sys.argv[1:] or ["MIMIC", "I0006"]
    if args[0] in ("agg", "raw"):                        # python3 s5_prepare.py agg|raw MIMIC I0006
        for s in args[1:]:
            (prepare_agg if args[0] == "agg" else prepare_raw)(s)
    else:
        for s in args:
            main(s)
