"""HIDDEN-DELTA Aşama 5 — kayıtlı analizler (ANALYSIS_PLAN v1.2; OSF https://osf.io/kvn28).

Site başına (`site` modu): H1 ve ikincil sonlanımlar, H2 testleri (5.1–5.4; ham tek yönlü p), 5.5, H4 (MIMIC), ayrık dalga
(koşullu lojistik), denetimler 7.1, 7.2 (MIMIC), 7.3 (Emory), 7.5 (MIMIC), 7.7, duyarlılık analizleri (§8: prospektif,
etiket varyantları, kovaryat düzeltmeli, yaş < 40). Birleştirme (`combine` modu): H2 ailesinde Holm, karar kuralları,
sonuç raporu.

Prova (`--rehearsal`): ~/data/hidden-delta/stage5/rehearsal/<site>/. Test setlerinde her eşleşmiş kümenin "vakası" rastgele
bir üyeyle değiştirilir (skorlar zaten rastgele modeldendir); MGH provası gerçek etiketle yapılır (geliştirme verisi).
Resmî çalıştırma bir kez yapılır; çıktı varsa çalışmaz.

Çalıştırma: python3 stage5/s5_analysis.py site MIMIC [--rehearsal];  python3 stage5/s5_analysis.py combine [--rehearsal]
"""
import argparse
import json
from datetime import date

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

import common5 as C
import statlib as S

TEST = {"MIMIC", "I0006"}
NB, NP = 2000, 10000


def base(site, rehearsal):
    return (C.P5 / "rehearsal" / site) if rehearsal else (C.P5 / site)


def load(site, rehearsal):
    src, b = C.P5 / site, base(site, rehearsal)
    sets = pd.read_parquet(src / "sets.parquet")
    idx = pd.read_parquet(src / "index.parquet")
    sc = pd.read_parquet(b / "scores.parquet")
    keep = ["k", "pid", "age", "male", "hr", "n_ecg", "shifted_year", "fyear", "filt", "strict_t0", "structural",
            "days_to_t2", "before_first_t2", "wave", "pqlow_fq_corr", "on_low_ms", "off_ms"]
    rows = sets.merge(idx[[c for c in keep if c in idx.columns]], on=["pid", "k"], how="left").merge(sc, on="k", how="left")
    if rehearsal and site in TEST:                       # prova: vaka kümede rastgele bir üyeyle değiştirilir
        rng = np.random.default_rng(C.SEED + 99)
        new = []
        for _, g in rows.groupby(["analysis", "set_id"], sort=False):
            g = g.copy()
            j = rng.integers(len(g))
            g["grp"] = np.where(np.arange(len(g)) == j, "case", "k1")
            new.append(g)
        rows = pd.concat(new, ignore_index=True)
    rows["y"] = rows.grp.isin(["case", "k2a"]).astype(int)
    return rows, idx, sc


def auc_block(d, col="logit_M1", seed=0):
    d = d.dropna(subset=[col])
    a, ci = S.auc_boot(d.y, d[col], d.set_id, NB, seed)
    return {"n_case": int(d.y.sum()), "n_ctrl": int((1 - d.y).sum()), "auroc": a, "ci": ci.tolist()}


def p_one_sided(vals):
    vals = np.asarray(vals)
    return float((np.sum(vals <= 0) + 1) / (len(vals) + 1))


def perm_agreement(A, B, rng):
    """Satır başına Pearson r (A_i, B_i); permütasyon: B satırları karıştırılır. Tek yönlü p."""
    def mean_r(Bm):
        a = A - A.mean(1, keepdims=True)
        b = Bm - Bm.mean(1, keepdims=True)
        return float(np.nanmean((a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))))
    obs = mean_r(B)
    perms = np.array([mean_r(B[rng.permutation(len(B))]) for _ in range(NP)])
    return obs, float(perms.mean()), float((np.sum(perms >= obs) + 1) / (NP + 1))


def site_analysis(site, rehearsal):
    b = base(site, rehearsal)
    out_json = b / "analysis.json"
    if not rehearsal and out_json.exists():
        raise SystemExit("analysis.json zaten var: analiz bir kez yapılır.")
    rows, idx, sc = load(site, rehearsal)
    fz_thr = json.loads((C.P4 / "m1" / "final_meta.json").read_text())["threshold_spec95"]
    rng = np.random.default_rng(C.SEED + 5)
    R = {"site": site, "rehearsal": rehearsal, "date": str(date.today())}
    d = rows[rows.analysis == "primary_k1"].copy()

    # ---- H1 ve ikincil
    R["H1"] = auc_block(d, seed=1)
    dd, dci, _ = S.auc_diff_boot(d.y, d.logit_M1, d.logit_demo, d.set_id, NB, 2)
    R["H1"]["minus_demo"], R["H1"]["minus_demo_ci"] = dd, dci.tolist()
    R["H1"]["met"] = bool(R["H1"]["ci"][0] > 0.60 and dci[0] > 0)
    R["demo"] = auc_block(d, "logit_demo", 3)
    yb, pr = d.y.to_numpy(bool), d.p_M1.to_numpy()
    sens, spec = float((pr[yb] > fz_thr).mean()), float((pr[~yb] <= fz_thr).mean())
    R["threshold"] = {"thr": fz_thr, "sens": sens, "spec": spec,
                      "ppv": {str(pv): (sens * pv / den if (den := sens * pv + (1 - spec) * (1 - pv)) > 0 else float("nan"))
                              for pv in (0.001, 0.01, 0.10)}}
    R["auprc"] = float(average_precision_score(d.y, d.p_M1))
    R["calibration"] = dict(zip(("intercept", "slope", "ece"), C.calibration(d.y.to_numpy(), pr)))
    m = d.dropna(subset=["logit_M0"])
    dm, dmci, _ = S.auc_diff_boot(m.y, m.logit_M1, m.logit_M0, m.set_id, NB, 4)
    R["M1_vs_M0"] = {"n": len(m), "diff": dm, "ci": dmci.tolist(), "delong_p": C.delong(m.y, m.logit_M1, m.logit_M0)[1]}
    R["M0"] = auc_block(d, "logit_M0", 5)
    R["M2"] = auc_block(d, "logit_M2", 6)

    # ---- H2 (5.1–5.4): ham tek yönlü p'ler
    h2 = {}
    for tag, tgt, rnd in (("5.1_onset", "occ_onset", "occ_rand_onset"), ("5.1_stt", "occ_stt", "occ_rand_stt")):
        y, st, sr = d.y.to_numpy(), d[tgt].to_numpy(), d[rnd].to_numpy()
        stat = lambda ix, y=y, st=st, sr=sr: S.auroc(y[ix], sr[ix]) - S.auroc(y[ix], st[ix]) if 0 < y[ix].sum() < len(ix) else np.nan
        ci, vals = S.cluster_boot(d.set_id, stat, NB, 7)
        h2[tag] = {"estimate": float(S.auroc(y, sr) - S.auroc(y, st)), "ci": ci.tolist(), "p": p_one_sided(vals),
                   "auroc_target_masked": float(S.auroc(y, st)), "auroc_random_masked": float(S.auroc(y, sr))}
    for tag, col in (("5.2_onset", "share_onset"), ("5.2_stt", "share_stt")):
        g = d.groupby("set_id")
        diff = (g.apply(lambda x: x.loc[x.y == 1, col].mean() - x.loc[x.y == 0, col].mean(), include_groups=False)).dropna()
        ci, vals = S.cluster_boot(diff.index.to_numpy(), lambda ix, v=diff.to_numpy(): v[ix].mean(), NB, 8)
        h2[tag] = {"estimate": float(diff.mean()), "ci": ci.tolist(), "p": p_one_sided(vals),
                   "share_case": float(d.loc[d.y == 1, col].mean()), "share_ctrl": float(d.loc[d.y == 0, col].mean())}
    cases = d[d.y == 1]
    if not (rehearsal and site in TEST):                 # 5.3–5.4 gerçek vaka gerektirir
        t2 = pd.read_parquet(C.P5 / site / "t2.parquet")
        mb_all = np.load(C.P5 / site / "index_mb.npy")
        pos = {k: i for i, k in enumerate(idx.k)}
        prof = t2.groupby("pid")[[f"slope20_{x}" for x in C.L8]].apply(
            lambda g: pd.Series(np.nanmean(g.to_numpy() / np.linalg.norm(g.to_numpy(), axis=1, keepdims=True), 0), index=g.columns))
        qint = t2.groupby("pid")[[f"qint_{x}" for x in C.L8]].mean()
        cc = cases[cases.pid.isin(prof.index)]
        A = np.array([C.initial_profile(mb_all[pos[k]], o) for k, o in zip(cc.k, cc.on_low_ms)])
        obs, null, p = perm_agreement(A, prof.loc[cc.pid].to_numpy(), rng)
        h2["5.3_delta"] = {"n": len(cc), "mean_r": obs, "perm_mean_r": null, "p": p}
        T = np.array([C.integral_profile(mb_all[pos[k]], f + 40, f + 400) for k, f in zip(cc.k, cc.off_ms)])
        obs, null, p = perm_agreement(T, qint.loc[cc.pid].to_numpy(), rng)
        h2["5.4_tmemory"] = {"n": len(cc), "mean_r": obs, "perm_mean_r": null, "p": p}
        strata = pd.cut(cases.days_to_t2, [-1, 30, 365, 1e9], labels=["0–30", "31–365", ">365"])
        R["5.4_time_strata"] = {str(s): {"n": int((strata == s).sum()), "median_p_M1": float(cases.p_M1[strata == s].median())}
                                for s in strata.cat.categories}
    R["H2"] = h2
    gib = lambda x: float((x.pqlow_fq_corr <= -0.5).mean())
    R["5.5_gibbs"] = {"case": gib(d[d.y == 1]), "ctrl": gib(d[d.y == 0])}

    # ---- H4 (MIMIC)
    if site == "MIMIC":
        k = rows[rows.analysis == "k2a_vs_k1"]
        g = k.groupby("set_id")
        diff = g.apply(lambda x: x.loc[x.y == 1, "p_M1"].mean() - x.loc[x.y == 0, "p_M1"].mean(), include_groups=False).dropna()
        ci, _ = S.cluster_boot(diff.index.to_numpy(), lambda ix, v=diff.to_numpy(): v[ix].mean(), NB, 9)
        R["H4"] = {"n_sets": len(diff), "diff": float(diff.mean()), "ci": ci.tolist(), "supported": bool(ci[0] <= 0)}

    # ---- ayrık dalga (koşullu lojistik)
    beta, se, n = C.clogit(d.set_id, d.y, d.wave)
    R["wave"] = {"n_sets": n, "or": float(np.exp(beta)), "ci": [float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se))],
                 "prev_case": float(d.wave[d.y == 1].mean()), "prev_ctrl": float(d.wave[d.y == 0].mean())}

    # ---- denetimler
    s = d[d.strict_t0.fillna(False)]
    s = s[s.groupby("set_id").y.transform("sum").eq(1) & s.groupby("set_id").y.transform("count").gt(1)]
    R["7.1_strict"] = auc_block(s, seed=10)
    R["7.1_strict"]["pass"] = bool(R["H1"]["ci"][0] <= R["7.1_strict"]["auroc"] <= R["H1"]["ci"][1] or R["7.1_strict"]["auroc"] >= 0.60)
    if site == "MIMIC":
        grp = {}
        for fg, g in d.groupby("filt"):
            if g.y.sum() >= 10 and (1 - g.y).sum() >= 10:
                grp[fg] = auc_block(g, seed=11)
        aucs = np.array([v["auroc"] for v in grp.values()])
        ses = np.array([(v["ci"][1] - v["ci"][0]) / 3.92 for v in grp.values()])
        if len(aucs) >= 2:
            w = 1 / ses ** 2
            q = float((w * (aucs - (w * aucs).sum() / w.sum()) ** 2).sum())
            from scipy.stats import chi2
            R["7.2_filter"] = {"groups": grp, "Q": q, "p": float(chi2.sf(q, len(aucs) - 1))}
        else:
            R["7.2_filter"] = {"groups": grp, "note": "heterojenlik için yeterli grup yok"}
        R["7.2_filter"]["pass"] = bool(R["7.2_filter"].get("p", 1) > 0.05)
        h3p = base(site, rehearsal) / "h3_scores.parquet"
        if h3p.exists():
            h3 = pd.read_parquet(h3p)
            win = pd.cut(h3.days_from_abl, [-1e9, -1, 90, 1e9], labels=["öncesi", "0–90 gün", ">90 gün"])
            R["7.5_ablation"] = {str(w): {"n_ecg": int((win == w).sum()), "n_pt": int(h3.pid[win == w].nunique()),
                                          "median_p_M1": float(h3.p_M1[win == w].median())} for w in win.cat.categories}
    if site == "I0006" and (C.P5 / site / "alt.parquet").exists():          # §7.3, betimsel
        alt = pd.read_parquet(C.P5 / site / "alt.parquet")
        narrow = np.where(alt.qrs_ms <= alt.c2_qrs_ms, alt.pqlow_mag, alt.c2_pqlow_mag).astype(float)
        k1_pq = idx.loc[idx.k.isin(set(d.k[d.y == 0])), "pqlow_mag"]
        R["7.3_altwpw"] = {"n_ecg": len(alt), "median_pqlow_mag_narrow_cluster": float(np.median(narrow)) if len(alt) else None,
                           "median_pqlow_mag_k1": float(k1_pq.median()), "note": "betimsel; iki kümeli 12SL 303 EKG'leri"}
    ns = d[~d.structural.fillna(False)]
    ns = ns[ns.groupby("set_id").y.transform("sum").eq(1) & ns.groupby("set_id").y.transform("count").gt(1)]
    R["7.7_structural"] = auc_block(ns, seed=12)
    R["7.7_structural"]["pass"] = bool(R["H1"]["ci"][0] <= R["7.7_structural"]["auroc"] <= R["H1"]["ci"][1])

    # ---- duyarlılık
    sens_ = {}
    for an in sorted(set(rows.analysis) - {"primary_k1", "primary_k2", "k2a_vs_k1"}):
        sens_[an] = auc_block(rows[rows.analysis == an], seed=13)
    young = d[d.set_id.isin(d.loc[(d.y == 1) & (d.age < 40), "set_id"])]
    sens_["age_lt40"] = auc_block(young, seed=14) if young.y.sum() >= 10 else {"n_case": int(young.y.sum()), "note": "az vaka"}
    R["K2_primary"] = auc_block(rows[rows.analysis == "primary_k2"], seed=15)
    from sklearn.linear_model import LogisticRegression
    cov = ["age", "hr", "male", "n_ecg", "shifted_year"] + (["fyear"] if site != "MIMIC" else [])
    a = d.dropna(subset=cov + ["logit_M1"]).copy()
    a["z"] = (a.logit_M1 - a.logit_M1.mean()) / a.logit_M1.std()
    a["log_n"] = np.log(a.n_ecg)
    X = np.c_[a.z, (a[[c for c in cov if c != "n_ecg"] + ["log_n"]] - a[[c for c in cov if c != "n_ecg"] + ["log_n"]].mean()) /
              a[[c for c in cov if c != "n_ecg"] + ["log_n"]].std()]
    coef = LogisticRegression(C=1e9, max_iter=2000).fit(X, a.y).coef_[0, 0]
    _, vals = S.cluster_boot(a.set_id, lambda ix: LogisticRegression(C=1e9, max_iter=2000).fit(X[ix], a.y.to_numpy()[ix]).coef_[0, 0]
                             if 0 < a.y.to_numpy()[ix].sum() < len(ix) else np.nan, 500, 16)
    sens_["adjusted_or_per_sd"] = {"or": float(np.exp(coef)), "ci": np.exp(np.percentile(vals, [2.5, 97.5])).tolist(),
                                   "covariates": cov, "n_boot": 500}
    if "logit_M1raw" in d.columns:                               # §8.7
        sens_["m1_raw10s"] = auc_block(d, "logit_M1raw", 18)
    ap_ = base(site, rehearsal) / "agg_scores.parquet"
    if ap_.exists():                                             # §8.4
        ag = pd.read_parquet(ap_)
        pt = ag.groupby(["set_id", "pid", "grp"]).logit_M1.agg(["mean", "max", "size"]).reset_index()
        if rehearsal and site in TEST:                           # prova: kümede rastgele vaka
            rr = np.random.default_rng(C.SEED + 98)
            pt["grp"] = pt.groupby("set_id").pid.transform(lambda x: np.where(np.arange(len(x)) == rr.integers(len(x)), "case", "k1"))
        pt["y"] = (pt.grp == "case").astype(int)
        pt = pt[pt.groupby("set_id").y.transform("sum").eq(1) & pt.groupby("set_id").y.transform("count").gt(1)]
        for stat in ("mean", "max"):
            a_, ci_ = S.auc_boot(pt.y, pt[stat], pt.set_id, NB, 17)
            sens_[f"agg_{stat}"] = {"n_case": int(pt.y.sum()), "n_ctrl": int((1 - pt.y).sum()), "auroc": a_, "ci": ci_.tolist(),
                                    "median_ecg_per_patient": float(pt["size"].median())}
    R["sensitivity"] = sens_
    out_json.write_text(json.dumps(R, indent=1, default=float))
    print(f"{site}{' (prova)' if rehearsal else ''}: analiz yazıldı → {out_json}")


def combine(rehearsal):
    sites = [s for s in ("MIMIC", "I0006") if (base(s, rehearsal) / "analysis.json").exists()]
    R = {s: json.loads((base(s, rehearsal) / "analysis.json").read_text()) for s in sites}
    tests = [(s, t, R[s]["H2"][t]["p"]) for s in sites for t in ("5.1_onset", "5.1_stt", "5.2_onset", "5.2_stt",
                                                                 "5.3_delta", "5.4_tmemory") if t in R[s]["H2"]]
    adj = C.holm([t[2] for t in tests]) if tests else []
    f = lambda x: f"{x:.3f}"
    md = [f"# HIDDEN-DELTA — kayıtlı test sonuçları{' (PROVA — anlamsız)' if rehearsal else ''}\n",
          "ANALYSIS_PLAN v1.2; OSF https://osf.io/kvn28. Yalnızca agregat.\n"]
    if "MIMIC" in R:
        h = R["MIMIC"]["H1"]
        md += ["## H1 (birincil, MIMIC)\n",
               f"AUROC {f(h['auroc'])} (%95 GA {f(h['ci'][0])}–{f(h['ci'][1])}); {h['n_case']} vaka, {h['n_ctrl']} kontrol. "
               f"M1 − demografik model {h['minus_demo']:+.3f} ({h['minus_demo_ci'][0]:+.3f} – {h['minus_demo_ci'][1]:+.3f}). "
               f"**H1 {'karşılandı' if h['met'] else 'karşılanmadı'}.**\n"]
    md += ["## H2 ailesi (Holm, tek yönlü α 0.05)\n", "| Site | Test | Tahmin / ortalama r | Ham p | Holm p |", "|---|---|---|---|---|"]
    for (s, t, p), pa in zip(tests, adj):
        e = R[s]["H2"][t]
        est = e.get("estimate", e.get("mean_r"))
        md.append(f"| {s} | {t} | {est:+.4f} | {p:.4f} | {pa:.4f} |")
    wave_ok = all(R[s]["wave"]["ci"][0] > 1 for s in sites) and len(sites) == 2
    md += ["", "## Ayrık dalga", *(f"- {s}: OR {f(R[s]['wave']['or'])} ({f(R[s]['wave']['ci'][0])}–{f(R[s]['wave']['ci'][1])}); "
                                    f"sıklık vaka {R[s]['wave']['prev_case']:.3f}, kontrol {R[s]['wave']['prev_ctrl']:.3f}" for s in sites),
           f"- Doğrulama (iki sitede de OR > 1, GA 1'i dışlıyor): {'evet' if wave_ok else 'hayır'}\n"]
    if "MIMIC" in R and "H4" in R["MIMIC"]:
        h4 = R["MIMIC"]["H4"]
        md.append(f"## H4 (MIMIC, AVNRT − K1, M1 olasılığı)\n\nFark {h4['diff']:+.4f} ({h4['ci'][0]:+.4f} – {h4['ci'][1]:+.4f}); "
                  f"{'destekleniyor' if h4['supported'] else 'desteklenmiyor'}.\n")
    md.append("Ayrıntılar: site başına `analysis.json` (ikincil sonlanımlar, denetimler, duyarlılık analizleri).")
    outd = (C.P5 / "rehearsal") if rehearsal else C.OUT
    outd.mkdir(parents=True, exist_ok=True)
    name = "results_rehearsal.md" if rehearsal else "results.md"
    if not rehearsal and (outd / name).exists():
        raise SystemExit("results.md zaten var.")
    (outd / name).write_text("\n".join(md) + "\n")
    (outd / name.replace(".md", "_holm.json")).write_text(json.dumps([{"site": t[0], "test": t[1], "p": t[2], "holm": float(a)}
                                                                      for t, a in zip(tests, adj)], indent=1))
    print("\n".join(md))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["site", "combine"])
    ap.add_argument("site", nargs="?", choices=["MIMIC", "I0006", "I0001"])
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    if a.mode == "site":
        if not a.rehearsal and a.site not in TEST:
            raise SystemExit("Resmî analiz yalnızca MIMIC ve I0006 için.")
        site_analysis(a.site, a.rehearsal)
    else:
        combine(a.rehearsal)


if __name__ == "__main__":
    main()
