"""HIDDEN-DELTA Aşama 5 — kayıt sonrası betimsel ekler (ANALYSIS_PLAN v1.2; OSF https://osf.io/kvn28).

Bu script resmî sonuçlar görüldükten sonra yazıldı (2026-11-08). Planın betimsel olarak istediği ama `s5_analysis.py`'de
eksik kalan çıktıları üretir. Test yok, karar kuralı yok; kayıt sonrası olarak raporlanır. Yeni skorlama yapılmaz:
dondurulmuş skorlar (`scores.parquet`) ve dondurulmuş indeks EKG'leri (`index_mb.npy`) kullanılır.

- 5.2: P dalgasına düşen atıf payı, vaka − kontrol (küme başına fark, küme bootstrap'ı).
- §6 ek kontroller (birincil K1 kümeleri):
  - dalganın polaritesi, aynı derivasyondaki ilk QRS sapmasına göre (0–15 ms ortalaması eksi başlangıçtaki değer;
    5.5 referansı (40–22 ms önce) dalga tepesini içerdiği için kullanılmadı);
  - dalganın 8 derivasyon profili ile hastanın kendi T2 ilk 20 ms eğim profilinin uyumu (5.3 ile aynı permütasyon);
  - zamanlama, genlik, genişlik;
  - MIMIC'te filtre grubuna göre sıklık.
- 7.2 ikinci koşul: "<not specified>" grubunun AUROC'u eksi ana grubun AUROC'u (grup içinde küme bootstrap'ı).
- 7.4 K2b (gizli yol, MIMIC): uygun hasta sayısı ve eşleşmiş kümelerdeki sayı.

Çalıştırma: python3 stage5/s5_supplement.py
"""
import json

import numpy as np
import pandas as pd

import common5 as C
import s5_analysis as A
import statlib as S

SEED = C.SEED + 500
LI = C.hdsig.LI
W60 = 30                                                # 60 ms, 500 Hz


def cluster_diff(d, col, seed):
    g = d.groupby("set_id")
    diff = g.apply(lambda x: x.loc[x.y == 1, col].mean() - x.loc[x.y == 0, col].mean(), include_groups=False).dropna()
    ci, _ = S.cluster_boot(diff.index.to_numpy(), lambda ix, v=diff.to_numpy(): v[ix].mean(), A.NB, seed)
    return {"estimate": float(diff.mean()), "ci": ci.tolist(),
            "share_case": float(d.loc[d.y == 1, col].mean()), "share_ctrl": float(d.loc[d.y == 0, col].mean())}


def wave_rows(d, idx, mb_all):
    """Her EKG için ayrık dalga ayrıntısı: polarite uyumu, 8 derivasyon profili, zamanlama, genlik, genişlik."""
    pos = {k: i for i, k in enumerate(idx.k)}
    out = []
    for k, o in zip(d.k, d.on_low_ms):
        mb = mb_all[pos[k]].astype(float)
        on = C.ms2s(o)
        ok, det = C.hdsig.discrete_wave(mb, C.FS, on)
        if not ok:
            out.append({"k": k, "wave_d": False})
            continue
        conc = []
        for lead, sgn in det["polarity"].items():
            x = mb[LI[lead]]
            fq = x[on:on + 8].mean() - x[on]               # QRS başlangıcındaki düzeye göre (dalga penceresinden bağımsız)
            conc.append(float(sgn == np.sign(fq)))
        W = mb[[LI[x] for x in C.L8], on - W60:on]
        Dv = W - np.median(W, axis=1, keepdims=True)
        pk = min(W60 - 1, max(0, W60 + int(round(det["peak_ms"] * C.FS / 1000))))
        out.append({"k": k, "wave_d": True, "concord": float(np.mean(conc)), "amp_uv": det["amp"] * 1000,
                    "fwhm_ms": det["fwhm"], "peak_ms": det["peak_ms"], **{f"wp_{x}": Dv[j, pk] for j, x in enumerate(C.L8)}})
    return d[["k", "pid", "set_id", "y", "wave", "filt"]].merge(pd.DataFrame(out), on="k", how="left")


def by_status(w, col, f=np.median):
    return {"case": float(f(w.loc[w.y == 1, col])), "ctrl": float(f(w.loc[w.y == 0, col]))}


def auc_group_diff(d, ga, gb, seed):
    """AUROC(gb) − AUROC(ga); her grupta kümeler ayrı yeniden örneklenir."""
    rng = np.random.default_rng(seed)
    parts = []
    for g in (ga, gb):
        x = d[d.filt == g].reset_index(drop=True)
        parts.append((x.y.to_numpy(), x.logit_M1.to_numpy(), list(x.groupby("set_id").indices.values())))
    est = S.auroc(parts[1][0], parts[1][1]) - S.auroc(parts[0][0], parts[0][1])
    vals = []
    for _ in range(A.NB):
        a = []
        for y, s, ind in parts:
            ix = np.concatenate([ind[j] for j in rng.integers(0, len(ind), len(ind))])
            a.append(S.auroc(y[ix], s[ix]) if 0 < y[ix].sum() < len(ix) else np.nan)
        vals.append(a[1] - a[0])
    return {"group_a": ga, "group_b": gb, "diff_b_minus_a": float(est),
            "ci": np.nanpercentile(vals, [2.5, 97.5]).tolist()}


def site(site):
    rows, idx, _ = A.load(site, False)
    d = rows[rows.analysis == "primary_k1"].copy()
    mb_all = np.load(C.P5 / site / "index_mb.npy")
    R = {"site": site, "note": "kayıt sonrası betimsel; resmî sonuçlardan sonra yazıldı"}
    R["5.2_p_share"] = cluster_diff(d, "share_p", SEED + 1)

    w = wave_rows(d, idx, mb_all)
    R["wave_recompute_agrees"] = bool((w.wave_d.astype(float) == w.wave).all())
    wv = w[w.wave_d.astype(bool)]
    R["wave_detail"] = {"n_case": int(wv.y.sum()), "n_ctrl": int((1 - wv.y).sum()),
                        "concord_mean": by_status(wv, "concord", np.mean),
                        "majority_opposite_share": by_status(wv.assign(opp=(wv.concord < 0.5).astype(float)), "opp", np.mean),
                        "amp_uv_median": by_status(wv, "amp_uv"), "fwhm_ms_median": by_status(wv, "fwhm_ms"),
                        "peak_ms_median": by_status(wv, "peak_ms")}

    t2 = pd.read_parquet(C.P5 / site / "t2.parquet")
    prof = t2.groupby("pid")[[f"slope20_{x}" for x in C.L8]].apply(
        lambda g: pd.Series(np.nanmean(g.to_numpy() / np.linalg.norm(g.to_numpy(), axis=1, keepdims=True), 0), index=g.columns))
    cw = wv[(wv.y == 1) & wv.pid.isin(prof.index)]
    obs, null, p = A.perm_agreement(cw[[f"wp_{x}" for x in C.L8]].to_numpy(float), prof.loc[cw.pid].to_numpy(),
                                    np.random.default_rng(SEED + 2))
    R["wave_own_delta"] = {"n": len(cw), "mean_r": obs, "perm_mean_r": null, "p_one_sided": p}

    if site == "MIMIC":
        R["wave_by_filter"] = {str(g): {"n_case": int(x.y.sum()), "n_ctrl": int((1 - x.y).sum()),
                                        "prev_all": float(x.wave.mean()), "prev_case": float(x.wave[x.y == 1].mean()),
                                        "prev_ctrl": float(x.wave[x.y == 0].mean())} for g, x in w.groupby("filt")}
        groups = [g for g, x in d.groupby("filt") if x.y.sum() >= 10 and (1 - x.y).sum() >= 10]
        ns = [g for g in groups if g.startswith("<not specified>")]
        main = [g for g in groups if g not in ns]
        if ns and main:
            R["7.2_second_condition"] = auc_group_diff(d, main[0], ns[0], SEED + 3)

        k2 = pd.read_parquet(C.P0 / "mimic_k2.parquet")
        b = lambda c: k2[c].astype("boolean").fillna(False).astype(bool)
        k2b = b("ap") & ~b("avnrt") & ~b("wpw")
        elig = k2b & (k2.n_t0 > 0) & ~b("any_spr") & ~b("any_abl")
        sets = pd.read_parquet(C.P5 / site / "sets.parquet")
        R["7.4_K2b"] = {"n_any": int(k2b.sum()), "n_eligible_t0": int(elig.sum()),
                        "n_in_matched_sets": int(sets[sets.pid.isin(set(k2.sid[k2b]))].pid.nunique()),
                        "note": "eşleşmiş kümede K2b yok; betimsel skor raporlanamaz"}
    (C.P5 / site / "supplement.json").write_text(json.dumps(R, indent=1, ensure_ascii=False))
    return R


def f3(x):
    return f"{x:.3f}"


def report(res):
    md = ["# HIDDEN-DELTA — kayıt sonrası betimsel ekler\n",
          "Resmî sonuçlardan sonra yazıldı (`stage5/s5_supplement.py`). Planın betimsel olarak istediği ama analiz kodunda eksik "
          "kalan çıktılar. Test ve karar kuralı yok. Yeni skorlama yok. Yalnızca agregat.\n"]
    for s, R in res.items():
        n = C.SITE_NAME[s]
        p = R["5.2_p_share"]
        wd, wo = R["wave_detail"], R["wave_own_delta"]
        md += [f"## {n}\n",
               f"- **5.2 P dalgası atıf payı:** vaka {f3(p['share_case'])}, kontrol {f3(p['share_ctrl'])}; küme farkı "
               f"{p['estimate']:+.4f} ({p['ci'][0]:+.4f} – {p['ci'][1]:+.4f}).",
               f"- **Ayrık dalga, ayrıntı** (dalgalı EKG: {wd['n_case']} vaka, {wd['n_ctrl']} kontrol; yeniden hesap resmî bayrakla "
               f"{'aynı' if R['wave_recompute_agrees'] else 'FARKLI'}):",
               f"  - İlk QRS sapmasıyla aynı polaritedeki derivasyon oranı: vaka {f3(wd['concord_mean']['case'])}, kontrol "
               f"{f3(wd['concord_mean']['ctrl'])}. Çoğunlukla ters polarite: vaka {f3(wd['majority_opposite_share']['case'])}, "
               f"kontrol {f3(wd['majority_opposite_share']['ctrl'])}.",
               f"  - Medyan genlik {wd['amp_uv_median']['case']:.1f} / {wd['amp_uv_median']['ctrl']:.1f} µV; genişlik "
               f"{wd['fwhm_ms_median']['case']:.1f} / {wd['fwhm_ms_median']['ctrl']:.1f} ms; tepe "
               f"{wd['peak_ms_median']['case']:.1f} / {wd['peak_ms_median']['ctrl']:.1f} ms (vaka / kontrol; QRS başlangıcına göre).",
               f"  - Kendi delta yönüyle uyum (dalgalı vakalar, n {wo['n']}): ortalama r {f3(wo['mean_r'])}, permütasyon "
               f"{f3(wo['perm_mean_r'])}; tek yönlü p {wo['p_one_sided']:.4f}."]
        if "wave_by_filter" in R:
            md.append("- **Dalga sıklığı, filtre grubuna göre:**")
            md += [f"  - {g}: {v['n_case']} vaka, {v['n_ctrl']} kontrol; genel {f3(v['prev_all'])}, vaka "
                   f"{f3(v['prev_case']) if v['n_case'] else '—'}, kontrol {f3(v['prev_ctrl'])}"
                   for g, v in R["wave_by_filter"].items()]
        if "7.2_second_condition" in R:
            c = R["7.2_second_condition"]
            md.append(f"- **7.2 ikinci koşul:** AUROC farkı ({c['group_b']} − {c['group_a']}) {c['diff_b_minus_a']:+.3f} "
                      f"({c['ci'][0]:+.3f} – {c['ci'][1]:+.3f}).")
        if "7.4_K2b" in R:
            k = R["7.4_K2b"]
            md.append(f"- **7.4 K2b (gizli yol):** notlarda {k['n_any']} hasta, uygun T0'ı olan {k['n_eligible_t0']}, eşleşmiş "
                      f"kümelerde {k['n_in_matched_sets']}. Betimsel skor raporlanamaz.")
        md.append("")
    (C.OUT / "supplement.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    report({s: site(s) for s in ("MIMIC", "I0006")})
