"""HIDDEN-DELTA G2 pilot polarite denetimi (concept.md §10 "G2 uygulaması"; ölçütler sonuçlardan önce yazıldı).

Birincil (karar): MGH, rastgele 30 vaka (birincil uzlaşı tanımı; ≥1 başarılı 500 Hz T0 ve T2) + 150 K1.
İkincil (betimsel): bütün uygun vakalar + bütün K1 havuzu. Hasta başına bir rastgele T0 EKG.

Lob: pqhigh_mag ≥ 0.020 mV (görünür başlangıçtan önce −20…−2 ms, 8 bağımsız derivasyon normu).
Karakter: pqhigh_fq_corr ≤ −0.5 Gibbs benzeri, ≥ +0.5 delta benzeri.
Delta vektörü uyumu: corr(T0 lob profili, hastanın T2 ilk 20 ms eğim profili), 8 derivasyon; tek yönlü permütasyon.
Çıktı yalnızca agregat: stage2/out/g2_pilot.{md,json}. Örneklem kimlikleri PRIVATE altında.
"""
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

P1 = Path("/home/fatih/data/hidden-delta/stage1")
P2 = Path("/home/fatih/data/hidden-delta/stage2/I0001")
P0 = Path("/home/fatih/data/hidden-delta/stage0")
OUT = Path(__file__).resolve().parent / "out"
SEED = 20261114
N_CASE, N_CTRL = 30, 150
LOBE_THR, GIBBS, DELTA = 0.020, -0.5, 0.5
L8 = ["I", "II", "V1", "V2", "V3", "V4", "V5", "V6"]
N_PERM, N_BOOT = 10000, 2000


def load():
    con = duckdb.connect()
    feat = con.execute(f"""SELECT * FROM read_parquet('{P2}/features_*.parquet', union_by_name=true)
                           WHERE status = 'ok'""").df()
    lab = con.execute(f"SELECT k, t2_cons, t0_lab FROM '{P1 / 'heedb_I0001_ecg_labels.parquet'}'").df()
    coh = con.execute(f"SELECT pid, case_cons FROM '{P1 / 'heedb_I0001_cohort.parquet'}'").df()
    age = con.execute(f"""SELECT k, age FROM '{P0 / 'heedb_I0001_ecg_flags.parquet'}'
                          WHERE k IN (SELECT k FROM read_parquet('{P2}/features_*.parquet'))""").df()
    feat = feat.merge(lab, on="k", how="left").merge(age, on="k", how="left")
    cases = set(coh.loc[coh.case_cons, "pid"])
    t0 = feat[(feat.role == "case_t0") & feat.t0_lab.fillna(False) & feat.pid.isin(cases)]
    t2 = feat[feat.t2_cons.fillna(False) & feat.pid.isin(cases)]
    ctrl = feat[feat.role == "k1"]
    eligible = sorted(set(t0.pid) & set(t2.pid))
    return t0, t2, ctrl, eligible


def unit(M):
    n = np.linalg.norm(M, axis=1, keepdims=True)
    return np.divide(M, n, out=np.zeros_like(M), where=n > 0)


def delta_profiles(t2):
    s = t2[["pid"] + [f"slope20_{lead}" for lead in L8]].copy()
    s[[f"slope20_{lead}" for lead in L8]] = unit(s[[f"slope20_{lead}" for lead in L8]].to_numpy())
    return s.groupby("pid").mean()


def corr_rows(A, B):
    A = A - A.mean(axis=1, keepdims=True)
    B = B - B.mean(axis=1, keepdims=True)
    den = np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1)
    return np.divide((A * B).sum(axis=1), den, out=np.full(len(A), np.nan), where=den > 0)


def boot_ci(x_case, x_ctrl, rng, stat=np.mean):
    d = [stat(rng.choice(x_case, len(x_case))) - stat(rng.choice(x_ctrl, len(x_ctrl))) for _ in range(N_BOOT)]
    return np.percentile(d, [2.5, 97.5])


def analyse(case_ecg, ctrl_ecg, dprof, rng):
    R = {"n_vaka": len(case_ecg), "n_kontrol": len(ctrl_ecg)}
    lc, lk = (case_ecg.pqhigh_mag >= LOBE_THR).to_numpy(), (ctrl_ecg.pqhigh_mag >= LOBE_THR).to_numpy()
    R["lob_vaka"], R["lob_kontrol"] = float(lc.mean()), float(lk.mean())
    R["lob_fark"], R["lob_fark_ci"] = float(lc.mean() - lk.mean()), boot_ci(lc, lk, rng).round(3).tolist()
    for tag, d, lob in (("vaka", case_ecg, lc), ("kontrol", ctrl_ecg, lk)):
        fc = d.pqhigh_fq_corr.to_numpy()[lob]
        R[f"lob_n_{tag}"] = int(lob.sum())
        R[f"gibbs_{tag}"] = float(np.mean(fc <= GIBBS)) if lob.any() else np.nan
        R[f"delta_{tag}"] = float(np.mean(fc >= DELTA)) if lob.any() else np.nan
        R[f"kararsiz_{tag}"] = float(np.mean((fc > GIBBS) & (fc < DELTA))) if lob.any() else np.nan
        R[f"pq_mag_medyan_{tag}"] = float(d.pqhigh_mag.median())
        R[f"slow_init_medyan_{tag}"] = float(d.slow_init_ms.median())
        R[f"septal_q_yok_{tag}"] = float(((d[["q_I", "q_aVL", "q_V5", "q_V6"]] == 0).all(axis=1)).mean())
        R[f"q_V6_var_{tag}"] = float((d.q_V6 > 0).mean())
        R[f"yas_medyan_{tag}"] = float(d.age.median())
        R[f"hr_medyan_{tag}"] = float(d.hr.median())
    R["slow_init_fark_ci"] = boot_ci(case_ecg.slow_init_ms.to_numpy(), ctrl_ecg.slow_init_ms.to_numpy(), rng,
                                     np.median).round(1).tolist()
    # delta vektörü uyumu
    pq = case_ecg[[f"pqhigh_{lead}" for lead in L8]].to_numpy()
    D = dprof.loc[case_ecg.pid, [f"slope20_{lead}" for lead in L8]].to_numpy()

    def stats(sel):
        own = corr_rows(pq[sel], D[sel])
        m_own = float(np.nanmean(own))
        idx = np.where(sel)[0]
        perms = []
        for _ in range(N_PERM):
            p = rng.permutation(len(D))
            perms.append(np.nanmean(corr_rows(pq[idx], D[p[: len(idx)]])))
        perms = np.array(perms)
        return m_own, float(perms.mean()), float((np.sum(perms >= m_own) + 1) / (N_PERM + 1))
    R["uyum_tum"] = stats(np.ones(len(pq), bool))
    R["uyum_loblu"] = stats(lc) if lc.sum() >= 5 else (np.nan, np.nan, np.nan)
    ctrl_pq = ctrl_ecg[[f"pqhigh_{lead}" for lead in L8]].to_numpy()
    Dr = D[rng.integers(0, len(D), len(ctrl_pq))]
    R["kontrol_rastgele_delta_corr"] = float(np.nanmean(corr_rows(ctrl_pq, Dr)))
    return R


def decide(R):
    cg, cd = R["gibbs_vaka"], R["delta_vaka"]
    m_own, m_null, p = R["uyum_tum"]
    artefakt = (cg >= 0.5 and cd < 0.25) or (R["lob_fark"] <= 0 and p >= 0.2)
    sinyal = (cd >= cg) and (p < 0.05)
    if artefakt:
        return "ARTEFAKT BASKIN"
    if sinyal:
        return "SİNYAL LEHİNE"
    return "KARARSIZ"


def table(R):
    f = lambda x: "—" if x is None or (isinstance(x, float) and np.isnan(x)) else (f"{x:.3f}" if isinstance(x, float) else str(x))
    rows = [("Hasta", R["n_vaka"], R["n_kontrol"]),
            ("Medyan yaş / KH", f"{f(R['yas_medyan_vaka'])} / {f(R['hr_medyan_vaka'])}",
             f"{f(R['yas_medyan_kontrol'])} / {f(R['hr_medyan_kontrol'])}"),
            (f"Lob (pqhigh_mag ≥ {LOBE_THR} mV)", f(R["lob_vaka"]), f(R["lob_kontrol"])),
            ("  loblu EKG", R["lob_n_vaka"], R["lob_n_kontrol"]),
            ("  Gibbs benzeri (fq_corr ≤ −0.5)", f(R["gibbs_vaka"]), f(R["gibbs_kontrol"])),
            ("  delta benzeri (fq_corr ≥ +0.5)", f(R["delta_vaka"]), f(R["delta_kontrol"])),
            ("  kararsız", f(R["kararsiz_vaka"]), f(R["kararsiz_kontrol"])),
            ("pqhigh_mag medyanı (mV)", f(R["pq_mag_medyan_vaka"]), f(R["pq_mag_medyan_kontrol"])),
            ("Yavaş başlangıç medyanı (ms)", f(R["slow_init_medyan_vaka"]), f(R["slow_init_medyan_kontrol"])),
            ("Septal q yok (I, aVL, V5, V6'nın hiçbirinde)", f(R["septal_q_yok_vaka"]), f(R["septal_q_yok_kontrol"])),
            ("V6'da q var", f(R["q_V6_var_vaka"]), f(R["q_V6_var_kontrol"]))]
    md = ["| | Vaka (T0) | K1 |", "|---|---|---|"] + [f"| {a} | {b} | {c} |" for a, b, c in rows]
    u, ul = R["uyum_tum"], R["uyum_loblu"]
    md += ["", f"- Lob sıklığı farkı (vaka − kontrol): {R['lob_fark']:.3f}, %95 GA {R['lob_fark_ci']}",
           f"- Yavaş başlangıç medyan farkı %95 GA: {R['slow_init_fark_ci']} ms",
           f"- Kendi delta vektörüyle uyum (tüm vakalar): ortalama r = {u[0]:.3f}; permütasyon ortalaması "
           f"{u[1]:.3f}; tek yönlü p = {u[2]:.4f}",
           f"- Aynısı yalnız loblu vakalarda: r = {f(ul[0])}, permütasyon {f(ul[1])}, p = {f(ul[2])}",
           f"- Kontrollerde rastgele bir vakanın delta vektörüyle ortalama r: {R['kontrol_rastgele_delta_corr']:.3f}"]
    return md


def main():
    rng = np.random.default_rng(SEED)
    t0, t2, ctrl, eligible = load()
    dprof = delta_profiles(t2)
    one = lambda d: d.sample(frac=1.0, random_state=int(rng.integers(2**31))).drop_duplicates("pid")
    # birincil pilot
    pilot_case_ids = rng.choice(eligible, N_CASE, replace=False)
    pilot_ctrl_ids = rng.choice(sorted(ctrl.pid.unique()), N_CTRL, replace=False)
    pc = one(t0[t0.pid.isin(pilot_case_ids)])
    pk = one(ctrl[ctrl.pid.isin(pilot_ctrl_ids)])
    Rp = analyse(pc, pk, dprof, rng)
    Rp["karar"] = decide(Rp)
    # ikincil: tüm geliştirme seti
    fc = one(t0[t0.pid.isin(eligible)])
    fk = one(ctrl)
    Rf = analyse(fc, fk, dprof, rng)
    pd.DataFrame({"pid": np.r_[pilot_case_ids, pilot_ctrl_ids],
                  "grup": ["vaka"] * N_CASE + ["kontrol"] * N_CTRL}).to_parquet(P2 / "g2_pilot_sample.parquet")
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# G2 pilot polarite denetimi (MGH)\n",
          "Ölçütler sonuçlardan önce concept.md §10'a yazıldı. Yalnızca agregat.\n",
          f"## Birincil pilot (30 vaka + 150 K1) — karar: **{Rp['karar']}**\n"] + table(Rp) + [
          "\n## İkincil, betimsel: bütün MGH geliştirme seti\n",
          f"Uygun vaka hastası: {len(eligible):,}; K1: {ctrl.pid.nunique():,}. Karar kuralına girmez.\n"] + table(Rf) + [
          "\n**Not.** Vakalar ve kontroller eşleştirilmedi (yaş, KH, cihaz). Farklar betimseldir; eşleştirme Aşama 3'te.\n"]
    (OUT / "g2_pilot.md").write_text("\n".join(md))
    (OUT / "g2_pilot.json").write_text(json.dumps({"pilot": Rp, "tum": Rf}, indent=1, default=float))
    print("\n".join(md))


if __name__ == "__main__":
    main()
