"""HIDDEN-DELTA Aşama 2 — pipeline doğrulaması (etiket kullanmaz).

1. Sentetik testler (gerçek veri yok): normal vuru, delta benzeri yavaş başlangıç (kısa PR, septal q yok),
   ve normal vurunun sıfır fazlı keskin FIR alçak geçirenden geçirilmiş hali (Gibbs pre-ringing).
   Beklenen: delta'da yavaş başlangıç uzun ve QRS öncesi lob ilk sapmayla aynı yönde (pqhigh_fq_corr > 0);
   Gibbs'te lob ilk sapmanın tersi (pqhigh_fq_corr < 0, pqhigh_opp_frac yüksek).
2. MIMIC makine ölçümleriyle uyum: preeksitasyon/kısa PR havuzlarının dışındaki rastgele T0 EKG'lerinde
   QRS süresi ve PR aralığı. Çıktı yalnızca agregat.

Çalıştırma: python3 stage2/s2_validate.py [--n 1000]
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import filtfilt, firwin

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage0"))
import hdsig  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
AMPS = [0.6, 1.0, 0.4, -0.8, 0.2, 0.7, -0.6, -0.2, 0.6, 1.2, 1.4, 1.1]


# ------------------------------------------------------------------ sentetik
def synth(kind, fs=500, hr=70, seed=0):
    rng = np.random.default_rng(seed)
    n = 10 * fs
    t = np.arange(n) / fs
    rr = 60 / hr
    X = np.zeros((12, n))
    starts = np.arange(0.5, 10 - 0.8, rr)
    for k, a in enumerate(AMPS):
        x = np.zeros(n)
        for t0 in starts:
            if kind == "delta":
                knots = [(-0.045, 0), (0.0, 0.18 * a), (0.030, 1.0 * a), (0.055, -0.3 * a), (0.080, 0)]
                p_c = t0 - 0.100
            else:
                knots = [(0.0, 0), (0.010, -0.12 * a), (0.035, 1.0 * a), (0.055, -0.3 * a), (0.075, 0)]
                p_c = t0 - 0.150
            kt = np.array([t0 + d for d, _ in knots])
            kv = np.array([v for _, v in knots])
            m = (t >= kt[0]) & (t <= kt[-1])
            x[m] += np.interp(t[m], kt, kv)
            x += 0.10 * np.sign(a + 0.3) * np.exp(-0.5 * ((t - p_c) / 0.022) ** 2)
            x += 0.25 * np.sign(a + 0.3) * np.exp(-0.5 * ((t - (t0 + 0.30)) / 0.05) ** 2)
        X[k] = x + 0.005 * rng.standard_normal(n)
    if kind == "gibbs":
        h = firwin(151, 30, fs=fs)                       # keskin, doğrusal fazlı; filtfilt ile sıfır fazlı
        X = filtfilt(h, [1.0], X, axis=1)
    return X


def synthetic_tests():
    keys = ["slow_init_ms", "qrs_ms", "pr_ms", "q_V6", "pqhigh_mag", "pqhigh_fq_corr", "pqhigh_opp_frac",
            "pqlow_fq_corr", "n_beats", "n_clusters"]
    rows = []
    for kind in ("normal", "delta", "gibbs"):
        rec, _ = hdsig.process_signal(synth(kind), 500)
        rows.append({"vaka": kind, **{k: rec.get(k) for k in keys}})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ MIMIC uyumu
def _one(path):
    try:
        rec, _ = hdsig.process_record(path)
        return {"ok": True, "pq_base_slope": rec.get("pq_base_slope"),
                **{k: rec.get(k) for k in ("qrs_ms", "qrs_high_ms", "pr_ms", "n_beats", "n_clusters",
                                           "dom_frac", "onset_capped", "slow_init_ms")}}
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "err": type(e).__name__ + ": " + str(e)[:40]}


def mimic_agreement(n, seed=20261108):
    import duckdb
    from common import MIMIC_ECG, PRIVATE
    con = duckdb.connect()
    df = con.execute(f"""
      WITH pool AS (SELECT f.sid, f.study FROM '{PRIVATE / 'mimic_ecg_flags.parquet'}' f
                    WHERE f.t0 AND f.sid NOT IN (SELECT sid FROM '{PRIVATE / 'mimic_ecg_flags.parquet'}' WHERE pre OR spr)),
           mm AS (SELECT study_id::BIGINT study, try_cast(p_onset AS INT) p_on, try_cast(qrs_onset AS INT) q_on,
                         try_cast(qrs_end AS INT) q_end
                  FROM read_csv('{MIMIC_ECG / 'machine_measurements.csv'}', all_varchar=true))
      SELECT pool.study, rl.path, mm.p_on, mm.q_on, mm.q_end FROM pool
      JOIN read_csv('{MIMIC_ECG / 'record_list.csv'}', all_varchar=true) rl ON rl.study_id::BIGINT = pool.study
      JOIN mm USING (study)""").df()
    df = df.sample(n=n, random_state=seed).reset_index(drop=True)
    with ProcessPoolExecutor(24) as ex:
        res = list(ex.map(_one, [str(MIMIC_ECG / p) for p in df.path], chunksize=8))
    r = pd.concat([df, pd.DataFrame(res)], axis=1)
    ok = r[r.ok]
    m_qrs = (ok.q_end - ok.q_on).where((ok.q_end < 29999) & (ok.q_on < 29999))
    m_pr = (ok.q_on - ok.p_on).where((ok.p_on < 29999) & (ok.q_on < 29999))

    def agree(ours, mach, name):
        d = pd.DataFrame({"o": ours, "m": mach}).dropna()
        diff = d.o - d.m
        return {"ölçüm": name, "n": len(d), "fark medyanı (ms)": round(float(diff.median()), 1),
                "mutlak fark medyanı (ms)": round(float(diff.abs().median()), 1),
                "|fark| ≤20 ms (%)": round(100 * float((diff.abs() <= 20).mean()), 1),
                "Spearman r": round(float(d.o.corr(d.m, method="spearman")), 3),
                "bizim medyan": round(float(d.o.median()), 1), "makine medyanı": round(float(d.m.median()), 1)}

    agr = pd.DataFrame([agree(ok.qrs_ms, m_qrs, "QRS süresi (düşük eşik)"),
                        agree(ok.qrs_high_ms, m_qrs, "QRS süresi (görünür başlangıç)"),
                        agree(ok.pr_ms, m_pr, "PR aralığı")])
    qc = {"n": len(r), "başarılı": int(r.ok.sum()),
          "hatalar": r.loc[~r.ok, "err"].str.split(":").str[0].value_counts().to_dict(),
          "baskın küme oranı medyanı": round(float(ok.dom_frac.median()), 3),
          "küme sayısı dağılımı": ok.n_clusters.value_counts().sort_index().to_dict(),
          "başlangıç sınıra dayanan (%)": round(100 * float(ok.onset_capped.mean()), 2),
          "PR-sonu taban penceresi eğimi medyanı (mV)": round(float(ok.pq_base_slope.median()), 4),
          "yavaş başlangıç medyanı (ms)": round(float(ok.slow_init_ms.median()), 1)}
    return agr, qc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--p-thr", type=float, default=None, help="hdsig.P_THR_FRAC denemesi")
    a = ap.parse_args()
    if a.p_thr is not None:                      # alt süreçlere ortam değişkeniyle geçer
        import os
        os.environ["HD_P_THR"] = str(a.p_thr)
        hdsig.P_THR_FRAC = a.p_thr
    OUT.mkdir(parents=True, exist_ok=True)
    syn = synthetic_tests()
    agr, qc = mimic_agreement(a.n)
    md = ["# Aşama 2 — pipeline doğrulaması\n", "## 1. Sentetik testler\n", syn.round(3).to_markdown(index=False),
          "\nBeklenen: delta'da `slow_init_ms` uzun ve `pqhigh_fq_corr` > 0; gibbs'te `pqhigh_fq_corr` < 0.\n",
          f"## 2. MIMIC makine ölçümleriyle uyum (rastgele {a.n} T0 EKG, preeksitasyon/kısa PR hastaları hariç)\n",
          agr.to_markdown(index=False), "\n**Kalite:**\n", "```", json.dumps(qc, ensure_ascii=False, indent=1), "```"]
    (OUT / "validation.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
