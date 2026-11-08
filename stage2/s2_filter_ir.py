"""HIDDEN-DELTA Aşama 2 — cihaz zincirinin dürtü yanıtının pace spike'larından ampirik tahmini (concept.md §6 adım 5).

Pace spike'ı (~0.5 ms) 500 Hz örneklemede yaklaşık bir dürtüdür. Kayıttaki şekli, anti-aliasing/desimasyon
filtresi dahil tüm cihaz zincirinin dürtü yanıtını taşır. Spike'tan ÖNCEKİ kısım (pre-ringing) temizdir:
spike'tan önce taban çizgisi vardır. Sonraki kısım ise ardından gelen QRS/P ile kirlenir, bu yüzden yalnızca
spike öncesi raporlanır.

Yöntem: taban çizgisi düzeltilmiş sinyalden 18 ms medyan filtre artığı → dar (≤6 ms) ve büyük (≥10×MAD)
sapmalar spike adayı. Spike, en büyük sapmanın olduğu derivasyonda işaret ve genlikle normalize edilir; ±30 ms
pencere kayıt içinde, sonra kayıtlar arasında ortalanır (her kayıt eşit ağırlık). Sentetik/yeniden
yapılandırılmış spike işareti için spike'ın derivasyonlar arası genlik profili ve genişliği raporlanır.

Gruplar: MIMIC'te filtering × bandwidth; HEEDB'de site. Çıktı yalnızca agregat (ortalama dürtü yanıtları).
Çalıştırma: python3 stage2/s2_filter_ir.py [--per-group 300]
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage0"))
import hdsig  # noqa: E402
from common import HEEDB, MIMIC_ECG, PRIVATE  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
FS = 500
HALF = 15                       # ±30 ms
SEED = 20261113


def spikes_in_record(rec):
    try:
        fs, X = hdsig.load_record(rec)
    except Exception:  # noqa: BLE001
        return None
    if fs != FS:
        return None
    X = hdsig.remove_baseline(X, fs)
    R = X - median_filter(X, size=(1, 9), mode="nearest")          # 18 ms medyan artığı
    mad = np.median(np.abs(R - np.median(R, axis=1, keepdims=True)), axis=1) + 1e-6
    Z = np.abs(R) / mad[:, None]
    score = Z.max(axis=0)
    cand = np.where(score > 10)[0]
    resp, widths, profiles, n_cand = [], [], [], 0
    last = -10
    for i in cand:
        if i - last < 50 or i < 30 or i > X.shape[1] - 30:           # en az 100 ms aralık
            continue
        n_cand += 1
        lead = int(np.argmax(Z[:, i]))
        p = i - 2 + int(np.argmax(np.abs(R[lead, i - 2:i + 3])))
        base = np.median(X[:, p - 12:p - 4], axis=1)                 # spike'tan 24–8 ms önce
        w = X[:, p - HALF:p + HALF + 1] - base[:, None]
        amp = w[lead, HALF]
        if not np.isfinite(w).all() or abs(amp) < 0.3:              # en az 0.3 mV
            continue
        pre_flat = X[lead, p - 12:p - 4].std() <= 0.05 * abs(amp)    # spike öncesi düz taban çizgisi
        narrow = (abs(w[lead, HALF - 3]) <= 0.5 * abs(amp)) and (abs(w[lead, HALF + 3]) <= 0.5 * abs(amp))
        multi = int((np.abs(w[:, HALF]) >= 0.3 * abs(amp)).sum()) >= 3
        if not (pre_flat and narrow and multi):
            continue
        above = np.abs(w[lead, HALF - 5:HALF + 6]) > 0.5 * abs(amp)
        resp.append(w[lead] / amp)
        widths.append(int(above.sum()))
        profiles.append(w[:, HALF] / amp)
        last = i
    if len(resp) < 2:
        return None
    return dict(h=np.median(resp, axis=0), n=len(resp), n_cand=n_cand, width=float(np.mean(widths)),
                prof=np.mean(profiles, axis=0), prof_sd=float(np.std(profiles, axis=0).mean()))


def mimic_groups(per_group, rng):
    con = duckdb.connect()
    df = con.execute(f"""SELECT f.study, coalesce(f.filtering,'NA') || ' | ' || coalesce(f.bandwidth,'NA') grp, rl.path
        FROM '{PRIVATE / 'mimic_ecg_flags.parquet'}' f
        JOIN read_csv('{MIMIC_ECG / 'record_list.csv'}', all_varchar=true) rl ON rl.study_id::BIGINT = f.study
        WHERE f.paced""").df()
    out = []
    for g, d in df.groupby("grp"):
        d = d.sample(n=min(per_group, len(d)), random_state=int(rng.integers(2**31)))
        out.append(pd.DataFrame({"dataset": "MIMIC", "grp": g, "rec": [str(MIMIC_ECG / p) for p in d.path]}))
    return pd.concat(out)


def heedb_groups(per_group, rng):
    out = []
    for site, name in (("I0001", "MGH"), ("I0006", "Emory")):
        con = duckdb.connect()
        d = con.execute(f"""SELECT k FROM '{PRIVATE / f'heedb_{site}_ecg_flags.parquet'}' WHERE paced
                            AND (hash(k || '{SEED}') % 20) = 0""").df()
        d = d.sample(n=min(4 * per_group, len(d)), random_state=int(rng.integers(2**31)))
        recs = [str(HEEDB / site / "WFDB" / k) for k in d.k]
        out.append(pd.DataFrame({"dataset": name, "grp": name, "rec": recs}))
    return pd.concat(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-group", type=int, default=300)
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)
    jobs = pd.concat([mimic_groups(a.per_group, rng), heedb_groups(a.per_group, rng)], ignore_index=True)
    with ProcessPoolExecutor(24) as ex:
        res = list(ex.map(spikes_in_record, jobs.rec, chunksize=4))
    jobs["res"] = res
    t_ms = (np.arange(-HALF, HALF + 1) * 1000 / FS).astype(int)
    rows, curves = [], []
    for (ds, g), d in jobs.groupby(["dataset", "grp"]):
        ok = [r for r in d.res if r is not None]
        if len(ok) < 5:
            rows.append(dict(veri=ds, grup=g, kayit=len(d), spike_li_kayit=len(ok)))
            continue
        H = np.stack([r["h"] for r in ok])                         # 500 Hz'de 250 Hz kayıtlar zaten elenir
        H = H[np.isfinite(H).all(axis=1)]
        h = np.median(H, axis=0)
        at = lambda ms: round(float(h[t_ms == ms][0]), 3)
        far = (t_ms >= -20) & (t_ms <= -10)
        rows.append({"veri": ds, "grup": g, "kayıt": len(d), "spike'lı kayıt": len(ok),
                     "spike": int(sum(r["n"] for r in ok)),
                     "genişlik (örnek)": round(float(np.mean([r["width"] for r in ok])), 2),
                     "h(−8)": at(-8), "h(−6)": at(-6), "h(−4)": at(-4), "h(−2)": at(-2), "h(+2)": at(2),
                     "max|h| −20…−10": round(float(np.abs(h[far]).max()), 3),
                     "kayıtlar arası SD (−8…−2)": round(float(H[:, (t_ms >= -8) & (t_ms <= -2)].std(axis=0).mean()), 3)})
        curves.append(pd.DataFrame({"veri": ds, "grup": g, "t_ms": t_ms, "h": h}))
    tab = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    if curves:
        pd.concat(curves).to_csv(OUT / "filter_ir_curves.csv", index=False)
    md = ["# Pace spike'larından cihaz dürtü yanıtı (yalnızca spike öncesi)\n",
          "Sıkı spike ölçütleri: ≥0.3 mV; spike'tan 24–8 ms önce düz taban çizgisi (SD ≤ genliğin %5'i); ±6 ms'de "
          "genliğin yarısının altına iniş; ≥3 derivasyonda görünürlük. h: spike tepesine (t = 0, ms) göre normalize "
          "medyan yanıt (önce kayıt içi, sonra kayıtlar arası medyan). Spike sonrası kısım ardından gelen P/QRS ve "
          "polarizasyon artığıyla kirlenebilir; pre-ringing yalnızca t < 0'dan okunur.\n",
          tab.to_markdown(index=False),
          "\n**Yorum.** MIMIC'te t = −4 ms'de spike'ın tersi işaretli ~%15'lik bir lob ve −6…−8 ms'de ~%3'lük bir lob "
          "var: doğrusal fazlı bir zincirin pre-ringing'i. Lob filtre ayarından bağımsız; kaynağı örnekleme/anti-aliasing "
          "zinciri. \"<not specified> | 0.0005-150 Hz\" grubunda lob daha küçük (farklı cart). HEEDB'de spike'lar tek örnek "
          "genişliğinde ve öncesi tam sıfır: büyük olasılıkla cihazın eklediği yapay pace işaretleri. HEEDB'nin dürtü yanıtı "
          "bu yolla tahmin edilemez.\n", "Eğriler: `stage2/out/filter_ir_curves.csv`.\n"]
    (OUT / "filter_ir.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
