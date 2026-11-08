"""HIDDEN-DELTA — gözlenen QRS öncesi ayrık dalga ekran filtresinin ürünü mü?

Gözlem (yazar): önceki IAB etiketleme çalışmasında (HEEDB + MIMIC) medyan vurularda Q'dan hemen önce,
QRS'ten ayrı, ~20 ms süren, ~10 µV'lik bir dalga. O çalışmanın medyan vuru kodu (LEAD-IAB/scripts/p_morphology.py)
ham sinyale önce 0.5–40 Hz, 2. derece Butterworth'u sıfır fazlı (sosfiltfilt) uyguluyor.

Test: (1) filtrenin basamak ve rampa kenarına yanıtı (analitik); (2) sentetik EKG; (3) rastgele MGH K1 kontrol
kayıtları: aynı kayıttan filtresiz (bizim zincir) ve IAB filtreli medyan vuru; QRS sınırları filtresiz vurudan.

Ölçüm (derivasyon başına, 8 bağımsız derivasyon): düşük eşikli QRS başlangıcından önceki −30…−4 ms'de, başlangıçtan
80–50 ms önceki düzeye göre en büyük sapma; filtreli − filtresiz fark filtrenin eklediği kısımdır. Kenar testinde de
sapma aynı yerel düzeye göre ölçülür (uzak sıfıra göre ölçmek 0.5 Hz yüksek geçirenin geniş kaymasını katar).
Çıktı yalnızca agregat: stage2/out/wave_check.md
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage0"))
import hdsig  # noqa: E402
from s2_validate import synth  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
P2 = Path("/home/fatih/data/hidden-delta/stage2/I0001")
FS = 500
BP = signal.butter(2, [0.5, 40], btype="band", fs=FS, output="sos")   # IAB çalışmasındaki filtre
SEED, N_REC = 20261115, 1000


def iab_filter(X):
    return signal.sosfiltfilt(BP, X, axis=1)


def edge_response():
    """Filtrenin QRS benzeri, taban çizgisine dönen geçici sinyallere yanıtı. Başlangıçtan önceki −40…−4 ms'de
    en büyük sapma; ilk QRS sapmasına göre işareti; yarı yükseklik genişliği ve QRS'ten ayrı olup olmadığı."""
    n, t0 = 4000, 2000
    t = (np.arange(n) - t0) * 1000 / FS                      # ms, başlangıç = 0
    tri = lambda c, w, a: a * np.clip(1 - np.abs(t - c) / w, 0, None)
    shapes = {"R 1 mV (üçgen, 40 ms)": tri(20, 20, 1.0),
              "qR: q −0.1 mV + R 1 mV": tri(5, 5, -0.1) + tri(30, 20, 1.0),
              "rS: r 0.2 mV + S −1 mV": tri(10, 10, 0.2) + tri(40, 20, -1.0),
              "geniş R 1 mV (üçgen, 80 ms)": tri(40, 40, 1.0)}
    rows = []
    for name, x in shapes.items():
        y = signal.sosfiltfilt(BP, x)
        y = y - np.median(y[(t >= -80) & (t <= -50)])          # yerel PR düzeyine göre (gerçek veriyle aynı)
        m = (t >= -40) & (t <= -4)
        i = int(np.argmax(np.abs(y[m])))
        amp = y[m][i]
        first = x[(t > 0) & (t <= 15)].mean()
        w = np.abs(y[m]) >= 0.5 * abs(amp)
        rows.append({"girdi": name, "QRS öncesi en büyük sapma (µV)": round(1000 * amp, 1),
                     "zamanı (ms)": round(float(t[m][i]), 1), "yarı yükseklik genişliği (ms)": round(w.sum() * 1000 / FS, 1),
                     "ilk QRS sapmasına göre": "ters" if np.sign(amp) != np.sign(first) else "aynı",
                     "−4 ms'deki değer (µV)": round(1000 * float(y[t == -4][0]), 1)})
    return pd.DataFrame(rows)


def pre_qrs_dev(mb, on, off):
    """Derivasyon başına (8 bağımsız): başlangıçtan önceki −30…−4 ms'de, başlangıçtan 80–50 ms önceki düzeye göre
    en büyük sapma (mV) ve QRS'in baskın sapması (başlangıç–bitiş arasında en büyük mutlak değer, işaretli)."""
    n = lambda ms: int(round(ms * FS / 1000))
    base = np.median(mb[:, on - n(80):on - n(50)], axis=1)
    Y = mb - base[:, None]
    W = Y[hdsig.IND8, on - n(30):on - n(4) + 1]
    e = W[np.arange(8), np.argmax(np.abs(W), axis=1)]
    Q = Y[hdsig.IND8, on:max(off, on + 1)]
    dom = Q[np.arange(8), np.argmax(np.abs(Q), axis=1)]
    return e, dom


def pair(X):
    """Aynı sinyalden filtresiz ve IAB filtreli medyan vuru; QRS sınırları filtresiz vurudan."""
    _, mb_raw = hdsig.process_signal(X, FS)
    _, mb_f = hdsig.process_signal(iab_filter(X), FS)
    raw, filt = mb_raw[0].astype(float), mb_f[0].astype(float)
    qb = hdsig.qrs_bounds(raw, FS)
    e_raw, dom = pre_qrs_dev(raw, qb["on_low"], qb["off"])
    e_f, _ = pre_qrs_dev(filt, qb["on_low"], qb["off"])
    return np.c_[e_raw, e_f, dom]


def one_record(rec):
    try:
        fs, X = hdsig.load_record(rec)
        if fs != FS:
            return None
        return pair(X)
    except Exception:  # noqa: BLE001
        return None


def summarise(res, label):
    A = np.vstack(res)                       # (kayıt × 8) × [ham, filtreli, baskın QRS]
    e_raw, e_f, dom = A[:, 0] * 1000, A[:, 1] * 1000, A[:, 2] * 1000
    added = e_f - e_raw
    big = np.abs(dom) >= 200
    return [{"veri": label, "derivasyon": len(A),
             "|sapma| medyanı, filtresiz (µV)": round(float(np.median(np.abs(e_raw))), 1),
             "|sapma| medyanı, IAB filtreli (µV)": round(float(np.median(np.abs(e_f))), 1),
             "≥10 µV, filtresiz (%)": round(100 * float(np.mean(np.abs(e_raw) >= 10)), 1),
             "≥10 µV, IAB filtreli (%)": round(100 * float(np.mean(np.abs(e_f) >= 10)), 1),
             "filtrenin eklediği, baskın QRS'e ters (%)": round(100 * float(np.mean(np.sign(added) != np.sign(dom))), 1),
             "aynısı, baskın QRS ≥0.2 mV olan derivasyonlarda (%)": round(100 * float(np.mean(np.sign(added[big]) != np.sign(dom[big]))), 1),
             "eklenen sapma / baskın QRS (medyan, %)": round(100 * float(np.median(np.abs(added) / np.maximum(np.abs(dom), 1))), 2),
             "Spearman r (|eklenen|, |baskın QRS|)": round(float(pd.Series(np.abs(added)).corr(pd.Series(np.abs(dom)), method="spearman")), 3)}]


def main():
    rng = np.random.default_rng(SEED)
    edge = edge_response()
    syn = [pair(synth("normal", seed=s)) for s in range(20)]
    jobs = pd.read_parquet(P2 / "jobs.parquet")
    feats = duckdb.connect().execute(f"""SELECT k FROM read_parquet('{P2}/features_*.parquet', union_by_name=true)
                                         WHERE role = 'k1' AND status = 'ok'""").df()
    recs = jobs[jobs.k.isin(feats.k)].sample(n=N_REC, random_state=int(rng.integers(2**31))).rec.tolist()
    with ProcessPoolExecutor(28) as ex:
        real = [r for r in ex.map(one_record, recs, chunksize=8) if r is not None]
    tab = pd.DataFrame(summarise(syn, "sentetik normal EKG (20)") + summarise(real, f"MGH K1 kontrol ({len(real)})"))
    OUT.mkdir(parents=True, exist_ok=True)
    md = ["# QRS öncesi ayrık dalga: ekran filtresi testi\n",
          "Gözlenen dalga: ~20 ms, ~10 µV, QRS'ten ayrı, Q'dan hemen önce (IAB etiketleme çalışması; medyan vurular "
          "0.5–40 Hz 2. derece Butterworth ile sıfır fazlı filtrelenmiş sinyalden).\n",
          "## 1. Filtrenin keskin kenara yanıtı\n", edge.to_markdown(index=False),
          "\n## 2. Aynı kayıtta filtresiz ve IAB filtreli medyan vuru: filtrenin QRS öncesine eklediği sapma\n",
          "Derivasyon başına (8 bağımsız derivasyon), düşük eşikli QRS başlangıcından önceki −30…−4 ms'de, başlangıçtan "
          "80–50 ms önceki düzeye göre en büyük sapma. P kuyruğu iki zincirde aynı olduğundan fark filtreye aittir. "
          "Baskın QRS: başlangıç–bitiş arasında en büyük mutlak değerli örnek. Kontroller preeksitasyonu olmayan "
          "hastalardır; etiket kullanılmaz.\n", tab.to_markdown(index=False), ""]
    (OUT / "wave_check.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
