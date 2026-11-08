"""HIDDEN-DELTA sinyal işleme kütüphanesi (concept.md §6, §6.1).

Zincir: yükle (mV, standart derivasyon sırası) → taban çizgisi (medyan filtre) → vuru tespiti (uzaysal hız)
→ çapraz korelasyonla hizalama → morfoloji kümeleme → küme başına medyan vuru → global QRS sınırları
(8 bağımsız derivasyonun uzaysal hızından) → Kors XYZ → §6.1 özellikleri.

Filtre ilkesi (§6): hiçbir doğrusal ya da sıfır fazlı alçak geçiren filtre uygulanmaz. Taban çizgisi
medyan filtreyle (0.2 s + 0.6 s) kestirilip çıkarılır; medyan filtre doğrusal değildir ve keskin kenarlarda
Gibbs salınımı üretmez. Uzaysal hızdaki kısa hareketli ortalamalar yalnızca sınır tespiti içindir; özellik
genlikleri ham (taban çizgisi düzeltilmiş) medyan vurudan ölçülür.

Zaman referansı: medyan vuru penceresi fiducial'ın 400 ms öncesinden 600 ms sonrasına uzanır; tüm *_ms
özellikleri fiducial'a göredir (negatif = önce), süreler ms cinsindendir. Genlikler mV.
"""
import os

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.signal import find_peaks
from scipy.spatial.distance import squareform

LEADS = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
LI = {name: i for i, name in enumerate(LEADS)}
IND8 = [0, 1, 6, 7, 8, 9, 10, 11]   # I, II, V1–V6: bağımsız derivasyonlar
# Kors regresyon matrisi (Kors JA ve ark., Eur Heart J 1990;11:1083–1092), sütunlar I, II, V1–V6.
# Frank eksenleri: X+ sol, Y+ aşağı, Z+ arka. Katsayılar Vondrak & Penhakert, IEEE JTEHM 2022;10:1900208,
# Tablo 1 ile birebir doğrulandı (2026-11-08).
KORS = np.array([[0.38, -0.07, -0.13, 0.05, -0.01, 0.14, 0.06, 0.54],
                 [-0.07, 0.93, 0.06, -0.02, -0.05, 0.06, -0.17, 0.13],
                 [0.11, -0.23, -0.43, -0.06, -0.14, -0.20, -0.11, 0.31]])
PRE_S, POST_S = 0.40, 0.60
CORR_MIN = 0.90          # kümeleme: QRS penceresinde korelasyon eşiği
P_THR_FRAC = float(os.environ.get("HD_P_THR", 0.20))   # P başlangıç/bitişi: P tepe uzaysal hızının bu oranı
MIN_BEATS = 3            # medyan vuru için en az vuru


class SignalError(Exception):
    pass


# ------------------------------------------------------------------ temel
def load_record(rec):
    import wfdb
    r = wfdb.rdrecord(str(rec))
    names = [n.upper() for n in r.sig_name]
    sig = np.nan_to_num(r.p_signal.astype(np.float64)).T
    scale = {"MV": 1.0, "UV": 1e-3, "V": 1e3}
    out = []
    for lead in LEADS:
        if lead.upper() not in names:
            raise SignalError(f"derivasyon yok: {lead}")
        j = names.index(lead.upper())
        out.append(sig[j] * scale.get(str(r.units[j]).upper(), 1.0))
    return int(round(r.fs)), np.vstack(out)


def remove_baseline(X, fs):
    w1, w2 = int(0.2 * fs) | 1, int(0.6 * fs) | 1
    base = median_filter(median_filter(X, size=(1, w1), mode="nearest"), size=(1, w2), mode="nearest")
    return X - base


def spatial_velocity(X, fs):
    d = np.diff(X[..., IND8, :], axis=-1) * fs         # mV/s; 12×N ya da vuru×12×L
    sv = np.sqrt((d ** 2).sum(axis=-2))
    return np.concatenate([sv[..., :1], sv], axis=-1)


def _ms(n, fs):
    return 1000.0 * n / fs


# ------------------------------------------------------------------ vurular
def detect_beats(X, fs):
    sv = uniform_filter1d(spatial_velocity(X, fs), size=max(1, int(0.04 * fs)))
    thr = 0.35 * np.percentile(sv, 99.5)
    pk, _ = find_peaks(sv, height=thr, distance=int(0.25 * fs))
    return pk


def segment(X, pk, fs):
    pre, post = int(PRE_S * fs), int(POST_S * fs)
    pk = pk[(pk >= pre) & (pk < X.shape[1] - post)]
    if len(pk) == 0:
        return pk, np.empty((0, X.shape[0], pre + post))
    return pk, np.stack([X[:, p - pre:p + post] for p in pk])


def align(beats, fs, maxlag_s=0.02):
    """Her vuruyu uzaysal hız şablonuna ±maxlag içinde hizalar; örnek cinsinden kaydırmaları döndürür."""
    sv = spatial_velocity(beats, fs)                    # nb × L
    pre, w, ml = int(PRE_S * fs), int(0.10 * fs), int(maxlag_s * fs)
    tmpl = np.median(sv[:, pre - w:pre + w], axis=0)
    tmpl = tmpl - tmpl.mean()
    lags = np.arange(-ml, ml + 1)
    out = []
    for s in sv:
        cc = [np.dot(s[pre - w + lag:pre + w + lag] - s[pre - w + lag:pre + w + lag].mean(), tmpl) for lag in lags]
        out.append(int(lags[int(np.argmax(cc))]))
    return np.array(out)


def cluster_beats(beats, fs, corr_min=CORR_MIN):
    """QRS penceresinde (fiducial −80…+100 ms, 8 bağımsız derivasyon) ortalama bağlantılı hiyerarşik kümeleme.
    Etiketler büyükten küçüğe 0, 1, … (0 = baskın küme). Ayrıca küme içi ortalama korelasyonu döndürür."""
    nb = len(beats)
    if nb == 1:
        return np.zeros(1, int), np.ones(1)
    pre = int(PRE_S * fs)
    F = beats[:, IND8, pre - int(0.08 * fs):pre + int(0.10 * fs)].reshape(nb, -1)
    C = np.corrcoef(F)
    D = np.clip(1 - C, 0, 2)
    np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), t=1 - corr_min, criterion="distance")
    sizes = {c: (lab == c).sum() for c in np.unique(lab)}
    order = {c: i for i, c in enumerate(sorted(sizes, key=lambda c: -sizes[c]))}
    lab = np.array([order[c] for c in lab])
    within = np.array([C[np.ix_(lab == c, lab == c)].mean() for c in range(len(order))])
    return lab, within


# ------------------------------------------------------------------ sınırlar
def _noise_floor(sv, lo, hi, win):
    if hi - lo <= win:
        seg = sv[max(0, lo):max(lo + win, hi)]
        return float(seg.mean()), float(seg.std())
    means = uniform_filter1d(sv[lo:hi], size=win, mode="nearest")
    j = int(np.argmin(means[win // 2: len(means) - win // 2])) + win // 2
    seg = sv[lo + j - win // 2: lo + j + win // 2]
    return float(seg.mean()), float(seg.std())


def qrs_bounds(mb, fs):
    """Global QRS başlangıcı/bitişi. on_low: gürültü tabanına göre düşük eşik (ilk elektriksel aktivite);
    on_high: tepe uzaysal hızın %15'i (görünür QRS başlangıcı). Aradaki fark yavaş başlangıç bileşenidir."""
    sv = uniform_filter1d(spatial_velocity(mb, fs), size=max(1, int(0.01 * fs)))
    pre = int(PRE_S * fs)
    a, b = pre - int(0.04 * fs), pre + int(0.04 * fs)
    pk = a + int(np.argmax(sv[a:b]))
    n_mu, n_sd = _noise_floor(sv, pk - int(0.30 * fs), pk - int(0.06 * fs), int(0.04 * fs))
    peak = float(sv[pk])
    thr_low = max(n_mu + 4 * n_sd, 0.04 * peak)
    thr_high = 0.15 * peak
    lim = pk - int(0.14 * fs)                           # en çok 140 ms geriye

    def back(thr):
        i = pk
        while i > lim and sv[i] > thr:
            i -= 1
        return i + 1, i <= lim

    on_low, capped = back(thr_low)
    on_high, _ = back(thr_high)
    # bitiş: son anlamlı QRS tepesinden (≥%25 tepe, tepeden en çok 120 ms sonra) ileri, eşiğin altına ilk düşüş
    rel, _ = find_peaks(sv[pk:pk + int(0.12 * fs)], height=0.25 * peak)
    last = pk + int(rel[-1]) if len(rel) else pk
    thr_off = max(n_mu + 4 * n_sd, 0.08 * peak)
    j, end = last, min(len(sv) - 1, last + int(0.12 * fs))
    while j < end and sv[j] > thr_off:
        j += 1
    off = j
    return dict(pk=pk, on_low=on_low, on_high=on_high, off=off, capped=bool(capped),
                sv_peak=peak, noise_mu=n_mu, noise_sd=n_sd)


def p_bounds(mb, fs, on):
    sv = uniform_filter1d(spatial_velocity(mb, fs), size=max(1, int(0.02 * fs)))
    a, b = max(1, on - int(0.30 * fs)), on - int(0.03 * fs)
    if b - a < int(0.05 * fs):
        return None
    pp = a + int(np.argmax(sv[a:b]))
    n_mu, n_sd = _noise_floor(sv, 0, max(int(0.06 * fs), pp - int(0.06 * fs)), int(0.04 * fs))
    if sv[pp] <= n_mu + 3 * n_sd:
        return None
    thr = max(n_mu + 3 * n_sd, P_THR_FRAC * float(sv[pp]))   # P sınırları: P tepe hızının belli oranı
    i = pp
    while i > 0 and sv[i] > thr:
        i -= 1
    j = pp
    while j < on - 1 and sv[j] > thr:
        j += 1
    return dict(p_peak=pp, p_on=i + 1, p_end=j, p_sv_peak=float(sv[pp]))


# ------------------------------------------------------------------ özellikler
def _xyz(v12):
    return KORS @ np.asarray(v12)[IND8]


def _angle(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return np.nan
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / (nu * nv), -1, 1))))


def beat_features(mb, fs):
    """Bir medyan vurudan §6.1 özellikleri."""
    f = {}
    pre = int(PRE_S * fs)
    n = lambda ms: int(round(ms * fs / 1000))
    qb = qrs_bounds(mb, fs)
    on, onh, off = qb["on_low"], qb["on_high"], qb["off"]
    f.update(on_low_ms=_ms(on - pre, fs), on_high_ms=_ms(onh - pre, fs), off_ms=_ms(off - pre, fs),
             slow_init_ms=_ms(onh - on, fs), qrs_ms=_ms(off - on, fs), qrs_high_ms=_ms(off - onh, fs),
             onset_capped=qb["capped"], sv_peak=qb["sv_peak"], sv_noise=qb["noise_mu"])
    pb = p_bounds(mb, fs, on)
    if pb:
        f.update(p_on_ms=_ms(pb["p_on"] - pre, fs), p_end_ms=_ms(pb["p_end"] - pre, fs),
                 pr_ms=_ms(on - pb["p_on"], fs), pr_seg_ms=_ms(on - pb["p_end"], fs))
    # İki referans: (1) morfoloji/vektör özellikleri kayıt düzeyi medyan filtre taban çizgisine göre (Y = mb);
    # (2) QRS öncesi lob özellikleri yerel PR-sonu düzeyine göre (Yl): düşük eşikli başlangıçtan 40–22 ms önce,
    # yani −20…−2 ms ölçüm penceresinin hemen öncesi. Kısa PR'de bu pencere P'ye taşabilir; eğimi kalite
    # özelliği olarak saklanır (pq_base_slope: pencerenin iki yarısı arasındaki fark, 8 derivasyon normu).
    Y = mb
    b0, bm, b1 = on - n(40), on - n(31), on - n(22)
    iso = np.median(mb[:, b0:b1], axis=1)
    Yl = mb - iso[:, None]
    half = np.median(mb[:, bm:b1], axis=1) - np.median(mb[:, b0:bm], axis=1)
    f["pq_base_slope"] = float(np.linalg.norm(half[IND8]))
    f["pq_base_level"] = float(np.linalg.norm(iso[IND8]))

    # başlangıç vektörleri ve eğimler
    for t in (10, 20, 40):
        v = _xyz(Y[:, on + n(t)])
        f.update({f"v{t}_x": v[0], f"v{t}_y": v[1], f"v{t}_z": v[2], f"v{t}_mag": float(np.linalg.norm(v))})
    for t in (20, 40):
        v = _xyz(Y[:, on:on + n(t)].mean(axis=1))
        f.update({f"i{t}_x": v[0], f"i{t}_y": v[1], f"i{t}_z": v[2]})
    slope = (Y[:, on + n(20)] - Y[:, on]) / 0.020   # mV/s, ilk 20 ms
    f.update({f"slope20_{lead}": float(slope[i]) for i, lead in enumerate(LEADS)})
    sv = uniform_filter1d(spatial_velocity(mb, fs), size=max(1, int(0.01 * fs)))
    for t in (5, 10, 20):
        f[f"svr{t}"] = float(sv[on + n(t)] / qb["sv_peak"])

    # septal q ve V1
    for lead in ("I", "aVL", "V5", "V6", "V1"):
        seg = Y[LI[lead], on:on + n(60)]
        pos = np.where(seg > 0.05)[0]
        head = seg[: pos[0]] if len(pos) else seg
        qa = float(max(0.0, -head.min())) if len(head) else 0.0
        f[f"q_{lead}"] = qa if qa >= 0.025 else 0.0
    v1 = Y[LI["V1"], on:max(off, on + 1)]
    r_amp, s_amp = float(max(0.0, v1.max())), float(max(0.0, -v1.min()))
    f.update(v1_r=r_amp, v1_s=s_amp, v1_rs=r_amp / s_amp if s_amp > 0.01 else np.nan)
    pk_pos, _ = find_peaks(v1, prominence=0.05)
    f["v1_rsr"] = bool(len(pk_pos) >= 2 and v1[pk_pos[1]] > 0.05)

    # QRS öncesi defleksiyon: hem düşük eşikli hem görünür (yüksek eşikli) başlangıca göre. Her biri kendi
    # başlangıcından sonraki ilk 15 ms'lik QRS sapmasıyla karşılaştırılır: filtre salınımında QRS öncesi lob
    # görünür ilk sapmanın tersidir (fq_corr < 0); delta dalgasında aynı yöndedir (fq_corr > 0).
    for tag, o in (("low", on), ("high", onh)):
        fq = Yl[:, o:o + n(15)].mean(axis=1)
        pq = Yl[:, o - n(20):o - n(2)].mean(axis=1)
        f.update({f"pq{tag}_{lead}": float(pq[i]) for i, lead in enumerate(LEADS)})
        v = _xyz(pq)
        f.update({f"pq{tag}_x": v[0], f"pq{tag}_y": v[1], f"pq{tag}_z": v[2],
                  f"pq{tag}_mag": float(np.linalg.norm(pq[IND8]))})
        sel = (np.abs(pq) > 0.010) & (np.abs(fq) > 0.020)
        f[f"pq{tag}_neval"] = int(sel.sum())
        f[f"pq{tag}_opp_frac"] = float(np.mean(np.sign(pq[sel]) == -np.sign(fq[sel]))) if sel.any() else np.nan
        f[f"pq{tag}_fq_corr"] = float(np.corrcoef(pq, fq)[0, 1]) if pq.std() > 0 and fq.std() > 0 else np.nan

    # QRS ve T vektörleri
    qrs_int = _xyz(Y[:, on:off].sum(axis=1) / fs * 1000)          # mV·ms
    t_end = min(Y.shape[1], off + n(400))
    t_a = min(off + n(40), t_end - 1)
    t_int = _xyz(Y[:, off:t_end].sum(axis=1) / fs * 1000)
    tx = KORS @ Y[IND8, t_a:t_end]
    tp = t_a + int(np.argmax(np.linalg.norm(tx, axis=0)))
    tv = _xyz(Y[:, tp])
    f.update(qrs_int_x=qrs_int[0], qrs_int_y=qrs_int[1], qrs_int_z=qrs_int[2],
             t_int_x=t_int[0], t_int_y=t_int[1], t_int_z=t_int[2],
             tpk_x=tv[0], tpk_y=tv[1], tpk_z=tv[2], tpk_ms=_ms(tp - pre, fs),
             qrst_angle=_angle(qrs_int, t_int), vg_mag=float(np.linalg.norm(qrs_int + t_int)))
    return f


# ------------------------------------------------------------------ kayıt düzeyi
def process_signal(X, fs):
    """Taban çizgisi düzeltilmemiş 12×N sinyalden (mV) özellikler ve medyan vurular.
    Döndürür: (özellik sözlüğü, {küme: medyan vuru 12×L float32})."""
    if fs != 500:
        raise SignalError(f"fs={fs}")
    flat = int((X.std(axis=1) < 0.005).sum())
    if flat > 2:
        raise SignalError(f"düz derivasyon: {flat}")
    Xb = remove_baseline(X, fs)
    pk = detect_beats(Xb, fs)
    pk, beats = segment(Xb, pk, fs)
    if len(pk) < MIN_BEATS:
        raise SignalError(f"az vuru: {len(pk)}")
    lags = align(beats, fs)
    pk, beats = segment(Xb, pk + lags, fs)
    if len(pk) < MIN_BEATS:
        raise SignalError(f"az vuru (hizalama sonrası): {len(pk)}")
    lab, within = cluster_beats(beats, fs)
    rr = np.diff(pk) / fs
    rec = dict(n_beats=len(pk), n_clusters=int(lab.max() + 1), n_flat=flat,
               rr_mean=float(rr.mean()) if len(rr) else np.nan, rr_sd=float(rr.std()) if len(rr) else np.nan,
               hr=float(60 / rr.mean()) if len(rr) else np.nan)
    mbs = {}
    for c in range(min(2, rec["n_clusters"])):
        idx = lab == c
        if idx.sum() < MIN_BEATS:
            continue
        mb = np.median(beats[idx], axis=0)
        mbs[c] = mb.astype(np.float32)
        try:
            bf = beat_features(mb, fs)
        except Exception as e:  # noqa: BLE001
            if c == 0:
                raise SignalError(f"özellik hatası: {type(e).__name__}: {e}")
            continue
        prefix = "" if c == 0 else "c2_"
        rec.update({prefix + k: v for k, v in bf.items()})
        rec[prefix + "n_in_cluster"] = int(idx.sum())
        rec[prefix + "within_corr"] = float(within[c])
    if 0 not in mbs:
        raise SignalError("baskın kümede az vuru")
    rec["dom_frac"] = rec["n_in_cluster"] / rec["n_beats"]
    return rec, mbs


def process_record(rec_path):
    fs, X = load_record(rec_path)
    return process_signal(X, fs)


# ------------------------------------------------------------------ ayrık QRS öncesi dalga (ANALYSIS_PLAN §6)
# Uygulama notu: plan "8 bağımsız derivasyon" ve "I–aVL, II–III–aVF, V1…V6 komşulukları" diyor; komşuluk şartı türetilmiş
# derivasyonlar olmadan anlamsız olduğundan dalga 12 derivasyonda aranır, komşuluk aşağıdaki çiftlerle tanımlanır.
WAVE_PAIRS = [("I", "aVL"), ("II", "III"), ("III", "aVF"), ("II", "aVF"),
              ("V1", "V2"), ("V2", "V3"), ("V3", "V4"), ("V4", "V5"), ("V5", "V6")]


def discrete_wave(mb, fs, on, amp_min=0.008, fwhm_ms=(6, 30), sync_ms=6, win_ms=60, lead_gap_ms=4):
    """QRS düşük eşikli başlangıcından (on, örnek) önceki 60 ms'de, pencere medyanına göre:
    tepe ≥ amp_min (mV), en geç başlangıçtan lead_gap_ms önce, yerel ekstremum; tepeden önce ve sonra (QRS'ten önce)
    tepenin yarısının altına iniş; yarı yükseklik genişliği fwhm_ms aralığında. EKG düzeyi: en az bir komşu çiftte,
    tepe zamanları sync_ms içinde. Döndürür: (bayrak, ayrıntı sözlüğü)."""
    n = lambda ms: int(round(ms * fs / 1000))
    a = on - n(win_ms)
    if a < 1:
        return False, {}
    W = mb[:, a:on]
    D = W - np.median(W, axis=1, keepdims=True)
    hits = {}
    for j, lead in enumerate(LEADS):
        d = D[j]
        last = len(d) - n(lead_gap_ms)
        best = None
        for i in range(1, last):
            if not ((d[i] > d[i - 1] and d[i] >= d[i + 1]) or (d[i] < d[i - 1] and d[i] <= d[i + 1])):
                continue
            amp = abs(d[i])
            if amp < amp_min or (best is not None and amp <= best[1]):
                continue
            half = 0.5 * amp
            if np.abs(d[:i]).min() >= half or np.abs(d[i + 1:]).min() >= half:
                continue
            lo = i
            while lo > 0 and abs(d[lo - 1]) >= half:
                lo -= 1
            hi = i
            while hi < len(d) - 1 and abs(d[hi + 1]) >= half:
                hi += 1
            w = (hi - lo + 1) * 1000 / fs
            if fwhm_ms[0] <= w <= fwhm_ms[1]:
                best = (i, amp, w, float(np.sign(d[i])))
        if best:
            hits[lead] = best
    pairs = [(x, y) for x, y in WAVE_PAIRS if x in hits and y in hits and abs(hits[x][0] - hits[y][0]) <= n(sync_ms)]
    if not pairs:
        return False, {"n_leads": len(hits)}
    leads = sorted({z for p in pairs for z in p})
    return True, {"n_leads": len(hits), "n_pair_leads": len(leads),
                  "amp": float(np.mean([hits[z][1] for z in leads])), "fwhm": float(np.mean([hits[z][2] for z in leads])),
                  "peak_ms": float(np.mean([(hits[z][0] - (on - a)) * 1000 / fs for z in leads])),
                  "polarity": {z: hits[z][3] for z in leads}}
