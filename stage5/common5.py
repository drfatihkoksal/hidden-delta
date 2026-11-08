"""HIDDEN-DELTA Aşama 5 — ortak araçlar (ANALYSIS_PLAN v1.2; OSF https://osf.io/kvn28).

Yollar, veri yükleme, türetilmiş özellikler, dondurulmuş modellerin uygulanması (M1, M0, M2, M-demo), oklüzyon,
integrated gradients, koşullu lojistik regresyon, Holm düzeltmesi.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
for sub in ("stage0", "stage2", "stage3", "stage4"):
    sys.path.insert(0, str(ROOT / sub))
import hdsig  # noqa: E402
import m0lib  # noqa: E402

D = Path("/home/fatih/data/hidden-delta")
P0, P1, P2, P3, P4 = (D / f"stage{i}" for i in range(5))
P5 = D / "stage5"
OUT = ROOT / "stage5" / "out"
FS, PRE, L = 500, 200, 500
SEED = 20261120
DEV = torch.device("cuda")
SITE_NAME = {"I0001": "MGH", "I0006": "Emory", "MIMIC": "MIMIC"}
L8 = ["I", "II", "V1", "V2", "V3", "V4", "V5", "V6"]


def ms2s(ms):
    return int(round(PRE + ms * FS / 1000))


# ------------------------------------------------------------------ veri
def load_features(site, keys):
    import duckdb
    con = duckdb.connect()
    con.register("want", pd.DataFrame({"k": list(keys)}))
    return con.execute(f"""SELECT f.* FROM read_parquet('{P2 / site}/features_*.parquet', union_by_name=true) f
                           JOIN want USING (k) WHERE f.status = 'ok'""").df()


def load_beats(site, keys, which="mb0"):
    want = {k: i for i, k in enumerate(keys)}
    out = np.full((len(keys), 12, L), np.nan, np.float32)
    for f in sorted((P2 / site).glob("beats_*.npz")):
        z = np.load(f)
        sel = [(j, want[k]) for j, k in enumerate(z["ids"]) if k in want]
        if sel:
            src, dst = zip(*sel)
            out[list(dst)] = z[which][list(src)]
    return out


def add_derived(df, mb):
    """M0 türetilmiş özellikleri ve ayrık dalga (plan §6). df ve mb aynı sırada."""
    df = df.copy()
    waves = [hdsig.discrete_wave(mb[i].astype(float), FS, ms2s(df.on_low_ms.iloc[i])) for i in range(len(df))]
    df["wave"] = [float(w[0]) for w in waves]
    df["septal_q_absent"] = (df[["q_I", "q_aVL", "q_V5", "q_V6"]] == 0).all(axis=1).astype(float)
    df["v1_rs_c"] = df.v1_r / np.maximum(df.v1_s, 0.01)
    df["male"] = (df.sex == "Male").astype(float)
    return df


def initial_profile(mb, on_low_ms):
    """İlk 20 ms'nin derivasyon başına ortalaması (8 bağımsız derivasyon)."""
    on = ms2s(on_low_ms)
    return mb[[hdsig.LI[x] for x in L8], on:on + 10].mean(axis=1)


def integral_profile(mb, a_ms, b_ms):
    a, b = ms2s(a_ms), min(L, ms2s(b_ms))
    return mb[[hdsig.LI[x] for x in L8], a:b].sum(axis=1) / FS * 1000


# ------------------------------------------------------------------ dondurulmuş modeller
class Frozen:
    def __init__(self):
        import s4_m1 as M1
        self.M1 = M1
        meta = json.loads((P4 / "m1" / "final_meta.json").read_text())
        assert meta["config"]["mode"] == "probe"
        self.depth = meta["config"]["depth"]
        self.enc = M1.Encoder(self.depth).to(DEV)
        self.enc.load_state_dict(torch.load(P4 / "m1" / "final_encoder.pt"))
        self.enc.eval()
        pr = np.load(P4 / "m1" / "final_probe.npz")
        self.w, self.b = pr["coef"].ravel().astype(np.float32), float(pr["intercept"][0])
        self.platt = meta["platt"]
        self.thr = meta["threshold_spec95"]
        self.jm = {n: json.loads((P4 / "final" / f"{n}.json").read_text()) for n in ("M0", "M-demo", "M2")}

    def embed(self, X):
        return self.M1.predict(self.enc, X.astype(np.float32), embed=True)

    def m1_logit(self, X):
        return self.embed(X) @ self.w + self.b

    def m1_prob(self, logit):
        a, c = self.platt
        return 1 / (1 + np.exp(-(a * logit + c)))

    def lin(self, name, df, Z=None):
        p = self.jm[name]
        ok = df[p["cont"] + p["bin"]].notna().all(axis=1).to_numpy()
        out = np.full(len(df), np.nan)
        X = m0lib.design(df[ok], p)
        if name == "M2":
            z = (Z[ok] - np.array(p["z_mean"])) / np.array(p["z_scale"])
            X = np.hstack([z, X])
        out[ok] = m0lib.lr_logit(X, p["lr"])
        return out

    def integrated_gradients(self, X, steps=64, batch=48):
        w = torch.tensor(self.w, device=DEV)
        out = np.empty_like(X, dtype=np.float32)
        alphas = (torch.arange(steps, device=DEV, dtype=torch.float32) + 0.5) / steps
        for i in range(0, len(X), batch):
            x = torch.tensor(X[i:i + batch], device=DEV, dtype=torch.float32)
            xs = (alphas[:, None, None, None] * x[None]).reshape(-1, *x.shape[1:]).requires_grad_(True)
            logit = self.enc.embed(xs) @ w + self.b
            g, = torch.autograd.grad(logit.sum(), xs)
            out[i:i + batch] = (x * g.reshape(steps, *x.shape).mean(0)).detach().cpu().numpy()
        return out


def windows(on_low_ms, off_ms):
    on, off = ms2s(on_low_ms), ms2s(off_ms)
    return (max(0, on - 20), min(L, on + 20)), (min(L, off + 20), min(L, off + 200))


def masked(X, wins):
    Y = X.copy()
    for i, (a, b) in enumerate(wins):
        Y[i, :, a:b] = 0.0
    return Y


def random_windows(widths, seed):
    rng = np.random.default_rng(seed)
    return [(s, s + w) for w in widths for s in [int(rng.integers(0, L - w + 1))]]


# ------------------------------------------------------------------ istatistik
def clogit(set_id, y, x, iters=50):
    """Tek ikili/sürekli değişkenli koşullu lojistik (küme başına bir vaka). Döndürür: beta, se."""
    df = pd.DataFrame({"s": set_id, "y": y, "x": x}).dropna()
    groups = [g for _, g in df.groupby("s") if g.y.sum() == 1 and len(g) > 1]
    beta = 0.0
    for _ in range(iters):
        grad = hess = 0.0
        for g in groups:
            xv = g.x.to_numpy(float)
            e = np.exp(beta * xv - (beta * xv).max())
            p = e / e.sum()
            m = (p * xv).sum()
            grad += g.x[g.y == 1].iloc[0] - m
            hess -= (p * xv ** 2).sum() - m ** 2
        if hess == 0:
            break
        step = grad / hess
        beta -= step
        if abs(step) < 1e-8:
            break
    return beta, float(np.sqrt(-1 / hess)) if hess < 0 else np.nan, len(groups)


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for r, i in enumerate(order):
        running = max(running, (len(p) - r) * p[i])
        adj[i] = min(1.0, running)
    return adj


def delong(y, s1, s2):
    """İki korelasyonlu AUROC için DeLong testi (iki yönlü p)."""
    from scipy.stats import norm
    y = np.asarray(y, bool)
    pos, neg = np.where(y)[0], np.where(~y)[0]

    def comps(s):
        sp, sn = s[pos][:, None], s[neg][None, :]
        psi = (sp > sn) + 0.5 * (sp == sn)
        return psi.mean(1), psi.mean(0), psi.mean()
    v10a, v01a, a1 = comps(np.asarray(s1, float))
    v10b, v01b, a2 = comps(np.asarray(s2, float))
    s10 = np.cov(np.vstack([v10a, v10b]))
    s01 = np.cov(np.vstack([v01a, v01b]))
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / len(pos) + (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / len(neg)
    z = (a1 - a2) / np.sqrt(var)
    return a1 - a2, float(2 * norm.sf(abs(z)))


def calibration(y, p):
    """Kalibrasyon kesişimi (eğim 1 sabitken), eğimi ve 10 eşit frekanslı kutuda ECE."""
    from sklearn.linear_model import LogisticRegression
    y = np.asarray(y, float)
    lp = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    slope = LogisticRegression(C=1e9, max_iter=1000).fit(lp[:, None], y).coef_[0, 0]
    a = 0.0
    for _ in range(50):                                     # ofsetli, yalnız kesişim
        mu = 1 / (1 + np.exp(-(lp + a)))
        a += (y - mu).sum() / (mu * (1 - mu)).sum()
    bins = np.array_split(np.argsort(p), 10)
    ece = sum(len(b) / len(p) * abs(y[b].mean() - p[b].mean()) for b in bins)
    return float(a), float(slope), float(ece)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
