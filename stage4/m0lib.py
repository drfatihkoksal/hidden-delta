"""M0 / M-demo / M2 için ortak tasarım: 3 düğümlü kısıtlı kübik spline (Harrell), standartlaştırma, L2 lojistik.

Eğitimde öğrenilen her şey (düğümler, ölçek, katsayılar) bir sözlükte saklanır; test verisine aynen uygulanır.
"""
import numpy as np
from sklearn.linear_model import LogisticRegressionCV

CONT = ["q_I", "q_aVL", "q_V5", "q_V6", "v1_rs_c", "qrs_ms", "slow_init_ms", "svr10",
        "v20_x", "v20_y", "v20_z", "v40_x", "v40_y", "v40_z", "i20_x", "i20_y", "i20_z",
        "pqlow_mag", "pqlow_fq_corr", "qrst_angle", "t_int_x", "t_int_y", "t_int_z", "pr_ms", "hr", "age"]
BIN = ["septal_q_absent", "v1_rsr", "wave", "male"]
MODELS = {"M0": (CONT, BIN), "M-demo": (["age", "hr"], ["male"])}
SEED = 20261118


def rcs(x, knots):
    k1, k2, k3 = knots
    p = lambda u: np.maximum(u, 0) ** 3
    nl = p(x - k1) - p(x - k2) * (k3 - k1) / (k3 - k2) + p(x - k3) * (k2 - k1) / (k3 - k2)
    return np.c_[x, nl / (k3 - k1) ** 2]


def fit_design(train, cont, binv):
    knots = {}
    for c in cont:
        kn = np.nanpercentile(train[c], [10, 50, 90])
        knots[c] = [float(v) for v in kn] if len(np.unique(kn)) == 3 else None
    params = {"cont": cont, "bin": binv, "knots": knots}
    raw = raw_design(train, params)
    params["mean"], params["scale"] = raw.mean(0).tolist(), np.where(raw.std(0) > 0, raw.std(0), 1.0).tolist()
    return params


def raw_design(df, params):
    cols = []
    for c in params["cont"]:
        kn = params["knots"][c]
        x = df[c].to_numpy(float)
        cols.append(x[:, None] if kn is None else rcs(x, kn))
    cols.append(df[params["bin"]].to_numpy(float))
    return np.hstack(cols)


def design(df, params):
    return (raw_design(df, params) - np.array(params["mean"])) / np.array(params["scale"])


def fit_lr(X, y):
    return LogisticRegressionCV(Cs=10, cv=5, solver="lbfgs", scoring="roc_auc", max_iter=5000,
                                random_state=SEED, n_jobs=-1).fit(X, y)          # n_jobs sonucu değiştirmez


def lr_params(m):
    return {"coef": m.coef_.ravel().tolist(), "intercept": float(m.intercept_[0]), "C": float(m.C_[0])}


def lr_logit(X, p):
    return X @ np.array(p["coef"]) + p["intercept"]
