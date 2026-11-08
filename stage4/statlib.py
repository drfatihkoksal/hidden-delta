"""Hızlı küme bootstrap'ı ve AUROC (eşleştirilmiş kümeler yerine koyarak çekilir; pandas birleştirmesi yok)."""
import numpy as np
from scipy.stats import rankdata


def auroc(y, s):
    """Mann-Whitney ile AUROC (bağlar ortalama sıra; vektörel)."""
    y = np.asarray(y, bool)
    ranks = rankdata(np.asarray(s, float))
    n1, n0 = y.sum(), (~y).sum()
    return (ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def cluster_index(set_ids):
    """Küme başına satır indeksleri."""
    set_ids = np.asarray(set_ids)
    uniq, inv = np.unique(set_ids, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    bounds = np.r_[0, np.cumsum(np.bincount(inv))]
    return [order[bounds[k]:bounds[k + 1]] for k in range(len(uniq))]


def cluster_boot(set_ids, stat, n_boot=2000, seed=0):
    """stat(idx) → değer; idx bir bootstrap örneğinin satır indeksleri. Yüzdelik %95 GA döndürür."""
    groups = cluster_index(set_ids)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        vals.append(stat(np.concatenate([groups[k] for k in pick])))
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    return np.percentile(vals, [2.5, 97.5]), vals


def auc_boot(y, s, set_ids, n_boot=2000, seed=0):
    y, s = np.asarray(y), np.asarray(s)

    def stat(idx):
        yy = y[idx]
        return auroc(yy, s[idx]) if 0 < yy.sum() < len(yy) else np.nan
    ci, _ = cluster_boot(set_ids, stat, n_boot, seed)
    return auroc(y, s), ci


def auc_diff_boot(y, s1, s2, set_ids, n_boot=2000, seed=0):
    y, s1, s2 = np.asarray(y), np.asarray(s1), np.asarray(s2)

    def stat(idx):
        yy = y[idx]
        return auroc(yy, s1[idx]) - auroc(yy, s2[idx]) if 0 < yy.sum() < len(yy) else np.nan
    ci, vals = cluster_boot(set_ids, stat, n_boot, seed)
    return auroc(y, s1) - auroc(y, s2), ci, vals
