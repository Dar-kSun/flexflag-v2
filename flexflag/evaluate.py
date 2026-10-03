"""Shared evaluation: cluster-grouped CV, cluster bootstrap, metrics.

Used by scripts/02_train_eval.py (whole-protein flag) and scripts/04_pocket_model.py
(pocket-conditioned flag) so both are scored identically.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from flexflag.config import ROOT
from flexflag.features import compute_features
from flexflag.splits import cluster, cluster_folds

RESULTS = ROOT / "results"
N_BOOT = 1000
N_SPLITS = 5
SEED = 0


def features_table(ds: pd.DataFrame) -> pd.DataFrame:
    """Whole-protein features per protein, cached in results/features.csv."""
    path = RESULTS / "features.csv"
    if path.exists():
        cached = pd.read_csv(path, index_col="uniprot")
        if set(ds.uniprot) <= set(cached.index):
            return cached.loc[ds.uniprot]
    with ThreadPoolExecutor(8) as pool:
        rows = list(pool.map(compute_features, ds.uniprot))
    feats = pd.DataFrame(rows, index=pd.Index(ds.uniprot, name="uniprot"))
    feats.to_csv(path)
    return feats


def clusters_table(ds: pd.DataFrame) -> pd.Series:
    """MMseqs2 cluster per protein, cached in results/clusters.csv."""
    path = RESULTS / "clusters.csv"
    if path.exists():
        cached = pd.read_csv(path, index_col="uniprot").cluster
        if set(ds.uniprot) <= set(cached.index):
            return cached.loc[ds.uniprot]
    c = cluster(list(ds.uniprot), list(ds.sequence))
    c.rename("cluster").rename_axis("uniprot").to_frame().to_csv(path)
    return c


def ece(y, p, bins=10) -> float:
    """Expected calibration error, equal-width bins, weighted by bin size."""
    idx = np.clip(np.digitize(p, np.linspace(0, 1, bins + 1)[1:-1]), 0, bins - 1)
    total = 0.0
    for b in range(bins):
        sel = idx == b
        if sel.any():
            total += sel.mean() * abs(y[sel].mean() - p[sel].mean())
    return float(total)


METRICS = {
    "auroc": roc_auc_score,
    "auprc": average_precision_score,
    "ece": ece,
    "brier": brier_score_loss,
}


def bootstrap(y, preds: dict, groups, rng, baseline: str, n_boot: int = N_BOOT):
    """Cluster bootstrap: resample whole clusters with replacement."""
    uniq = np.unique(groups)
    members = {g: np.flatnonzero(groups == g) for g in uniq}
    draws = {m: {k: [] for k in preds} for m in METRICS}
    diffs = {m: {k: [] for k in preds if k != baseline} for m in ("auroc", "auprc")}
    for _ in range(n_boot):
        idx = np.concatenate([members[g] for g in rng.choice(uniq, size=len(uniq))])
        if y[idx].min() == y[idx].max():
            continue
        vals = {m: {k: f(y[idx], p[idx]) for k, p in preds.items()} for m, f in METRICS.items()}
        for m in METRICS:
            for k in preds:
                draws[m][k].append(vals[m][k])
        for m in diffs:
            for k in diffs[m]:
                diffs[m][k].append(vals[m][k] - vals[m][baseline])
    return draws, diffs


def ci(xs) -> list[float]:
    return [float(np.percentile(xs, 2.5)), float(np.percentile(xs, 97.5))]


def evaluate(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray, models: dict) -> dict:
    """Out-of-fold predicted probabilities for every model, same folds for all."""
    preds = {k: np.zeros(len(y)) for k in models}
    for train, test in cluster_folds(y, groups, N_SPLITS, SEED):
        for name, (make, cols) in models.items():
            m = make()
            m.fit(X.iloc[train][cols], y[train])
            preds[name][test] = m.predict_proba(X.iloc[test][cols])[:, 1]
    return preds
