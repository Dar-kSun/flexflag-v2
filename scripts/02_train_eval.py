"""Cluster-aware cross-validated comparison: pLDDT-only baselines vs the full model.

    python scripts/02_train_eval.py

Reads results/dataset.csv. Writes results/features.csv, results/clusters.csv,
results/predictions.csv, results/metrics.json, results/metrics.md and three figures.
"""

import json
from concurrent.futures import ThreadPoolExecutor

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from flexflag.config import LARGE_CHANGE_A, LARGE_CHANGE_ALTERNATIVES_A, ROOT
from flexflag.features import compute_features
from flexflag.model import gradient_boosting, logistic
from flexflag.splits import CLUSTER_IDENTITY, cluster, cluster_folds

matplotlib.use("Agg")

RESULTS = ROOT / "results"
N_BOOT = 1000
N_SPLITS = 5
MIN_POSITIVES = 25
SEED = 0
BASELINE = "plddt_mean"

# Validated default chart palette (dataviz reference instance), light surface.
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e8e7e3"
BLUE, ORANGE = "#2a78d6", "#eb6834"

LABELS = {
    "plddt_mean": "Baseline: mean pLDDT",
    "plddt_frac_below70": "Baseline: fraction pLDDT < 70",
    "logistic_all": "Logistic regression, all features",
    "full_gb": "Full model (gradient boosting)",
}


def features_table(ds: pd.DataFrame) -> pd.DataFrame:
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


def bootstrap(y, preds: dict, groups, rng):
    """Cluster bootstrap: resample whole clusters with replacement."""
    uniq = np.unique(groups)
    members = {g: np.flatnonzero(groups == g) for g in uniq}
    draws = {m: {k: [] for k in preds} for m in METRICS}
    diffs = {m: {k: [] for k in preds if k != BASELINE} for m in ("auroc", "auprc")}
    for _ in range(N_BOOT):
        idx = np.concatenate([members[g] for g in rng.choice(uniq, size=len(uniq))])
        if y[idx].min() == y[idx].max():
            continue
        vals = {m: {k: f(y[idx], p[idx]) for k, p in preds.items()} for m, f in METRICS.items()}
        for m in METRICS:
            for k in preds:
                draws[m][k].append(vals[m][k])
        for m in diffs:
            for k in diffs[m]:
                diffs[m][k].append(vals[m][k] - vals[m][BASELINE])
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


def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def reliability_figure(y, preds: dict, colors: dict, thr: float, path, n_bins: int = 5):
    """Equal-count bins per model, so a model that predicts a narrow range stays visible."""
    fig, ax = plt.subplots(figsize=(5.6, 4.8), facecolor=SURFACE)
    _style(ax)
    curves = {}
    for name, p in preds.items():
        chunks = np.array_split(np.argsort(p), n_bins)
        curves[name] = ([p[c].mean() for c in chunks], [y[c].mean() for c in chunks])
    hi = max(0.2, 1.25 * max(max(xs + ys) for xs, ys in curves.values()))
    ax.plot([0, hi], [0, hi], color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.text(hi * 0.97, hi * 0.9, "perfect calibration", color=MUTED, fontsize=8, ha="right")
    for name, p in preds.items():
        xs, ys = curves[name]
        span = f"predicts {p.min():.2f}–{p.max():.2f}"
        ax.plot(xs, ys, color=colors[name], lw=2, zorder=2)
        ax.scatter(xs, ys, s=64, color=colors[name], edgecolor=SURFACE, linewidth=2,
                   zorder=3, label=f"{LABELS[name]} ({span})")  # fmt: skip
    ax.axhline(y.mean(), color=GRID, lw=1, zorder=0)
    ax.text(hi * 0.99, y.mean() + 0.005, f"base rate {y.mean():.1%}", color=MUTED,
            fontsize=8, ha="right", va="bottom")  # fmt: skip
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xlabel("Predicted probability of large change", color=INK2)
    ax.set_ylabel("Observed fraction with large change", color=INK2)
    ax.set_title(f"Calibration, site RMSD > {thr:g} Å (out-of-fold)", color=INK, loc="left")
    ax.legend(frameon=False, labelcolor=INK2, loc="upper left", fontsize=8)
    fig.text(0.01, 0.01, f"{n_bins} equal-count bins per model ({len(y) // n_bins} proteins each)",
             color=MUTED, fontsize=7)  # fmt: skip
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def comparison_figure(out: dict, path):
    """AUROC with 95% cluster-bootstrap CI, baseline vs full model, at every threshold."""
    rows = [(t, r) for t, r in out["thresholds"].items() if "models" in r]
    fig, ax = plt.subplots(figsize=(7.2, 0.9 + 0.75 * len(rows)), facecolor=SURFACE)
    _style(ax)
    ax.grid(axis="y", visible=False)
    ax.axvline(0.5, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax.text(0.497, len(rows) - 0.45, "chance", color=MUTED, fontsize=8, ha="right")
    for i, (_thr, res) in enumerate(reversed(rows)):
        for name, color, dy in ((BASELINE, ORANGE, 0.14), ("full_gb", BLUE, -0.14)):
            m = res["models"][name]["auroc"]
            ax.plot(m["ci95"], [i + dy] * 2, color=color, lw=2, solid_capstyle="round")
            ax.scatter([m["value"]], [i + dy], s=64, color=color, edgecolor=SURFACE,
                       linewidth=2, zorder=3)  # fmt: skip
            ax.text(m["ci95"][1] + 0.008, i + dy, f"{m['value']:.2f}", color=INK2,
                    fontsize=8, va="center")  # fmt: skip
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(
        [f"> {t} Å  ({r['base_rate']:.0%} positive)" for t, r in reversed(rows)], color=INK2
    )
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_xlim(0.35, 0.8)
    ax.set_xlabel("AUROC, out-of-fold (95% cluster-bootstrap CI)", color=INK2)
    ax.set_title(f"pLDDT alone vs full model (n = {out['n_proteins']} proteins)",
                 color=INK, loc="left")  # fmt: skip
    ax.scatter([], [], s=64, color=ORANGE, label=LABELS[BASELINE])
    ax.scatter([], [], s=64, color=BLUE, label=LABELS["full_gb"])
    ax.legend(frameon=False, labelcolor=INK2, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def distribution_figure(rmsd: np.ndarray, thresholds, path):
    fig, ax = plt.subplots(figsize=(7.2, 3.8), facecolor=SURFACE)
    _style(ax)
    hi = 8.0
    bins = np.arange(0, hi + 0.25, 0.25)
    ax.hist(np.clip(rmsd, 0, hi - 1e-9), bins=bins, color=BLUE, edgecolor=SURFACE, linewidth=1)
    top = ax.get_ylim()[1]
    for t in thresholds:
        main = t == LARGE_CHANGE_A
        ax.axvline(t, color=INK if main else MUTED, lw=1, ls="-" if main else (0, (3, 3)))
        ax.text(t + 0.07, top * 0.97, f"> {t:g} Å: {np.mean(rmsd > t):.0%}",
                color=INK if main else INK2, fontsize=8, va="top", rotation=90)  # fmt: skip
    ax.set_xlabel(f"Binding-site Cα RMSD, apo vs holo (Å; last bin is ≥ {hi - 0.25:g})",
                  color=INK2)  # fmt: skip
    ax.set_ylabel("Proteins", color=INK2)
    ax.set_title(f"How much the binding site changes shape (n = {len(rmsd)} proteins)",
                 color=INK, loc="left")  # fmt: skip
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def univariate_table(X: pd.DataFrame, rmsd: np.ndarray) -> pd.DataFrame:
    """In-sample AUROC of each feature on its own (0.5 = no signal, < 0.5 = inverse).

    Descriptive only: shows which way each feature points, not out-of-fold skill.
    """
    rows = {}
    for col in X.columns:
        v = X[col].to_numpy()
        ok = ~np.isnan(v)
        rows[col] = {"non_missing": ok.mean()}
        for thr in (1.0, LARGE_CHANGE_A, 3.0):
            y = rmsd[ok] > thr
            rows[col][f"auroc_{thr:g}A"] = roc_auc_score(y, v[ok])
    df = pd.DataFrame(rows).T
    return df.reindex(
        (df[f"auroc_{LARGE_CHANGE_A:g}A"] - 0.5).abs().sort_values(ascending=False).index
    )


def markdown(out: dict) -> str:
    def fmt(d):
        return f"{d['value']:.3f} [{d['ci95'][0]:.3f}, {d['ci95'][1]:.3f}]"

    lines = [
        f"n = {out['n_proteins']} proteins in {out['n_clusters']} sequence clusters "
        f"(MMseqs2, {out['cluster_identity']:.0%} identity). {out['cv_folds']}-fold "
        f"cluster-grouped CV; 95% CIs from {out['bootstrap']} cluster bootstraps.",
        "",
    ]
    for thr, res in out["thresholds"].items():
        lines.append(f"### Site RMSD > {thr} Å")
        if "skipped" in res:
            lines += [f"Skipped: {res['skipped']}.", ""]
            continue
        lines += [
            f"{res['n_positive']} positives ({res['base_rate']:.1%} base rate).",
            "",
            "| Model | AUROC | AUPRC | ECE | Brier |",
            "|---|---|---|---|---|",
        ]
        for k, r in res["models"].items():
            lines.append(f"| {LABELS[k]} | " + " | ".join(fmt(r[m]) for m in METRICS) + " |")
        lines += ["", "Difference from the mean-pLDDT baseline (paired cluster bootstrap):", ""]
        for k, r in res["models"].items():
            if "minus_baseline" in r:
                d = r["minus_baseline"]
                lines.append(
                    f"- {LABELS[k]}: AUROC {d['auroc']['mean']:+.3f} "
                    f"[{d['auroc']['ci95'][0]:+.3f}, {d['auroc']['ci95'][1]:+.3f}], "
                    f"AUPRC {d['auprc']['mean']:+.3f} "
                    f"[{d['auprc']['ci95'][0]:+.3f}, {d['auprc']['ci95'][1]:+.3f}]"
                )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    ds = pd.read_csv(RESULTS / "dataset.csv")
    X = features_table(ds).reset_index(drop=True)
    groups = clusters_table(ds).to_numpy()
    rmsd = ds.site_rmsd.to_numpy()
    cols = list(X.columns)
    models = {
        "plddt_mean": (logistic, ["plddt_mean"]),
        "plddt_frac_below70": (logistic, ["plddt_frac_below70"]),
        "logistic_all": (logistic, cols),
        "full_gb": (gradient_boosting, cols),
    }
    thresholds = sorted({LARGE_CHANGE_A, *LARGE_CHANGE_ALTERNATIVES_A})
    rng = np.random.default_rng(SEED)
    out = {
        "n_proteins": len(ds),
        "n_clusters": int(len(np.unique(groups))),
        "cluster_identity": CLUSTER_IDENTITY,
        "cv_folds": N_SPLITS,
        "bootstrap": N_BOOT,
        "seed": SEED,
        "features": cols,
        "thresholds": {},
    }
    pred_table = ds[["uniprot", "site_rmsd"]].assign(cluster=groups)
    for thr in thresholds:
        y = (rmsd > thr).astype(int)
        if min(y.sum(), (1 - y).sum()) < MIN_POSITIVES:
            out["thresholds"][f"{thr:g}"] = {"skipped": f"fewer than {MIN_POSITIVES} per class"}
            continue
        preds = evaluate(X, y, groups, models)
        draws, diffs = bootstrap(y, preds, groups, rng)
        res = {"n_positive": int(y.sum()), "base_rate": float(y.mean()), "models": {}}
        for k, p in preds.items():
            res["models"][k] = {
                m: {"value": float(f(y, p)), "ci95": ci(draws[m][k])} for m, f in METRICS.items()
            }
            if k != BASELINE:
                res["models"][k]["minus_baseline"] = {
                    m: {
                        "mean": float(np.mean(diffs[m][k])),
                        "ci95": ci(diffs[m][k]),
                        "frac_draws_le_0": float(np.mean(np.array(diffs[m][k]) <= 0)),
                    }
                    for m in diffs
                }
            pred_table[f"p_{k}_{thr:g}A"] = p
        out["thresholds"][f"{thr:g}"] = res
        if thr == LARGE_CHANGE_A:
            reliability_figure(
                y,
                {k: preds[k] for k in (BASELINE, "full_gb")},
                {BASELINE: ORANGE, "full_gb": BLUE},
                thr,
                RESULTS / "calibration.png",
            )
        print(f"threshold {thr:g} A: {y.sum()} positives / {len(y)}", flush=True)

    distribution_figure(rmsd, thresholds, RESULTS / "site_rmsd_distribution.png")
    comparison_figure(out, RESULTS / "baseline_vs_model.png")
    pred_table.to_csv(RESULTS / "predictions.csv", index=False)
    univariate_table(X, rmsd).to_csv(RESULTS / "univariate_auroc.csv", float_format="%.3f")
    (RESULTS / "metrics.json").write_text(json.dumps(out, indent=2))
    (RESULTS / "metrics.md").write_text(markdown(out), encoding="utf-8")
    print(markdown(out))


if __name__ == "__main__":
    main()
