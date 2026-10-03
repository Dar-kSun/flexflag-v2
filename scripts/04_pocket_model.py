"""Pocket-conditioned flag, part B of docs/plan-v0.2-pocket.md.

    python scripts/04_pocket_model.py

Needs results/pocket.csv (scripts/03_pocket.py) and the v0.1 features and clusters
(scripts/02_train_eval.py). Writes results/pocket_metrics.json and
results/pocket_metrics.md.
"""

import json

import numpy as np
import pandas as pd

from flexflag.config import LARGE_CHANGE_A, LARGE_CHANGE_ALTERNATIVES_A
from flexflag.evaluate import (
    METRICS,
    N_BOOT,
    N_SPLITS,
    RESULTS,
    SEED,
    bootstrap,
    ci,
    clusters_table,
    evaluate,
    features_table,
)
from flexflag.features.pocket import POCKET_FEATURES
from flexflag.model import gradient_boosting, logistic

BASELINE = "site_plddt"
LABELS = {
    "site_plddt": "Pocket baseline: site mean pLDDT",
    "pocket_gb": "Pocket model (site + whole-protein features, boosting)",
    "site_size": "Confound reference: site size alone (not a usable model)",
}


def main() -> None:
    ds = pd.read_csv(RESULTS / "dataset.csv")
    pocket = pd.read_csv(RESULTS / "pocket.csv").set_index("uniprot")
    ds = ds[ds.uniprot.isin(pocket.dropna(subset=POCKET_FEATURES[:4]).index)].reset_index(drop=True)
    protein = features_table(ds).reset_index(drop=True)
    X = pd.concat(
        [protein, pocket.loc[ds.uniprot, POCKET_FEATURES].reset_index(drop=True)], axis=1
    ).assign(n_site=ds.n_site.to_numpy())
    groups = clusters_table(ds).to_numpy()
    rmsd = ds.site_rmsd.to_numpy()
    feature_cols = list(protein.columns) + POCKET_FEATURES  # n_site deliberately excluded
    models = {
        "site_plddt": (logistic, ["site_plddt_mean"]),
        "pocket_gb": (gradient_boosting, feature_cols),
        "site_size": (logistic, ["n_site"]),
    }
    rng = np.random.default_rng(SEED)
    out = {"n_proteins": len(ds), "n_clusters": int(len(np.unique(groups))), "thresholds": {}}
    lines = [
        f"Pocket-conditioned mode. n = {len(ds)} proteins, {out['n_clusters']} clusters, "
        f"{N_SPLITS}-fold cluster CV, {N_BOOT} cluster bootstraps. Pre-declared in "
        "docs/plan-v0.2-pocket.md.",
        "",
    ]
    for thr in sorted({LARGE_CHANGE_A, *LARGE_CHANGE_ALTERNATIVES_A}):
        y = (rmsd > thr).astype(int)
        preds = evaluate(X, y, groups, models)
        draws, diffs = bootstrap(y, preds, groups, rng, BASELINE)
        res = {"n_positive": int(y.sum()), "base_rate": float(y.mean()), "models": {}}
        lines += [
            f"### Site RMSD > {thr:g} Å ({y.sum()} positives, {y.mean():.1%})",
            "",
            "| Model | AUROC | AUPRC | ECE | Brier |",
            "|---|---|---|---|---|",
        ]
        for k, p in preds.items():
            r = {m: {"value": float(f(y, p)), "ci95": ci(draws[m][k])} for m, f in METRICS.items()}
            if k != BASELINE:
                r["minus_baseline"] = {
                    m: {"mean": float(np.mean(diffs[m][k])), "ci95": ci(diffs[m][k])} for m in diffs
                }
            res["models"][k] = r
            cells = [f"{r[m]['value']:.3f} [{r[m]['ci95'][0]:.3f}, {r[m]['ci95'][1]:.3f}]"
                     for m in METRICS]  # fmt: skip
            lines.append(f"| {LABELS[k]} | " + " | ".join(cells) + " |")
        d = res["models"]["pocket_gb"]["minus_baseline"]["auroc"]
        lines += ["", f"Pocket model − pocket baseline, AUROC: {d['mean']:+.3f} "
                  f"[{d['ci95'][0]:+.3f}, {d['ci95'][1]:+.3f}]", ""]  # fmt: skip
        out["thresholds"][f"{thr:g}"] = res
        print(f"threshold {thr:g} A done", flush=True)

    (RESULTS / "pocket_metrics.json").write_text(json.dumps(out, indent=2))
    (RESULTS / "pocket_metrics.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
