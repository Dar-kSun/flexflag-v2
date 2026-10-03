"""E1 of docs/plan-v0.5-overnight.md, step 2: do ESM-2 pocket features add to pocket pLDDT?

    python scripts/12_esm_eval.py

Needs scripts/11_esm_embed.py to have run. Writes results/esm/metrics.{json,md}.
"""

import json

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from flexflag.config import CACHE_DIR, LARGE_CHANGE_A
from flexflag.data import alphafold
from flexflag.evaluate import METRICS, RESULTS, SEED, bootstrap, ci, evaluate

OUT = RESULTS / "esm"
EMB = [f"esm_{i}" for i in range(1280)]
LABELS = {
    "plddt": "Pocket pLDDT (logistic)",
    "esm": "ESM-2 pocket embedding (PCA 32, logistic)",
    "esm_plddt": "ESM-2 pocket embedding + pocket pLDDT",
}


def esm_model():
    return make_pipeline(
        StandardScaler(), PCA(32, random_state=SEED), LogisticRegression(C=0.1, max_iter=5000)
    )


def esm_plddt_model():
    pre = ColumnTransformer(
        [
            ("esm", make_pipeline(StandardScaler(), PCA(32, random_state=SEED)), EMB),
            ("plddt", StandardScaler(), ["pocket_plddt"]),
        ]
    )
    return make_pipeline(pre, LogisticRegression(C=0.1, max_iter=5000))


def plddt_model():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))


MODELS = {
    "plddt": (plddt_model, ["pocket_plddt"]),
    "esm": (esm_model, EMB),
    "esm_plddt": (esm_plddt_model, EMB + ["pocket_plddt"]),
}


def table(name: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    path = RESULTS / ("dataset.csv" if name == "apobind" else "external/dataset.csv")
    cl = RESULTS / ("clusters.csv" if name == "apobind" else "external/clusters.csv")
    ds = pd.read_csv(path)
    pockets = json.loads((RESULTS / "robustness" / "pockets.json").read_text())
    rows = []
    for u, rmsd in zip(ds.uniprot, ds.site_rmsd, strict=True):
        pocket = pockets.get(f"{name}|{u}")
        f = CACHE_DIR / "esm" / f"{u}.npy"
        if not pocket or len(pocket) < 3 or not f.exists():
            continue
        idx = np.array(pocket) - 1
        emb = np.load(f).astype(np.float32)[idx].mean(axis=0)
        rows.append(
            {
                "uniprot": u,
                "site_rmsd": rmsd,
                "pocket_plddt": float(alphafold.plddt(u)[idx].mean()),
                **dict(zip(EMB, emb, strict=True)),
            }
        )
    X = pd.DataFrame(rows)
    groups = pd.read_csv(cl, index_col="uniprot").cluster.loc[X.uniprot].to_numpy()
    y = (X.site_rmsd.to_numpy() > LARGE_CHANGE_A).astype(int)
    return X, y, groups


def score(y, preds, groups, rng) -> dict:
    draws, diffs = bootstrap(y, preds, groups, rng, "plddt")
    out = {}
    for k, p in preds.items():
        out[k] = {m: {"value": float(f(y, p)), "ci95": ci(draws[m][k])} for m, f in METRICS.items()}
        if k != "plddt":
            out[k]["minus_plddt"] = {
                m: {"mean": float(np.mean(diffs[m][k])), "ci95": ci(diffs[m][k])} for m in diffs
            }
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    Xa, ya, ga = table("apobind")
    Xe, ye, ge = table("external")
    res = {"n_apobind": len(Xa), "n_external": len(Xe)}
    res["apobind_cv"] = score(ya, evaluate(Xa, ya, ga, MODELS), ga, rng)
    frozen = {
        k: make().fit(Xa[cols], ya).predict_proba(Xe[cols])[:, 1]
        for k, (make, cols) in MODELS.items()
    }
    res["external_frozen"] = score(ye, frozen, ge, rng)
    lines = [
        f"ESM-2 pocket features (E1). APObind n = {len(Xa)} (5-fold cluster CV); "
        f"external n = {len(Xe)} (frozen). > {LARGE_CHANGE_A:g} Å.",
        "",
    ]
    for part in ("apobind_cv", "external_frozen"):
        lines += [f"### {part}", "", "| Model | AUROC | AUPRC |", "|---|---|---|"]
        for k, r in res[part].items():
            a, p = r["auroc"], r["auprc"]
            lines.append(
                f"| {LABELS[k]} | {a['value']:.3f} [{a['ci95'][0]:.3f}, {a['ci95'][1]:.3f}] "
                f"| {p['value']:.3f} |"
            )
        for k, r in res[part].items():
            if "minus_plddt" in r:
                d = r["minus_plddt"]["auroc"]
                lines.append(
                    f"- {LABELS[k]} − pocket pLDDT, AUROC {d['mean']:+.3f} "
                    f"[{d['ci95'][0]:+.3f}, {d['ci95'][1]:+.3f}]"
                )
        lines.append("")
    (OUT / "metrics.json").write_text(json.dumps(res, indent=2))
    (OUT / "metrics.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
