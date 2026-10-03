"""Export the frozen pocket-pLDDT rule used by `flexflag check` to flexflag/rule.json.

    python scripts/07_export_rule.py

The rule is the logistic regression on site mean pLDDT fit on all APObind
proteins (identical to the frozen model scored in scripts/06_external_eval.py),
plus the APObind band edges and the observed rates in both datasets.
"""

import json

import numpy as np
import pandas as pd

from flexflag.config import LARGE_CHANGE_A, ROOT
from flexflag.evaluate import RESULTS
from flexflag.model import logistic

OUT = ROOT / "flexflag" / "rule.json"


def main() -> None:
    ds = pd.read_csv(RESULTS / "dataset.csv")
    pk = pd.read_csv(RESULTS / "pocket.csv").set_index("uniprot")
    ds = ds[ds.uniprot.isin(pk.dropna(subset=["site_plddt_mean"]).index)]
    x = pk.loc[ds.uniprot, ["site_plddt_mean"]].reset_index(drop=True)
    y = (ds.site_rmsd.to_numpy() > LARGE_CHANGE_A).astype(int)
    pipe = logistic().fit(x, y)
    scaler, lr = pipe[-2], pipe[-1]
    # Fold the standardisation into the coefficients: p = sigmoid(a + b * site_plddt).
    b = float(lr.coef_[0, 0] / scaler.scale_[0])
    a = float(lr.intercept_[0] - lr.coef_[0, 0] * scaler.mean_[0] / scaler.scale_[0])
    check = 1 / (1 + np.exp(-(a + b * x.site_plddt_mean.to_numpy())))
    assert np.allclose(check, pipe.predict_proba(x)[:, 1])

    ext = json.loads((RESULTS / "external" / "metrics.json").read_text(encoding="utf-8"))
    ext2 = ext["thresholds"][f"{LARGE_CHANGE_A:g}"]["all"]["models"]
    rule = {
        "description": "P(binding-site C-alpha RMSD > 2 A between apo and holo) from the mean "
        "AlphaFold pLDDT over the pocket residues. Fit on APObind; validated on PDB-2019+.",
        "threshold_A": LARGE_CHANGE_A,
        "intercept": a,
        "coef_site_plddt": b,
        "band_edges": ext["band_edges"],
        "bands": ext["bands"],
        "n_apobind": ext["n_train"],
        "n_external": ext["n_external"],
        "external_auroc_site_plddt": ext2["site_plddt"]["auroc"],
        "external_auroc_protein_plddt": ext2["protein_plddt"]["auroc"],
        "alphafold_state_moving_external": ext["alphafold_state"]["external"]["moving"],
        "alphafold_state_moving_apobind": ext["alphafold_state"]["apobind"]["moving"],
    }
    OUT.write_text(json.dumps(rule, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: p = sigmoid({a:.3f} + {b:.4f} * site_plddt)")


if __name__ == "__main__":
    main()
