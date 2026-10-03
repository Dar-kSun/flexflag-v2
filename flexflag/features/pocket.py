"""Pocket-conditioned features: AlphaFold DB confidence at given pocket residues.

This is a separate mode from the whole-protein flag. The pocket's residue
*positions* are an input (as they are at docking time); everything computed from
them comes from AlphaFold DB only. See docs/plan-v0.2-pocket.md.
"""

import numpy as np

from flexflag.config import LOW_PLDDT
from flexflag.features.pae import domains

POCKET_FEATURES = [
    "site_plddt_mean",
    "site_plddt_min",
    "site_pae_within",
    "site_pae_to_rest",
    "site_n_domains",
    "site_frac_low_plddt",
]


def pocket_features(site: np.ndarray, plddt: np.ndarray, pae: np.ndarray) -> dict[str, float]:
    """`site`: 0-based indices into the AlphaFold sequence."""
    site = np.asarray(sorted(set(site)))
    sym = (pae + pae.T) / 2
    rest = np.setdiff1d(np.flatnonzero(plddt >= LOW_PLDDT), site)
    within = sym[np.ix_(site, site)][np.triu_indices(len(site), k=1)]
    in_domain = [np.isin(site, d).sum() for d in domains(pae, plddt)]
    return {
        "site_plddt_mean": float(plddt[site].mean()),
        "site_plddt_min": float(plddt[site].min()),
        "site_pae_within": float(within.mean()),
        "site_pae_to_rest": float(sym[np.ix_(site, rest)].mean()) if len(rest) else np.nan,
        "site_n_domains": float(sum(n >= 2 for n in in_domain)),
        "site_frac_low_plddt": float(np.mean(plddt[site] < LOW_PLDDT)),
    }
