"""PAE (predicted aligned error) features, centred on inter-domain uncertainty.

A pair of rigid domains joined by a flexible hinge shows up in the PAE matrix as
low-error blocks on the diagonal (each domain is confident internally) and
high-error off-diagonal blocks (their relative placement is not). We find the
blocks by clustering confident residues on PAE, then summarise the contrast.
"""

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from flexflag.config import LOW_PLDDT

DOMAIN_PAE_CUT = 10.0  # residues closer than this (average linkage) share a domain
MIN_DOMAIN_SIZE = 30


def domains(pae: np.ndarray, plddt: np.ndarray) -> list[np.ndarray]:
    """Residue-index arrays of PAE domains (confident residues only)."""
    idx = np.flatnonzero(plddt >= LOW_PLDDT)
    if len(idx) < MIN_DOMAIN_SIZE:
        return []
    d = pae[np.ix_(idx, idx)]
    d = (d + d.T) / 2
    np.fill_diagonal(d, 0.0)
    labels = fcluster(linkage(squareform(d, checks=False), "average"), DOMAIN_PAE_CUT, "distance")
    out = [idx[labels == k] for k in np.unique(labels)]
    return [dom for dom in out if len(dom) >= MIN_DOMAIN_SIZE]


def pae_features(pae: np.ndarray, plddt: np.ndarray) -> dict[str, float]:
    sym = (pae + pae.T) / 2
    doms = domains(pae, plddt)
    feats = {
        "pae_mean": float(pae.mean()),
        "pae_frac_above15": float(np.mean(pae > 15)),
        "pae_n_domains": float(len(doms)),
        "pae_intra_mean": np.nan,
        "pae_inter_mean": np.nan,
        "pae_inter_max": np.nan,
        "pae_inter_intra_ratio": np.nan,
    }
    if doms:
        intra = np.concatenate([sym[np.ix_(d, d)].ravel() for d in doms])
        feats["pae_intra_mean"] = float(intra.mean())
    if len(doms) >= 2:
        inter = [sym[np.ix_(a, b)].mean() for i, a in enumerate(doms) for b in doms[i + 1 :]]
        feats["pae_inter_mean"] = float(np.mean(inter))
        feats["pae_inter_max"] = float(np.max(inter))
        feats["pae_inter_intra_ratio"] = feats["pae_inter_mean"] / feats["pae_intra_mean"]
    return feats
