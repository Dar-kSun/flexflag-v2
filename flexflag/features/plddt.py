"""pLDDT summary features (AlphaFold DB per-residue confidence)."""

import numpy as np

from flexflag.config import LOW_PLDDT


def mean_plddt(scores: np.ndarray) -> float:
    return float(np.mean(scores))


def frac_low_plddt(scores: np.ndarray, threshold: float = LOW_PLDDT) -> float:
    return float(np.mean(scores < threshold))


def longest_run(mask: np.ndarray) -> int:
    best = run = 0
    for m in mask:
        run = run + 1 if m else 0
        best = max(best, run)
    return best


def plddt_features(scores: np.ndarray) -> dict[str, float]:
    low = scores < LOW_PLDDT
    return {
        "plddt_mean": mean_plddt(scores),
        "plddt_min": float(scores.min()),
        "plddt_std": float(scores.std()),
        "plddt_frac_below70": frac_low_plddt(scores),
        # pLDDT < 50 is AlphaFold's own proxy for intrinsic disorder.
        "plddt_frac_below50": float(np.mean(scores < 50)),
        "plddt_longest_low_run": float(longest_run(low)),
        "plddt_longest_low_run_frac": longest_run(low) / len(scores),
    }
