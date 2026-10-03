"""pLDDT summary features."""

import numpy as np

from flexflag.config import LOW_PLDDT


def mean_plddt(scores: np.ndarray) -> float:
    return float(np.mean(scores))


def frac_low_plddt(scores: np.ndarray, threshold: float = LOW_PLDDT) -> float:
    return float(np.mean(scores < threshold))
