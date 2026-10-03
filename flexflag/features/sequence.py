"""Sequence-only features."""

import math
from collections import Counter

GROUPS = {
    "hydrophobic": set("AVILMFWC"),
    "charged": set("DEKR"),
    "polar": set("STNQHY"),
    "gly_pro": set("GP"),
    "aromatic": set("FWY"),
}
LC_WINDOW, LC_ENTROPY = 12, 2.2  # SEG-like low-complexity: window entropy below this (bits)


def low_complexity_frac(seq: str) -> float:
    if len(seq) < LC_WINDOW:
        return 0.0
    flagged = [False] * len(seq)
    for i in range(len(seq) - LC_WINDOW + 1):
        counts = Counter(seq[i : i + LC_WINDOW])
        h = -sum(c / LC_WINDOW * math.log2(c / LC_WINDOW) for c in counts.values())
        if h < LC_ENTROPY:
            for j in range(i, i + LC_WINDOW):
                flagged[j] = True
    return sum(flagged) / len(seq)


def sequence_features(seq: str) -> dict[str, float]:
    n = len(seq)
    feats = {"seq_length": float(n), "seq_log_length": math.log(n)}
    for name, aas in GROUPS.items():
        feats[f"seq_frac_{name}"] = sum(a in aas for a in seq) / n
    feats["seq_frac_cys"] = seq.count("C") / n
    feats["seq_low_complexity_frac"] = low_complexity_frac(seq)
    return feats
