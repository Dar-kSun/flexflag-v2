"""Cluster-aware cross-validation splits, with a hard leakage check.

Proteins are clustered with MMseqs2 (easy-cluster) at CLUSTER_IDENTITY sequence
identity; every member of a cluster lands in the same fold.
"""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

CLUSTER_IDENTITY = 0.3
CLUSTER_COVERAGE = 0.5
# On Windows the MMseqs2 build needs a path without spaces and its BusyBox helpers
# installed (see README). Override with FLEXFLAG_MMSEQS.
MMSEQS_WINDOWS = Path(os.environ.get("LOCALAPPDATA", "")) / "flexflag/mmseqs/bin/mmseqs.exe"


class LeakageError(AssertionError):
    """A sequence cluster appears in both train and test."""


def mmseqs_binary() -> str:
    if os.environ.get("FLEXFLAG_MMSEQS"):
        return os.environ["FLEXFLAG_MMSEQS"]
    found = shutil.which("mmseqs")
    if found:
        return found
    if MMSEQS_WINDOWS.exists():
        return str(MMSEQS_WINDOWS)
    raise FileNotFoundError("MMseqs2 not found; see README for install")


def cluster(ids: list[str], seqs: list[str]) -> pd.Series:
    """Cluster representative accession for each id (indexed by id)."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / "in.fasta"
        # Unix newlines and relative paths: the Windows (Cygwin) build of MMseqs2
        # crashes on CRLF FASTA and mishandles C:\ paths.
        text = "".join(f">{i}\n{s}\n" for i, s in zip(ids, seqs, strict=True))
        fasta.write_text(text, newline="\n")
        subprocess.run(
            [mmseqs_binary(), "easy-cluster", "in.fasta", "out", "work",
             "--min-seq-id", str(CLUSTER_IDENTITY), "-c", str(CLUSTER_COVERAGE),
             "--cov-mode", "0", "-s", "7.5", "--threads", "4", "-v", "1"],
            check=True,
            cwd=tmp,
        )  # fmt: skip
        tsv = pd.read_csv(tmp / "out_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    clusters = tsv.set_index("member").rep
    missing = set(ids) - set(clusters.index)
    if missing:
        raise RuntimeError(f"MMseqs2 dropped {len(missing)} sequences")
    return clusters.loc[ids]


def assert_no_leakage(groups: np.ndarray, train_idx: np.ndarray, test_idx: np.ndarray) -> None:
    shared = set(groups[train_idx]) & set(groups[test_idx])
    if shared:
        raise LeakageError(
            f"{len(shared)} clusters in both train and test, e.g. {next(iter(shared))}"
        )


def cluster_folds(y: np.ndarray, groups: np.ndarray, n_splits: int = 5, seed: int = 0):
    """Yield (train_idx, test_idx), stratified on y, grouped by cluster, leakage-checked."""
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train_idx, test_idx in cv.split(np.zeros(len(y)), y, groups):
        assert_no_leakage(groups, train_idx, test_idx)
        yield train_idx, test_idx
