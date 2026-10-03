import numpy as np
import pytest

from flexflag.splits import LeakageError, assert_no_leakage, cluster_folds


def test_folds_never_share_a_cluster():
    rng = np.random.default_rng(0)
    groups = rng.integers(0, 40, size=300)
    y = rng.integers(0, 2, size=300)
    seen = []
    for train, test in cluster_folds(y, groups):
        assert not set(groups[train]) & set(groups[test])
        seen.extend(test)
    assert sorted(seen) == list(range(300))  # every protein tested exactly once


def test_leakage_assertion_fires():
    groups = np.array(["a", "a", "b", "c"])
    with pytest.raises(LeakageError):
        assert_no_leakage(groups, np.array([0, 2]), np.array([1, 3]))
