"""Every tunable number in one place, so results can quote them exactly."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
CACHE_DIR = DATA_DIR / "cache"

# Pair list derived from APObind's apobind_all.csv (see docs/data-choice.md).
PAIRS_CSV = DATA_DIR / "apobind_pairs.csv"

# Binding site = holo-chain residues with any heavy atom within this distance
# of any ligand heavy atom.
SITE_CUTOFF_A = 5.0

# Site RMSD above this counts as "large change". Sensitivity is reported at
# the alternatives too.
LARGE_CHANGE_A = 2.0
LARGE_CHANGE_ALTERNATIVES_A = (1.0, 1.5, 3.0)

# pLDDT below this is "low confidence" (AlphaFold's own convention).
LOW_PLDDT = 70.0

HTTP_TIMEOUT_S = 60
