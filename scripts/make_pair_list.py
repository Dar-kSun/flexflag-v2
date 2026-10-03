"""Write data/apobind_pairs.csv (committed) from APObind's apobind_all.csv (not committed).

apobind_all.csv comes from the APObind data folder linked at
https://github.com/devalab/Apobind (manual browser download; see docs/data-choice.md).
"""

import pandas as pd

from flexflag.config import DATA_DIR, PAIRS_CSV
from flexflag.data.apoholo import PAIR_COLUMNS

full = pd.read_csv(DATA_DIR / "apobind_all.csv", index_col=0)
full[PAIR_COLUMNS].to_csv(PAIRS_CSV, index=False)
print(f"wrote {len(full)} pairs to {PAIRS_CSV.relative_to(DATA_DIR.parent)}")
