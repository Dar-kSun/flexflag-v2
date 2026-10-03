"""Pocket-level analyses pre-declared in docs/plan-v0.2-pocket.md (parts A and C,
plus the pocket features for part B).

    python scripts/03_pocket.py

For every protein in results/dataset.csv: re-derive the binding site (and check the
label reproduces), map it onto the AlphaFold DB model, and compute
  - site C-alpha RMSD of the AlphaFold model vs holo and vs apo (part A),
  - pocket features from AlphaFold confidence at those residues (part B).
Writes results/pocket.csv and prints the part A and C summaries.
"""

from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from flexflag.config import LARGE_CHANGE_A, ROOT
from flexflag.pocket_analysis import AF2_PDB_CUTOFF, analyse

RESULTS = ROOT / "results"


def main() -> None:
    ds = pd.read_csv(RESULTS / "dataset.csv")
    with ThreadPoolExecutor(16) as pool:
        rows = list(pool.map(analyse, ds.itertuples()))
    pk = pd.DataFrame(rows)
    pk.to_csv(RESULTS / "pocket.csv", index=False, float_format="%.4f", lineterminator="\n")

    d = ds.merge(pk, on="uniprot")
    ok = d.dropna(subset=["rmsd_af_holo"])
    moving = ok[ok.site_rmsd > LARGE_CHANGE_A]
    print(f"site mapped onto AlphaFold for {len(ok)}/{len(d)} proteins")
    post = ok.holo_release > AF2_PDB_CUTOFF
    both = post & (ok.apo_release > AF2_PDB_CUTOFF)
    groups = [
        ("all", ok),
        (f"moving (> {LARGE_CHANGE_A:g} A)", moving),
        (f"moving, holo released after {AF2_PDB_CUTOFF}", moving[post[moving.index]]),
        (f"moving, holo and apo released after {AF2_PDB_CUTOFF}", moving[both[moving.index]]),
        (f"all, holo released after {AF2_PDB_CUTOFF}", ok[post]),
    ]
    for name, sub in groups:
        holo_like = (sub.rmsd_af_holo < sub.rmsd_af_apo).mean()
        print(f"\n[A] {name}, n = {len(sub)}")
        print(f"  AlphaFold closer to holo than apo: {holo_like:.1%}")
        print(f"  median site RMSD  AF-holo {sub.rmsd_af_holo.median():.2f} A,"
              f"  AF-apo {sub.rmsd_af_apo.median():.2f} A,"
              f"  apo-holo {sub.rmsd_apo_holo_common.median():.2f} A")  # fmt: skip
    print(f"\n[C] moving pockets, n = {len(moving)}")
    for t in (90, 70):
        print(f"  site pLDDT > {t}: {(moving.site_plddt_mean > t).mean():.1%}")


if __name__ == "__main__":
    main()
