"""Fetch 10 apo-holo pairs end to end: one label (site RMSD) and one feature (mean pLDDT).

    python scripts/00_feasibility.py [--n 10] [--seed 0]

Prints per-pair timings and outcomes. If data/apobind_all.csv is present it also
prints `in_apobind`: the fraction of our site residues that are in APObind's own
site list (theirs uses 6 A and, it appears, hydrogens, so theirs is larger; a
value near 1 means we picked the same pocket). Exits non-zero if
fewer than half the pairs make it through.
"""

import argparse
import sys
import time

import pandas as pd

from flexflag.config import DATA_DIR, SITE_CUTOFF_A
from flexflag.data import alphafold
from flexflag.data.apoholo import load_pairs, row_to_pair, structure_path, uniprot_for_chain
from flexflag.features.plddt import mean_plddt
from flexflag.labels import site_rmsd


def apobind_site(holo_id: str):
    path = DATA_DIR / "apobind_all.csv"
    if not path.exists():
        return None
    full = pd.read_csv(path, index_col=0, usecols=[0, 1, 11])
    row = full[full.holo_id == holo_id]
    return set(map(int, row.iloc[0].holo_bind_indices.split())) if len(row) else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    pairs = load_pairs().sample(n=args.n, random_state=args.seed)
    rows, t_all = [], time.perf_counter()
    for _, r in pairs.iterrows():
        p = row_to_pair(r)
        out = {"holo": p.holo_id, "apo": p.apo_id}
        t0 = time.perf_counter()
        try:
            hp, apath = structure_path(p.holo_id), structure_path(p.apo_id)
            out["t_pdb"] = time.perf_counter() - t0

            t1 = time.perf_counter()
            lab = site_rmsd(hp, p.holo_chains, apath, p.apo_chains)
            out["t_label"] = time.perf_counter() - t1
            out.update(
                site_rmsd=round(lab.rmsd, 2),
                site=f"{lab.n_matched}/{lab.n_site}",
                ligand=lab.ligand,
                chains=f"{lab.holo_chain}->{lab.apo_chain}",
            )
            theirs = apobind_site(p.holo_id)
            if theirs:
                ours = set(lab.site_indices)  # APObind lists 0-based chain positions
                out["in_apobind"] = round(len(ours & theirs) / len(ours), 2)

            t2 = time.perf_counter()
            acc = uniprot_for_chain(p.holo_id, lab.holo_chain)
            out["uniprot"] = acc
            out["mean_plddt"] = round(mean_plddt(alphafold.plddt(acc)), 1) if acc else None
            out["t_afdb"] = time.perf_counter() - t2
            out["status"] = "ok" if acc else "no uniprot"
        except Exception as e:  # report every failure mode, don't stop
            out["status"] = f"{type(e).__name__}: {e}"[:70]
        out["t_total"] = time.perf_counter() - t0
        rows.append(out)
        print(f"{p.holo_id}-{p.apo_id}: {out['status']} ({out['t_total']:.1f}s)", flush=True)

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print(f"\nsite cutoff {SITE_CUTOFF_A} A, seed {args.seed}\n")
    print(df.round(1).to_string(index=False))
    n_ok = (df.status == "ok").sum()
    print(f"\n{n_ok}/{len(df)} pairs fully processed in {time.perf_counter() - t_all:.0f}s")
    return 0 if n_ok >= len(df) / 2 else 1


if __name__ == "__main__":
    sys.exit(main())
