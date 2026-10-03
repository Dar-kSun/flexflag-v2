"""Build the per-protein dataset: one labelled apo-holo pair per UniProt accession.

    python scripts/01_build_dataset.py [--workers 16] [--limit N]

Writes results/dataset.csv (one row per protein) and results/dropped.csv (every
pair or protein that was filtered out, with the reason). Downloads are cached in
data/cache/, so re-running resumes rather than starting over.
"""

import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from flexflag.config import ROOT
from flexflag.data.apoholo import load_pairs, sifts_chain_map
from flexflag.dataset import label_protein

RESULTS = ROOT / "results"

# Pair-level filters on APObind's own metrics, applied before any download.
MIN_TMSCORE = 0.5  # below this the apo is probably not the same fold
MIN_IDENTITY = 0.9
MIN_COVERAGE = 0.9


def prefilter(pairs: pd.DataFrame, chain_map) -> tuple[pd.DataFrame, list[dict]]:
    dropped = []

    def drop(mask, reason):
        for _, r in pairs[mask].iterrows():
            dropped.append({"holo_id": r.holo_id, "apo_id": r.apo_id, "reason": reason})
        return pairs[~mask]

    pairs = drop(pairs.apo_resolution <= 0, "apo has no resolution (NMR or unknown)")
    pairs = drop(pairs.tmscore < MIN_TMSCORE, f"APObind TM-score < {MIN_TMSCORE}")
    pairs = drop(pairs.sequence_identity < MIN_IDENTITY, f"identity < {MIN_IDENTITY}")
    pairs = drop(pairs.sequence_coverage < MIN_COVERAGE, f"coverage < {MIN_COVERAGE}")

    def acc(pdb_id, chains):
        accs = {chain_map.get((pdb_id.lower(), c)) for c in chains.split()}
        accs.discard(None)
        return accs.pop() if len(accs) == 1 else None

    pairs = pairs.assign(
        holo_uniprot=[acc(h, c) for h, c in zip(pairs.holo_id, pairs.holo_chains, strict=False)],
        apo_uniprot=[acc(a, c) for a, c in zip(pairs.apo_id, pairs.apo_chains, strict=False)],
    )
    pairs = drop(pairs.holo_uniprot.isna(), "holo chains do not map to one UniProt")
    pairs = drop(pairs.holo_uniprot != pairs.apo_uniprot, "apo and holo map to different UniProt")
    return pairs, dropped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit", type=int, help="only the first N proteins (for testing)")
    args = ap.parse_args()

    t0 = time.perf_counter()
    pairs = load_pairs()
    n_pairs = len(pairs)
    pairs, dropped = prefilter(pairs, sifts_chain_map())
    # Best candidates first: closest sequence match, then best apo resolution.
    pairs = pairs.sort_values(
        ["holo_uniprot", "sequence_identity", "sequence_coverage", "apo_resolution", "holo_id"],
        ascending=[True, False, False, True, True],
    )
    groups = list(pairs.groupby("holo_uniprot", sort=True))
    if args.limit:
        groups = groups[: args.limit]
    print(f"{n_pairs} pairs -> {len(pairs)} after filters -> {len(groups)} proteins", flush=True)

    rows = []
    with ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(label_protein, acc, g): acc for acc, g in groups}
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                row, drops = fut.result()
            except Exception as e:  # network trouble etc.: record and carry on
                row, drops = None, [{"uniprot": futures[fut], "reason": f"error: {e}"[:200]}]
            dropped += drops
            if row:
                rows.append(row)
            if i % 100 == 0:
                print(
                    f"{i}/{len(groups)} proteins, {len(rows)} labelled, "
                    f"{time.perf_counter() - t0:.0f}s",
                    flush=True,
                )

    RESULTS.mkdir(exist_ok=True)
    pd.DataFrame(rows).sort_values("uniprot").to_csv(RESULTS / "dataset.csv", index=False)
    dropped = pd.DataFrame(dropped)
    sort_cols = [c for c in ("uniprot", "holo_id", "apo_id", "reason") if c in dropped]
    dropped.sort_values(sort_cols, na_position="first").to_csv(RESULTS / "dropped.csv", index=False)
    print(f"done: {len(rows)} proteins labelled in {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
