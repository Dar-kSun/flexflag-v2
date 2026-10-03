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
from flexflag.data import alphafold
from flexflag.data.apoholo import load_pairs, sifts_chain_map, structure_path
from flexflag.data.cache import NotFound
from flexflag.labels import LabelError, site_rmsd

RESULTS = ROOT / "results"

# Pair-level filters on APObind's own metrics, applied before any download.
MIN_TMSCORE = 0.5  # below this the apo is probably not the same fold
MIN_IDENTITY = 0.9
MIN_COVERAGE = 0.9
MAX_CANDIDATES = 3  # pairs tried per protein before giving up
MAX_LENGTH = 1500  # AFDB PAE files get very large; also keeps the set to single chains


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


def label_protein(uniprot: str, candidates: pd.DataFrame) -> tuple[dict | None, list[dict]]:
    """Try candidate pairs in order; return the first that labels cleanly."""
    dropped = []
    try:
        info = alphafold.entry_info(uniprot)
        if len(info["sequence"]) > MAX_LENGTH:
            raise LookupError(f"AFDB model longer than {MAX_LENGTH}")
        alphafold.plddt(uniprot)
        alphafold.pae(uniprot)
    except (NotFound, LookupError) as e:
        reason = "no AFDB F1 model" if isinstance(e, NotFound) else str(e)
        return None, [{"uniprot": uniprot, "reason": reason}]

    for _, r in candidates.head(MAX_CANDIDATES).iterrows():
        try:
            lab = site_rmsd(
                structure_path(r.holo_id),
                r.holo_chains.split(),
                structure_path(r.apo_id),
                r.apo_chains.split(),
            )
        except (LabelError, NotFound, ValueError, RuntimeError) as e:
            dropped.append(
                {"uniprot": uniprot, "holo_id": r.holo_id, "apo_id": r.apo_id, "reason": str(e)}
            )
            continue
        return {
            "uniprot": uniprot,
            "holo_id": r.holo_id,
            "holo_chain": lab.holo_chain,
            "apo_id": r.apo_id,
            "apo_chain": lab.apo_chain,
            "ligand": lab.ligand,
            "n_site": lab.n_site,
            "n_matched": lab.n_matched,
            "site_rmsd": round(lab.rmsd, 3),
            "length": len(info["sequence"]),
            "organism": info["organism"],
            "description": info["description"],
            "n_apobind_pairs": len(candidates),
            "sequence": info["sequence"],
        }, dropped
    dropped.append({"uniprot": uniprot, "reason": "no candidate pair could be labelled"})
    return None, dropped


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
    pd.DataFrame(dropped).to_csv(RESULTS / "dropped.csv", index=False)
    print(f"done: {len(rows)} proteins labelled in {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
