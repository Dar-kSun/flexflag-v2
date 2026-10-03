"""Build the external PDB-2019+ validation set, as fixed in docs/plan-v0.3-external.md.

    python scripts/05_external_dataset.py [--workers 24] [--dry-run]

Independent of APObind: holo entries released on or after 2019-01-01 (after
AlphaFold2's training cutoff), paired with strict-apo entries of the same UniProt.
Labels use the same code as v0.1. Writes results/external/dataset.csv and
results/external/dropped.csv.
"""

import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from flexflag.config import ROOT
from flexflag.data.apoholo import SIFTS_TSV
from flexflag.data.cache import fetch
from flexflag.data.rcsb import entry_details, search_entries
from flexflag.dataset import MAX_CANDIDATES, label_protein
from flexflag.labels import NOT_LIGANDS

OUT = ROOT / "results" / "external"
HOLO_FROM = "2019-01-01"
MIN_LIGAND_WEIGHT = 150.0
MIN_RANGE_OVERLAP = 0.9

COMMON = [
    ("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
    ("rcsb_entry_info.resolution_combined", "less_or_equal", 2.5),
    ("rcsb_entry_info.polymer_entity_count_protein", "equals", 1),
    ("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
]
HOLO_TERMS = COMMON + [
    ("rcsb_accession_info.initial_release_date", "greater_or_equal", HOLO_FROM),
    ("rcsb_entry_info.nonpolymer_entity_count", "greater_or_equal", 1),
]
APO_TERMS = COMMON + [("rcsb_entry_info.nonpolymer_entity_count", "equals", 0)]


def sifts_entries() -> pd.DataFrame:
    """One row per PDB entry mapping to exactly one UniProt: chains and UniProt range."""
    path = fetch(SIFTS_TSV, "sifts", "pdb_chain_uniprot.tsv.gz")
    # keep_default_na=False: chain IDs such as "NA" are real names, not missing values.
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=["PDB", "CHAIN", "SP_PRIMARY", "SP_BEG", "SP_END"],
        dtype={"PDB": str, "CHAIN": str, "SP_PRIMARY": str},
        keep_default_na=False,
    )
    df["CHAIN"] = df.CHAIN.astype(str)
    g = df.groupby("PDB")
    one = g.SP_PRIMARY.nunique() == 1
    df = df[df.PDB.isin(one[one].index)]
    return df.groupby("PDB").agg(
        uniprot=("SP_PRIMARY", "first"),
        chains=("CHAIN", lambda c: " ".join(sorted(set(c)))),
        beg=("SP_BEG", "min"),
        end=("SP_END", "max"),
    )


def overlap_ok(a, b) -> bool:
    inter = min(a.end, b.end) - max(a.beg, b.beg) + 1
    return inter >= MIN_RANGE_OVERLAP * (a.end - a.beg + 1) and inter >= MIN_RANGE_OVERLAP * (
        b.end - b.beg + 1
    )


def has_real_ligand(ligands) -> bool:
    return any(cid not in NOT_LIGANDS and (fw or 0) >= MIN_LIGAND_WEIGHT for cid, fw in ligands)


def candidates(sifts: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict]:
    holo = [i for i in search_entries("holo_2019", HOLO_TERMS) if i in sifts.index]
    apo = [i for i in search_entries("apo_strict", APO_TERMS) if i in sifts.index]
    holo_df, apo_df = sifts.loc[holo], sifts.loc[apo]
    shared = set(holo_df.uniprot) & set(apo_df.uniprot)
    holo_df = holo_df[holo_df.uniprot.isin(shared)]
    apo_df = apo_df[apo_df.uniprot.isin(shared)]
    info = entry_details(sorted(set(holo_df.index) | set(apo_df.index)))
    for df in (holo_df, apo_df):
        df["resolution"] = [info.get(i, {}).get("resolution") for i in df.index]
        df["release"] = [info.get(i, {}).get("release") for i in df.index]
    holo_df = holo_df[[has_real_ligand(info.get(i, {}).get("ligands", [])) for i in holo_df.index]]
    stats = {
        "holo_search": len(holo),
        "apo_search": len(apo),
        "uniprot_with_both": len(shared),
        "holo_with_real_ligand": len(holo_df),
    }

    groups = {}
    apo_by_acc = {acc: g.sort_values(["resolution"]) for acc, g in apo_df.groupby("uniprot")}
    for acc, hs in holo_df.sort_values(["resolution"]).groupby("uniprot", sort=True):
        rows = []
        apos = list(apo_by_acc[acc].reset_index(names="pdb").itertuples())
        for h in hs.reset_index(names="pdb").sort_values(["resolution", "pdb"]).itertuples():
            match = next((a for a in apos if overlap_ok(h, a)), None)
            if match is not None:
                rows.append(
                    {
                        "holo_id": h.pdb,
                        "holo_chains": h.chains,
                        "apo_id": match.pdb,
                        "apo_chains": match.chains,
                        "holo_release": h.release,
                        "apo_release": match.release,
                    }  # fmt: skip
                )
            if len(rows) == MAX_CANDIDATES:
                break
        if rows:
            groups[acc] = pd.DataFrame(rows)
    stats["uniprot_with_pair"] = len(groups)
    return groups, stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--dry-run", action="store_true", help="print counts, download nothing")
    args = ap.parse_args()

    t0 = time.perf_counter()
    groups, stats = candidates(sifts_entries())
    print(stats, flush=True)
    if args.dry_run:
        return

    rows, dropped = [], []
    with ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(label_protein, acc, g): (acc, g) for acc, g in groups.items()}
        for i, fut in enumerate(as_completed(futures), 1):
            acc, g = futures[fut]
            try:
                row, drops = fut.result()
            except Exception as e:  # network trouble etc.: record and carry on
                row, drops = None, [{"uniprot": acc, "reason": f"error: {e}"[:200]}]
            dropped += drops
            if row:
                pair = g[g.holo_id == row["holo_id"]].iloc[0]
                rows.append({**row, "holo_release": pair.holo_release,
                             "apo_release": pair.apo_release})  # fmt: skip
            if i % 200 == 0:
                print(f"{i}/{len(groups)} proteins, {len(rows)} labelled, "
                      f"{time.perf_counter() - t0:.0f}s", flush=True)  # fmt: skip

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).sort_values("uniprot").to_csv(OUT / "dataset.csv", index=False)
    dropped = pd.DataFrame(dropped)
    dropped.sort_values([c for c in ("uniprot", "holo_id", "reason") if c in dropped]).to_csv(
        OUT / "dropped.csv", index=False
    )
    (OUT / "build_stats.txt").write_text(
        "\n".join(f"{k}: {v}" for k, v in {**stats, "labelled": len(rows)}.items()) + "\n"
    )
    print(f"done: {len(rows)} proteins labelled in {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
