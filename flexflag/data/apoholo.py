"""Apo-holo pair list (from APObind) and PDB structure retrieval."""

import json
from dataclasses import dataclass

import pandas as pd

from flexflag.config import PAIRS_CSV
from flexflag.data.cache import fetch

# Columns kept from APObind's apobind_all.csv. Only identifiers and APObind's
# own quality metrics; their RMSDs are kept for cross-checking, never as labels.
PAIR_COLUMNS = [
    "holo_id",
    "holo_chains",
    "apo_id",
    "apo_chains",
    "apo_resolution",
    "sequence_identity",
    "sequence_coverage",
    "tmscore",
    "backbone_rmsd",
    "side_chain_rmsd",
]


@dataclass(frozen=True)
class Pair:
    holo_id: str
    holo_chains: tuple[str, ...]
    apo_id: str
    apo_chains: tuple[str, ...]


def load_pairs(path=PAIRS_CSV) -> pd.DataFrame:
    return pd.read_csv(path)


def row_to_pair(row) -> Pair:
    return Pair(
        holo_id=row.holo_id.lower(),
        holo_chains=tuple(row.holo_chains.split()),
        apo_id=row.apo_id.lower(),
        apo_chains=tuple(row.apo_chains.split()),
    )


def structure_path(pdb_id: str):
    pdb_id = pdb_id.lower()
    return fetch(f"https://files.rcsb.org/download/{pdb_id}.cif.gz", "pdb", f"{pdb_id}.cif.gz")


def uniprot_for_chain(pdb_id: str, chain: str) -> str | None:
    """UniProt accession for an author chain ID, via PDBe SIFTS."""
    pdb_id = pdb_id.lower()
    path = fetch(
        f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}", "sifts", f"{pdb_id}.json"
    )
    data = json.loads(path.read_text())
    for acc, entry in data.get(pdb_id, {}).get("UniProt", {}).items():
        if any(m["chain_id"] == chain for m in entry["mappings"]):
            return acc
    return None
