"""Label one protein from a list of candidate apo-holo pairs (shared by dataset builds)."""

import pandas as pd

from flexflag.data import alphafold
from flexflag.data.apoholo import structure_path
from flexflag.data.cache import NotFound
from flexflag.labels import LabelError, site_rmsd

MAX_CANDIDATES = 3  # pairs tried per protein before giving up
MAX_LENGTH = 1500  # AFDB PAE files get very large; also keeps the set to single chains


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
