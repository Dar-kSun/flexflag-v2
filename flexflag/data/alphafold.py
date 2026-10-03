"""AlphaFold DB confidence retrieval (pLDDT and PAE). Prediction-only: no PDB data."""

import json

import numpy as np

from flexflag.data.cache import fetch

AFDB = "https://alphafold.ebi.ac.uk"


def _entry(uniprot: str) -> dict:
    path = fetch(f"{AFDB}/api/prediction/{uniprot}", "afdb", f"{uniprot}.api.json")
    entries = json.loads(path.read_text())
    # F1 is the full-length model; longer proteins are split into fragments we skip.
    for e in entries:
        if e["entryId"] == f"AF-{uniprot}-F1":
            return e
    raise LookupError(f"no AF-{uniprot}-F1 entry")


def plddt(uniprot: str) -> np.ndarray:
    url = _entry(uniprot)["plddtDocUrl"]
    path = fetch(url, "afdb", url.rsplit("/", 1)[-1])
    return np.asarray(json.loads(path.read_text())["confidenceScore"], dtype=float)


def pae(uniprot: str) -> np.ndarray:
    url = _entry(uniprot)["paeDocUrl"]
    path = fetch(url, "afdb", url.rsplit("/", 1)[-1])
    doc = json.loads(path.read_text())
    doc = doc[0] if isinstance(doc, list) else doc
    return np.asarray(doc["predicted_aligned_error"], dtype=float)


def entry_info(uniprot: str) -> dict:
    """Sequence and annotation for the F1 model (sequence is what AlphaFold predicted)."""
    e = _entry(uniprot)
    return {
        "sequence": e["sequence"],
        "organism": e.get("organismScientificName"),
        "description": e.get("uniprotDescription"),
        "af_version": e.get("latestVersion"),
    }


def model_path(uniprot: str):
    """Local path to the AlphaFold DB model (mmCIF) for the F1 entry."""
    url = _entry(uniprot)["cifUrl"]
    return fetch(url, "afdb_models", url.rsplit("/", 1)[-1])
