"""Feature pipeline. Inputs are a UniProt accession and AlphaFold DB files only.

Nothing here may read a PDB entry, the holo structure or the ligand: that would
leak the label. tests/test_no_holo_leakage.py enforces it.
"""

from flexflag.data import alphafold
from flexflag.features.pae import pae_features
from flexflag.features.plddt import plddt_features
from flexflag.features.sequence import sequence_features

BASELINE_FEATURES = ["plddt_mean"]


def compute_features(uniprot: str) -> dict[str, float]:
    plddt = alphafold.plddt(uniprot)
    pae = alphafold.pae(uniprot)
    seq = alphafold.entry_info(uniprot)["sequence"]
    if not (len(plddt) == len(seq) == pae.shape[0]):
        raise ValueError(f"{uniprot}: AFDB file lengths disagree")
    return {**plddt_features(plddt), **pae_features(pae, plddt), **sequence_features(seq)}
