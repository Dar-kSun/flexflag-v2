"""Per-protein pocket analysis: AlphaFold model vs apo and holo at the site, and
pocket features. Shared by scripts/03_pocket.py and scripts/06_external_eval.py."""

import gemmi
import numpy as np

from flexflag.data import alphafold
from flexflag.data.apoholo import structure_path
from flexflag.features.pocket import pocket_features
from flexflag.labels import align_indices, one_letter, protein_residues, read_model, site_rmsd

# AlphaFold2 (which AlphaFold DB models come from) was trained on PDB entries
# released up to this date. Entries released later cannot have been memorised.
AF2_PDB_CUTOFF = "2018-04-30"


def release_date(path) -> str:
    """Initial release date: the first entry in the mmCIF revision history."""
    block = gemmi.cif.read(str(path)).sole_block()
    return min(block.find_values("_pdbx_audit_revision_history.revision_date"))


def ca(res):
    return res.find_atom("CA", "*").pos


def rmsd(a, b) -> float:
    return float(gemmi.superpose_positions(a, b).rmsd)


def analyse(row) -> dict:
    holo_p, apo_p = structure_path(row.holo_id), structure_path(row.apo_id)
    lab = site_rmsd(holo_p, [row.holo_chain], apo_p, [row.apo_chain])
    if abs(lab.rmsd - row.site_rmsd) > 1e-3:
        raise AssertionError(f"{row.uniprot}: label did not reproduce")

    holo_res = protein_residues(read_model(holo_p).find_chain(row.holo_chain))
    apo_res = protein_residues(read_model(apo_p).find_chain(row.apo_chain))
    af_res = protein_residues(read_model(alphafold.model_path(row.uniprot))[0])
    holo_seq = one_letter(holo_res)
    to_apo = align_indices(holo_seq, one_letter(apo_res))
    to_af = align_indices(holo_seq, one_letter(af_res))

    site_af = [to_af[i] for i in lab.site_indices if i in to_af]
    common = [i for i in lab.site_indices if i in to_apo and i in to_af]
    out = {
        "uniprot": row.uniprot,
        "holo_release": release_date(holo_p),
        "apo_release": release_date(apo_p),
        "n_site_af": len(site_af),
        "n_common": len(common),
    }
    if len(common) >= 3:
        h = [ca(holo_res[i]) for i in common]
        a = [ca(apo_res[to_apo[i]]) for i in common]
        f = [ca(af_res[to_af[i]]) for i in common]
        out.update(rmsd_apo_holo_common=rmsd(h, a), rmsd_af_holo=rmsd(h, f),
                   rmsd_af_apo=rmsd(a, f))  # fmt: skip
    if len(site_af) >= 3:
        out.update(pocket_features(np.array(site_af), alphafold.plddt(row.uniprot),
                                   alphafold.pae(row.uniprot)))  # fmt: skip
    return out
