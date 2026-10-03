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

import gemmi
import numpy as np
import pandas as pd

from flexflag.config import LARGE_CHANGE_A, ROOT
from flexflag.data import alphafold
from flexflag.data.apoholo import structure_path
from flexflag.features.pocket import pocket_features
from flexflag.labels import align_indices, one_letter, protein_residues, read_model, site_rmsd

RESULTS = ROOT / "results"
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
