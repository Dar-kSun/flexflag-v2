"""Binding-site RMSD between apo and holo structures.

Fixed method, stated once:
  1. First model only, hydrogens and waters removed, first altloc kept.
  2. Ligand: the non-polymer, non-additive residue (>= MIN_LIGAND_ATOMS heavy atoms)
     with the most heavy atoms within SITE_CUTOFF_A of a listed holo chain.
  3. Site: residues of that holo chain with any heavy atom within SITE_CUTOFF_A
     of the ligand.
  4. Holo and apo chains are matched residue-by-residue by global sequence alignment.
  5. Site C-alpha atoms are superposed (Kabsch, gemmi.superpose_positions) and the
     C-alpha RMSD of the site is the label. Superposing on the site itself measures
     how much the pocket deforms, not how far it moves with the rest of the protein.
"""

import re
from dataclasses import dataclass

import gemmi
import numpy as np

from flexflag.config import LARGE_CHANGE_A, SITE_CUTOFF_A

MIN_LIGAND_ATOMS = 6

# Crystallisation additives, buffers and glycans that are never the ligand of interest.
NOT_LIGANDS = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "P6G", "PE4", "MPD", "DMS",
    "SO4", "PO4", "ACT", "FMT", "CIT", "TRS", "MES", "EPE", "BME", "IMD", "MLI", "TAR",
    "NAG", "NDG", "MAN", "BMA", "FUC", "GAL", "GLC", "BOG", "LDA", "CL", "NA", "K",
}  # fmt: skip


class LabelError(Exception):
    """This pair cannot be labelled; the message says why."""


@dataclass
class SiteRMSD:
    rmsd: float
    n_site: int  # site residues in the holo chain
    n_matched: int  # of those, present in the apo chain and used for the RMSD
    holo_chain: str
    apo_chain: str
    ligand: str  # e.g. "STI A/301"
    site_indices: list[int]  # 0-based positions in the holo chain's amino-acid residues
    site_auth_numbers: list[int]

    @property
    def large_change(self) -> bool:
        return self.rmsd > LARGE_CHANGE_A


def read_model(path) -> gemmi.Model:
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_waters()
    return st[0]


def protein_residues(chain: gemmi.Chain) -> list[gemmi.Residue]:
    """Amino-acid residues of a chain that have a C-alpha atom."""
    out = []
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        if info is not None and info.is_amino_acid() and res.find_atom("CA", "*") is not None:
            out.append(res)
    return out


def _coords(residues) -> np.ndarray:
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for r in residues for a in r])


def _near(a: np.ndarray, b: np.ndarray, cutoff: float) -> int:
    """Number of points in `a` within `cutoff` of any point in `b`."""
    d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)
    return int((d2.min(axis=1) <= cutoff**2).sum())


def find_ligand(model: gemmi.Model, chains, cutoff: float = SITE_CUTOFF_A):
    """Return (ligand residue, its chain name, the holo chain it contacts most)."""
    chain_xyz = {}
    for name in chains:
        ch = model.find_chain(name)
        if ch is not None and protein_residues(ch):
            chain_xyz[name] = _coords(protein_residues(ch))
    if not chain_xyz:
        raise LabelError(f"none of holo chains {chains} found")

    best, best_score = None, 0
    for ch in model:
        for res in ch:
            if res.het_flag != "H" or res.name in NOT_LIGANDS or len(res) < MIN_LIGAND_ATOMS:
                continue
            if gemmi.find_tabulated_residue(res.name) is not None and (
                gemmi.find_tabulated_residue(res.name).is_amino_acid()
            ):
                continue
            lig = _coords([res])
            for name, xyz in chain_xyz.items():
                score = _near(lig, xyz, cutoff)
                if score > best_score:
                    best, best_score = (res, ch.name, name), score
    if best is None:
        raise LabelError("no small-molecule ligand (peptide ligands are out of scope)")
    return best


def site_residues(chain_res, ligand: gemmi.Residue, cutoff: float = SITE_CUTOFF_A):
    lig = _coords([ligand])
    return [i for i, r in enumerate(chain_res) if _near(_coords([r]), lig, cutoff) > 0]


def one_letter(residues) -> list[str]:
    return [gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in residues]


def align_indices(holo_seq: list[str], apo_seq: list[str]) -> dict[int, int]:
    """Map holo index -> apo index for aligned positions with identical residues."""
    result = gemmi.align_string_sequences(holo_seq, apo_seq, [])
    mapping, i, j = {}, 0, 0
    for n, op in re.findall(r"(\d+)([MID])", result.cigar_str()):
        n = int(n)
        if op == "M":
            for k in range(n):
                if holo_seq[i + k] == apo_seq[j + k]:
                    mapping[i + k] = j + k
            i, j = i + n, j + n
        elif op == "I":  # in holo only
            i += n
        else:  # "D": in apo only
            j += n
    return mapping


def _ca(res: gemmi.Residue) -> gemmi.Position:
    return res.find_atom("CA", "*").pos


def site_rmsd(
    holo_path, holo_chains, apo_path, apo_chains, cutoff: float = SITE_CUTOFF_A
) -> SiteRMSD:
    holo = read_model(holo_path)
    apo = read_model(apo_path)

    ligand, lig_chain, holo_chain = find_ligand(holo, holo_chains, cutoff)
    holo_res = protein_residues(holo.find_chain(holo_chain))
    site = site_residues(holo_res, ligand, cutoff)
    if len(site) < 3:
        raise LabelError(f"binding site has {len(site)} residues")
    holo_seq = one_letter(holo_res)

    # Of the listed apo chains, use the one covering the most site residues.
    best = None
    for name in apo_chains:
        ch = apo.find_chain(name)
        if ch is None:
            continue
        apo_res = protein_residues(ch)
        mapping = align_indices(holo_seq, one_letter(apo_res))
        matched = [i for i in site if i in mapping]
        if best is None or len(matched) > len(best[2]):
            best = (name, apo_res, matched, mapping)
    if best is None:
        raise LabelError(f"none of apo chains {apo_chains} found")
    apo_chain, apo_res, matched, mapping = best
    if len(matched) < max(3, 0.8 * len(site)):
        raise LabelError(f"only {len(matched)}/{len(site)} site residues present in apo")

    fixed = [_ca(holo_res[i]) for i in matched]
    moving = [_ca(apo_res[mapping[i]]) for i in matched]
    sup = gemmi.superpose_positions(fixed, moving)
    return SiteRMSD(
        rmsd=float(sup.rmsd),
        n_site=len(site),
        n_matched=len(matched),
        holo_chain=holo_chain,
        apo_chain=apo_chain,
        ligand=f"{ligand.name} {lig_chain}/{ligand.seqid.num}",
        site_indices=site,
        site_auth_numbers=[holo_res[i].seqid.num for i in site],
    )
