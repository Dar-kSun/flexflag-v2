"""`flexflag check`: pocket-change risk for one protein, from AlphaFold DB only.

The pocket is an input (UniProt residue numbers, or residues near a ligand in a PDB
entry you supply). The risk comes from the rule frozen on APObind and validated on
PDB-2019+ (flexflag/rule.json, written by scripts/07_export_rule.py).
"""

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from flexflag.config import LOW_PLDDT, SITE_CUTOFF_A
from flexflag.data import alphafold

RULE = json.loads((Path(__file__).parent / "rule.json").read_text(encoding="utf-8"))


@dataclass
class Report:
    uniprot: str
    description: str | None
    length: int
    protein_plddt: float
    pocket: list[int] = field(default_factory=list)  # 1-based UniProt positions
    pocket_source: str | None = None
    pocket_plddt: float | None = None
    pocket_low: int | None = None
    band: int | None = None  # 0 = lowest pocket pLDDT (highest risk) ... 4
    probability: float | None = None


def parse_residues(spec: str) -> list[int]:
    """'13,31,35-38' -> [13, 31, 35, 36, 37, 38]."""
    out = set()
    for part in spec.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            out.update(range(int(lo), int(hi) + 1))
        else:
            out.add(int(part))
    return sorted(out)


def pocket_from_pdb(uniprot: str, pdb: str) -> tuple[list[int], str]:
    """UniProt positions within SITE_CUTOFF_A of the main ligand in a PDB entry or file."""
    # Imported here: the whole-protein feature pipeline must never touch structure code.
    from flexflag.data.apoholo import structure_path
    from flexflag.labels import (
        align_indices,
        find_ligand,
        one_letter,
        protein_residues,
        read_model,
        site_residues,
    )

    path = Path(pdb) if Path(pdb).exists() else structure_path(pdb)
    model = read_model(path)
    chains = [ch.name for ch in model if protein_residues(ch)]
    ligand, _, chain = find_ligand(model, chains)
    res = protein_residues(model.find_chain(chain))
    site = site_residues(res, ligand, SITE_CUTOFF_A)
    seq = alphafold.entry_info(uniprot)["sequence"]
    mapping = align_indices(one_letter(res), list(seq))
    positions = sorted(mapping[i] + 1 for i in site if i in mapping)
    if len(positions) < 3:
        raise ValueError(f"only {len(positions)} pocket residues map onto {uniprot}")
    return positions, f"{ligand.name} in {path.name.split('.')[0]} chain {chain}"


def check(uniprot: str, pocket: list[int] | None = None, source: str | None = None) -> Report:
    info = alphafold.entry_info(uniprot)
    plddt = alphafold.plddt(uniprot)
    rep = Report(
        uniprot=uniprot,
        description=info["description"],
        length=len(plddt),
        protein_plddt=float(plddt.mean()),
    )
    if not pocket:
        return rep
    bad = [p for p in pocket if not 1 <= p <= len(plddt)]
    if bad:
        raise ValueError(f"residues outside 1..{len(plddt)}: {bad[:5]}")
    scores = plddt[np.array(pocket) - 1]
    rep.pocket, rep.pocket_source = pocket, source
    rep.pocket_plddt = float(scores.mean())
    rep.pocket_low = int((scores < LOW_PLDDT).sum())
    rep.band = int(np.digitize(rep.pocket_plddt, RULE["band_edges"]))
    z = RULE["intercept"] + RULE["coef_site_plddt"] * rep.pocket_plddt
    rep.probability = 1 / (1 + math.exp(-z))
    return rep


RISK_WORDS = ["highest", "elevated", "moderate", "low", "lowest"]


def format_report(rep: Report) -> str:
    thr = RULE["threshold_A"]
    lines = [f"{rep.uniprot}  {rep.description or ''}  ({rep.length} residues)", ""]
    if rep.pocket_plddt is None:
        ext = RULE["external_auroc_protein_plddt"]["value"]
        lines += [
            f"Whole-protein mean pLDDT: {rep.protein_plddt:.1f}",
            "",
            "No pocket given, so no flag. Whole-protein pLDDT does not predict binding-site",
            f"change (AUROC {ext:.2f} on {RULE['n_external']} unseen proteins). Pass the pocket:",
            f"  flexflag check {rep.uniprot} --residues 45,46,49-52",
            f"  flexflag check {rep.uniprot} --from-pdb <PDB ID or file with a ligand>",
        ]
        return "\n".join(lines)

    a = RULE["bands"]["apobind"][rep.band]
    e = RULE["bands"]["external"][rep.band]
    lines += [
        f"Pocket: {len(rep.pocket)} residues"
        + (f" ({rep.pocket_source})" if rep.pocket_source else ""),
        f"Pocket mean pLDDT: {rep.pocket_plddt:.1f}   (whole protein: {rep.protein_plddt:.1f})",
        "",
        f"Pocket-change risk (> {thr:g} Å, empty vs bound): {RISK_WORDS[rep.band].upper()}",
        f"  band {rep.band + 1} of 5 by pocket pLDDT ({a['band']})",
        f"  observed in this band: {a['frac_gt_2A']:.0%} of {a['n']} APObind pockets,",
        f"                         {e['frac_gt_2A']:.0%} of {e['n']} unseen PDB-2019+ pockets",
        f"  smoothed estimate: {rep.probability:.0%}  (base rate about 10%)",
        "",
        "Why:",
    ]
    if rep.band <= 1:
        lines += [
            "- Pocket pLDDT is in the lowest 40% of pockets studied. Even above 90 ('very",
            "  high' in AlphaFold's usual reading), lower pocket pLDDT marks more change.",
        ]
    elif rep.band >= 3:
        lines.append("- Pocket pLDDT is in the top 40% of pockets studied; these rarely reshape.")
    else:
        lines.append("- Pocket pLDDT is mid-range for the pockets studied.")
    if rep.pocket_low:
        lines.append(f"- {rep.pocket_low} pocket residue(s) have pLDDT < {LOW_PLDDT:g}.")
    st = RULE["alphafold_state_moving_external"]
    auroc = RULE["external_auroc_site_plddt"]["value"]
    lines += [
        "- When pockets do move, the AlphaFold model resembled the bound (holo) state in",
        f"  {st['frac_af_closer_to_holo']:.0%} of unseen cases. If you need the empty state, one",
        "  AlphaFold structure is likely the wrong one.",
        "",
        f"Limits: flags risk only, does not predict the other conformation. AUROC {auroc:.2f}",
        "on unseen proteins; it can miss rigid-domain closures over a confident pocket",
        "(adenylate kinase scores 'moderate' but moves 6 Å).",
    ]
    return "\n".join(lines)
