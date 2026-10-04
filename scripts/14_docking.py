"""Docking validation (D1 of docs/plan-v0.6-docking.md): does a flagged pocket make
docking into the AlphaFold model fail?

    python scripts/14_docking.py [--workers 18] [--limit N] [--only UNIPROT ...]

For each protein, the crystal ligand (bond orders from the RCSB ModelServer) is
re-docked with AutoDock Vina into three receptors, all superposed on the holo pocket:
the holo crystal chain (control), the apo crystal chain, and the AlphaFold DB model.
Success = top-ranked pose within 2 Å heavy-atom RMSD (symmetry-aware) of the crystal
pose. This is a validation experiment only; flexflag itself does not dock.

Needs: pip install meeko rdkit openbabel-wheel, and the Vina 1.2.7 executable
(FLEXFLAG_VINA, else %LOCALAPPDATA%/flexflag/vina_1.2.7_win.exe, else `vina` on PATH).
Writes results/docking/runs.csv (resumable: finished proteins are skipped) and
results/docking/summary.md via --summarise.
"""

import argparse
import os
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import requests

from flexflag.config import CACHE_DIR, ROOT
from flexflag.data import alphafold
from flexflag.data.apoholo import structure_path
from flexflag.evaluate import RESULTS
from flexflag.labels import (
    NOT_LIGANDS,
    align_indices,
    find_ligand,
    one_letter,
    protein_residues,
    read_model,
    site_residues,
)

OUT = RESULTS / "docking"
LIG_CACHE = CACHE_DIR / "ligands"
LOCAL = Path(os.environ.get("LOCALAPPDATA", "")) / "flexflag"

# Pre-declared limits (docs/plan-v0.6-docking.md).
MAX_HEAVY_ATOMS = 60
MAX_ROTATABLE = 15
COVALENT_CUTOFF = 2.0  # Å: a ligand heavy atom this close to protein -> covalent, skip
COFACTOR_RADIUS = 8.0  # Å: holo cofactors and metals this close to the ligand are kept
BOX_PAD = 10.0  # Å added to the ligand extent on each axis
BOX_MIN = 20.0
EXHAUSTIVENESS = 8
SUCCESS_RMSD = 2.0
SEED = 0


def vina_exe() -> str:
    if os.environ.get("FLEXFLAG_VINA"):
        return os.environ["FLEXFLAG_VINA"]
    local = LOCAL / "vina_1.2.7_win.exe"
    return str(local) if local.exists() else (shutil.which("vina") or "vina")


def obabel_exe() -> str:
    here = Path(os.sys.executable).parent
    for name in ("obabel.exe", "obabel"):
        if (here / name).exists():
            return str(here / name)
    return shutil.which("obabel") or "obabel"


def ligand_sdf(pdb: str, chain: str, num: int, comp: str) -> str:
    """One ligand instance with CCD bond orders, from the RCSB ModelServer. Cached."""
    path = LIG_CACHE / f"{pdb}_{chain}_{num}_{comp}.sdf"
    if path.exists():
        return path.read_text()
    url = (
        f"https://models.rcsb.org/v1/{pdb}/ligand?auth_asym_id={chain}&auth_seq_id={num}"
        f"&label_comp_id={comp}&encoding=sdf"
    )
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(r.text)
    return r.text


def crystal_ligand(sdf_text: str):
    from rdkit import Chem

    sup = Chem.SDMolSupplier()
    sup.SetData(sdf_text, removeHs=False, sanitize=True)
    mols = [m for m in sup if m is not None]
    if not mols:
        raise ValueError("ligand SDF could not be parsed")
    return Chem.RemoveHs(mols[0])


def chain_pdb(residues, path: Path, transform=None, hetero=()) -> None:
    """Write protein residues (optionally transformed) plus hetero groups (already in the
    holo frame, never transformed) as a minimal PDB."""
    st = gemmi.Structure()
    model = gemmi.Model("1")
    ch = gemmi.Chain("A")
    for r in residues:
        nr = gemmi.Residue()
        nr.name, nr.seqid, nr.het_flag = r.name, r.seqid, "A"
        for a in r:
            na = gemmi.Atom()
            na.name, na.element, na.b_iso, na.occ = a.name, a.element, a.b_iso, 1.0
            na.pos = gemmi.Position(transform.apply(a.pos)) if transform else a.pos
            nr.add_atom(na)
        ch.add_residue(nr)
    model.add_chain(ch)
    if hetero:
        hch = gemmi.Chain("Z")
        for k, r in enumerate(hetero, 1):
            nr = gemmi.Residue()
            nr.name, nr.seqid, nr.het_flag = r.name, gemmi.SeqId(k, " "), "H"
            for a in r:
                nr.add_atom(a)
            hch.add_residue(nr)
        model.add_chain(hch)
    st.add_model(model)
    st.setup_entities()
    st.write_pdb(str(path))


def cofactors(model, ligand, lig_xyz: np.ndarray) -> list:
    """Non-polymer groups (cofactors, metal ions) within COFACTOR_RADIUS of the ligand,
    excluding the docked ligand itself, waters and crystallisation additives."""
    keep = []
    for ch in model:
        for r in ch:
            if r.het_flag != "H" or r.name in NOT_LIGANDS or r.is_water():
                continue
            if r.entity_type == gemmi.EntityType.Polymer:
                continue
            if r.name == ligand.name and r.seqid == ligand.seqid:
                continue
            xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in r])
            if np.min(np.linalg.norm(xyz[:, None] - lig_xyz[None], axis=-1)) <= COFACTOR_RADIUS:
                keep.append(r)
    return keep


def prepare_receptor(pdb: Path, pdbqt: Path) -> None:
    subprocess.run(
        [obabel_exe(), str(pdb), "-xr", "-p", "7.4", "-O", str(pdbqt)],
        check=True,
        capture_output=True,
    )


def prepare_ligand(mol, pdbqt: Path):
    """Fresh ETKDG conformer (no memory of the crystal pose), protonated as given."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.AddHs(Chem.Mol(mol))
    if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
        raise ValueError("RDKit could not embed the ligand")
    AllChem.MMFFOptimizeMolecule(m)
    setups = MoleculePreparation().prepare(m)
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"meeko: {err}")
    pdbqt.write_text(text)


def docked_rmsds(out_pdbqt: Path, ref) -> tuple[list[float], float]:
    """Heavy-atom, symmetry-aware RMSD of every pose (in Vina's rank order) to the crystal
    pose, without re-alignment, and the top pose's score."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    pm = PDBQTMolecule(out_pdbqt.read_text(), is_dlg=False, skip_typing=True)
    probe = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
    rmsds = [float(rdMolAlign.CalcRMS(probe, ref, prbId=c.GetId())) for c in probe.GetConformers()]
    score = np.nan
    for line in out_pdbqt.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            score = float(line.split()[3])
            break
    return rmsds, score


def dock_one(row: dict) -> list[dict]:
    from rdkit.Chem import Descriptors

    out = []
    base = {
        "uniprot": row["uniprot"],
        "dataset": row["dataset"],
        "holo_id": row["holo_id"],
        "apo_id": row["apo_id"],
        "site_rmsd": row["site_rmsd"],
    }
    try:
        holo = read_model(structure_path(row["holo_id"]))
        ligand, lig_chain, chain = find_ligand(holo, [row["holo_chain"]])
        holo_res = protein_residues(holo.find_chain(chain))
        site = site_residues(holo_res, ligand)
        ref = crystal_ligand(ligand_sdf(row["holo_id"], lig_chain, ligand.seqid.num, ligand.name))
        n_heavy = ref.GetNumHeavyAtoms()
        n_rot = Descriptors.NumRotatableBonds(ref)
        base.update(
            ligand=f"{ligand.name} {lig_chain}/{ligand.seqid.num}",
            n_heavy=n_heavy,
            n_rotatable=n_rot,
        )
        if n_heavy > MAX_HEAVY_ATOMS or n_rot > MAX_ROTATABLE:
            return [{**base, "skipped": "ligand too large or flexible"}]
        lig_xyz = ref.GetConformer().GetPositions()
        prot_xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in holo_res for a in r])
        if np.min(np.linalg.norm(lig_xyz[:, None] - prot_xyz[None], axis=-1)) < COVALENT_CUTOFF:
            return [{**base, "skipped": "covalent ligand"}]

        # Receptors, all in the holo frame, superposed on the pocket C-alpha atoms.
        hetero = cofactors(holo, ligand, lig_xyz)
        base["cofactors"] = " ".join(r.name for r in hetero)
        receptors = {"holo": (holo_res, None)}
        holo_seq = one_letter(holo_res)
        apo_res = protein_residues(
            read_model(structure_path(row["apo_id"])).find_chain(row["apo_chain"])
        )
        af_res = protein_residues(read_model(alphafold.model_path(row["uniprot"]))[0])
        for name, res in (("apo", apo_res), ("alphafold", af_res)):
            mp = align_indices(holo_seq, one_letter(res))
            common = [i for i in site if i in mp]
            sup = gemmi.superpose_positions(
                [holo_res[i].find_atom("CA", "*").pos for i in common],
                [res[mp[i]].find_atom("CA", "*").pos for i in common],
            )
            receptors[name] = (res, sup.transform)
            base[f"pocket_ca_rmsd_{name}"] = float(sup.rmsd)

        lo, hi = lig_xyz.min(axis=0), lig_xyz.max(axis=0)
        center = (lo + hi) / 2
        size = np.maximum(hi - lo + BOX_PAD, BOX_MIN)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            prepare_ligand(ref, tmp / "lig.pdbqt")
            for name, (res, tf) in receptors.items():
                t0 = time.perf_counter()
                rec = {**base, "receptor": name}
                try:
                    chain_pdb(res, tmp / f"{name}.pdb", tf, hetero)
                    prepare_receptor(tmp / f"{name}.pdb", tmp / f"{name}.pdbqt")
                    subprocess.run(
                        [
                            vina_exe(),
                            "--receptor",
                            str(tmp / f"{name}.pdbqt"),
                            "--ligand",
                            str(tmp / "lig.pdbqt"),
                            "--center_x",
                            f"{center[0]:.3f}",
                            "--center_y",
                            f"{center[1]:.3f}",
                            "--center_z",
                            f"{center[2]:.3f}",
                            "--size_x",
                            f"{size[0]:.1f}",
                            "--size_y",
                            f"{size[1]:.1f}",
                            "--size_z",
                            f"{size[2]:.1f}",
                            "--exhaustiveness",
                            str(EXHAUSTIVENESS),
                            "--num_modes",
                            "9",
                            "--seed",
                            str(SEED),
                            "--cpu",
                            "1",
                            "--out",
                            str(tmp / f"{name}_out.pdbqt"),
                        ],
                        check=True,
                        capture_output=True,
                        timeout=3600,
                    )
                    rmsds, score = docked_rmsds(tmp / f"{name}_out.pdbqt", ref)
                    rec.update(
                        pose_rmsd=rmsds[0],
                        best_top3_rmsd=min(rmsds[:3]),
                        vina_score=score,
                        success=rmsds[0] <= SUCCESS_RMSD,
                        success_top3=min(rmsds[:3]) <= SUCCESS_RMSD,
                    )
                except Exception as e:  # one receptor failing does not lose the others
                    rec["error"] = f"{type(e).__name__}: {e}"[:150]
                rec["seconds"] = round(time.perf_counter() - t0, 1)
                out.append(rec)
    except Exception as e:
        out.append({**base, "error": f"{type(e).__name__}: {e}"[:150]})
    return out


def proteins() -> pd.DataFrame:
    a = pd.read_csv(RESULTS / "dataset.csv").assign(dataset="apobind")
    e = pd.read_csv(RESULTS / "external" / "dataset.csv").assign(dataset="external")
    return pd.concat([a, e], ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=18)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--only", nargs="*", help="UniProt accessions (testing)")
    ap.add_argument("--summarise", action="store_true", help="only rewrite summary.md")
    ap.add_argument("--runs", default=str(OUT / "runs.csv"), help="per-docking output CSV")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    runs_path = Path(args.runs)
    if args.summarise:
        summarise(pd.read_csv(runs_path))
        return
    ds = proteins()
    if args.only:
        ds = ds[ds.uniprot.isin(args.only)]
    done = set()
    if runs_path.exists():
        prev = pd.read_csv(runs_path)
        done = set(zip(prev.dataset, prev.uniprot, strict=True))
    todo = [r for r in ds.to_dict("records") if (r["dataset"], r["uniprot"]) not in done]
    todo.sort(key=lambda r: r["dataset"] != "external")  # primary endpoint first
    if args.limit:
        todo = todo[: args.limit]
    print(f"{len(todo)} proteins to dock ({len(done)} already done)", flush=True)
    t0 = time.perf_counter()
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(dock_one, r) for r in todo]
        for k, fut in enumerate(as_completed(futures), 1):
            rows = pd.DataFrame(fut.result())
            rows.to_csv(
                runs_path,
                mode="a",
                header=not runs_path.exists(),
                index=False,
                float_format="%.3f",
                lineterminator="\n",
            )
            if k % 20 == 0 or k == len(todo):
                print(f"{k}/{len(todo)} proteins, {time.perf_counter() - t0:.0f}s", flush=True)


def _auroc(y, score, groups, rng, n_boot=1000):
    from sklearn.metrics import roc_auc_score

    ok = ~np.isnan(score)
    y, score, groups = y[ok], score[ok], groups[ok]
    if y.min() == y.max():
        return None
    uniq = np.unique(groups)
    members = {g: np.flatnonzero(groups == g) for g in uniq}
    draws = []
    for _ in range(n_boot):
        idx = np.concatenate([members[g] for g in rng.choice(uniq, size=len(uniq))])
        if 0 < y[idx].sum() < len(idx):
            draws.append(roc_auc_score(y[idx], score[idx]))
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return f"{roc_auc_score(y, score):.3f} [{lo:.3f}, {hi:.3f}] ({int(y.sum())}/{len(y)})"


def summarise(runs: pd.DataFrame) -> None:
    """Endpoints 1-5 of docs/plan-v0.6-docking.md."""
    import json

    rule = json.loads((ROOT / "flexflag" / "rule.json").read_text(encoding="utf-8"))
    edges = rule["band_edges"]
    rng = np.random.default_rng(SEED)
    lines = [
        "Docking validation (docs/plan-v0.6-docking.md). AutoDock Vina 1.2.7, "
        "success = top pose within 2 Å of the crystal pose.",
        "",
    ]
    for name, sub in (("apobind", ""), ("external", "external/")):
        r = runs[runs.dataset == name]
        meta = r.drop_duplicates("uniprot").set_index("uniprot")
        skipped = meta.skipped.value_counts().to_dict() if "skipped" in meta else {}
        ok = r[r.receptor.notna() & r.pose_rmsd.notna()]
        wide = ok.pivot_table(
            index="uniprot", columns="receptor", values="success", aggfunc="first"
        )
        wide3 = ok.pivot_table(
            index="uniprot", columns="receptor", values="success_top3", aggfunc="first"
        )
        wide = wide.dropna(subset=[c for c in ("holo", "apo", "alphafold") if c in wide])
        if wide.empty:
            continue
        wide = wide.astype(bool)
        pk = pd.read_csv(RESULTS / f"{sub}pocket.csv").set_index("uniprot")
        feats = pd.read_csv(RESULTS / f"{sub}features.csv", index_col="uniprot")
        cl = pd.read_csv(RESULTS / f"{sub}clusters.csv", index_col="uniprot").cluster
        d = wide.join(meta[["site_rmsd"]]).join(pk.site_plddt_mean).join(feats.plddt_mean)
        d = d.dropna(subset=["site_plddt_mean"])
        groups = cl.loc[d.index].to_numpy()
        moving = d.site_rmsd > 2.0
        af_fail = ~d.alphafold
        struct_fail = af_fail & d.holo
        lines += [
            f"### {name}: {len(d)} proteins docked in all three receptors",
            f"(excluded before docking: {skipped})",
            "",
            "| Receptor | Success (top 1) | Success (any of top 3) |",
            "|---|---|---|",
        ]
        for rec in ("holo", "apo", "alphafold"):
            t3 = wide3.loc[d.index, rec].astype(bool).mean()
            lines.append(f"| {rec} | {d[rec].mean():.1%} | {t3:.1%} |")
        lines += [
            "",
            "| Pocket | n | AlphaFold fails | AlphaFold fails while holo succeeds | apo fails |",
            "|---|---|---|---|---|",
        ]
        for label, m in (("moving (> 2 Å)", moving), ("not moving", ~moving)):
            lines.append(
                f"| {label} | {m.sum()} | {af_fail[m].mean():.1%} | "
                f"{struct_fail[m].mean():.1%} | {(~d.apo)[m].mean():.1%} |"
            )
        band = np.digitize(d.site_plddt_mean, edges)
        lines += [
            "",
            "| Pocket-pLDDT band | n | AlphaFold fails | AlphaFold fails while holo succeeds |",
            "|---|---|---|---|",
        ]
        for k, bname in enumerate(b["band"] for b in rule["bands"]["apobind"]):
            m = band == k
            if m.any():
                lines.append(
                    f"| {bname} | {m.sum()} | {af_fail[m].mean():.1%} | "
                    f"{struct_fail[m].mean():.1%} |"
                )
        hs = d.holo.to_numpy()
        fail = af_fail.to_numpy().astype(int)
        site_score = -d.site_plddt_mean.to_numpy()
        auc_all = _auroc(fail, site_score, groups, rng)
        auc_primary = _auroc(fail[hs], site_score[hs], groups[hs], rng)
        auc_protein = _auroc(fail[hs], -d.plddt_mean.to_numpy()[hs], groups[hs], rng)
        auc_oracle = _auroc(fail[hs], d.site_rmsd.to_numpy()[hs], groups[hs], rng)
        lines += [
            "",
            "AUROC for AlphaFold docking failure (95% cluster bootstrap; positives/n):",
            f"- pocket pLDDT, all proteins: {auc_all}",
            f"- **pocket pLDDT, where holo docking succeeded (primary):** {auc_primary}",
            f"- whole-protein pLDDT, where holo docking succeeded: {auc_protein}",
            f"- measured apo–holo pocket RMSD, where holo succeeded (oracle): {auc_oracle}",
            "",
        ]
    text = "\n".join(lines)
    (OUT / "summary.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
