"""Robustness checks R1-R4, as fixed in docs/plan-v0.4-robustness.md.

    python scripts/08_robustness.py

For each protein in both datasets: decoy-pocket pLDDT (R2), pocket pLDDT at other
pocket cutoffs (R3), and two alternative labels (R4). Then AUROCs of (-pLDDT) with
cluster-bootstrap CIs, and a permuted-label control (R1). Writes
results/robustness/{per_protein.csv, summary.md, summary.json}.
"""

import json
import zlib
from concurrent.futures import ThreadPoolExecutor

import gemmi
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from flexflag.config import LARGE_CHANGE_A
from flexflag.data import alphafold
from flexflag.data.apoholo import structure_path
from flexflag.evaluate import RESULTS, SEED
from flexflag.labels import (
    align_indices,
    find_ligand,
    one_letter,
    protein_residues,
    read_model,
    site_residues,
    site_rmsd,
)
from flexflag.splits import cluster

OUT = RESULTS / "robustness"
CUTOFFS = (4.0, 6.0, 8.0)
N_DECOYS = 20
DECOY_MIN_DIST = 15.0
N_BOOT = 1000


def ca(res) -> np.ndarray:
    p = res.find_atom("CA", "*").pos
    return np.array([p.x, p.y, p.z])


def heavy_atoms(res) -> dict[str, gemmi.Position]:
    return {a.name: a.pos for a in res if a.element.name != "H"}


def analyse(row) -> dict:
    holo_p, apo_p = structure_path(row.holo_id), structure_path(row.apo_id)
    lab = site_rmsd(holo_p, [row.holo_chain], apo_p, [row.apo_chain])
    holo = read_model(holo_p)
    ligand, _, _ = find_ligand(holo, [row.holo_chain])
    holo_res = protein_residues(holo.find_chain(row.holo_chain))
    apo_res = protein_residues(read_model(apo_p).find_chain(row.apo_chain))
    af_res = protein_residues(read_model(alphafold.model_path(row.uniprot))[0])
    plddt = alphafold.plddt(row.uniprot)
    holo_seq = one_letter(holo_res)
    to_apo = align_indices(holo_seq, one_letter(apo_res))
    to_af = align_indices(holo_seq, one_letter(af_res))
    out = {"uniprot": row.uniprot, "site_rmsd": lab.rmsd}

    # R3: pocket pLDDT with the pocket defined at other cutoffs (label unchanged).
    for c in CUTOFFS:
        idx = [to_af[i] for i in site_residues(holo_res, ligand, c) if i in to_af]
        out[f"pocket_plddt_{c:g}A"] = float(plddt[idx].mean()) if len(idx) >= 3 else np.nan
    site_af = sorted({to_af[i] for i in lab.site_indices if i in to_af})
    out["pocket_plddt_5A"] = float(plddt[site_af].mean()) if len(site_af) >= 3 else np.nan

    # R2: decoy pockets of the same size, far from the real pocket, on the AF model.
    xyz = np.array([ca(r) for r in af_res])
    if len(site_af) >= 3:
        d_to_site = np.linalg.norm(xyz[:, None] - xyz[site_af][None], axis=-1).min(axis=1)
        seeds = np.flatnonzero(d_to_site > DECOY_MIN_DIST)
        if len(seeds):
            rng = np.random.default_rng(zlib.crc32(row.uniprot.encode()) + SEED)
            vals = []
            for s in rng.choice(seeds, size=min(N_DECOYS, len(seeds)), replace=False):
                patch = np.argsort(np.linalg.norm(xyz - xyz[s], axis=1))[: len(site_af)]
                vals.append(plddt[patch].mean())
            out["decoy_plddt"] = float(np.mean(vals))
            out["n_decoy_seeds"] = int(len(seeds))

    # R4: alternative labels on the same pocket residues.
    common = [i for i in lab.site_indices if i in to_apo]
    fixed = [holo_res[i].find_atom("CA", "*").pos for i in common]
    moving = [apo_res[to_apo[i]].find_atom("CA", "*").pos for i in common]
    sup = gemmi.superpose_positions(fixed, moving)
    sq = []
    for i in common:  # (a) all heavy atoms, site-superposed on C-alpha
        h, a = heavy_atoms(holo_res[i]), heavy_atoms(apo_res[to_apo[i]])
        for name in h.keys() & a.keys():
            sq.append(h[name].dist(gemmi.Position(sup.transform.apply(a[name]))) ** 2)
    out["label_heavy_atom"] = float(np.sqrt(np.mean(sq)))
    all_common = sorted(to_apo)  # (b) whole-chain superposition, RMSD on the site
    sup_all = gemmi.superpose_positions(
        [holo_res[i].find_atom("CA", "*").pos for i in all_common],
        [apo_res[to_apo[i]].find_atom("CA", "*").pos for i in all_common],
    )
    d2 = [
        holo_res[i]
        .find_atom("CA", "*")
        .pos.dist(
            gemmi.Position(sup_all.transform.apply(apo_res[to_apo[i]].find_atom("CA", "*").pos))
        )
        ** 2
        for i in common
    ]
    out["label_global_sup"] = float(np.sqrt(np.mean(d2)))
    return out


def per_protein(ds: pd.DataFrame, name: str) -> pd.DataFrame:
    def safe(row):
        try:
            return analyse(row)
        except Exception as e:  # keep going; failures are counted in the summary
            return {"uniprot": row.uniprot, "error": f"{type(e).__name__}: {e}"[:120]}

    with ThreadPoolExecutor(16) as pool:
        rows = list(pool.map(safe, ds.itertuples()))
    df = pd.DataFrame(rows).assign(dataset=name)
    bad = df.site_rmsd.notna() & (abs(df.site_rmsd - ds.site_rmsd.to_numpy()) > 1e-3)
    if bad.any():
        raise AssertionError(f"{name}: {bad.sum()} labels did not reproduce")
    return df


def auroc_ci(y, score, groups, rng) -> dict:
    ok = ~np.isnan(score)
    y, score, groups = y[ok], score[ok], groups[ok]
    uniq = np.unique(groups)
    members = {g: np.flatnonzero(groups == g) for g in uniq}
    draws = []
    for _ in range(N_BOOT):
        idx = np.concatenate([members[g] for g in rng.choice(uniq, size=len(uniq))])
        if 0 < y[idx].sum() < len(idx):
            draws.append(roc_auc_score(y[idx], score[idx]))
    return {
        "auroc": float(roc_auc_score(y, score)),
        "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
        "n": int(len(y)),
        "n_positive": int(y.sum()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = {
        "apobind": (RESULTS / "dataset.csv", RESULTS / "clusters.csv"),
        "external": (RESULTS / "external" / "dataset.csv", RESULTS / "external" / "clusters.csv"),
    }
    path = OUT / "per_protein.csv"
    if path.exists():
        pp = pd.read_csv(path)
    else:
        pp = pd.concat([per_protein(pd.read_csv(d), n) for n, (d, _) in sets.items()])
        pp.to_csv(path, index=False, float_format="%.4f", lineterminator="\n")

    rng = np.random.default_rng(SEED)
    summary, lines = (
        {},
        [
            "Robustness checks (docs/plan-v0.4-robustness.md). AUROC of "
            "(-pLDDT), 95% cluster-bootstrap CI.",
            "",
        ],
    )
    for name, (ds_path, cl_path) in sets.items():
        ds = pd.read_csv(ds_path)
        if not cl_path.exists():  # external clusters (30%, external set alone)
            c = cluster(list(ds.uniprot), list(ds.sequence))
            c.rename("cluster").rename_axis("uniprot").to_frame().to_csv(cl_path)
        clusters = pd.read_csv(cl_path, index_col="uniprot").cluster
        d = pp[pp.dataset == name].dropna(subset=["pocket_plddt_5A"]).copy()
        groups = clusters.loc[d.uniprot].to_numpy()
        y = (d.site_rmsd.to_numpy() > LARGE_CHANGE_A).astype(int)
        res = {"n_failed": int(pp[pp.dataset == name].get("error", pd.Series()).notna().sum())}
        res["pocket_5A"] = auroc_ci(y, -d.pocket_plddt_5A.to_numpy(), groups, rng)
        # R1: permuted labels.
        perm = [roc_auc_score(rng.permutation(y), -d.pocket_plddt_5A) for _ in range(200)]
        res["R1_permuted"] = {
            "mean": float(np.mean(perm)),
            "range95": [float(np.percentile(perm, 2.5)), float(np.percentile(perm, 97.5))],
        }
        # R2: decoy pockets.
        res["R2_decoy"] = auroc_ci(y, -d.decoy_plddt.to_numpy(), groups, rng)
        # R3: other pocket cutoffs.
        for c in CUTOFFS:
            res[f"R3_pocket_{c:g}A"] = auroc_ci(
                y, -d[f"pocket_plddt_{c:g}A"].to_numpy(), groups, rng
            )
        # R4: alternative labels.
        for lab in ("label_heavy_atom", "label_global_sup"):
            y2 = (d[lab].to_numpy() > LARGE_CHANGE_A).astype(int)
            res[f"R4_{lab}_2A"] = auroc_ci(y2, -d.pocket_plddt_5A.to_numpy(), groups, rng)
        # (a) at the cutoff matching the primary label's positive rate.
        cut = float(np.quantile(d.label_heavy_atom, 1 - y.mean()))
        y3 = (d.label_heavy_atom.to_numpy() > cut).astype(int)
        res["R4_label_heavy_atom_matched_rate"] = {
            **auroc_ci(y3, -d.pocket_plddt_5A.to_numpy(), groups, rng),
            "cutoff_A": cut,
        }
        res["spearman_primary_vs_heavy_atom"] = float(
            d.site_rmsd.corr(d.label_heavy_atom, method="spearman")
        )
        res["spearman_primary_vs_global_sup"] = float(
            d.site_rmsd.corr(d.label_global_sup, method="spearman")
        )
        summary[name] = res

        def row(label, r):
            return (
                f"| {label} | {r['auroc']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}] | "
                f"{r['n_positive']}/{r['n']} |"
            )

        lines += [
            f"### {name} (failed: {res['n_failed']})",
            "",
            "| Check | AUROC | positives |",
            "|---|---|---|",
            row("Pocket pLDDT, 5 Å pocket (reference)", res["pocket_5A"]),
            f"| R1 permuted labels (200×) | mean {res['R1_permuted']['mean']:.3f}, 95% range "
            f"[{res['R1_permuted']['range95'][0]:.3f}, {res['R1_permuted']['range95'][1]:.3f}] | |",
            row("R2 decoy patch, same size, > 15 Å away", res["R2_decoy"]),
        ]
        lines += [row(f"R3 pocket defined at {c:g} Å", res[f"R3_pocket_{c:g}A"]) for c in CUTOFFS]
        lines += [
            row("R4a label: all-heavy-atom RMSD > 2 Å", res["R4_label_heavy_atom_2A"]),
            row(
                f"R4a label: all-heavy-atom RMSD > {cut:.2f} Å (matched rate)",
                res["R4_label_heavy_atom_matched_rate"],
            ),
            row(
                "R4b label: whole-chain superposition, site RMSD > 2 Å",
                res["R4_label_global_sup_2A"],
            ),
            "",
            "Spearman with the primary label: "
            f"heavy-atom {res['spearman_primary_vs_heavy_atom']:.2f}, "
            f"whole-chain superposition {res['spearman_primary_vs_global_sup']:.2f}.",
            "",
        ]
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
