"""Ligand-free mode, R6 of docs/plan-v0.4-robustness.md: pockets predicted by P2Rank
on the AlphaFold DB model, then the frozen pocket-pLDDT rule.

    python scripts/10_p2rank.py [--threads 16]

Needs Java 17+ and P2Rank 2.5 (see README). Looks for them via FLEXFLAG_JAVA and
FLEXFLAG_P2RANK, else under %LOCALAPPDATA%/flexflag/. Writes
results/robustness/p2rank.csv and p2rank.md.
"""

import argparse
import importlib.util
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from flexflag.config import CACHE_DIR, LARGE_CHANGE_A, ROOT
from flexflag.data import alphafold
from flexflag.data.apoholo import structure_path
from flexflag.evaluate import RESULTS, SEED
from flexflag.labels import align_indices, one_letter, protein_residues, read_model, site_rmsd

OUT = RESULTS / "robustness"
P2_OUT = CACHE_DIR / "p2rank"
LOCAL = Path(os.environ.get("LOCALAPPDATA", "")) / "flexflag"
TOP_K = 3
MIN_COVERAGE = 0.5


def tool_paths() -> tuple[str, Path]:
    java = os.environ.get("FLEXFLAG_JAVA") or next(
        (str(p / "bin" / "java.exe") for p in sorted(LOCAL.glob("jdk-2*"))), "java"
    )
    p2rank = Path(os.environ.get("FLEXFLAG_P2RANK") or LOCAL / "p2rank_2.5.1")
    return java, p2rank


def run_p2rank(models: list[Path], threads: int) -> None:
    """One batch run over all models; skipped if every prediction already exists."""
    todo = [m for m in models if not (P2_OUT / f"{m.name}_predictions.csv").exists()]
    if not todo:
        return
    java, p2rank = tool_paths()
    model_dir = models[0].parent
    # File names relative to the .ds file: P2Rank splits dataset lines on whitespace,
    # and the repo path contains a space.
    ds = model_dir / "p2rank_todo.ds"
    ds.write_text("\n".join(m.name for m in todo) + "\n", newline="\n")
    classpath = f"{p2rank / 'bin' / 'p2rank.jar'};{p2rank / 'bin' / 'lib' / '*'}"
    subprocess.run(
        [
            java,
            "-Xmx8g",
            "-cp",
            classpath,
            "cz.siret.prank.program.Main",
            "predict",
            "-c",
            "alphafold",
            ds.name,
            "-o",
            str(P2_OUT),
            "-threads",
            str(threads),
            "-visualizations",
            "0",
        ],
        check=True,
        cwd=model_dir,
    )


def predicted_pockets(model: Path) -> list[set[int]]:
    df = pd.read_csv(P2_OUT / f"{model.name}_predictions.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    return [
        {int(r.split("_")[1]) for r in str(ids).split()}
        for ids in df.sort_values("rank").residue_ids.head(TOP_K)
    ]


def true_pocket(row) -> set[int]:
    """The 5 Å pocket mapped onto AlphaFold residue numbers."""
    lab = site_rmsd(
        structure_path(row.holo_id), [row.holo_chain], structure_path(row.apo_id), [row.apo_chain]
    )
    holo_res = protein_residues(read_model(structure_path(row.holo_id)).find_chain(row.holo_chain))
    af_res = protein_residues(read_model(alphafold.model_path(row.uniprot))[0])
    to_af = align_indices(one_letter(holo_res), one_letter(af_res))
    return {af_res[to_af[i]].seqid.num for i in lab.site_indices if i in to_af}


def _robustness():
    spec = importlib.util.spec_from_file_location("rob", ROOT / "scripts" / "08_robustness.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sets = {
        "apobind": (RESULTS / "dataset.csv", RESULTS / "clusters.csv"),
        "external": (RESULTS / "external" / "dataset.csv", RESULTS / "external" / "clusters.csv"),
    }
    data = {n: pd.read_csv(d) for n, (d, _) in sets.items()}
    uniprots = sorted(set().union(*(set(d.uniprot) for d in data.values())))
    with ThreadPoolExecutor(16) as pool:
        models = dict(zip(uniprots, pool.map(alphafold.model_path, uniprots), strict=True))
    run_p2rank(sorted(set(models.values())), args.threads)

    def per_protein(row):
        try:
            truth = true_pocket(row)
            preds = predicted_pockets(models[row.uniprot])
        except Exception as e:  # keep going; counted below
            return {"uniprot": row.uniprot, "error": f"{type(e).__name__}: {e}"[:120]}
        plddt = alphafold.plddt(row.uniprot)
        cover = [len(p & truth) / len(truth) for p in preds]

        def mean_plddt(res):
            return float(np.mean([plddt[n - 1] for n in res])) if res else np.nan

        best = int(np.argmax(cover)) if cover else None
        return {
            "uniprot": row.uniprot,
            "site_rmsd": row.site_rmsd,
            "n_true": len(truth),
            "n_predicted": len(preds),
            "top1_coverage": cover[0] if cover else 0.0,
            "best_top3_coverage": max(cover) if cover else 0.0,
            "true_pocket_plddt": mean_plddt(truth),
            "top1_pocket_plddt": mean_plddt(preds[0]) if preds else np.nan,
            "best_top3_pocket_plddt": mean_plddt(preds[best]) if preds else np.nan,
        }

    rob = _robustness()
    rng = np.random.default_rng(SEED)
    lines = ["Ligand-free mode (R6): P2Rank 2.5.1 (`-c alphafold`) on AlphaFold DB models.", ""]
    frames = []
    for name, (_, cl_path) in sets.items():
        with ThreadPoolExecutor(16) as pool:
            df = pd.DataFrame(list(pool.map(per_protein, data[name].itertuples())))
        df["dataset"] = name
        frames.append(df)
        ok = df.dropna(subset=["true_pocket_plddt"])
        groups = pd.read_csv(cl_path, index_col="uniprot").cluster.loc[ok.uniprot].to_numpy()
        y = (ok.site_rmsd.to_numpy() > LARGE_CHANGE_A).astype(int)

        def row(label, col, ok=ok, y=y, groups=groups):
            r = rob.auroc_ci(y, -ok[col].to_numpy(), groups, rng)
            return (
                f"| {label} | {r['auroc']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}] "
                f"| {r['n_positive']}/{r['n']} |"
            )

        recovered = (ok.best_top3_coverage >= MIN_COVERAGE).mean()
        top1 = (ok.top1_coverage >= MIN_COVERAGE).mean()
        lines += [
            f"### {name}: {len(ok)} proteins ({df.get('error', pd.Series()).notna().sum()} failed)",
            "",
            f"- True pocket recovered (≥ {MIN_COVERAGE:.0%} of its residues) by a top-{TOP_K} "
            f"P2Rank pocket: {recovered:.1%}; by the top-1 pocket: {top1:.1%}",
            "",
            "| Pocket used for the rule | AUROC (> 2 Å) | positives |",
            "|---|---|---|",
            row("True pocket (ligand-derived), reference", "true_pocket_plddt"),
            row(f"Best-overlapping top-{TOP_K} P2Rank pocket", "best_top3_pocket_plddt"),
            row("Top-1 P2Rank pocket, no ligand information at all", "top1_pocket_plddt"),
            "",
        ]
    pd.concat(frames).to_csv(
        OUT / "p2rank.csv", index=False, float_format="%.4f", lineterminator="\n"
    )
    (OUT / "p2rank.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
