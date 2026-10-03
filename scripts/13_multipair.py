"""E2 of docs/plan-v0.5-overnight.md: many ligands per protein, apo held fixed.

    python scripts/13_multipair.py [--workers 24] [--limit N]

Labels up to MAX_EXTRA further holo entries per protein against the first pair's
apo structure. Writes results/multipair/pairs.csv (one row per labelled pair),
proteins.csv (per-protein summary) and summary.md. Re-running resumes from the
download cache.
"""

import argparse
import importlib.util
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from flexflag.config import LARGE_CHANGE_A, ROOT
from flexflag.data.apoholo import load_pairs, sifts_chain_map, structure_path
from flexflag.evaluate import RESULTS
from flexflag.labels import site_rmsd

OUT = RESULTS / "multipair"
MAX_EXTRA = 10


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def apobind_tasks(first: pd.DataFrame) -> list[dict]:
    pairs, _ = _script("01_build_dataset").prefilter(load_pairs(), sifts_chain_map())
    pairs = pairs.sort_values(
        ["holo_uniprot", "sequence_identity", "sequence_coverage", "apo_resolution", "holo_id"],
        ascending=[True, False, False, True, True],
    )
    tasks = []
    for u in first.itertuples():
        g = pairs[
            (pairs.holo_uniprot == u.uniprot)
            & (pairs.apo_id == u.apo_id)
            & (pairs.holo_id != u.holo_id)
        ].drop_duplicates("holo_id")
        for r in g.head(MAX_EXTRA).itertuples():
            tasks.append(
                {
                    "dataset": "apobind",
                    "uniprot": u.uniprot,
                    "holo_id": r.holo_id,
                    "holo_chains": r.holo_chains,
                    "apo_id": u.apo_id,
                    "apo_chain": u.apo_chain,
                }
            )
    return tasks


def external_tasks(first: pd.DataFrame) -> list[dict]:
    ext = _script("05_external_dataset")
    sifts = ext.sifts_entries()
    holo = [i for i in ext.search_entries("holo_2019", ext.HOLO_TERMS) if i in sifts.index]
    holo_df = sifts.loc[holo]
    holo_df = holo_df[holo_df.uniprot.isin(set(first.uniprot))].copy()
    info = ext.entry_details(sorted(holo_df.index))
    holo_df["resolution"] = [info.get(i, {}).get("resolution") for i in holo_df.index]
    holo_df = holo_df[
        [ext.has_real_ligand(info.get(i, {}).get("ligands", [])) for i in holo_df.index]
    ]
    tasks = []
    for u in first.itertuples():
        apo = sifts.loc[u.apo_id]
        hs = holo_df[(holo_df.uniprot == u.uniprot) & (holo_df.index != u.holo_id)]
        hs = hs.reset_index(names="pdb").sort_values(["resolution", "pdb"])
        n = 0
        for h in hs.itertuples():
            if ext.overlap_ok(h, apo):
                tasks.append(
                    {
                        "dataset": "external",
                        "uniprot": u.uniprot,
                        "holo_id": h.pdb,
                        "holo_chains": h.chains,
                        "apo_id": u.apo_id,
                        "apo_chain": u.apo_chain,
                    }
                )
                n += 1
            if n == MAX_EXTRA:
                break
    return tasks


def label(task: dict) -> dict:
    try:
        lab = site_rmsd(
            structure_path(task["holo_id"]),
            task["holo_chains"].split(),
            structure_path(task["apo_id"]),
            [task["apo_chain"]],
        )
        return {**task, "site_rmsd": lab.rmsd, "ligand": lab.ligand, "n_site": lab.n_site}
    except Exception as e:  # recorded, not fatal
        return {**task, "error": f"{type(e).__name__}: {e}"[:150]}


def summarise(pairs: pd.DataFrame) -> list[str]:
    lines = ["Many ligands per protein (E2): apo fixed, up to 10 further holo entries.", ""]
    proteins = []
    for name in ("apobind", "external"):
        first = pd.read_csv(
            RESULTS / ("dataset.csv" if name == "apobind" else "external/dataset.csv")
        )
        pk = pd.read_csv(RESULTS / ("pocket.csv" if name == "apobind" else "external/pocket.csv"))
        p = pairs[(pairs.dataset == name) & pairs.site_rmsd.notna()]
        allp = pd.concat([first[["uniprot", "site_rmsd"]], p[["uniprot", "site_rmsd"]]])
        g = allp.groupby("uniprot").site_rmsd
        s = pd.DataFrame(
            {
                "n_ligands": g.size(),
                "frac_moving": g.apply(lambda x: (x > LARGE_CHANGE_A).mean()),
                "max_rmsd": g.max(),
                "median_rmsd": g.median(),
            }
        )
        s = s.join(pk.set_index("uniprot").site_plddt_mean).dropna(subset=["site_plddt_mean"])
        s["dataset"] = name
        proteins.append(s.reset_index())
        multi = s[s.n_ligands >= 3]
        consistent = ((multi.frac_moving == 0) | (multi.frac_moving == 1)).mean()
        any_move = (s.frac_moving > 0).astype(int)
        most_move = (s.frac_moving > 0.5).astype(int)
        rho = spearmanr(-s.site_plddt_mean, s.frac_moving).statistic
        lines += [
            f"### {name}: {len(s)} proteins, {len(p)} extra labelled pairs, "
            f"{(s.n_ligands >= 3).sum()} proteins with ≥ 3 ligands",
            "",
            "- proteins with ≥ 3 ligands that are consistent (all or none move > 2 Å): "
            f"{consistent:.1%}",
            f"- AUROC of (−first-pocket pLDDT) for 'any ligand moves': "
            f"{roc_auc_score(any_move, -s.site_plddt_mean):.3f} ({any_move.sum()} positive)",
            f"- AUROC for 'most ligands move': {roc_auc_score(most_move, -s.site_plddt_mean):.3f} "
            f"({most_move.sum()} positive)",
            f"- Spearman(−pocket pLDDT, fraction of ligands moving): {rho:.2f}",
            "",
        ]
    pd.concat(proteins).to_csv(
        OUT / "proteins.csv", index=False, float_format="%.4f", lineterminator="\n"
    )
    return lines


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--limit", type=int, help="only this many proteins per dataset (testing)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    tasks = []
    for name, make in (("apobind", apobind_tasks), ("external", external_tasks)):
        first = pd.read_csv(
            RESULTS / ("dataset.csv" if name == "apobind" else "external/dataset.csv")
        )
        if args.limit:
            first = first.head(args.limit)
        tasks += make(first)
    print(f"{len(tasks)} pairs to label ({time.perf_counter() - t0:.0f}s)", flush=True)
    rows = []
    with ThreadPoolExecutor(args.workers) as pool:
        for k, fut in enumerate(as_completed([pool.submit(label, t) for t in tasks]), 1):
            rows.append(fut.result())
            if k % 500 == 0:
                print(f"{k}/{len(tasks)} labelled, {time.perf_counter() - t0:.0f}s", flush=True)
    pairs = pd.DataFrame(rows).sort_values(["dataset", "uniprot", "holo_id"])
    if "error" not in pairs:
        pairs["error"] = np.nan
    pairs.to_csv(OUT / "pairs.csv", index=False, float_format="%.4f", lineterminator="\n")
    lines = summarise(pairs)
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    (OUT / "summary.json").write_text(
        json.dumps({"n_tasks": len(tasks), "n_ok": int(pairs.site_rmsd.notna().sum())})
    )
    print("\n".join(lines))
    print(f"done in {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
