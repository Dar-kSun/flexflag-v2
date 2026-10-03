"""Label test-retest, R5 of docs/plan-v0.4-robustness.md.

    python scripts/09_retest.py [--workers 24]

For each protein, label a second, independent pair: a different holo entry *and* a
different apo entry from the first pair, chosen in the same order as the original
builds. Compares the two labels and scores pocket pLDDT (from the first pair's
pocket) against the retest label. Writes results/robustness/retest.csv and
retest.md.
"""

import argparse
import importlib.util
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score, roc_auc_score

from flexflag.config import LARGE_CHANGE_A, ROOT
from flexflag.data.apoholo import load_pairs, sifts_chain_map, structure_path
from flexflag.dataset import MAX_CANDIDATES
from flexflag.evaluate import RESULTS
from flexflag.labels import site_rmsd

OUT = RESULTS / "robustness"


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def apobind_candidates(first: pd.DataFrame) -> dict[str, list[dict]]:
    build = _script("01_build_dataset")
    pairs, _ = build.prefilter(load_pairs(), sifts_chain_map())
    pairs = pairs.sort_values(
        ["holo_uniprot", "sequence_identity", "sequence_coverage", "apo_resolution", "holo_id"],
        ascending=[True, False, False, True, True],
    )
    used = first.set_index("uniprot")
    out = {}
    for acc, g in pairs.groupby("holo_uniprot"):
        if acc not in used.index:
            continue
        u = used.loc[acc]
        g = g[(g.holo_id != u.holo_id) & (g.apo_id != u.apo_id)]
        out[acc] = g[["holo_id", "holo_chains", "apo_id", "apo_chains"]].to_dict("records")
    return out


def external_candidates(first: pd.DataFrame) -> dict[str, list[dict]]:
    ext = _script("05_external_dataset")
    sifts = ext.sifts_entries()
    holo = [i for i in ext.search_entries("holo_2019", ext.HOLO_TERMS) if i in sifts.index]
    apo = [i for i in ext.search_entries("apo_strict", ext.APO_TERMS) if i in sifts.index]
    holo_df, apo_df = sifts.loc[holo].copy(), sifts.loc[apo].copy()
    wanted = set(first.uniprot)
    holo_df = holo_df[holo_df.uniprot.isin(wanted)]
    apo_df = apo_df[apo_df.uniprot.isin(wanted)]
    info = ext.entry_details(sorted(set(holo_df.index) | set(apo_df.index)))
    for df in (holo_df, apo_df):
        df["resolution"] = [info.get(i, {}).get("resolution") for i in df.index]
    holo_df = holo_df[
        [ext.has_real_ligand(info.get(i, {}).get("ligands", [])) for i in holo_df.index]
    ]
    used = first.set_index("uniprot")
    out = {}
    for acc, u in used.iterrows():
        hs = holo_df[(holo_df.uniprot == acc) & (holo_df.index != u.holo_id)]
        aps = apo_df[(apo_df.uniprot == acc) & (apo_df.index != u.apo_id)]
        aps = list(aps.reset_index(names="pdb").sort_values(["resolution", "pdb"]).itertuples())
        rows = []
        for h in hs.reset_index(names="pdb").sort_values(["resolution", "pdb"]).itertuples():
            a = next((a for a in aps if ext.overlap_ok(h, a)), None)
            if a is not None:
                rows.append(
                    {
                        "holo_id": h.pdb,
                        "holo_chains": h.chains,
                        "apo_id": a.pdb,
                        "apo_chains": a.chains,
                    }
                )
        out[acc] = rows
    return out


def retest_one(acc: str, cands: list[dict]) -> dict:
    for c in cands[:MAX_CANDIDATES]:
        try:
            lab = site_rmsd(
                structure_path(c["holo_id"]),
                c["holo_chains"].split(),
                structure_path(c["apo_id"]),
                c["apo_chains"].split(),
            )
        except Exception:  # any failure: try the next candidate
            continue
        return {
            "uniprot": acc,
            "retest_holo": c["holo_id"],
            "retest_apo": c["apo_id"],
            "retest_rmsd": lab.rmsd,
        }
    return {"uniprot": acc, "retest_rmsd": np.nan, "n_candidates": len(cands)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["Label test-retest (R5). Retest pair: different holo and different apo.", ""]
    rows = []
    for name, ds_path, pk_path, make in (
        ("apobind", RESULTS / "dataset.csv", RESULTS / "pocket.csv", apobind_candidates),
        (
            "external",
            RESULTS / "external" / "dataset.csv",
            RESULTS / "external" / "pocket.csv",
            external_candidates,
        ),
    ):
        first = pd.read_csv(ds_path)
        cands = make(first)
        with ThreadPoolExecutor(args.workers) as pool:
            res = list(
                pool.map(lambda kv: retest_one(*kv), [(a, c) for a, c in cands.items() if c])
            )
        r = pd.DataFrame(res).merge(first[["uniprot", "site_rmsd"]], on="uniprot")
        r = r.merge(pd.read_csv(pk_path)[["uniprot", "site_plddt_mean"]], on="uniprot")
        r = r.dropna(subset=["retest_rmsd", "site_plddt_mean"]).assign(dataset=name)
        rows.append(r)
        y1 = (r.site_rmsd > LARGE_CHANGE_A).astype(int)
        y2 = (r.retest_rmsd > LARGE_CHANGE_A).astype(int)
        rho = spearmanr(r.site_rmsd, r.retest_rmsd).statistic
        kappa = cohen_kappa_score(y1, y2)
        auc_first = roc_auc_score(y1, -r.site_plddt_mean)
        auc_retest = roc_auc_score(y2, -r.site_plddt_mean)
        ceiling = roc_auc_score(y2, r.site_rmsd)
        lines += [
            f"### {name}: {len(r)} proteins with a retest pair",
            "",
            f"- Spearman, first vs retest site RMSD: {rho:.2f}",
            f"- 2 Å label agreement: {np.mean(y1 == y2):.1%}, Cohen's κ = {kappa:.2f}",
            f"- moving in first pair: {y1.sum()}, in retest: {y2.sum()}, "
            f"in both: {(y1 & y2).sum()}",
            f"- AUROC of pocket pLDDT vs first labels (this subset): {auc_first:.3f}",
            f"- AUROC of pocket pLDDT vs **retest** labels: {auc_retest:.3f}",
            f"- AUROC of the first pair's RMSD as a predictor of the retest label "
            f"(label ceiling): {ceiling:.3f}",
            "",
        ]
    pd.concat(rows).to_csv(
        OUT / "retest.csv", index=False, float_format="%.4f", lineterminator="\n"
    )
    (OUT / "retest.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
