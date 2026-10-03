"""External validation, exactly as fixed in docs/plan-v0.3-external.md.

    python scripts/06_external_eval.py

Fits the frozen models once on all APObind proteins, applies them unchanged to the
PDB-2019+ set (results/external/dataset.csv), and reports discrimination, calibration,
risk bands, subsets and the AlphaFold-state analysis. Writes results/external/*.
"""

import json
from concurrent.futures import ThreadPoolExecutor

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from flexflag.config import LARGE_CHANGE_A, LARGE_CHANGE_ALTERNATIVES_A
from flexflag.evaluate import METRICS, RESULTS, SEED, bootstrap, ci, features_table
from flexflag.features.pocket import POCKET_FEATURES
from flexflag.model import gradient_boosting, logistic
from flexflag.pocket_analysis import AF2_PDB_CUTOFF, analyse
from flexflag.splits import cluster

matplotlib.use("Agg")

EXT = RESULTS / "external"
BASELINE = "protein_plddt"
LABELS = {
    "protein_plddt": "Whole-protein mean pLDDT (frozen logistic)",
    "site_plddt": "Pocket mean pLDDT (frozen logistic)",
    "pocket_gb": "Pocket model (frozen boosting)",
}
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e8e7e3"
BLUE, ORANGE = "#2a78d6", "#eb6834"


def load(ds_path, pocket_path, feat_path, compute_pocket=False):
    ds = pd.read_csv(ds_path)
    if compute_pocket:
        if pocket_path.exists():
            pocket = pd.read_csv(pocket_path)
        else:
            with ThreadPoolExecutor(16) as pool:
                pocket = pd.DataFrame(list(pool.map(analyse, ds.itertuples())))
            pocket.to_csv(pocket_path, index=False, float_format="%.4f", lineterminator="\n")
    else:
        pocket = pd.read_csv(pocket_path)
    pocket = pocket.set_index("uniprot")
    ds = ds[ds.uniprot.isin(pocket.dropna(subset=POCKET_FEATURES[:4]).index)].reset_index(drop=True)
    protein = features_table(ds, feat_path).reset_index(drop=True)
    X = pd.concat([protein, pocket.loc[ds.uniprot, POCKET_FEATURES].reset_index(drop=True)], axis=1)
    return ds, X, pocket.loc[ds.uniprot].reset_index()


def band_edges(site_plddt: pd.Series) -> np.ndarray:
    return np.quantile(site_plddt, [0.2, 0.4, 0.6, 0.8])


def band_table(site_plddt, rmsd, edges) -> pd.DataFrame:
    idx = np.digitize(site_plddt, edges)
    names = [f"< {edges[0]:.1f}"] + [
        f"{a:.1f}–{b:.1f}" for a, b in zip(edges, edges[1:], strict=False)
    ]
    names.append(f"> {edges[-1]:.1f}")
    rows = []
    for b, name in enumerate(names):
        sel = idx == b
        rows.append(
            {
                "band": name,
                "n": int(sel.sum()),
                "frac_gt_2A": float(np.mean(rmsd[sel] > LARGE_CHANGE_A)) if sel.any() else np.nan,
                "n_gt_2A": int(np.sum(rmsd[sel] > LARGE_CHANGE_A)),
            }
        )
    return pd.DataFrame(rows)


def band_figure(train: pd.DataFrame, ext: pd.DataFrame, n_train: int, n_ext: int, path):
    fig, ax = plt.subplots(figsize=(7.4, 4.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    x = np.arange(len(train))
    w = 0.38
    for off, df, color, label in (
        (-w / 2 - 0.01, train, BLUE, f"APObind (discovery, n = {n_train})"),
        (w / 2 + 0.01, ext, ORANGE, f"PDB-2019+ (external, n = {n_ext})"),
    ):
        bars = ax.bar(x + off, df.frac_gt_2A * 100, width=w, color=color, label=label)
        for bar, k, n in zip(bars, df.n_gt_2A, df.n, strict=True):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.4,
                f"{k}/{n}",
                ha="center",
                va="bottom",
                fontsize=7,
                color=INK2,
            )
    ax.set_xticks(x)
    ax.set_xticklabels(train.band, color=INK2)
    ax.set_xlabel("Mean pLDDT over the binding-site residues (bands fixed on APObind)", color=INK2)
    ax.set_ylabel("Binding sites moving > 2 Å (%)", color=INK2)
    ax.set_title(
        "Lower pocket pLDDT, more pocket change — replicated on unseen structures",
        color=INK,
        loc="left",
        fontsize=11,
    )
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK2)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, labelcolor=INK2, fontsize=8)
    fig.text(
        0.01,
        0.01,
        "All bands above the first lie above 90, AlphaFold's 'very high confidence' cutoff. "
        "Labels: moving / total.",
        color=MUTED,
        fontsize=7,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def score(y, preds, groups, rng) -> dict:
    draws, diffs = bootstrap(y, preds, groups, rng, BASELINE)
    out = {"n": int(len(y)), "n_positive": int(y.sum()), "models": {}}
    for k, p in preds.items():
        r = {m: {"value": float(f(y, p)), "ci95": ci(draws[m][k])} for m, f in METRICS.items()}
        if k != BASELINE:
            r["minus_protein_plddt"] = {
                m: {"mean": float(np.mean(diffs[m][k])), "ci95": ci(diffs[m][k])} for m in diffs
            }
        out["models"][k] = r
    return out


def agreement(pocket: pd.DataFrame, rmsd: np.ndarray) -> dict:
    ok = pocket.assign(site_rmsd=rmsd).dropna(subset=["rmsd_af_holo"])
    res = {}
    for name, sub in {
        "all": ok,
        "moving": ok[ok.site_rmsd > LARGE_CHANGE_A],
        "moving_apo_before_cutoff": ok[
            (ok.site_rmsd > LARGE_CHANGE_A) & (ok.apo_release <= AF2_PDB_CUTOFF)
        ],
        "moving_apo_after_cutoff": ok[
            (ok.site_rmsd > LARGE_CHANGE_A) & (ok.apo_release > AF2_PDB_CUTOFF)
        ],
    }.items():
        res[name] = {
            "n": len(sub),
            "frac_af_closer_to_holo": float((sub.rmsd_af_holo < sub.rmsd_af_apo).mean())
            if len(sub)
            else None,
            "median_af_holo": float(sub.rmsd_af_holo.median()) if len(sub) else None,
            "median_af_apo": float(sub.rmsd_af_apo.median()) if len(sub) else None,
            "median_apo_holo": float(sub.rmsd_apo_holo_common.median()) if len(sub) else None,
            "frac_site_plddt_gt90": float((sub.site_plddt_mean > 90).mean()) if len(sub) else None,
        }
    return res


def main() -> None:
    EXT.mkdir(parents=True, exist_ok=True)
    tr_ds, tr_X, tr_pk = load(
        RESULTS / "dataset.csv", RESULTS / "pocket.csv", RESULTS / "features.csv"
    )
    ex_ds, ex_X, ex_pk = load(
        EXT / "dataset.csv", EXT / "pocket.csv", EXT / "features.csv", compute_pocket=True
    )
    for name, d in (("APObind", tr_ds), ("external", ex_ds)):
        print(f"{name}: {len(d)} proteins, {(d.site_rmsd > LARGE_CHANGE_A).sum()} move > 2 A")

    # Joint clustering: bootstrap groups for the external set, and subset (c).
    ids = [f"A|{u}" for u in tr_ds.uniprot] + [f"E|{u}" for u in ex_ds.uniprot]
    joint = cluster(ids, list(tr_ds.sequence) + list(ex_ds.sequence))
    ex_groups = joint.loc[[f"E|{u}" for u in ex_ds.uniprot]].to_numpy()
    train_clusters = set(joint.loc[[f"A|{u}" for u in tr_ds.uniprot]])
    subsets = {
        "all": np.ones(len(ex_ds), bool),
        "new_uniprot": ~ex_ds.uniprot.isin(set(tr_ds.uniprot)).to_numpy(),
        "new_cluster": np.array([g not in train_clusters for g in ex_groups]),
    }

    feat_cols = [c for c in tr_X.columns]
    frozen = {
        "protein_plddt": (logistic, ["plddt_mean"]),
        "site_plddt": (logistic, ["site_plddt_mean"]),
        "pocket_gb": (lambda: gradient_boosting(SEED), feat_cols),
    }
    out = {
        "n_train": len(tr_ds),
        "n_external": len(ex_ds),
        "subset_sizes": {k: int(v.sum()) for k, v in subsets.items()},
        "thresholds": {},
    }
    rng = np.random.default_rng(SEED)
    for thr in sorted({LARGE_CHANGE_A, *LARGE_CHANGE_ALTERNATIVES_A}):
        y_tr = (tr_ds.site_rmsd.to_numpy() > thr).astype(int)
        y_ex = (ex_ds.site_rmsd.to_numpy() > thr).astype(int)
        preds = {}
        for name, (make, cols) in frozen.items():
            m = make().fit(tr_X[cols], y_tr)
            preds[name] = m.predict_proba(ex_X[cols])[:, 1]
        out["thresholds"][f"{thr:g}"] = {}
        for sname, mask in subsets.items():
            if min(y_ex[mask].sum(), (1 - y_ex[mask]).sum()) < 5:
                out["thresholds"][f"{thr:g}"][sname] = {"skipped": "fewer than 5 per class"}
                continue
            out["thresholds"][f"{thr:g}"][sname] = score(
                y_ex[mask], {k: p[mask] for k, p in preds.items()}, ex_groups[mask], rng
            )
        if thr == LARGE_CHANGE_A:
            pd.DataFrame(
                {
                    "uniprot": ex_ds.uniprot,
                    "site_rmsd": ex_ds.site_rmsd,
                    **{f"p_{k}": v for k, v in preds.items()},
                }
            ).to_csv(EXT / "predictions.csv", index=False, float_format="%.4f")

    edges = band_edges(tr_pk.site_plddt_mean)
    tr_bands = band_table(tr_pk.site_plddt_mean.to_numpy(), tr_ds.site_rmsd.to_numpy(), edges)
    ex_bands = band_table(ex_pk.site_plddt_mean.to_numpy(), ex_ds.site_rmsd.to_numpy(), edges)
    out["band_edges"] = [float(e) for e in edges]
    out["bands"] = {"apobind": tr_bands.to_dict("records"), "external": ex_bands.to_dict("records")}
    out["alphafold_state"] = {
        "apobind": agreement(tr_pk, tr_ds.site_rmsd.to_numpy()),
        "external": agreement(ex_pk, ex_ds.site_rmsd.to_numpy()),
    }
    band_figure(tr_bands, ex_bands, len(tr_ds), len(ex_ds), RESULTS / "pocket_plddt_bands.png")
    (EXT / "metrics.json").write_text(json.dumps(out, indent=2))
    (EXT / "metrics.md").write_text(markdown(out), encoding="utf-8")
    print(markdown(out))


def markdown(out: dict) -> str:
    def cell(r):
        return f"{r['value']:.3f} [{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]"

    lines = [
        f"Frozen models fit on {out['n_train']} APObind proteins, applied to "
        f"{out['n_external']} PDB-2019+ proteins. Subset sizes: {out['subset_sizes']}. "
        "95% CIs: cluster bootstrap on the external set.",
        "",
    ]
    for thr, subs in out["thresholds"].items():
        for sname, res in subs.items():
            lines.append(f"### > {thr} Å, subset: {sname}")
            if "skipped" in res:
                lines += [res["skipped"], ""]
                continue
            lines += [
                f"{res['n_positive']}/{res['n']} positive.",
                "",
                "| Model | AUROC | AUPRC | ECE | Brier |",
                "|---|---|---|---|---|",
            ]
            for k, r in res["models"].items():
                lines.append(f"| {LABELS[k]} | " + " | ".join(cell(r[m]) for m in METRICS) + " |")
            for k, r in res["models"].items():
                if "minus_protein_plddt" in r:
                    d = r["minus_protein_plddt"]["auroc"]
                    lines.append(
                        f"- {LABELS[k]} − whole-protein pLDDT, AUROC: {d['mean']:+.3f} "
                        f"[{d['ci95'][0]:+.3f}, {d['ci95'][1]:+.3f}]"
                    )
            lines.append("")
    lines += [
        "### Risk by pocket-pLDDT band (edges fixed on APObind)",
        "",
        "| Band | APObind > 2 Å | External > 2 Å |",
        "|---|---|---|",
    ]
    for a, e in zip(out["bands"]["apobind"], out["bands"]["external"], strict=True):
        lines.append(
            f"| {a['band']} | {a['frac_gt_2A']:.1%} ({a['n_gt_2A']}/{a['n']}) | "
            f"{e['frac_gt_2A']:.1%} ({e['n_gt_2A']}/{e['n']}) |"
        )
    lines += [
        "",
        "### Which conformation does AlphaFold return?",
        "",
        "| Set | Subset | n | AF closer to holo | median AF–holo | median AF–apo | "
        "median apo–holo | site pLDDT > 90 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for setname, res in out["alphafold_state"].items():
        for sname, r in res.items():
            if not r["n"]:
                continue
            lines.append(
                f"| {setname} | {sname} | {r['n']} | {r['frac_af_closer_to_holo']:.1%} | "
                f"{r['median_af_holo']:.2f} | {r['median_af_apo']:.2f} | "
                f"{r['median_apo_holo']:.2f} | {r['frac_site_plddt_gt90']:.1%} |"
            )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
