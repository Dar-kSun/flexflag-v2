# flexflag

**Before you dock into an AlphaFold structure, check whether the pocket holds still.**
On 833 protein structures released *after* AlphaFold2's training cutoff, pockets in the
lowest fifth by **pocket pLDDT** reshaped by more than 2 Å between empty and
ligand-bound **13× as often** as pockets in the top fifth (24% vs 1.8%). The rule was
fixed on a separate dataset beforehand. **Whole-protein pLDDT, the usual check, is at
chance** (AUROC 0.44 vs 0.76).

![Risk of pocket change by pocket pLDDT band, discovery and external sets](results/pocket_plddt_bands.png)

> **Status: v0.1.0.** Dataset, labels, cluster-aware evaluation, the pLDDT baseline
> comparison, pre-declared external validation and the `flexflag check` CLI are done
> and reproducible. See [Status](#status) for what is scoped out.

## What this answers, and why

Docking screens molecules against one protein structure, increasingly an AlphaFold
prediction. Proteins flex, and binding pockets can reshape when a ligand binds.
AlphaFold returns one frame, and its confidence score (pLDDT) says how sure it is of
*that* frame, not whether the protein has others. Ensemble methods (MD,
AlphaFold2-RAVE, MSA subsampling) can sample other states but are expensive.
**flexflag is the cheap check in front of them:** is this pocket likely to change shape?
It flags risk. It does not predict the other conformation.

## Quickstart

```bash
pip install -e .            # Python ≥ 3.11
flexflag check Q9HWI0 --from-pdb 8evw        # pocket = residues near the ligand in 8EVW
flexflag check P69441 --residues 13,31,35-38  # or give pocket residues (UniProt numbering)
```

```
Q9HWI0  D-alanine--D-alanine ligase A  (346 residues)

Pocket: 30 residues (ATP in 8evw chain A)
Pocket mean pLDDT: 84.5   (whole protein: 89.3)

Pocket-change risk (> 2 Å, empty vs bound): HIGHEST
  band 1 of 5 by pocket pLDDT (< 91.8)
  observed in this band: 21% of 186 APObind pockets,
                         24% of 139 unseen PDB-2019+ pockets
  smoothed estimate: 22%  (base rate about 10%)

Why:
- Pocket pLDDT is in the lowest 40% of pockets studied. Even above 90 ('very
  high' in AlphaFold's usual reading), lower pocket pLDDT marks more change.
- 5 pocket residue(s) have pLDDT < 70.
- When pockets do move, the AlphaFold model resembled the bound (holo) state in
  69% of unseen cases. If you need the empty state, one
  AlphaFold structure is likely the wrong one.
...
```

This example's pocket does move (5.7 Å, 8EVW vs 8EVV). It comes from the evaluation
set, so it illustrates the method rather than being a new prediction. The tool also
misses cases: **adenylate kinase**, the textbook induced-fit enzyme, moves 6 Å but has a
confident pocket (pLDDT 94.7) and scores "moderate". Pocket pLDDT does not see rigid
domains closing over a well-predicted pocket.

Without a pocket, `flexflag check <UniProt>` reports whole-protein pLDDT and declines to
flag, because whole-protein pLDDT carries no signal (below).

## Results

Label: binding-site Cα RMSD between apo and holo after superposing the site, > 2 Å =
"large change" (other thresholds below). Two independent datasets:

| | Source | Proteins | Moving > 2 Å | Evaluation |
|---|---|---|---|---|
| **Discovery** | APObind (PDBbind v2019) | 929 | 101 (10.9%) | 5-fold CV grouped by 30%-identity clusters |
| **External** | PDB-2019+ (built from RCSB + SIFTS) | 833 | 84 (10.1%) | Models frozen on discovery, applied unchanged |

Every holo structure in the external set was released after 2019-01-01, so after
AlphaFold2's training cutoff (2018-04-30) and after PDBbind v2019. All CIs are 95%
cluster bootstraps (clusters resampled, not proteins).

### 1. Whole-protein pLDDT does not predict pocket change

| | Discovery (CV) | External (frozen) |
|---|---|---|
| Whole-protein mean pLDDT, AUROC | 0.52 [0.46, 0.58] | **0.44** [0.38, 0.50] |

A high-confidence AlphaFold model is no less likely to have a pocket that reshapes.
This is the field's default heuristic, and it fails at every threshold tested.
v0.1 tried richer whole-protein features (PAE inter-domain blocks, sequence) and got
AUROC 0.59 in CV. The gain was modest and mostly a protein-family proxy (see
[findings](docs/findings.md)).

### 2. pLDDT at the pocket does, and the result replicates on unseen structures

| Model (2 Å) | Discovery AUROC (CV) | External AUROC (frozen) | External AUPRC |
|---|---|---|---|
| Whole-protein mean pLDDT | 0.521 [0.461, 0.578] | 0.440 [0.379, 0.503] | 0.086 |
| **Pocket mean pLDDT** (one number) | **0.686** [0.630, 0.739] | **0.763** [0.710, 0.811] | 0.237 |
| Pocket model (pocket + whole-protein features, boosting) | 0.716 [0.655, 0.767] | 0.791 [0.736, 0.840] | 0.344 |
| *Confound reference: pocket size alone* | *0.614* | — | — |

- **External gain of pocket pLDDT over whole-protein pLDDT: +0.32 AUROC
  [+0.22, +0.42].**
- **The simple heuristic does most of the work.** The 29-feature pocket model adds
  ~0.03 AUROC over pocket pLDDT alone, with overlapping CIs on both datasets. So the
  CLI ships the one-number rule.
- **The hardest test:** 480 external proteins with **no** discovery-set protein in the
  same 30% sequence cluster. Pocket pLDDT scores AUROC **0.714** [0.640, 0.785];
  whole-protein pLDDT scores 0.481.
- **It is not a pocket-size effect.** Within small, medium and large pockets
  separately, pocket pLDDT still scores AUROC 0.65–0.74 (discovery). Pocket size is
  never a feature. (`results/pocket_metrics.md`)
- **Calibration:** the frozen rule's external ECE is 0.047. Its band rates are
  reported directly rather than relying on the fitted curve.

Threshold sensitivity, pocket pLDDT (frozen), external set:

| Threshold | Positives | Whole-protein pLDDT | Pocket pLDDT |
|---|---|---|---|
| > 1 Å | 202 | 0.569 | 0.733 [0.696, 0.769] |
| > 1.5 Å | 122 | 0.434 | 0.755 [0.713, 0.796] |
| **> 2 Å** | 84 | 0.440 | **0.763** [0.710, 0.811] |
| > 3 Å | 35 | 0.449 | 0.707 [0.618, 0.790] |

Risk by pocket-pLDDT band (edges fixed on the discovery set; the figure above):

| Pocket pLDDT | Discovery: moving > 2 Å | External: moving > 2 Å |
|---|---|---|
| < 91.8 | 21.0% (39/186) | 23.7% (33/139) |
| 91.8–94.6 | 13.4% (25/186) | 18.5% (27/146) |
| 94.6–96.4 | 10.3% (19/185) | 9.2% (15/163) |
| 96.4–97.9 | 7.0% (13/186) | 2.8% (6/215) |
| > 97.9 | 2.7% (5/186) | 1.8% (3/170) |

Nearly all these pockets are "very high confidence" by the usual reading of pLDDT (> 90).
**Within that range, the exact value still matters.**

### 3. When a pocket moves, AlphaFold usually hands you the bound shape, confidently

For each protein, the AlphaFold DB model's pocket was compared with both crystal
structures (same residues, site-superposed Cα RMSD):

| Moving pockets (> 2 Å) | n | AlphaFold closer to holo | Median AF–holo | Median AF–apo | Pocket pLDDT > 90 |
|---|---|---|---|---|---|
| Discovery | 101 | 80% | 0.60 Å | 2.47 Å | 71% |
| **External** (holo never seen by AlphaFold2) | 84 | **69%** | 0.92 Å | 2.54 Å | 69% |

AlphaFold's preference for the bound-like pocket holds on holo structures it cannot
have memorised. It is weaker there (69% vs 80%), so memorisation may explain part of
the discovery-set figure. It holds whether the apo structure is old (70%, n = 40) or
new (68%, n = 44). In practice: docking hits found in an AlphaFold pocket may look
better than the empty protein supports. If you need the empty or cryptic state, one
AlphaFold model is likely the wrong structure, and AlphaFold's confidence won't warn you.

### How these analyses were kept honest

The pocket analyses and the external validation were **written down and committed
before they were run**: [`docs/plan-v0.2-pocket.md`](docs/plan-v0.2-pocket.md)
(commit `a301a49`) and [`docs/plan-v0.3-external.md`](docs/plan-v0.3-external.md)
(commit `ad7767c`). The git history shows each plan before its results. One check was
added afterwards and is labelled as such: splitting the AlphaFold-state analysis by
release date. The v0.1 whole-protein result, a mostly negative one, is kept as it was.

## How the label is defined

For each apo–holo pair (`flexflag/labels.py`, tested in `tests/test_labels.py`):

1. Take the first model, drop hydrogens and waters, and keep the first altloc.
2. **Ligand:** the non-polymer, non-additive residue (≥ 6 heavy atoms) with the most
   heavy atoms within 5 Å of the protein. Buffers, glycans and ions are excluded;
   peptide ligands are out of scope.
3. **Pocket:** residues with any heavy atom within **5 Å** of the ligand.
4. Match apo and holo residues by sequence alignment (handles renumbering and chain
   renaming, which is tested).
5. **Label:** Cα RMSD of the pocket after superposing the pocket itself, so it measures
   deformation rather than whole-protein motion. Large change = > 2 Å.

Checks: adenylate kinase scores 6.1 Å and trypsin + benzamidine 0.16 Å. On APObind, our
pockets fall 96–100% inside APObind's own site lists.

**What the pocket features may use:** only the pocket's residue *positions*, as at
docking time. At docking you know where you are docking. In this evaluation those
positions come from the holo ligand. No holo coordinates, ligand identity, ligand size
or pocket size enter any feature. Pocket pLDDT is read from AlphaFold DB.
`tests/test_no_holo_leakage.py` checks that the whole-protein pipeline never reads a PDB
structure.

## Data

- **Discovery: [APObind](https://github.com/devalab/Apobind).** 12,267 pairs →
  quality filters (APObind TM-score ≥ 0.5, identity and coverage ≥ 0.9, crystal apo,
  same UniProt on both sides) → 1,183 proteins → 933 labelled → 929 with the pocket
  mapped onto AlphaFold. One pair per protein, so heavily studied proteins don't
  dominate.
- **External: PDB-2019+** (`scripts/05_external_dataset.py`). Holo: X-ray ≤ 2.5 Å, one
  protein entity, a ligand of ≥ 150 Da, released ≥ 2019-01-01 (43,016 search hits;
  40,469 map to one UniProt in SIFTS). Apo: X-ray ≤ 2.5 Å, one protein entity, no
  non-polymer entities at all (18,190 hits; 16,821 mapped). Pairs
  need the same single UniProt and ≥ 90% overlapping UniProt ranges → 971 proteins →
  833 labelled.
- Every dropped pair and protein is listed with its reason: `results/dropped.csv`,
  `results/external/dropped.csv`.
- Structures come from RCSB, UniProt mapping from SIFTS, and predictions and confidence
  from AlphaFold DB (v6 files).

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"   # Windows: .venv\Scripts\
pytest                                  # 18 tests, offline, a few seconds

python scripts/00_feasibility.py        # 10 APObind pairs end to end, with timings
python scripts/01_build_dataset.py      # discovery set, ~15 min, ~1.3 GB cache in data/cache/
python scripts/02_train_eval.py         # v0.1 whole-protein comparison (needs MMseqs2)
python scripts/03_pocket.py             # AlphaFold vs apo/holo at the pocket; pocket features
python scripts/04_pocket_model.py       # pocket-level comparison, discovery CV
python scripts/05_external_dataset.py   # external PDB-2019+ set, ~10 min
python scripts/06_external_eval.py      # frozen rule on the external set + band figure
python scripts/07_export_rule.py        # writes flexflag/rule.json for the CLI
```

`data/apobind_pairs.csv` (the APObind pair list) is committed. Regenerating it from
APObind's `apobind_all.csv` needs one manual browser download (see
[data-choice](docs/data-choice.md)). Results are deterministic (seed 0), and the
committed `results/` match the numbers above.

**MMseqs2 on Windows:** unzip the `mmseqs-win64` release to a path without spaces
(default looked up: `%LOCALAPPDATA%\flexflag\mmseqs`). Then run
`mmseqs\bin\busybox.exe --install <that bin folder>` once; it installs the helpers as
hard links, which needs no admin rights. Or set `FLEXFLAG_MMSEQS`.

## Limitations

- **Moderate discrimination.** AUROC 0.76 means a useful prior, not a verdict. Even in
  the riskiest band, three in four pockets do *not* move more than 2 Å.
- **Misses confident domain closures.** Adenylate kinase is the example. Hinge motions
  that close rigid, well-predicted domains over a pocket are invisible to pLDDT.
- **The pocket must be supplied.** The flag is per pocket, not per protein. A pocket
  that only forms on binding (cryptic) may be hard to specify.
- **Datasets.** Both lean toward well-studied, crystallisable drug targets (30–40%
  human). Viral polyproteins (no AlphaFold DB model), proteins over 1,500 residues and
  peptide ligands are excluded. The external set's "strict apo" rule (no ligands, ions
  or buffers at all) is conservative, and its holo ligands include cofactors.
- **One pair per protein.** A different ligand may move the same pocket differently.
- **The label is backbone-only (Cα).** Side-chain rearrangements, which also break
  docking, are not captured. Larger pockets score somewhat larger RMSDs; the band
  results hold within pocket-size tertiles.
- **Crystal artefacts.** Crystal packing, different constructs and resolution
  differences can create or hide apparent motion.
- **Clustering at 30% identity** is the standard but does not remove all remote
  homology.
- This repo flags risk. **It does not predict conformations, run docking, or make any
  claim about drug efficacy or clinical outcomes.**

## Status

**v0.1.0.** Done: both datasets, labels, cluster-aware CV with a leakage assertion,
pLDDT baselines, the pocket-level analysis, pre-declared external validation,
calibration, threshold sensitivity, and the CLI.

Next:

- Pockets predicted on the AlphaFold model itself (e.g. fpocket or P2Rank), so the tool
  needs no ligand at all.
- Side-chain-aware labels, and labels aggregated over several ligands per protein.
- SHAP explanations for the pocket model, and family-level analysis.
- Features aimed at hinge closures (the adenylate kinase failure mode).

## Citations

- Jumper et al. (2021). Highly accurate protein structure prediction with AlphaFold.
  *Nature* 596, 583–589.
- Varadi et al. (2024). AlphaFold Protein Structure Database in 2024. *Nucleic Acids
  Research* 52, D368–D375.
- Aggarwal et al. (2021). APObind: a dataset of ligand unbound protein conformations
  for machine learning applications in de novo drug design. ICML Workshop on
  Computational Biology; arXiv:2108.09926.
- Liu et al. (2017). Forging the basis for developing protein–ligand interaction
  scoring functions (PDBbind). *Accounts of Chemical Research* 50, 302–309.
- Dana et al. (2019). SIFTS: updated Structure Integration with Function, Taxonomy and
  Sequences resource. *Nucleic Acids Research* 47, D482–D489.
- Burley et al. (2023). RCSB Protein Data Bank. *Nucleic Acids Research* 51,
  D488–D508.
- Steinegger & Söding (2017). MMseqs2 enables sensitive protein sequence searching for
  the analysis of massive data sets. *Nature Biotechnology* 35, 1026–1028.
- Vani, Aranganathan, Wang & Tiwary (2023). AlphaFold2-RAVE: from sequence to Boltzmann
  ranking. *Journal of Chemical Theory and Computation* 19, 4351–4354.
- del Alamo et al. (2022). Sampling alternative conformational states of transporters
  and receptors with AlphaFold2. *eLife* 11, e75751.
- Buttenschoen, Morris & Deane (2024). PoseBusters: AI-based docking methods fail to
  generate physically valid poses or generalise to novel sequences. *Chemical Science*
  15, 3130–3139.

MIT licence. Built by Aryan Sinha (IIT Kharagpur).
