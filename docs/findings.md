# Findings

Everything here comes from the committed `results/` files, produced by
`scripts/01_build_dataset.py` and `scripts/02_train_eval.py` (seed 0).

## M3: dataset

### From pairs to proteins

| Step | Count |
|---|---|
| APObind pairs | 12,267 |
| After pair-level filters | 7,687 |
| Distinct proteins (holo UniProt) | 1,183 |
| **Labelled proteins** | **933** |
| Sequence clusters (MMseqs2, 30% identity, 50% coverage) | 680 |

Pair-level filters, applied in order using APObind's own metrics and SIFTS (each pair
counted once, at the first filter it fails):

| Reason | Pairs |
|---|---|
| APObind TM-score < 0.5 (likely wrong chain or protein) | 1,957 |
| Apo and holo map to different UniProt accessions | 1,074 |
| Sequence identity < 0.9 | 557 |
| Holo chains don't map to exactly one UniProt | 456 |
| Apo has no resolution (NMR or unknown) | 438 |
| Sequence coverage < 0.9 | 98 |

Protein-level drops (1,183 → 933):

| Reason | Proteins |
|---|---|
| No candidate pair could be labelled (up to 3 tried per protein) | 151 |
| No AlphaFold DB F1 model (mostly viral polyproteins) | 67 |
| AlphaFold model longer than 1,500 residues (skipped) | 32 |

Why individual pairs failed to label (several can belong to one protein):

| Reason | Pairs |
|---|---|
| Ligand is a peptide, not a small molecule (out of scope) | 251 |
| Fewer than 80% of site residues resolved in the apo structure | 25 |
| Binding site under 3 residues | 1 |

The set leans heavily toward drug-discovery targets: 40% human (370), then *E. coli*
(59), mouse (31), rat (23), *P. aeruginosa* (21) and *M. tuberculosis* (21).

### Label distribution

Binding-site Cα RMSD (site-superposed), 933 proteins:

| Quantile | 10% | 25% | 50% | 75% | 90% | 95% | 99% |
|---|---|---|---|---|---|---|---|
| Å | 0.15 | 0.23 | 0.43 | 0.94 | 2.08 | 2.72 | 4.93 |

| Threshold | Positives |
|---|---|
| > 1 Å | 218 (23.4%) |
| > 1.5 Å | 145 (15.5%) |
| **> 2 Å** | **102 (10.9%)** |
| > 3 Å | 38 (4.1%) |

Most binding sites barely move. The top of the distribution is plausible:

- adenylate kinase (5.7 Å, found independently of the test fixture),
- metallo-β-lactamase (8.6 Å; its active-site loop is known to be flexible),
- BTK kinase (5.1 Å),
- EPSP synthase (7.0 Å; known domain closure on substrate binding).

Kinases and transferases make up the largest named groups above 2 Å.

**Not yet checked by hand:** the single largest value, sirtuin-1 (5BTR vs 4ZZH,
13.9 Å). The holo ligand is a Fluor-de-Lys peptide–fluorophore substrate (`FDL`).
That could be real domain closure or an artefact of the pairing. One protein cannot
move the results much, but it should be inspected.

## M5–M6: baseline vs model

Headline numbers and figures are in the README. Full tables are in
`results/metrics.md` and `results/metrics.json`.

### Which way each feature points

In-sample AUROC of each feature on its own at 2 Å, from
`results/univariate_auroc.csv`. Below 0.5 means a higher value goes with *less*
change. This shows direction, not out-of-fold skill.

| Feature | AUROC @ 2 Å |
|---|---|
| Cysteine fraction | 0.39 |
| Polar fraction | 0.42 |
| pLDDT fraction < 50 | 0.43 |
| Longest low-pLDDT run | 0.43 |
| pLDDT std | 0.44 |
| pLDDT fraction < 70 | 0.44 |
| PAE fraction > 15 Å | 0.44 |
| Sequence length | 0.45 |
| Inter-domain PAE (mean / max) | 0.46 / 0.45 (only 27% of proteins have ≥ 2 PAE domains) |
| PAE domain count | 0.46 |
| Mean pLDDT | 0.53 |
| Minimum pLDDT | 0.57 |
| Hydrophobic fraction | 0.58 |

Three things stand out:

1. **AlphaFold's "flexibility" signals point the wrong way, weakly.** More
   low-confidence residues, more PAE domains and higher inter-domain PAE all go with
   slightly *less* binding-site change. A plausible reading: low pLDDT and high PAE
   mark disordered tails and floppy linkers, which are mostly *not* where ligands
   bind. The pocket itself sits in a well-predicted core, and whether that core
   reshapes on binding is invisible to the confidence scores. This is the gap the
   project set out to probe: **confidence in a shape is not evidence of a single
   shape.**
2. **The PAE hinge features did not deliver.** Only 27% of proteins split into ≥ 2
   confident PAE domains at the 10 Å cut, and among those, inter-domain PAE carries no
   useful signal.
3. **Composition carries most of the signal, and it is probably a family proxy.**
   Cysteine-rich and polar proteins (often secreted, disulfide-stabilised enzymes such
   as serine proteases) are rigid. More hydrophobic proteins (often cytoplasmic enzymes
   and kinases) move more. This is a real, cluster-robust pattern, but it says "what
   kind of protein is this," not "this protein is flexible." It may not transfer to
   families absent from PDBbind.

### Label caveat: site size

Larger binding sites score larger RMSDs (site residue count alone gives AUROC 0.61 on
the 2 Å label). Site size comes from the holo ligand, so it is **not** a feature (it
would be leakage). But it means the label partly measures ligand size, and the next
label version should consider a size-normalised RMSD.

## What would most likely improve this

- **Pocket-level features from the AlphaFold model alone:** predict pockets on the
  AlphaFold structure (e.g. fpocket), then take local pLDDT and PAE around them.
  Today's whole-protein summaries dilute any local signal.
- **Side-chain-aware labels:** many docking failures are side-chain rearrangements,
  which a Cα label misses.
- **More pairs per protein:** label each protein by its maximum or median change
  across several ligands, rather than one pair.

---

## v0.2: pocket level (pre-declared in `docs/plan-v0.2-pocket.md`, commit a301a49)

All numbers from `scripts/03_pocket.py` and `scripts/04_pocket_model.py`.
929 of 933 proteins had their pocket map onto the AlphaFold model.

**A. Which conformation does AlphaFold return?** For moving pockets (> 2 Å, n = 101),
the AlphaFold model is closer to holo than to apo in 80% of cases (median site RMSD
0.60 Å to holo, 2.47 Å to apo). Across all proteins it is 57%, as expected when most
pockets barely move.

**C. Confident but moving:** 71% of moving pockets have pocket pLDDT > 90, and 97%
are > 70. For comparison, 89% of non-moving pockets are > 90.

**B. Pocket-conditioned flag**, 5-fold cluster CV:

| > 2 Å | AUROC |
|---|---|
| Pocket mean pLDDT | 0.686 [0.630, 0.739] |
| Pocket model (6 pocket + 23 whole-protein features) | 0.716 [0.655, 0.767] |
| Pocket size alone (confound reference) | 0.614 [0.553, 0.671] |

The pocket model minus pocket pLDDT is +0.029 [−0.034, +0.091]: no reliable gain.
Pocket pLDDT holds within pocket-size tertiles (AUROC 0.685 / 0.738 / 0.652 for
small / medium / large), and correlates only weakly with pocket size (Spearman
−0.19).

**Added after the plan (labelled as such):** the AlphaFold-state analysis split by
AlphaFold2's training cutoff. Only 9 discovery-set moving pockets have a post-cutoff
holo, too few, which motivated v0.3.

## v0.3: external validation (pre-declared in `docs/plan-v0.3-external.md`, commit ad7767c)

All numbers from `scripts/05_external_dataset.py` and `scripts/06_external_eval.py`.

**Set:** 971 proteins with a valid pair → 833 labelled; 84 (10.1%) move > 2 Å.
Main drops: too few pocket residues resolved in the apo (85 pairs), no AlphaFold DB
model (51 proteins), peptide ligands (36 pairs). 30% human. Holo structures all
released ≥ 2019-01-01. 658 proteins are UniProt accessions not in the discovery set;
480 have no discovery protein in their 30% cluster.

**Frozen rules on the external set, > 2 Å:**

| Subset | n (moving) | Whole-protein pLDDT | Pocket pLDDT | Pocket model |
|---|---|---|---|---|
| All | 833 (84) | 0.440 [0.379, 0.503] | 0.763 [0.710, 0.811] | 0.791 [0.736, 0.840] |
| New UniProt | 658 (58) | 0.462 | 0.743 [0.684, 0.801] | 0.786 |
| New 30% cluster | 480 (46) | 0.481 | 0.714 [0.640, 0.785] | 0.749 |

The pocket result replicates, and is no weaker than in discovery CV. It is somewhat
weaker on the most distant proteins (0.71) but still well clear of chance.
Whole-protein pLDDT stays at or below chance throughout.

**AlphaFold state, external:** for moving pockets, AlphaFold is closer to holo in 69%
of cases (n = 84; median 0.92 Å to holo, 2.54 Å to apo). With apo released before the
cutoff it is 70% (n = 40); after, 68% (n = 44). No holo here can have been in
AlphaFold2's training set. The holo preference is therefore not mainly memorisation
of the bound entry, though it is lower than on the discovery set (80%).

**What did not go as expected:**

- The PAE hinge features (v0.1), expected to be the strongest signal, contributed
  nothing.
- Adenylate kinase, the textbook flexible enzyme, is missed by the pocket rule
  ("moderate", pocket pLDDT 94.7). Confident rigid-domain closures are a blind spot.
- The richer pocket model never reliably beat the one-number rule. The CLI therefore
  uses the rule.
