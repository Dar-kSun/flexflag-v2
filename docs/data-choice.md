# Data choice

## Source: APObind (pair list only)

**Decision:** apo–holo pairs come from APObind (Aggarwal et al., 2021,
arXiv:2108.09926; https://github.com/devalab/Apobind). APObind pairs each
PDBbind v2019 complex with an apo structure of the same protein: BLAST at >80%
sequence identity, rejected if no ligand sits within 4 Å of the crystal pose,
backbone RMSD above 15 Å, or ligand above 1000 Da.

**Only the pair list is used.** APObind's data lives in a University of Pittsburgh
SharePoint folder that returns HTTP 403 to scripted requests, so
`apobind_all.csv` was downloaded once by hand in a browser (2026-10-03). Every
structure is then fetched fresh from RCSB by PDB ID, so every structure in the
dataset traces back to a public PDB entry.

The upstream repo has no licence file. So the full CSV is not committed
(`.gitignore`). `scripts/make_pair_list.py` writes `data/apobind_pairs.csv`,
which keeps only the identifiers and APObind's own quality metrics.

**Alternatives considered:**

| Source | Why not |
|---|---|
| PLINDER (apo links) | Public GCS bucket returned 403 on every object tried |
| AHoJ-DB v2c (apoholo.cz) | 38 GB archive. Its per-entry REST API works and is the fallback if APObind proves unusable |
| Build from SIFTS ourselves | Multi-day job; CLAUDE.md says prefer a packaged set |

## Counts (from `apobind_all.csv`)

| | |
|---|---|
| Pairs (one per holo entry) | 12,267 |
| Distinct holo entries | 12,267 |
| Distinct apo entries | **2,873** |
| Apo resolution recorded as 0 (likely NMR) | 438 |
| Holo with more than one listed chain | 4,619 |
| Sequence identity < 0.95 | 2,264 |
| Sequence coverage < 0.90 | 620 |

APObind's paper reports 10,599 pairs. The CSV has 12,267, and we use the CSV.

**Many holo entries share one apo**, e.g. a single empty protease paired with
dozens of ligand complexes. The number of distinct *proteins* is therefore well
below 12k. The distinct UniProt count is determined in M3, and the dataset size
quoted in results will be that count, not the pair count.

APObind's own backbone RMSD has a long tail (90th percentile 21.6 Å, some pairs
with TM-score 0.19). Those look like wrong-chain or wrong-protein matches rather
than motion, so M3 will need a quality filter (planned: APObind TM-score).

## Feasibility run (`scripts/00_feasibility.py`, 10 pairs, seed 0)

End to end: RCSB mmCIF for both structures → binding-site Cα RMSD → SIFTS
UniProt mapping → AlphaFold DB pLDDT → mean pLDDT.

- **7/10 pairs fully processed.** Cold run 106 s, about 11 s per pair, almost
  all of it network time (each RCSB file takes 4–8 s). Labelling itself takes
  ~0.1 s. The full set will need parallel downloads (M3).
- Failures, each a real category to count in M3:
  - **1tkz** (HIV-1 protease, P04585): no AlphaFold DB entry. Viral
    polyproteins are not in AFDB, so a known protein family drops out.
  - **3eyu** (antibody Fab): the ligand is a 6-residue peptide. **Peptide
    ligands are scoped out**; we label small molecules only.
  - **4x60**: only 15 of 22 site residues are resolved in the apo structure,
    below the 80% coverage we require.
- **Ligand-picking check:** APObind lists its binding-site residues (6 Å,
  0-based chain positions). For the 8 pairs where we found a ligand, 96–100% of
  our site residues are in APObind's list. Theirs are larger, which fits their
  measuring with hydrogens included. We pick the same pocket.
- **All 7 site Cα RMSDs were small (0.1–1.6 Å).** None would count as "large
  change" at 2 Å. APObind's side-chain RMSD (a different measure) puts 46% of
  pairs above 2 Å. Our Cα-only site RMSD is more conservative, so the label
  threshold must be set from the full M3 distribution and checked at
  alternatives, not fixed now.

## Label method (as implemented in `flexflag/labels.py`)

First model, no hydrogens, no waters, first altloc. Ligand = the non-polymer,
non-additive residue (≥ 6 heavy atoms) with the most heavy atoms within 5 Å of
the listed holo chains. Site = holo-chain residues within 5 Å of it. Holo and
apo chains are matched by global sequence alignment, keeping identical residues
only. The label is the site Cα RMSD after Kabsch superposition on the site Cα
atoms.
