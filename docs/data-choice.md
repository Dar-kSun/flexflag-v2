# Where the data comes from, and why

To learn whether a pocket changes shape, you need proteins that have been solved
twice: once empty (the **apo** structure) and once with a molecule bound (the
**holo** structure). This note explains where those pairs come from, what was
filtered out, and what the first small test run showed.

## The discovery set: APObind

The pairs come from [APObind](https://github.com/devalab/Apobind) (Aggarwal, Gupta &
Priyakumar, 2021, arXiv:2108.09926). For each ligand-bound complex in PDBbind v2019,
APObind found an empty structure of the same protein (BLAST, over 80% sequence
identity). It dropped pairs where a ligand sits within 4 Å of the crystal pose in the
"empty" structure, where the backbone RMSD is above 15 Å, or where the ligand is
heavier than 1,000 Da.

We use only APObind's **list of pairs**, not its structure files. That list was
downloaded once by hand (2026-10-03), because APObind's data folder (a University of
Pittsburgh SharePoint) refuses scripted downloads. Every structure is then fetched
fresh from the PDB by its ID, so everything in the dataset traces back to a public PDB
entry.

The APObind repository has no licence file, so its full CSV is not committed (it is in
`.gitignore`). `scripts/make_pair_list.py` writes `data/apobind_pairs.csv`, which keeps
only the PDB identifiers and APObind's own quality scores.

### Other sources we looked at

| Source | Why we did not use it |
|---|---|
| PLINDER (apo links) | Its public storage bucket returned "access denied" for every file tried |
| AHoJ-DB v2c (apoholo.cz) | A 38 GB archive. Its per-entry web API works, and was the fallback |
| Building pairs from SIFTS ourselves | A multi-day job; a packaged set was preferred |

We later did build an independent set ourselves, for validation (below).

## What the APObind list contains

| | |
|---|---|
| Pairs (one per holo entry) | 12,267 |
| Distinct apo entries | **2,873** |
| Apo entries with no resolution recorded (probably NMR) | 438 |
| Holo entries listing more than one chain | 4,619 |
| Sequence identity below 0.95 | 2,264 |
| Sequence coverage below 0.90 | 620 |

APObind's paper reports 10,599 pairs; the CSV has 12,267, and we use the CSV.

Two things in this table shaped the whole project:

- **Many holo entries share one apo.** A single empty protease can be paired with
  dozens of drug complexes. So 12,267 pairs cover only about 1,200 distinct proteins.
  We keep one pair per protein, and every dataset size quoted in the results counts
  proteins, not pairs.
- **Some pairs are not real matches.** APObind's own backbone RMSD has a long tail
  (90th percentile 21.6 Å, and some pairs with a TM-score of 0.19). Those look like
  the wrong chain or the wrong protein rather than real motion. We therefore require
  an APObind TM-score of at least 0.5 and the same UniProt accession on both sides,
  among other filters. The exact counts removed at each step are in
  [findings](findings.md), and every dropped pair is listed in `results/dropped.csv`.

After all filters, the discovery set has **916 labelled proteins**.

## The external set: PDB-2019+

To check that the results hold up on structures nobody tuned anything on, we built a
second set directly from the PDB (`scripts/05_external_dataset.py`). It uses ligand-bound
structures released on or after 2019-01-01, which is after both PDBbind v2019 and
AlphaFold2's training cutoff (April 2018). They are paired with strictly empty
structures (no ligands, ions or buffers at all) of the same protein. That gives **833
labelled proteins**. How it was built was fixed in advance, in
[plan-v0.3-external.md](plan-v0.3-external.md).

## The first test run (10 pairs)

Before building anything at scale, `scripts/00_feasibility.py` took 10 random APObind
pairs all the way through: download both structures, measure how much the pocket
moves, look up the protein in UniProt and AlphaFold DB, and read its pLDDT.

- **7 of 10 went through.** The first run took 106 s, about 11 s per pair, almost all
  of it download time. The calculation itself takes about 0.1 s. The full build
  therefore downloads in parallel and caches everything.
- **The three failures were each a real category**, later counted at scale:
  - 1TKZ (HIV-1 protease): viral polyproteins have no AlphaFold DB model.
  - 3EYU (an antibody fragment): the ligand is a short peptide. Peptide ligands are
    out of scope; we only label small molecules.
  - 4X60: too few pocket residues are resolved in the empty structure (15 of 22;
    we require 80%).
- **Our pocket matches APObind's.** APObind lists its own binding-site residues. For
  the 8 pairs where we found a ligand, 96–100% of our pocket residues are on APObind's
  list. Theirs are a little larger, which fits their measuring with hydrogens included.
- **Every pocket in the test moved only a little** (0.1–1.6 Å). That warned us that
  "large change" would be rare. So the 2 Å threshold was checked against the full
  distribution and reported alongside 1, 1.5 and 3 Å, rather than taken on trust.

## How "pocket change" is measured

The method is in `flexflag/labels.py` and tested in `tests/test_labels.py`.

1. Take the first model of each structure, without hydrogens or waters, and keep one
   alternative conformation per atom.
2. Find the ligand: the bound molecule (not part of a protein or peptide chain, not a
   buffer or additive, at least 6 heavy atoms) with the most atoms near the protein.
3. The pocket is every residue with an atom within 5 Å of that ligand.
4. Line up the empty and bound structures residue by residue, using their sequences.
   This copes with different numbering and chain names.
5. Superpose the two pockets on their Cα atoms and measure the Cα RMSD. A pocket that
   moves more than 2 Å counts as "large change".

One correction came later. The original code also accepted unusual residues *inside*
a peptide chain as ligands, for example the fluorescent substrate peptide in 5BTR.
That contradicted step 2. It was fixed, tested, and every result was rebuilt; see the
deviations section of [findings](findings.md).
