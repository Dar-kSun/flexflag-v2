# Pre-declared plan: docking validation (written before the full run)

CLAUDE.md lists "running docking" as out of scope for the *tool*. The project owner
approved this one-off **validation experiment**: it tests whether the flag predicts the
outcome it exists to warn about. flexflag itself still does not dock.

## Question

When the AlphaFold model of a protein is used for docking, does the pocket-pLDDT flag
(and the measured apo–holo pocket change) predict docking failure, beyond failures
that happen even in the holo crystal structure?

## Method (`scripts/14_docking.py`)

- Proteins: both datasets (APObind 916, PDB-2019+ 833), each with its labelled pair.
- Ligand: the crystal ligand that defines the label, with bond orders from the RCSB
  ModelServer (CCD chemistry), protonated as given, re-embedded from scratch with RDKit
  ETKDG (seed 0) so the input carries no memory of the crystal pose.
- **Excluded before docking:** ligands with > 60 heavy atoms or > 15 rotatable bonds
  (outside Vina's reliable range), and covalent ligands (any ligand heavy atom within
  2.0 Å of a protein heavy atom). Exclusions are counted and reported.
- Receptors, all single protein chains with no waters, superposed on the holo pocket Cα
  atoms. The holo structure's cofactors and metal ions within 8 Å of the ligand
  (non-polymer, not the docked ligand, not an additive) are kept and placed identically
  in all three receptors, so only the protein conformation differs. The receptors are: (1) **holo** crystal chain (control), (2) **apo** crystal
  chain, (3) **AlphaFold DB model**. Hydrogens added by Open Babel at pH 7.4.
- AutoDock Vina 1.2.7: box centred on the crystal ligand, ligand extent + 10 Å per
  axis (minimum 20 Å), exhaustiveness 8, seed 0, 9 modes.
- **Success:** top-ranked pose within 2.0 Å heavy-atom RMSD (RDKit, symmetry-aware, no
  re-alignment) of the crystal pose. Also recorded: success within the top 3 poses, to
  separate sampling from scoring failures.
- Order: the external set first (it carries the primary endpoint), then APObind.

## Pilot (before this plan was frozen)

The pipeline was piloted on 8 random proteins (seed 1) and on two classic re-docking
cases, 3PTB benzamidine (0.42 Å) and 1STP biotin (0.78 Å), which succeed as expected.
The pilot showed that stripping all hetero groups removes cofactors that form the
pocket, so cofactors and metals are now kept (above), and top-3 success was added.
Holo success in the pilot was 2/8. The 8 pilot proteins stay in the full run.

## Endpoints (fixed now)

Per dataset, over proteins where all three receptors docked:

1. Success rate per receptor.
2. **Structure-caused AlphaFold failure** = AlphaFold fails *and* holo succeeds. Its
   rate in moving (> 2 Å) vs non-moving pockets, and per pocket-pLDDT band (the frozen
   APObind band edges).
3. AUROC of (−pocket pLDDT) for AlphaFold failure: (a) among all proteins, and
   (b) among proteins where holo docking succeeded (the primary endpoint). 95% cluster
   bootstrap CIs.
4. The same for whole-protein pLDDT, as the baseline.
5. Apo-docking failure rate in moving vs non-moving pockets (sanity check: if the
   label means anything, apo docking should fail more where pockets move).

Interpretation fixed now: the flag is useful for docking if (3b) has a lower CI bound
above 0.5 on the external set and exceeds whole-protein pLDDT. If not, the README says
so in the results.

## Amendment (2026-10-04, before any docking outcome was examined)

The full run needs about 10 hours, which did not fit the time available. The run was
stopped after 11 proteins, whose outcomes were not looked at. Those rows were then
discarded, because of a CSV bug (rows with different fields were appended under one
header, so values were misaligned; fixed by writing a fixed column set). Changes:

- **External set only** (it carries the primary endpoint). APObind is not docked.
- Proteins are processed in a **random order (seed 0)**, and no new protein is started
  after a **2.5-hour budget**. The docked proteins are therefore a random sample of the
  external set. The sample size is reported with every number.
- Docking settings and endpoints are unchanged.
- Vina results were not bit-reproducible between runs on this machine: tiny numerical
  differences, and occasionally two near-tied top poses swap order. This affects all
  three receptors alike.
