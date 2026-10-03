# Pre-declared plan: robustness and end-to-end checks (written before running them)

The v0.3 claim: mean AlphaFold pLDDT over the pocket residues flags pockets that
reshape between apo and holo (external AUROC 0.76 at 2 Å), while whole-protein
pLDDT does not. These checks try to break that claim. Each is scored on **both**
datasets (APObind = discovery, PDB-2019+ = external) with the in-sample AUROC of
(−pocket pLDDT): no fitting, so no CV is needed. 95% CIs are cluster bootstraps.
All results go in `docs/findings.md` whichever way they go.

## R1. Negative control: permuted labels

Shuffle the 2 Å labels 200 times. The AUROC of pocket pLDDT must centre on 0.5. A
check of the scoring code, not of the science.

## R2. Pocket-specificity control: a same-size patch elsewhere

For each protein, 20 random "decoy pockets" on the AlphaFold model. Each is a random
seed residue whose Cα is > 15 Å from every pocket Cα, plus its nearest residues (by
Cα), giving the same residue count as the real pocket. Score = mean pLDDT over the
decoy, averaged over the 20 decoys. **If the decoy AUROC comes close to the pocket
AUROC (within 0.05), the signal is not pocket-specific.** Seed 0.

## R3. Pocket-definition sensitivity

Pocket residues for the *feature* at 4, 6 and 8 Å from the ligand, with the label
fixed at the 5 Å definition. Pass: AUROC stays within 0.05 of the 5 Å value.

## R4. Label variants

Same pairs and the same pocket feature, with the label replaced by:

- (a) all-heavy-atom pocket RMSD (atoms matched by name; site-superposed on Cα);
- (b) pocket Cα RMSD after **whole-chain** Cα superposition, which includes rigid
  movement of the pocket relative to the rest.

Each is binarised at 2 Å. (a) has many more positives, so it is also reported at
the cutoff that gives the same positive rate as the primary label. Pass: AUROC well
above 0.5 for each.

## R5. Label test–retest

For proteins with a second usable pair (a different holo entry, with the best
available apo other than the first pair's), recompute the label. Report the Spearman
correlation of the two site RMSDs, agreement of the 2 Å labels (Cohen's κ), and the
AUROC of pocket pLDDT against the retest labels. This sets the ceiling: a label that
does not reproduce cannot be predicted well. APObind only, plus the external set
where a second holo exists.

## R6. Ligand-free mode (P2Rank on the AlphaFold model)

Run P2Rank (default model) on each AlphaFold DB model, with no ligand information.

- Pocket recovery: fraction of proteins where any of the top-3 predicted pockets
  covers ≥ 50% of the true pocket residues.
- AUROC of the rule when the pocket is (i) the best-overlapping top-3 P2Rank pocket,
  and (ii) simply the top-1 P2Rank pocket, whether or not it is the true one. (ii) is
  the realistic no-ligand use.

## R7. Tool checks

- Fresh `pip install git+https://github.com/Dar-kSun/flexflag` in a clean venv:
  `flexflag check` runs, and `rule.json` ships in the package.
- Bad inputs give a one-line error, not a traceback: unknown UniProt, a UniProt
  missing from AlphaFold DB, a PDB entry with no ligand, a PDB entry of a different
  protein, a bad residue range.
