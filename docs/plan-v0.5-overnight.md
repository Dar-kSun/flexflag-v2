# Pre-declared plan: overnight analyses (written before running them)

## E1. Does a protein language model add anything at the pocket?

ESM-2 650M (`esm2_t33_650M_UR50D`, fair-esm) per-residue embeddings for the
AlphaFold DB sequence of every protein in both datasets (layer 33). Pocket feature =
mean embedding over the 5 Å pocket residues (the same residues as the pocket pLDDT
feature). Models, scored like v0.3:

- **pLDDT:** logistic on pocket pLDDT (the frozen rule's form).
- **ESM:** PCA to 32 components (fit inside each training fold), then logistic (C = 0.1).
- **ESM + pLDDT:** the same PCA components plus pocket pLDDT, logistic (C = 0.1).

APObind: 5-fold cluster CV, seed 0. External: fit on all of APObind, applied
unchanged. 95% cluster-bootstrap CIs. Question: does ESM + pLDDT beat pLDDT alone on
the external set (paired bootstrap of the AUROC difference, CI excluding 0)? If not,
the README says so and the rule stays as it is.

## E2. Many ligands per protein

For each protein, keep the first pair's **apo** fixed and label up to 10 further holo
entries of the same protein (APObind: its other pairs with that apo; external: other
post-2019 holo entries with a real ligand and ≥ 90% range overlap), in the build
order. Per protein: number of labelled ligands, fraction moving > 2 Å, and max and
median site RMSD (each ligand defines its own pocket).

Report:

- how often proteins with ≥ 3 ligands are consistent (all move or none move);
- the AUROC of the first pocket's pLDDT for "any ligand moves > 2 Å" and for "most
  ligands move", and the Spearman correlation with the fraction moving.
