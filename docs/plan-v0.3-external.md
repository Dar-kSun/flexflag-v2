# Pre-declared plan: external, post-AlphaFold validation (written before building it)

Everything up to v0.2 uses one dataset (APObind, built on PDBbind v2019). This plan
fixes, **before** the external set is built or labelled, how the v0.2 pocket result
will be tested on an independent set that AlphaFold2 cannot have memorised.

## External set ("PDB-2019+")

Built directly from the RCSB search API and SIFTS, independently of APObind:

- **Holo candidates:** X-ray, resolution ≤ 2.5 Å, exactly one protein entity, no
  nucleic acid, ≥ 1 non-polymer entity, **initially released on or after
  2019-01-01**. That is after PDBbind v2019 and after AlphaFold2's PDB training
  cutoff (2018-04-30). (43,016 entries at query time.)
- **Apo candidates:** X-ray, resolution ≤ 2.5 Å, exactly one protein entity, no
  nucleic acid, **no non-polymer entities at all** (strict apo: no ligands, ions or
  buffers), any release date. (18,190 entries at query time.)
- **Pairing:** every protein chain of both entries maps to the same single UniProt
  accession (SIFTS). The UniProt ranges covered by the two entries overlap by ≥ 90%
  of each entry's range, so the constructs are similar.
- **Holo must contain a real ligand:** at least one non-polymer component not in the
  v0.1 additive list and with formula weight ≥ 150 Da (RCSB chem_comp). The v0.1
  ligand picker still runs on the structure.
- **One protein, one pair:** per UniProt, holo candidates are ordered by resolution
  (ties broken by PDB ID), and the apo is the best-resolution apo entry. Up to 3 holo
  candidates are tried. Labels, site definition and filters are exactly v0.1's
  (`flexflag/labels.py`, 5 Å, ≥ 80% site residues resolved in the apo).
- AlphaFold DB F1 model required, length ≤ 1,500 (as v0.1).

## What is frozen from APObind (fit once on all 929 APObind proteins)

1. **Pocket-pLDDT rule:** logistic regression on `site_plddt_mean`.
2. **Pocket-pLDDT bands:** the APObind quintile edges of site pLDDT (91.8 / 94.6 /
   96.4 / 97.9).
3. **Pocket model:** v0.2 gradient boosting on site + whole-protein features.
4. **Whole-protein baseline:** logistic on whole-protein mean pLDDT.

Nothing is refit on the external set.

## What will be reported on the external set

- AUROC and AUPRC of 1–4 at > 2 Å (primary) and at 1, 1.5 and 3 Å, with 95% CIs from
  a cluster bootstrap (MMseqs2 30% clusters of the external set).
- Fraction moving > 2 Å in each frozen pLDDT band.
- Calibration (ECE) of the frozen rule.
- Subsets: (a) all; (b) UniProt accessions not in the APObind set; (c) proteins with
  no APObind protein in the same 30% cluster (clustered jointly).
- **AlphaFold state (v0.2 part A) on the external set**, where no holo structure can
  have been in AlphaFold2's training data: the fraction of moving pockets where the
  AlphaFold model is closer to holo than to apo. Also split by whether the apo was
  released before or after the cutoff.

Whatever these show goes into the README, including if the APObind result does not
replicate.
