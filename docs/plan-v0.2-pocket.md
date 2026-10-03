# Pre-declared plan: pocket-level analysis (written before running it)

v0.1 found that whole-protein AlphaFold confidence does not predict binding-site
change. This file fixes the follow-up analyses, features and comparisons **before**
any of them are run, so the results cannot be tuned to look good. Any deviation will
be reported in `docs/findings.md`.

## A. Which conformation does AlphaFold return?

For each of the 933 proteins, superpose the AlphaFold DB model on the binding site
(the same site residues as the label, Cα, site-superposed) and compute:

- `rmsd_af_holo`: AlphaFold model vs holo crystal structure
- `rmsd_af_apo`: AlphaFold model vs apo crystal structure

AlphaFold is "holo-like" if `rmsd_af_holo < rmsd_af_apo`. Reported for all proteins
and separately for the moving set (apo–holo site RMSD > 2 Å). Same residues in all
three structures. Questions: how often does AlphaFold return the bound-like pocket,
and how far is it from the empty one?

## C. Confident but moving

Fraction of moving pockets (> 2 Å) whose mean site pLDDT is > 90 and > 70.

## B. Pocket-conditioned flag (separate mode)

**Input change:** the pocket's residue positions are given, as they are at docking
time. In this dataset they come from the holo ligand (5 Å). Only the residue
**positions** are used, mapped onto the AlphaFold sequence. No holo coordinates,
ligand identity or ligand properties are used. **Site size is excluded** because it
partly encodes ligand size, and is reported only as a confound reference.

Features, computed on AlphaFold DB data at the site residues (fixed now):

1. `site_plddt_mean`
2. `site_plddt_min`
3. `site_pae_within`: mean PAE between pairs of site residues
4. `site_pae_to_rest`: mean PAE between site residues and all other residues
   with pLDDT ≥ 70
5. `site_n_domains`: number of distinct PAE domains (v0.1 definition) containing
   at least 2 site residues
6. `site_frac_low_plddt`: fraction of site residues with pLDDT < 70

Models, all with the same 5-fold cluster CV, seed 0 and cluster bootstrap as v0.1:

- pocket baseline: `site_plddt_mean` alone (logistic)
- pocket model: the 6 site features + the v0.1 whole-protein features (gradient
  boosting, same hyperparameters)
- confound reference (not a usable model): site size alone

The primary comparison is pocket model vs pocket baseline at 2 Å, with 1, 1.5 and 3 Å
reported as in v0.1. The v0.1 whole-protein result stays as it is.
