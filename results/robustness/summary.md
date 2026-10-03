Robustness checks (docs/plan-v0.4-robustness.md). AUROC of (-pLDDT), 95% cluster-bootstrap CI.

### apobind (failed: 0)

| Check | AUROC | positives |
|---|---|---|
| Pocket pLDDT, 5 Å pocket (reference) | 0.693 [0.636, 0.746] | 101/929 |
| R1 permuted labels (200×) | mean 0.497, 95% range [0.437, 0.568] | |
| R2 decoy patch, same size, > 15 Å away | 0.461 [0.406, 0.524] | 101/929 |
| R3 pocket defined at 4 Å | 0.700 [0.645, 0.754] | 101/920 |
| R3 pocket defined at 6 Å | 0.683 [0.628, 0.734] | 101/929 |
| R3 pocket defined at 8 Å | 0.674 [0.619, 0.729] | 101/929 |
| R4a label: all-heavy-atom RMSD > 2 Å | 0.690 [0.644, 0.738] | 166/929 |
| R4a label: all-heavy-atom RMSD > 2.56 Å (matched rate) | 0.686 [0.623, 0.745] | 101/929 |
| R4b label: whole-chain superposition, site RMSD > 2 Å | 0.691 [0.642, 0.742] | 141/929 |

Spearman with the primary label: heavy-atom 0.86, whole-chain superposition 0.94.

### external (failed: 0)

| Check | AUROC | positives |
|---|---|---|
| Pocket pLDDT, 5 Å pocket (reference) | 0.763 [0.713, 0.809] | 84/833 |
| R1 permuted labels (200×) | mean 0.504, 95% range [0.446, 0.566] | |
| R2 decoy patch, same size, > 15 Å away | 0.520 [0.458, 0.584] | 84/830 |
| R3 pocket defined at 4 Å | 0.763 [0.710, 0.812] | 83/826 |
| R3 pocket defined at 6 Å | 0.752 [0.699, 0.799] | 84/833 |
| R3 pocket defined at 8 Å | 0.743 [0.691, 0.791] | 84/833 |
| R4a label: all-heavy-atom RMSD > 2 Å | 0.721 [0.680, 0.758] | 157/833 |
| R4a label: all-heavy-atom RMSD > 2.63 Å (matched rate) | 0.727 [0.673, 0.778] | 84/833 |
| R4b label: whole-chain superposition, site RMSD > 2 Å | 0.761 [0.717, 0.806] | 116/833 |

Spearman with the primary label: heavy-atom 0.82, whole-chain superposition 0.92.
