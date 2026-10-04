Robustness checks (docs/plan-v0.4-robustness.md). AUROC of (-pLDDT), 95% cluster-bootstrap CI.

### apobind (failed: 0)

| Check | AUROC | positives |
|---|---|---|
| Pocket pLDDT, 5 Å pocket (reference) | 0.693 [0.636, 0.749] | 99/912 |
| R1 permuted labels (200×) | mean 0.496, 95% range [0.443, 0.557] | |
| R2 decoy patch, same size, > 15 Å away | 0.459 [0.405, 0.518] | 99/912 |
| R3 pocket defined at 4 Å | 0.700 [0.642, 0.756] | 99/903 |
| R3 pocket defined at 6 Å | 0.683 [0.628, 0.740] | 99/912 |
| R3 pocket defined at 8 Å | 0.674 [0.615, 0.727] | 99/912 |
| R4a label: all-heavy-atom RMSD > 2 Å | 0.689 [0.639, 0.737] | 160/912 |
| R4a label: all-heavy-atom RMSD > 2.56 Å (matched rate) | 0.691 [0.631, 0.752] | 99/912 |
| R4b label: whole-chain superposition, site RMSD > 2 Å | 0.690 [0.641, 0.738] | 137/912 |

Spearman with the primary label: heavy-atom 0.87, whole-chain superposition 0.94.

### external (failed: 0)

| Check | AUROC | positives |
|---|---|---|
| Pocket pLDDT, 5 Å pocket (reference) | 0.762 [0.711, 0.810] | 84/833 |
| R1 permuted labels (200×) | mean 0.505, 95% range [0.433, 0.566] | |
| R2 decoy patch, same size, > 15 Å away | 0.520 [0.449, 0.581] | 84/830 |
| R3 pocket defined at 4 Å | 0.763 [0.712, 0.810] | 83/826 |
| R3 pocket defined at 6 Å | 0.751 [0.697, 0.802] | 84/833 |
| R3 pocket defined at 8 Å | 0.743 [0.693, 0.793] | 84/833 |
| R4a label: all-heavy-atom RMSD > 2 Å | 0.721 [0.678, 0.760] | 157/833 |
| R4a label: all-heavy-atom RMSD > 2.63 Å (matched rate) | 0.727 [0.676, 0.774] | 84/833 |
| R4b label: whole-chain superposition, site RMSD > 2 Å | 0.761 [0.716, 0.807] | 116/833 |

Spearman with the primary label: heavy-atom 0.82, whole-chain superposition 0.92.
