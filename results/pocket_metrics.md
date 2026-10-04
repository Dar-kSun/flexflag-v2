Pocket-conditioned mode. n = 912 proteins, 669 clusters, 5-fold cluster CV, 1000 cluster bootstraps. Pre-declared in docs/plan-v0.2-pocket.md.

### Site RMSD > 1 Å (212 positives, 23.2%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.698 [0.653, 0.735] | 0.384 [0.320, 0.450] | 0.057 [0.034, 0.085] | 0.168 [0.151, 0.184] |
| Pocket model (site + whole-protein features, boosting) | 0.719 [0.676, 0.754] | 0.403 [0.338, 0.476] | 0.020 [0.013, 0.050] | 0.161 [0.146, 0.176] |
| Confound reference: site size alone (not a usable model) | 0.607 [0.564, 0.649] | 0.326 [0.272, 0.389] | 0.007 [0.009, 0.044] | 0.174 [0.157, 0.190] |

Pocket model − pocket baseline, AUROC: +0.021 [-0.014, +0.059]

### Site RMSD > 1.5 Å (142 positives, 15.6%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.683 [0.634, 0.734] | 0.276 [0.224, 0.352] | 0.016 [0.012, 0.045] | 0.127 [0.109, 0.145] |
| Pocket model (site + whole-protein features, boosting) | 0.720 [0.673, 0.763] | 0.306 [0.244, 0.384] | 0.025 [0.013, 0.048] | 0.121 [0.106, 0.137] |
| Confound reference: site size alone (not a usable model) | 0.608 [0.555, 0.665] | 0.242 [0.192, 0.311] | 0.026 [0.012, 0.051] | 0.129 [0.111, 0.147] |

Pocket model − pocket baseline, AUROC: +0.036 [-0.012, +0.088]

### Site RMSD > 2 Å (99 positives, 10.9%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.691 [0.633, 0.743] | 0.229 [0.170, 0.316] | 0.022 [0.013, 0.043] | 0.093 [0.078, 0.109] |
| Pocket model (site + whole-protein features, boosting) | 0.746 [0.686, 0.799] | 0.302 [0.221, 0.406] | 0.022 [0.012, 0.041] | 0.088 [0.074, 0.102] |
| Confound reference: site size alone (not a usable model) | 0.607 [0.543, 0.671] | 0.187 [0.140, 0.262] | 0.007 [0.004, 0.031] | 0.095 [0.079, 0.111] |

Pocket model − pocket baseline, AUROC: +0.054 [-0.007, +0.114]

### Site RMSD > 3 Å (38 positives, 4.2%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.678 [0.587, 0.765] | 0.097 [0.059, 0.184] | 0.004 [0.002, 0.017] | 0.040 [0.028, 0.051] |
| Pocket model (site + whole-protein features, boosting) | 0.715 [0.626, 0.801] | 0.149 [0.082, 0.277] | 0.007 [0.002, 0.020] | 0.038 [0.027, 0.050] |
| Confound reference: site size alone (not a usable model) | 0.633 [0.521, 0.739] | 0.118 [0.059, 0.215] | 0.004 [0.002, 0.017] | 0.039 [0.028, 0.051] |

Pocket model − pocket baseline, AUROC: +0.034 [-0.075, +0.153]

### Pocket pLDDT within pocket-size tertiles (> 2 Å, in-sample)

- small pockets (3–14 residues, n = 306, 24 moving): AUROC 0.691
- medium pockets (15–19 residues, n = 304, 29 moving): AUROC 0.737
- large pockets (20–47 residues, n = 302, 46 moving): AUROC 0.646