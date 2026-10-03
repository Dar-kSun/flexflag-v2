Pocket-conditioned mode. n = 929 proteins, 678 clusters, 5-fold cluster CV, 1000 cluster bootstraps. Pre-declared in docs/plan-v0.2-pocket.md.

### Site RMSD > 1 Å (217 positives, 23.4%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.699 [0.656, 0.737] | 0.385 [0.325, 0.448] | 0.056 [0.035, 0.083] | 0.169 [0.153, 0.186] |
| Pocket model (site + whole-protein features, boosting) | 0.720 [0.684, 0.758] | 0.394 [0.335, 0.467] | 0.031 [0.017, 0.057] | 0.162 [0.147, 0.177] |
| Confound reference: site size alone (not a usable model) | 0.612 [0.566, 0.652] | 0.326 [0.273, 0.386] | 0.004 [0.010, 0.044] | 0.174 [0.158, 0.191] |

Pocket model − pocket baseline, AUROC: +0.022 [-0.016, +0.062]

### Site RMSD > 1.5 Å (144 positives, 15.5%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.683 [0.633, 0.729] | 0.271 [0.215, 0.343] | 0.012 [0.010, 0.041] | 0.127 [0.111, 0.145] |
| Pocket model (site + whole-protein features, boosting) | 0.713 [0.664, 0.756] | 0.326 [0.258, 0.405] | 0.017 [0.011, 0.040] | 0.120 [0.106, 0.135] |
| Confound reference: site size alone (not a usable model) | 0.614 [0.556, 0.665] | 0.244 [0.197, 0.318] | 0.016 [0.009, 0.043] | 0.128 [0.111, 0.146] |

Pocket model − pocket baseline, AUROC: +0.029 [-0.021, +0.083]

### Site RMSD > 2 Å (101 positives, 10.9%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.686 [0.630, 0.739] | 0.214 [0.162, 0.299] | 0.023 [0.011, 0.044] | 0.094 [0.079, 0.110] |
| Pocket model (site + whole-protein features, boosting) | 0.716 [0.655, 0.767] | 0.292 [0.217, 0.391] | 0.010 [0.008, 0.034] | 0.089 [0.075, 0.103] |
| Confound reference: site size alone (not a usable model) | 0.614 [0.553, 0.671] | 0.192 [0.143, 0.269] | 0.016 [0.007, 0.038] | 0.095 [0.079, 0.110] |

Pocket model − pocket baseline, AUROC: +0.029 [-0.034, +0.091]

### Site RMSD > 3 Å (38 positives, 4.1%)

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Pocket baseline: site mean pLDDT | 0.679 [0.580, 0.766] | 0.092 [0.054, 0.182] | 0.005 [0.003, 0.019] | 0.039 [0.028, 0.051] |
| Pocket model (site + whole-protein features, boosting) | 0.763 [0.686, 0.843] | 0.184 [0.102, 0.329] | 0.015 [0.005, 0.027] | 0.037 [0.026, 0.048] |
| Confound reference: site size alone (not a usable model) | 0.659 [0.561, 0.745] | 0.104 [0.060, 0.196] | 0.003 [0.002, 0.018] | 0.039 [0.028, 0.050] |

Pocket model − pocket baseline, AUROC: +0.086 [-0.008, +0.191]

### Pocket pLDDT within pocket-size tertiles (> 2 Å, in-sample)

- small pockets (3–14 residues, n = 323, 25 moving): AUROC 0.685
- medium pockets (15–19 residues, n = 303, 29 moving): AUROC 0.738
- large pockets (20–47 residues, n = 303, 47 moving): AUROC 0.652