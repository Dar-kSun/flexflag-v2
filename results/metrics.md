n = 916 proteins in 671 sequence clusters (MMseqs2, 30% identity). 5-fold cluster-grouped CV; 95% CIs from 1000 cluster bootstraps.

### Site RMSD > 1 Å
213 positives (23.3% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.418 [0.367, 0.466] | 0.200 [0.169, 0.242] | 0.008 [0.004, 0.040] | 0.180 [0.164, 0.197] |
| Baseline: fraction pLDDT < 70 | 0.424 [0.372, 0.477] | 0.208 [0.176, 0.256] | 0.015 [0.005, 0.044] | 0.180 [0.164, 0.197] |
| Logistic regression, all features | 0.574 [0.528, 0.617] | 0.275 [0.231, 0.324] | 0.039 [0.029, 0.074] | 0.181 [0.165, 0.198] |
| Full model (gradient boosting) | 0.562 [0.517, 0.605] | 0.273 [0.232, 0.325] | 0.013 [0.009, 0.046] | 0.177 [0.162, 0.194] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.007 [-0.018, +0.030], AUPRC +0.010 [-0.004, +0.026]
- Logistic regression, all features: AUROC +0.156 [+0.094, +0.219], AUPRC +0.075 [+0.035, +0.115]
- Full model (gradient boosting): AUROC +0.143 [+0.081, +0.206], AUPRC +0.075 [+0.035, +0.116]

### Site RMSD > 1.5 Å
143 positives (15.6% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.486 [0.435, 0.536] | 0.148 [0.122, 0.181] | 0.001 [0.001, 0.030] | 0.132 [0.115, 0.148] |
| Baseline: fraction pLDDT < 70 | 0.512 [0.463, 0.563] | 0.160 [0.132, 0.199] | 0.001 [0.001, 0.030] | 0.132 [0.115, 0.148] |
| Logistic regression, all features | 0.584 [0.532, 0.634] | 0.193 [0.155, 0.235] | 0.028 [0.021, 0.058] | 0.135 [0.118, 0.151] |
| Full model (gradient boosting) | 0.555 [0.501, 0.605] | 0.184 [0.147, 0.231] | 0.001 [0.004, 0.033] | 0.131 [0.114, 0.148] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.025 [+0.006, +0.045], AUPRC +0.013 [+0.003, +0.025]
- Logistic regression, all features: AUROC +0.097 [+0.036, +0.160], AUPRC +0.046 [+0.013, +0.079]
- Full model (gradient boosting): AUROC +0.068 [+0.000, +0.128], AUPRC +0.038 [+0.005, +0.072]

### Site RMSD > 2 Å
100 positives (10.9% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.528 [0.472, 0.582] | 0.107 [0.085, 0.135] | 0.000 [0.002, 0.027] | 0.097 [0.082, 0.113] |
| Baseline: fraction pLDDT < 70 | 0.552 [0.492, 0.609] | 0.119 [0.092, 0.151] | 0.004 [0.002, 0.028] | 0.097 [0.082, 0.112] |
| Logistic regression, all features | 0.596 [0.538, 0.654] | 0.140 [0.108, 0.185] | 0.022 [0.013, 0.045] | 0.099 [0.084, 0.115] |
| Full model (gradient boosting) | 0.574 [0.507, 0.643] | 0.159 [0.118, 0.226] | 0.011 [0.004, 0.031] | 0.096 [0.081, 0.112] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.024 [+0.001, +0.049], AUPRC +0.013 [+0.003, +0.024]
- Logistic regression, all features: AUROC +0.070 [+0.011, +0.130], AUPRC +0.035 [+0.011, +0.065]
- Full model (gradient boosting): AUROC +0.046 [-0.027, +0.124], AUPRC +0.057 [+0.017, +0.115]

### Site RMSD > 3 Å
38 positives (4.1% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.470 [0.387, 0.553] | 0.038 [0.026, 0.056] | 0.000 [0.000, 0.015] | 0.040 [0.028, 0.052] |
| Baseline: fraction pLDDT < 70 | 0.522 [0.435, 0.612] | 0.043 [0.029, 0.064] | 0.000 [0.000, 0.015] | 0.040 [0.028, 0.052] |
| Logistic regression, all features | 0.537 [0.441, 0.644] | 0.046 [0.031, 0.072] | 0.013 [0.007, 0.026] | 0.041 [0.029, 0.053] |
| Full model (gradient boosting) | 0.541 [0.445, 0.652] | 0.057 [0.036, 0.120] | 0.001 [0.001, 0.015] | 0.040 [0.028, 0.052] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.052 [+0.002, +0.105], AUPRC +0.005 [+0.000, +0.013]
- Logistic regression, all features: AUROC +0.068 [-0.034, +0.169], AUPRC +0.009 [-0.002, +0.025]
- Full model (gradient boosting): AUROC +0.072 [-0.030, +0.168], AUPRC +0.026 [+0.003, +0.076]
