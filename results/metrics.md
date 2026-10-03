n = 933 proteins in 680 sequence clusters (MMseqs2, 30% identity). 5-fold cluster-grouped CV; 95% CIs from 1000 cluster bootstraps.

### Site RMSD > 1 Å
218 positives (23.4% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.461 [0.413, 0.508] | 0.214 [0.183, 0.255] | 0.002 [0.001, 0.037] | 0.180 [0.164, 0.197] |
| Baseline: fraction pLDDT < 70 | 0.462 [0.418, 0.508] | 0.221 [0.189, 0.269] | 0.012 [0.004, 0.045] | 0.180 [0.164, 0.197] |
| Logistic regression, all features | 0.591 [0.546, 0.632] | 0.281 [0.237, 0.331] | 0.038 [0.026, 0.069] | 0.180 [0.165, 0.197] |
| Full model (gradient boosting) | 0.567 [0.520, 0.611] | 0.270 [0.228, 0.328] | 0.011 [0.007, 0.045] | 0.178 [0.162, 0.195] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.002 [-0.047, +0.046], AUPRC +0.008 [-0.018, +0.037]
- Logistic regression, all features: AUROC +0.130 [+0.063, +0.195], AUPRC +0.067 [+0.028, +0.105]
- Full model (gradient boosting): AUROC +0.106 [+0.042, +0.170], AUPRC +0.058 [+0.014, +0.106]

### Site RMSD > 1.5 Å
145 positives (15.5% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.467 [0.411, 0.523] | 0.142 [0.116, 0.176] | 0.005 [0.002, 0.031] | 0.132 [0.115, 0.150] |
| Baseline: fraction pLDDT < 70 | 0.510 [0.455, 0.563] | 0.161 [0.130, 0.202] | 0.009 [0.003, 0.034] | 0.132 [0.115, 0.149] |
| Logistic regression, all features | 0.594 [0.543, 0.644] | 0.199 [0.160, 0.253] | 0.025 [0.017, 0.053] | 0.132 [0.116, 0.149] |
| Full model (gradient boosting) | 0.580 [0.529, 0.627] | 0.198 [0.154, 0.254] | 0.002 [0.003, 0.031] | 0.130 [0.113, 0.147] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.042 [+0.009, +0.077], AUPRC +0.019 [+0.007, +0.036]
- Logistic regression, all features: AUROC +0.127 [+0.062, +0.191], AUPRC +0.059 [+0.023, +0.100]
- Full model (gradient boosting): AUROC +0.113 [+0.049, +0.173], AUPRC +0.058 [+0.020, +0.104]

### Site RMSD > 2 Å
102 positives (10.9% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.521 [0.461, 0.578] | 0.108 [0.085, 0.139] | 0.005 [0.002, 0.027] | 0.097 [0.082, 0.113] |
| Baseline: fraction pLDDT < 70 | 0.551 [0.492, 0.611] | 0.126 [0.097, 0.172] | 0.004 [0.002, 0.026] | 0.097 [0.082, 0.113] |
| Logistic regression, all features | 0.618 [0.564, 0.673] | 0.145 [0.114, 0.190] | 0.020 [0.011, 0.042] | 0.098 [0.083, 0.113] |
| Full model (gradient boosting) | 0.592 [0.529, 0.655] | 0.169 [0.123, 0.250] | 0.005 [0.003, 0.029] | 0.096 [0.081, 0.112] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.030 [+0.009, +0.051], AUPRC +0.020 [+0.007, +0.039]
- Logistic regression, all features: AUROC +0.098 [+0.036, +0.157], AUPRC +0.039 [+0.014, +0.067]
- Full model (gradient boosting): AUROC +0.073 [+0.002, +0.142], AUPRC +0.067 [+0.020, +0.133]

### Site RMSD > 3 Å
38 positives (4.1% base rate).

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Baseline: mean pLDDT | 0.498 [0.408, 0.585] | 0.041 [0.029, 0.065] | 0.000 [0.000, 0.014] | 0.039 [0.028, 0.051] |
| Baseline: fraction pLDDT < 70 | 0.524 [0.438, 0.615] | 0.044 [0.031, 0.083] | 0.000 [0.000, 0.014] | 0.039 [0.028, 0.051] |
| Logistic regression, all features | 0.563 [0.459, 0.658] | 0.052 [0.034, 0.090] | 0.008 [0.003, 0.021] | 0.040 [0.029, 0.051] |
| Full model (gradient boosting) | 0.612 [0.502, 0.715] | 0.079 [0.044, 0.152] | 0.002 [0.001, 0.016] | 0.039 [0.028, 0.051] |

Difference from the mean-pLDDT baseline (paired cluster bootstrap):

- Baseline: fraction pLDDT < 70: AUROC +0.027 [-0.004, +0.059], AUPRC +0.005 [-0.003, +0.028]
- Logistic regression, all features: AUROC +0.066 [-0.028, +0.168], AUPRC +0.013 [-0.009, +0.044]
- Full model (gradient boosting): AUROC +0.114 [+0.012, +0.221], AUPRC +0.044 [+0.007, +0.103]
