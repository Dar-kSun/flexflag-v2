Frozen models fit on 929 APObind proteins, applied to 833 PDB-2019+ proteins. Subset sizes: {'all': 833, 'new_uniprot': 658, 'new_cluster': 480}. 95% CIs: cluster bootstrap on the external set.

### > 1 Å, subset: all
202/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.569 [0.527, 0.612] | 0.282 [0.240, 0.338] | 0.010 [0.001, 0.041] | 0.183 [0.168, 0.200] |
| Pocket mean pLDDT (frozen logistic) | 0.733 [0.696, 0.769] | 0.425 [0.370, 0.501] | 0.079 [0.054, 0.108] | 0.172 [0.156, 0.189] |
| Pocket model (frozen boosting) | 0.746 [0.705, 0.783] | 0.462 [0.391, 0.546] | 0.030 [0.021, 0.061] | 0.159 [0.145, 0.174] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.164 [+0.123, +0.203]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.177 [+0.117, +0.234]

### > 1 Å, subset: new_uniprot
151/658 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.541 [0.490, 0.590] | 0.248 [0.206, 0.304] | 0.002 [0.001, 0.039] | 0.177 [0.159, 0.195] |
| Pocket mean pLDDT (frozen logistic) | 0.703 [0.657, 0.749] | 0.380 [0.313, 0.463] | 0.061 [0.035, 0.094] | 0.168 [0.150, 0.187] |
| Pocket model (frozen boosting) | 0.721 [0.671, 0.768] | 0.415 [0.342, 0.519] | 0.018 [0.017, 0.057] | 0.158 [0.141, 0.176] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.162 [+0.110, +0.213]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.179 [+0.111, +0.252]

### > 1 Å, subset: new_cluster
115/480 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.505 [0.449, 0.562] | 0.244 [0.197, 0.310] | 0.008 [0.001, 0.047] | 0.182 [0.161, 0.203] |
| Pocket mean pLDDT (frozen logistic) | 0.683 [0.630, 0.735] | 0.370 [0.303, 0.462] | 0.055 [0.031, 0.098] | 0.175 [0.154, 0.195] |
| Pocket model (frozen boosting) | 0.719 [0.668, 0.769] | 0.419 [0.347, 0.512] | 0.018 [0.018, 0.067] | 0.163 [0.143, 0.180] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.177 [+0.124, +0.229]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.213 [+0.145, +0.287]

### > 1.5 Å, subset: all
122/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.434 [0.382, 0.492] | 0.123 [0.102, 0.156] | 0.010 [0.001, 0.033] | 0.125 [0.110, 0.144] |
| Pocket mean pLDDT (frozen logistic) | 0.755 [0.713, 0.796] | 0.311 [0.253, 0.404] | 0.027 [0.012, 0.052] | 0.119 [0.103, 0.136] |
| Pocket model (frozen boosting) | 0.767 [0.718, 0.810] | 0.397 [0.312, 0.500] | 0.038 [0.025, 0.061] | 0.110 [0.097, 0.126] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.320 [+0.231, +0.404]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.332 [+0.255, +0.399]

### > 1.5 Å, subset: new_uniprot
92/658 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.452 [0.395, 0.509] | 0.122 [0.096, 0.159] | 0.018 [0.001, 0.046] | 0.121 [0.101, 0.140] |
| Pocket mean pLDDT (frozen logistic) | 0.731 [0.684, 0.779] | 0.279 [0.216, 0.370] | 0.024 [0.011, 0.051] | 0.116 [0.097, 0.135] |
| Pocket model (frozen boosting) | 0.739 [0.680, 0.793] | 0.356 [0.268, 0.463] | 0.038 [0.024, 0.066] | 0.109 [0.092, 0.125] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.279 [+0.185, +0.366]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.287 [+0.203, +0.364]

### > 1.5 Å, subset: new_cluster
74/480 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.473 [0.407, 0.540] | 0.141 [0.111, 0.191] | 0.003 [0.000, 0.036] | 0.130 [0.108, 0.153] |
| Pocket mean pLDDT (frozen logistic) | 0.705 [0.639, 0.760] | 0.268 [0.199, 0.367] | 0.014 [0.008, 0.048] | 0.127 [0.106, 0.151] |
| Pocket model (frozen boosting) | 0.719 [0.655, 0.776] | 0.342 [0.256, 0.444] | 0.037 [0.019, 0.070] | 0.119 [0.101, 0.140] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.230 [+0.115, +0.340]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.245 [+0.143, +0.338]

### > 2 Å, subset: all
84/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.440 [0.379, 0.503] | 0.086 [0.068, 0.113] | 0.026 [0.008, 0.046] | 0.091 [0.076, 0.107] |
| Pocket mean pLDDT (frozen logistic) | 0.763 [0.710, 0.811] | 0.237 [0.173, 0.321] | 0.047 [0.029, 0.068] | 0.088 [0.072, 0.103] |
| Pocket model (frozen boosting) | 0.791 [0.736, 0.840] | 0.344 [0.248, 0.457] | 0.035 [0.019, 0.055] | 0.079 [0.066, 0.093] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.323 [+0.222, +0.419]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.351 [+0.272, +0.421]

### > 2 Å, subset: new_uniprot
58/658 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.462 [0.393, 0.531] | 0.080 [0.060, 0.117] | 0.026 [0.007, 0.048] | 0.081 [0.064, 0.100] |
| Pocket mean pLDDT (frozen logistic) | 0.743 [0.684, 0.801] | 0.195 [0.135, 0.291] | 0.043 [0.025, 0.065] | 0.079 [0.062, 0.098] |
| Pocket model (frozen boosting) | 0.786 [0.719, 0.845] | 0.324 [0.211, 0.466] | 0.038 [0.024, 0.060] | 0.072 [0.057, 0.087] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.279 [+0.165, +0.398]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.324 [+0.220, +0.415]

### > 2 Å, subset: new_cluster
46/480 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.481 [0.399, 0.567] | 0.091 [0.066, 0.145] | 0.017 [0.003, 0.044] | 0.087 [0.067, 0.108] |
| Pocket mean pLDDT (frozen logistic) | 0.714 [0.640, 0.785] | 0.177 [0.123, 0.276] | 0.042 [0.020, 0.069] | 0.086 [0.066, 0.108] |
| Pocket model (frozen boosting) | 0.749 [0.671, 0.823] | 0.262 [0.171, 0.405] | 0.029 [0.017, 0.057] | 0.080 [0.062, 0.098] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.233 [+0.087, +0.371]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.269 [+0.149, +0.381]

### > 3 Å, subset: all
35/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.449 [0.359, 0.535] | 0.036 [0.025, 0.056] | 0.000 [0.000, 0.017] | 0.040 [0.028, 0.054] |
| Pocket mean pLDDT (frozen logistic) | 0.707 [0.618, 0.790] | 0.091 [0.054, 0.179] | 0.006 [0.002, 0.022] | 0.040 [0.029, 0.055] |
| Pocket model (frozen boosting) | 0.767 [0.666, 0.856] | 0.255 [0.131, 0.412] | 0.014 [0.004, 0.026] | 0.037 [0.026, 0.050] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.257 [+0.111, +0.409]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.318 [+0.185, +0.445]

### > 3 Å, subset: new_uniprot
26/658 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.455 [0.359, 0.562] | 0.035 [0.022, 0.058] | 0.002 [0.000, 0.019] | 0.038 [0.024, 0.052] |
| Pocket mean pLDDT (frozen logistic) | 0.702 [0.600, 0.793] | 0.088 [0.048, 0.181] | 0.005 [0.001, 0.020] | 0.038 [0.025, 0.052] |
| Pocket model (frozen boosting) | 0.715 [0.587, 0.832] | 0.209 [0.077, 0.420] | 0.011 [0.004, 0.026] | 0.036 [0.023, 0.048] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.251 [+0.058, +0.417]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.261 [+0.085, +0.416]

### > 3 Å, subset: new_cluster
18/480 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.460 [0.352, 0.569] | 0.034 [0.021, 0.059] | 0.004 [0.000, 0.021] | 0.036 [0.021, 0.052] |
| Pocket mean pLDDT (frozen logistic) | 0.661 [0.538, 0.775] | 0.070 [0.037, 0.158] | 0.003 [0.001, 0.022] | 0.037 [0.022, 0.053] |
| Pocket model (frozen boosting) | 0.628 [0.476, 0.782] | 0.132 [0.044, 0.328] | 0.009 [0.003, 0.026] | 0.035 [0.021, 0.050] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.198 [-0.012, +0.396]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.166 [-0.050, +0.352]

### Risk by pocket-pLDDT band (edges fixed on APObind)

| Band | APObind > 2 Å | External > 2 Å |
|---|---|---|
| < 91.8 | 21.0% (39/186) | 23.7% (33/139) |
| 91.8–94.6 | 13.4% (25/186) | 18.5% (27/146) |
| 94.6–96.4 | 10.3% (19/185) | 9.2% (15/163) |
| 96.4–97.9 | 7.0% (13/186) | 2.8% (6/215) |
| > 97.9 | 2.7% (5/186) | 1.8% (3/170) |

### Which conformation does AlphaFold return?

| Set | Subset | n | AF closer to holo | median AF–holo | median AF–apo | median apo–holo | site pLDDT > 90 |
|---|---|---|---|---|---|---|---|
| apobind | all | 929 | 57.4% | 0.39 | 0.46 | 0.42 | 87.1% |
| apobind | moving | 101 | 80.2% | 0.60 | 2.47 | 2.60 | 71.3% |
| apobind | moving_apo_before_cutoff | 92 | 81.5% | 0.61 | 2.49 | 2.60 | 73.9% |
| apobind | moving_apo_after_cutoff | 9 | 66.7% | 0.41 | 2.17 | 2.64 | 44.4% |
| external | all | 833 | 56.5% | 0.41 | 0.46 | 0.41 | 89.2% |
| external | moving | 84 | 69.0% | 0.92 | 2.54 | 2.80 | 69.0% |
| external | moving_apo_before_cutoff | 40 | 70.0% | 0.82 | 2.65 | 2.76 | 65.0% |
| external | moving_apo_after_cutoff | 44 | 68.2% | 0.97 | 2.51 | 2.85 | 72.7% |
