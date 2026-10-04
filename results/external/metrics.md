Frozen models fit on 912 APObind proteins, applied to 833 PDB-2019+ proteins. Subset sizes: {'all': 833, 'new_uniprot': 661, 'new_cluster': 485}. 95% CIs: cluster bootstrap on the external set.

### > 1 Å, subset: all
202/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.569 [0.525, 0.615] | 0.282 [0.241, 0.341] | 0.011 [0.001, 0.041] | 0.184 [0.168, 0.200] |
| Pocket mean pLDDT (frozen logistic) | 0.732 [0.695, 0.767] | 0.424 [0.368, 0.503] | 0.073 [0.049, 0.101] | 0.172 [0.156, 0.188] |
| Pocket model (frozen boosting) | 0.746 [0.709, 0.783] | 0.472 [0.401, 0.556] | 0.023 [0.018, 0.056] | 0.158 [0.144, 0.173] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.163 [+0.123, +0.204]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.177 [+0.120, +0.234]

### > 1 Å, subset: new_uniprot
152/661 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.537 [0.486, 0.586] | 0.247 [0.205, 0.307] | 0.001 [0.000, 0.036] | 0.177 [0.160, 0.195] |
| Pocket mean pLDDT (frozen logistic) | 0.699 [0.655, 0.742] | 0.375 [0.313, 0.455] | 0.051 [0.030, 0.084] | 0.169 [0.151, 0.188] |
| Pocket model (frozen boosting) | 0.717 [0.663, 0.763] | 0.417 [0.336, 0.520] | 0.011 [0.016, 0.054] | 0.158 [0.142, 0.176] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.162 [+0.111, +0.214]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.178 [+0.106, +0.253]

### > 1 Å, subset: new_cluster
117/485 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.504 [0.446, 0.558] | 0.245 [0.200, 0.306] | 0.010 [0.001, 0.051] | 0.183 [0.161, 0.205] |
| Pocket mean pLDDT (frozen logistic) | 0.679 [0.626, 0.732] | 0.370 [0.301, 0.462] | 0.044 [0.024, 0.084] | 0.176 [0.153, 0.198] |
| Pocket model (frozen boosting) | 0.722 [0.669, 0.773] | 0.435 [0.354, 0.531] | 0.018 [0.019, 0.066] | 0.162 [0.142, 0.182] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.174 [+0.121, +0.230]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.217 [+0.150, +0.290]

### > 1.5 Å, subset: all
122/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.434 [0.382, 0.484] | 0.123 [0.102, 0.153] | 0.011 [0.001, 0.034] | 0.125 [0.110, 0.143] |
| Pocket mean pLDDT (frozen logistic) | 0.754 [0.714, 0.798] | 0.311 [0.253, 0.393] | 0.027 [0.014, 0.052] | 0.119 [0.103, 0.135] |
| Pocket model (frozen boosting) | 0.766 [0.721, 0.807] | 0.380 [0.298, 0.479] | 0.032 [0.023, 0.057] | 0.111 [0.097, 0.126] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.321 [+0.241, +0.407]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.332 [+0.261, +0.400]

### > 1.5 Å, subset: new_uniprot
92/661 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.452 [0.396, 0.510] | 0.121 [0.097, 0.158] | 0.019 [0.001, 0.046] | 0.120 [0.102, 0.139] |
| Pocket mean pLDDT (frozen logistic) | 0.730 [0.678, 0.783] | 0.276 [0.210, 0.369] | 0.024 [0.011, 0.053] | 0.115 [0.097, 0.134] |
| Pocket model (frozen boosting) | 0.740 [0.688, 0.790] | 0.338 [0.251, 0.446] | 0.040 [0.024, 0.066] | 0.109 [0.093, 0.125] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.279 [+0.178, +0.376]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.288 [+0.205, +0.368]

### > 1.5 Å, subset: new_cluster
76/485 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.475 [0.411, 0.539] | 0.143 [0.114, 0.192] | 0.001 [0.001, 0.037] | 0.132 [0.110, 0.155] |
| Pocket mean pLDDT (frozen logistic) | 0.699 [0.635, 0.756] | 0.269 [0.206, 0.363] | 0.017 [0.008, 0.052] | 0.129 [0.107, 0.152] |
| Pocket model (frozen boosting) | 0.718 [0.656, 0.775] | 0.327 [0.244, 0.444] | 0.038 [0.020, 0.067] | 0.121 [0.101, 0.143] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.223 [+0.114, +0.328]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.243 [+0.153, +0.335]

### > 2 Å, subset: all
84/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.440 [0.377, 0.498] | 0.086 [0.067, 0.114] | 0.026 [0.007, 0.047] | 0.092 [0.076, 0.109] |
| Pocket mean pLDDT (frozen logistic) | 0.762 [0.713, 0.810] | 0.237 [0.181, 0.328] | 0.047 [0.029, 0.068] | 0.088 [0.072, 0.104] |
| Pocket model (frozen boosting) | 0.781 [0.721, 0.833] | 0.343 [0.247, 0.452] | 0.035 [0.021, 0.054] | 0.079 [0.065, 0.093] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.323 [+0.229, +0.423]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.340 [+0.257, +0.409]

### > 2 Å, subset: new_uniprot
58/661 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.462 [0.394, 0.534] | 0.079 [0.060, 0.115] | 0.026 [0.008, 0.048] | 0.081 [0.066, 0.098] |
| Pocket mean pLDDT (frozen logistic) | 0.742 [0.679, 0.802] | 0.193 [0.137, 0.282] | 0.043 [0.025, 0.065] | 0.079 [0.063, 0.095] |
| Pocket model (frozen boosting) | 0.771 [0.705, 0.831] | 0.308 [0.206, 0.440] | 0.045 [0.031, 0.065] | 0.072 [0.058, 0.086] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.282 [+0.155, +0.394]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.310 [+0.206, +0.404]

### > 2 Å, subset: new_cluster
48/485 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.483 [0.396, 0.564] | 0.094 [0.067, 0.139] | 0.014 [0.003, 0.043] | 0.090 [0.068, 0.109] |
| Pocket mean pLDDT (frozen logistic) | 0.704 [0.634, 0.773] | 0.178 [0.126, 0.263] | 0.041 [0.020, 0.067] | 0.089 [0.068, 0.109] |
| Pocket model (frozen boosting) | 0.738 [0.664, 0.812] | 0.257 [0.178, 0.375] | 0.031 [0.017, 0.057] | 0.082 [0.064, 0.100] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.223 [+0.087, +0.355]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.257 [+0.143, +0.370]

### > 3 Å, subset: all
35/833 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.449 [0.366, 0.536] | 0.036 [0.025, 0.056] | 0.000 [0.000, 0.017] | 0.040 [0.028, 0.055] |
| Pocket mean pLDDT (frozen logistic) | 0.706 [0.611, 0.780] | 0.091 [0.051, 0.170] | 0.005 [0.002, 0.022] | 0.040 [0.028, 0.055] |
| Pocket model (frozen boosting) | 0.760 [0.656, 0.852] | 0.240 [0.118, 0.404] | 0.017 [0.007, 0.030] | 0.037 [0.026, 0.050] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.253 [+0.098, +0.395]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.309 [+0.174, +0.439]

### > 3 Å, subset: new_uniprot
26/661 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.455 [0.355, 0.557] | 0.035 [0.023, 0.057] | 0.003 [0.000, 0.019] | 0.038 [0.025, 0.053] |
| Pocket mean pLDDT (frozen logistic) | 0.702 [0.595, 0.796] | 0.088 [0.049, 0.184] | 0.004 [0.001, 0.020] | 0.038 [0.025, 0.053] |
| Pocket model (frozen boosting) | 0.706 [0.578, 0.826] | 0.182 [0.072, 0.354] | 0.014 [0.005, 0.029] | 0.036 [0.024, 0.049] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.246 [+0.053, +0.412]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.250 [+0.075, +0.413]

### > 3 Å, subset: new_cluster
19/485 positive.

| Model | AUROC | AUPRC | ECE | Brier |
|---|---|---|---|---|
| Whole-protein mean pLDDT (frozen logistic) | 0.470 [0.370, 0.586] | 0.036 [0.023, 0.061] | 0.003 [0.000, 0.020] | 0.038 [0.023, 0.054] |
| Pocket mean pLDDT (frozen logistic) | 0.631 [0.508, 0.746] | 0.068 [0.034, 0.150] | 0.004 [0.001, 0.024] | 0.039 [0.024, 0.055] |
| Pocket model (frozen boosting) | 0.617 [0.450, 0.757] | 0.102 [0.043, 0.250] | 0.009 [0.004, 0.026] | 0.037 [0.022, 0.053] |
- Pocket mean pLDDT (frozen logistic) − whole-protein pLDDT, AUROC: +0.159 [-0.052, +0.350]
- Pocket model (frozen boosting) − whole-protein pLDDT, AUROC: +0.142 [-0.068, +0.323]

### Risk by pocket-pLDDT band (edges fixed on APObind)

| Band | APObind > 2 Å | External > 2 Å |
|---|---|---|
| < 91.8 | 21.3% (39/183) | 23.6% (33/140) |
| 91.8–94.6 | 12.6% (23/182) | 18.6% (27/145) |
| 94.6–96.4 | 10.4% (19/182) | 9.0% (15/167) |
| 96.4–97.9 | 7.1% (13/182) | 2.9% (6/209) |
| > 97.9 | 2.7% (5/183) | 1.7% (3/172) |

### Which conformation does AlphaFold return?

| Set | Subset | n | AF closer to holo | median AF–holo | median AF–apo | median apo–holo | site pLDDT > 90 |
|---|---|---|---|---|---|---|---|
| apobind | all | 912 | 57.5% | 0.40 | 0.46 | 0.42 | 87.2% |
| apobind | moving | 99 | 79.8% | 0.60 | 2.48 | 2.64 | 70.7% |
| apobind | moving_apo_before_cutoff | 90 | 81.1% | 0.61 | 2.51 | 2.63 | 73.3% |
| apobind | moving_apo_after_cutoff | 9 | 66.7% | 0.41 | 2.17 | 2.64 | 44.4% |
| external | all | 833 | 56.5% | 0.41 | 0.46 | 0.41 | 89.2% |
| external | moving | 84 | 69.0% | 0.92 | 2.54 | 2.80 | 69.0% |
| external | moving_apo_before_cutoff | 40 | 70.0% | 0.82 | 2.65 | 2.76 | 65.0% |
| external | moving_apo_after_cutoff | 44 | 68.2% | 0.97 | 2.51 | 2.85 | 72.7% |
