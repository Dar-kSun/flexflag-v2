ESM-2 pocket features (E1). APObind n = 912 (5-fold cluster CV); external n = 833 (frozen). > 2 Å.

### apobind_cv

| Model | AUROC | AUPRC |
|---|---|---|
| Pocket pLDDT (logistic) | 0.691 [0.632, 0.747] | 0.229 |
| ESM-2 pocket embedding (PCA 32, logistic) | 0.613 [0.551, 0.676] | 0.168 |
| ESM-2 pocket embedding + pocket pLDDT | 0.697 [0.643, 0.753] | 0.242 |
- ESM-2 pocket embedding (PCA 32, logistic) − pocket pLDDT, AUROC -0.076 [-0.165, +0.019]
- ESM-2 pocket embedding + pocket pLDDT − pocket pLDDT, AUROC +0.007 [-0.067, +0.084]

### external_frozen

| Model | AUROC | AUPRC |
|---|---|---|
| Pocket pLDDT (logistic) | 0.762 [0.717, 0.808] | 0.237 |
| ESM-2 pocket embedding (PCA 32, logistic) | 0.599 [0.535, 0.665] | 0.189 |
| ESM-2 pocket embedding + pocket pLDDT | 0.698 [0.632, 0.759] | 0.239 |
- ESM-2 pocket embedding (PCA 32, logistic) − pocket pLDDT, AUROC -0.163 [-0.243, -0.085]
- ESM-2 pocket embedding + pocket pLDDT − pocket pLDDT, AUROC -0.065 [-0.131, -0.002]
