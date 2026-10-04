Docking validation (docs/plan-v0.6-docking.md). AutoDock Vina 1.2.7, success = top pose within 2 Å of the crystal pose.

### external: 145 proteins docked in all three receptors
(excluded before docking: {'covalent ligand': 25, 'ligand too large or flexible': 9})

| Receptor | Success (top 1) | Success (any of top 3) |
|---|---|---|
| holo | 27.6% | 47.6% |
| apo | 3.4% | 8.3% |
| alphafold | 9.7% | 15.2% |

| Pocket | n | AlphaFold fails | AlphaFold fails while holo succeeds | apo fails |
|---|---|---|---|---|
| moving (> 2 Å) | 14 | 92.9% | 35.7% | 100.0% |
| not moving | 131 | 90.1% | 19.1% | 96.2% |

| Pocket-pLDDT band | n | AlphaFold fails | AlphaFold fails while holo succeeds |
|---|---|---|---|
| < 91.8 | 20 | 90.0% | 25.0% |
| 91.8–94.6 | 29 | 89.7% | 24.1% |
| 94.6–96.4 | 28 | 92.9% | 25.0% |
| 96.4–97.9 | 45 | 91.1% | 20.0% |
| > 97.9 | 23 | 87.0% | 8.7% |

AUROC for AlphaFold docking failure (95% cluster bootstrap; positives/n):
- pocket pLDDT, all proteins: 0.495 [0.309, 0.663] (131/145)
- **pocket pLDDT, where holo docking succeeded (primary):** 0.593 [0.379, 0.793] (30/40)
- whole-protein pLDDT, where holo docking succeeded: 0.587 [0.359, 0.809] (30/40)
- measured apo–holo pocket RMSD, where holo succeeded (oracle): 0.627 [0.429, 0.813] (30/40)
