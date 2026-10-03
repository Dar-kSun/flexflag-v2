# flexflag

**Status: v0.1, in development.** Nothing here is a result yet.

A cheap triage flag that says "this protein target changes shape, so a single
predicted structure may mislead you," trained on proteins the PDB has solved
both empty (apo) and bound (holo).

flexflag flags risk. It does not predict alternative conformations, run docking,
or make any claim about drug efficacy.

The plan is to report, side by side, a trained model and the simple pLDDT-only
baseline. If pLDDT alone does as well, this README will say so in the results.
