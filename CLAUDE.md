# CLAUDE.md — `flexflag`

> Before you dock into a predicted structure, find out whether the protein holds still.
> **Mission: a cheap triage flag that says "this target changes shape — one structure will mislead you here," trained on proteins the PDB has solved both empty and bound.**

---

## 0. Context for you, the agent

Built by Aryan Sinha (3rd-year B.Tech Biotechnology, IIT Kharagpur), to be read by professional research reviewers.

- **Honest beats impressive.** For this project in particular: if the simple existing heuristic (pLDDT alone) predicts flexibility as well as the trained model, **say so prominently**. A clean negative result on the field's default heuristic is a real contribution and reads as judgment.
- Every README number reproducible from a repo command.
- This is the heaviest of the three projects — structure data wrangling eats time. **Scope aggressively and say what was scoped out.**

Small commits, real messages, never backdated.

---

## 1. The problem, stated precisely

Proteins are the cell's machines; function follows 3D shape. AlphaFold predicts that shape from sequence, and hundreds of millions of predicted structures are now freely available.

The dominant downstream use is **docking**: take the structure, computationally screen millions of candidate molecules against its binding pocket, shortlist the ones that fit. It's the first step of most computational drug discovery.

**The problem.** Proteins aren't statues. They flex and often reorganise substantially on binding — pockets open that weren't there, flaps close. AlphaFold returns **one** structure: a single frozen frame.

And it isn't a neutral frame. Training data (the PDB) is skewed toward structures solved *with* a ligand bound. The documented consequence: for E3 ubiquitin ligases — genuinely open when empty, closed when occupied — AlphaFold predicts the **closed** form for both. It hands you, confidently, the wrong shape for the empty protein, which is exactly the shape you need when searching for something to put in it.

The usual safety net doesn't catch this. **pLDDT says "I'm confident about this shape." It does not say "and this protein only has one shape."** A structure can be high-confidence and still be the wrong conformational state for your purpose. Separately, AF3 shows a **4.4% chirality violation rate** on the PoseBusters benchmark and can hallucinate order in genuinely disordered regions.

So teams dock into one confident-looking structure, rank hits, order compounds, and discover months and lakhs of rupees later that they were docking into a shape the empty protein never adopts.

**The gap.** There is good research on sampling multiple conformations (AlphaFold2-RAVE, MD, ensemble methods). All of it is expensive and expert-only. What doesn't exist is the **cheap triage step at the front**: before spending anything, a flag saying *this target is shape-shifty*. That's a far easier question than predicting the full ensemble — and it's the one that changes what someone does on Monday morning.

**What this repo is.** That flag.

---

## 2. Scope

### v0.1 — must ship
1. A curated set of **apo–holo pairs** (same protein solved empty and bound) with a measured conformational-change label.
2. Features computed from **public, free sources** — AlphaFold DB confidence outputs and sequence-derived properties.
3. A calibrated gradient-boosting classifier predicting "large conformational change."
4. **An explicit, prominent comparison against the pLDDT-only baseline.**
5. CLI: UniProt ID or sequence in → risk band + reasons out. README with real numbers.

### v0.2 — if time allows
6. SHAP-based per-prediction explanations in the CLI output.
7. Family-level analysis: which protein families are most shape-shifty.
8. Pocket-specific rather than whole-protein flags.

### Explicitly out of scope
- Predicting the actual alternative conformation. **This repo flags risk; it does not model ensembles.**
- Running docking.
- Any claim about drug efficacy or clinical relevance.
- Molecular dynamics.

---

## 3. Feasibility first — do this before anything else

Structure data is where this project can quietly consume two days. Before building:

1. Pick **one** well-curated apo–holo source and commit to it. Candidates: **APObind** (packaged apo conformations for PDBbind complexes), PDBFlex (conformational variability), or a SIFTS/UniProt-mapped PDB subset built yourself. **Prefer a pre-packaged set** — building pairs from scratch is a multi-day job on its own.
2. Verify you can download it and that it has enough pairs (target: ≥500 usable; if fewer, scope down the claims, not the rigour).
3. Verify AlphaFold DB entries (PDB + confidence files) are retrievable for those UniProt IDs.
4. Write the decision into `docs/data-choice.md` with counts and any filtering.

Write `scripts/00_feasibility.py` that fetches 10 pairs end to end, computes one label and one feature, prints timings, and exits. **Commit this before writing the model.**

---

## 4. Labels — the part that must be right

The label is "how much does this protein change shape between empty and bound."

- **Primary label**: RMSD between apo and holo structures after superposition, computed on **binding-site residues** (not whole-chain — whole-chain RMSD is dominated by loops and domain motions that may be irrelevant to docking).
- Define the binding site as residues within a fixed cutoff (e.g. 5 Å) of the ligand in the holo structure. Put the cutoff in config.
- **Secondary label**: pocket volume change, if a pocket-detection tool is available cheaply. Optional for v0.1.
- Binarise at a stated threshold (e.g. > 2 Å site RMSD = "large change"). **Report the threshold everywhere, and show how results change at a couple of alternatives** — a result that only holds at one arbitrary cutoff is not a result.

Pitfalls to handle explicitly, and to document:
- Superposition method changes RMSD. Fix one, state it.
- Multiple chains, altlocs, missing residues, differing residue numbering between apo and holo entries. Write the alignment carefully and test it.
- Crystal-packing artefacts can inflate apparent motion. Note as a limitation.

---

## 5. Features

All from free sources, all computed in `flexflag/features/`:

**From AlphaFold DB:**
- pLDDT statistics: mean, min, std, fraction below 70, longest contiguous low-confidence run.
- PAE matrix structure: inter-domain PAE blocks. A rigid pair of domains joined by a flexible hinge leaves a characteristic block signature — this is the most promising feature and worth implementing carefully.

**From sequence:**
- Length, predicted disorder fraction, low-complexity content, composition.
- MSA depth, if cheaply available.

**Annotations:**
- Domain count, Pfam family, whether multi-domain.

**Feature hygiene:** no feature may encode the answer. Specifically, nothing derived from the holo structure or the ligand may enter the feature set — that's leakage, and a reviewer will look for it. Write a test asserting the feature pipeline never touches holo files.

---

## 6. Evaluation — leakage is the main risk

Homologous proteins are near-duplicates. A random split will produce a beautiful, meaningless AUROC.

- **Cluster by sequence identity** (MMseqs2 or CD-HIT at e.g. 30%) and split **by cluster**. If clustering tools are unavailable, split by Pfam family and say so.
- Assert in code that no cluster appears on both sides; raise on violation.
- Report AUROC **and** AUPRC with bootstrap CIs over clusters.
- **Calibration matters more than discrimination here** — the output is a risk band someone acts on. Report expected calibration error and a reliability diagram.

### The comparison that defines this repo

Report, side by side:
1. pLDDT-only baseline (e.g. mean pLDDT, or fraction below 70, as a single feature).
2. The full model.

If the full model's gain over pLDDT alone is small or within CI, **the README says that in the results section, not buried in limitations.** That honesty is worth more than a fake win, and experienced reviewers will notice which way you went.

---

## 7. Repo layout

```
flexflag/
  __init__.py
  data/
    apoholo.py       # fetch + pair curation
    alphafold.py     # AFDB structure + confidence retrieval
  labels.py          # site RMSD, superposition, binarisation
  features/
    plddt.py
    pae.py
    sequence.py
  splits.py          # cluster-aware splitting + leakage assertion
  model.py           # gradient boosting + calibration
  explain.py         # SHAP (v0.2)
  cli.py             # typer: flexflag check <uniprot_id>
scripts/
  00_feasibility.py
  01_build_dataset.py
  02_train_eval.py
tests/
  fixtures/          # 2-3 small structures
  test_no_holo_leakage.py
docs/
  data-choice.md
  findings.md
results/             # committed metrics + figures
README.md
pyproject.toml
```

---

## 8. Build order

- **M0** — Skeleton, pyproject, ruff, pytest, CI, MIT licence, `.gitignore` (exclude structure downloads), README stub "v0.1, in development."
- **M1** — `00_feasibility.py` passing on 10 pairs. Commit `docs/data-choice.md`.
- **M2** — Label pipeline + tests on committed fixtures. Verify a known flexible protein scores high and a known rigid one scores low; write that as a test.
- **M3** — Dataset build across the full pair set. Report pair count, label distribution, and what was filtered out, in `docs/findings.md`.
- **M4** — `splits.py` with cluster-aware splitting and leakage assertion. **Before any training.**
- **M5** — pLDDT-only baseline. Get this number first, so the model has a real bar.
- **M6** — Feature pipeline + full model + calibration. Compare against M5.
- **M7** — CLI + README with real numbers and the honest comparison. Tag `v0.1.0`.
- **M8+** (v0.2) — SHAP explanations, family analysis, pocket-level flags.

---

## 9. README requirements

1. One sentence on what and why, with the headline number (e.g. "X% of targets in this set shift more than 2 Å between empty and bound").
2. Signature figure: distribution of conformational change across targets, or the calibration curve.
3. Quickstart: install, then `flexflag check P12345` producing a risk band and reasons.
4. How the label is defined, including the threshold and sensitivity to it.
5. Results: model vs pLDDT baseline, AUROC/AUPRC with CIs, calibration error. **The comparison is the result, whichever way it goes.**
6. **Limitations**: dataset size and source, crystal-packing caveats, what "large change" does and doesn't mean for docking, families not represented, and the fact that this flags risk rather than predicting conformations.
7. Status: "v0.1 — active development," with what's scoped out and why.
8. Citations for AlphaFold limitations, APObind, and the ensemble-method prior work.

---

## 10. Honesty checklist

- [ ] No holo-derived information in features (test enforces it).
- [ ] Splits cluster-aware; leakage assertion active.
- [ ] pLDDT baseline reported prominently, whatever the outcome.
- [ ] Calibration reported, not just discrimination.
- [ ] Label threshold stated, with sensitivity shown.
- [ ] Prior ensemble work credited; this is positioned as cheap triage *in front of* those methods, not a replacement.
- [ ] No efficacy, clinical, or drug-discovery outcome claims.
- [ ] Dataset size and composition stated wherever results appear.

---

## 11. Definition of done for v0.1

A reviewer clones, runs `flexflag check` on a UniProt ID, and gets a calibrated risk band with a human-readable reason — backed by a README that honestly reports whether the model beats simply looking at pLDDT.
