# Crossing changes in alternating diagrams

This directory supports Sections 3–4 of *How to Compute the Unknotting Number: Theory and Computations*. It studies the knots obtained by changing each crossing of an alternating diagram, and it supplies McCoy tests and the conditional Bernhard–Jablan analysis.

The electrical description in Section 4 uses a positive-definite Goeritz matrix `G`. Changing an edge with incidence vector `x` gives `G - 2xxᵀ`; effective resistance is `xᵀG⁻¹x`. Determinant and linking-pairing updates are evaluated exactly. `sig.py` implements the Gordon–Litherland formula with exact rational inertia, including its type-II correction, so a floating-point eigenvalue tolerance cannot change a signature bound.

## Identification and inference

Canonical diagram codes, numerical SnapPy isometries against named KnotInfo diagrams, and identified factors in a diagrammatic connected-sum decomposition provide separate identification routes. Matching a determinant and a Jones polynomial is only a hypothesis. `candidates.py` keeps that hypothesis in a separate classification and prevents it from becoming a verified upper-bound witness. Likewise, a decomposition with only one nontrivial factor does not justify the lower bound for a nontrivial connected sum.

The Appendix F candidate classification is conditional on the Bernhard–Jablan conjecture and the recorded identification and bound premises. It does not establish `u = 3` unconditionally. A failure to find an unknotting crossing in a general diagram does not supply a lower bound. The special alternating case is treated separately by McCoy's theorem.

## McCoy verification

`mccoy_u1_check.py` checks the alternating-diagram hypothesis, changes every crossing, and decides determinant-one neighbors by simplification or knot Floer genus detection. An exception or unfinished recognition is recorded as undecided, never as an obstruction. The default targets are alternating knots with an open range beginning at 1 in the frozen `results/paper_v1_1_snapshot.json`.

```sh
python crossing_changes/mccoy_u1_check.py --out generated/mccoy_alternating_u1.json
python crossing_changes/mccoy_u1_check.py --names 13a_1069 --out generated/mccoy_subset.json
python crossing_changes/mccoy_u1_check.py --all /path/to/archived-wheel --out generated/mccoy_alternating_all.json
python -m unittest discover -s tests -p test_upper_bounds.py -v
```

The optional `--all` command reads the named archived April 2026 database wheel and is a separate historical comparison; it is not needed to reproduce the manuscript snapshot. Choose names actually present among the frozen targets for a bounded subset. New McCoy output goes to `generated/` unless an output path is explicitly supplied.

## Pipeline files

| File | Purpose |
|---|---|
| `crossing_dataset.py`, `augment.py` | Generate crossing-neighbor records and exact signature/resistance values |
| `identify_rest.py`, `relabel.py` | Identify endpoints and attach justified labels |
| `apply833.py`, `rank833.py` | Classify neighbors of open alternating knots and order candidate searches |
| `candidates.py` | Separate verified bounds, unresolved identifications, and polynomial hypotheses |
| `flype727.py` | Examine the flype orbit for the minimal-diagram analysis |
| `merge_tiers.py`, `resolve_unres.py` | Assemble conditional candidate tiers and revisit unresolved records |
| `lickorish_electrical.py`, `residual_breakdown.py` | Evaluate the linking-pairing obstruction and organize recorded exclusions |
| `analyze_v2.py` | Analyze the deposited data and fit the exploratory crossing-selection model |
| `mccoy_u1_check.py` | Apply the alternating-diagram theorem with explicit undecided outcomes |

`run_all.sh` is an exploratory rebuild that writes work files and depends on the supplied reference tables and installed database; it is not the default frozen-manuscript verification command. The retained script names reflect historical campaign sizes. The regression suite checks exact invariants, diagram moves, certificate rejection, and the distinction between verified identifications and polynomial hypotheses. It is not a fresh topological revalidation of every Appendix F entry.
