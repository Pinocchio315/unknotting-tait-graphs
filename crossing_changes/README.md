# Crossing changes in alternating diagrams

This directory supports Sections 3–5 of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*, manuscript v1.8. It studies the knots obtained by changing each crossing of an alternating diagram, and it supplies McCoy tests and the minimal-diagram analysis in Section 5.2. These crossing-change calculations retain their frozen `v1.3` inputs; the [repository README](../README.md) explains its input manifests and the separate current KnotInfo comparison.

The electrical description in Section 4 uses a positive-definite Goeritz matrix `G`. Changing an edge with incidence vector `x` gives `G - 2xxᵀ`; effective resistance is `xᵀG⁻¹x`. Determinant and linking-pairing updates are evaluated exactly. `sig.py` implements the Gordon–Litherland formula with exact rational inertia, including its type-II correction, so a floating-point eigenvalue tolerance cannot change a signature bound.

## Identification and inference

Canonical diagram codes, numerical SnapPy isometries against named KnotInfo diagrams, and identified factors in a diagrammatic connected-sum decomposition provide separate identification routes. Matching a determinant and a Jones polynomial is only a hypothesis. `candidates.py` keeps that hypothesis in a separate classification and prevents it from becoming a verified upper-bound witness. Likewise, a decomposition with only one nontrivial factor does not justify the lower bound for a nontrivial connected sum.

For the 959 alternating knots in Appendix F, the recorded invariant tests and
verified identifications exclude a neighbour of unknotting number one after
every crossing change in a minimal diagram. Thus `u_BJ^s >= 3`, and any of
these knots with `u = 2` would fail the strong Bernhard–Jablan equality. Their
unknotting numbers remain in `[2,3]`; the analysis does not establish an
alternating counterexample. The 1,027 open alternating parents divide into
these 959 cases, 22 parents with 11 named possible children, and 46 parents
whose exclusions depend on unverified identifications. The latter 46 are
excluded from Appendix F.

A failure to find an unknotting crossing in a general diagram does not supply
a lower bound. McCoy's theorem supplies the additional conclusion for a
reduced alternating diagram.

## McCoy verification

`mccoy_u1_check.py` checks the alternating-diagram hypothesis, changes every crossing, and decides determinant-one neighbors by simplification or knot Floer genus detection. An exception or unfinished recognition is recorded as undecided, never as an obstruction. The default targets are alternating knots with an open range beginning at 1 in the frozen `results/paper_v1_1_snapshot.json`.

```sh
python crossing_changes/mccoy_u1_check.py --out generated/mccoy_alternating_u1.json
python crossing_changes/mccoy_u1_check.py --names 13a_1069 --out generated/mccoy_subset.json
python crossing_changes/mccoy_u1_check.py --all /path/to/archived-wheel --out generated/mccoy_alternating_all.json
python -m unittest discover -s tests -p test_upper_bounds.py -v
```

The deposited full computation covers all 3,627 alternating knots whose
April 2026 package lower bound was one. Exactly 710 have a recorded
unknotting crossing; the other 2,917 have none, agreeing with the lower
bounds two already recorded in June. The 27 changes in the historical
baseline are a subset of these 2,917 verifications and are not new lower
bounds relative to June. The separate crossing dataset covers 3,105 knots
in the cohort and agrees with every stored unknotting-crossing list in
that overlap; the other 522 rely on the deposited full McCoy computation.

The optional `--all` command reruns the computation from the named April wheel.
Choose names from the frozen targets for a bounded subset. New output goes to
`generated/` unless an explicit output path is supplied.

## Crossing-formula audit

From the repository root, run `python audit_crossing_formulas.py` to compare
the determinant, signature, dual-resistance, and self-linking formulas with
all 69,632 deposited crossing changes from 5,546 parent knots. The audit uses
exact rational arithmetic and reports zero mismatches for the frozen inputs.
It records hashes identifying the compared inputs and checks stored values;
it does not recompute diagrams or establish their knot identifications.
Output is printed unless `--out PATH` is explicitly supplied. The deposited
summary is `results/comparison/crossing_formula_verification.json`.

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
| `mccoy_u1_check.py` | Apply the alternating-diagram theorem with explicit undecided outcomes |

`run_all.sh` and `run_open23.sh` are exploratory rebuilds that overwrite files
under the deposited data paths and depend on historical reference tables and
the installed database. Run them only in a separate copy intended for new
computations; they are not frozen-manuscript verification commands. The retained
script names reflect historical campaign sizes. The regression suite checks
exact invariants, diagram moves, certificate rejection, and the distinction
between verified identifications and polynomial hypotheses. It is not a fresh
topological revalidation of every Appendix F entry.
