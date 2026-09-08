# Crossing changes in alternating diagrams

This directory supports Sections 3–5 of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*. It studies the knots obtained by changing each crossing of an alternating diagram, and it supplies McCoy tests and the conditional Bernhard–Jablan analysis in Section 5.2.

The electrical description in Section 4 uses a positive-definite Goeritz matrix `G`. Changing an edge with incidence vector `x` gives `G - 2xxᵀ`; effective resistance is `xᵀG⁻¹x`. Determinant and linking-pairing updates are evaluated exactly. `sig.py` implements the Gordon–Litherland formula with exact rational inertia, including its type-II correction, so a floating-point eigenvalue tolerance cannot change a signature bound.

## Identification and inference

Canonical diagram codes, numerical SnapPy isometries against named KnotInfo diagrams, and identified factors in a diagrammatic connected-sum decomposition provide separate identification routes. Matching a determinant and a Jones polynomial is only a hypothesis. `candidates.py` keeps that hypothesis in a separate classification and prevents it from becoming a verified upper-bound witness. Likewise, a decomposition with only one nontrivial factor does not justify the lower bound for a nontrivial connected sum.

The Appendix F candidate classification is conditional on the Bernhard–Jablan conjecture and the recorded identification and bound premises. It does not establish `u = 3` unconditionally. A failure to find an unknotting crossing in a general diagram does not supply a lower bound. The special alternating case is treated separately by McCoy's theorem.

The published `11a_14` figure can be checked with
`python audit_11a14_figure.py` from the repository root. This read-only
command replays three numerical exterior comparisons from deposited PD
codes: the left diagram represents `11a_14`, its boxed crossing change
leads to `8_8`, and the right diagram represents `10_129`, allowing mirrors.
SnapPy uses floating-point canonical isometries, not interval certification.
A separate exact Jones state sum on the small reference diagrams verifies
`V_8_8(t^-1) = V_10_129(t)`; it is not used for identification. The dated
reference values and the scope of this single-figure diagnostic are stored
in `results/comparison/11a14_figure_verification.json`.

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

`python audit_knotinfo_releases.py --verify` checks this comparison from
frozen public range extracts and deposited records, without restarting a
search. The optional `--all` command above instead reruns the computation
from the named April wheel. Choose names actually present among the frozen
targets for a bounded subset. New McCoy output goes to `generated/` unless
an output path is explicitly supplied.

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

`u1_minimal_diagram_scan.py` tests every crossing of one tabulated minimal
diagram per knot. The archived positive cohort has 1,516 knots and the open
cohort has 546. A negative result applies only to the chosen diagram; the
program does not enumerate other minimal diagrams of a non-alternating knot.
Determinant-one changes are decided by simplification or by knot Floer
homology's Seifert-genus detection, with failures left undecided.

From the repository root:

```sh
python crossing_changes/u1_minimal_diagram_scan.py 3_1 11n_116
python crossing_changes/u1_minimal_diagram_scan.py --cohort open --out /tmp/open_diagrams.json
```

The program prints results unless an explicit `--out` path is supplied. New
records include the PD and per-crossing evidence; the existing archives are
not replaced by a named diagnostic run. The independent verification in
`results/review_v1_3/u1_scanner_audit.json` checks one positive witness for
every one of the 1,516 known knots and all 7,052 crossing changes of the 546
chosen open-range diagrams.

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
| `u1_minimal_diagram_scan.py` | Inspect one chosen minimal diagram per knot, with explicit determinant and unknot-recognition evidence |

`run_all.sh` is an exploratory rebuild that writes work files and depends on the supplied reference tables and installed database; it is not the default frozen-manuscript verification command. The retained script names reflect historical campaign sizes. The regression suite checks exact invariants, diagram moves, certificate rejection, and the distinction between verified identifications and polynomial hypotheses. It is not a fresh topological revalidation of every Appendix F entry.
