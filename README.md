# [Computation of Unknotting Numbers: Which Knot Breaks the Bernhard-Jablan Conjecture?][paper]

This repository contains the code, input data, and computational results for
all calculations reported in the paper by **Seong-Jin Lee**. It includes
lower-bound obstructions, searches for crossing-change constructions,
verification of the resulting diagrams, and the analysis of minimal diagrams.
The paper's arXiv link will be added after announcement.

## Results

For prime knots with at most 13 crossings, the computations determine **2,525
unknotting numbers that were undetermined in the official KnotInfo snapshot
retrieved on 9 September 2026**. Relative to the archived starting ranges,
there are 2,719 exact determinations and 390 further improved ranges.

| Unknotting number | Exact determinations from the archived starting ranges |
| --- | ---: |
| 2 | 2,149 |
| 3 | 373 |
| 4 | 172 |
| 5 | 25 |
| **Total** | **2,719** |

Of these 2,719 values, 194 are already recorded in the dated official snapshot
and 2,525 are still open. Of the 390 improved ranges, 333 are narrower than
that snapshot and 57 are already recorded. Lower-bound obstructions account
for 2,711 exact values; eight explicit diagrams give the matching upper bounds
for the other eight. Greene's model supplies 1,492 additional exact values in
the enlarged sweep, or 1,494 including `12n_491` and `13n_3370`.

The calculations establish

```text
u(12n_288) = u(12n_491) = u(12n_501) = u(13n_3370) = 2,
u_BJ^s(13n_3370) = 3.
```

Together with [Brittenham-Hermiller's minimal-diagram enumeration, Theorem 1.3](https://arxiv.org/abs/1705.05985v2),
this identifies **13n3370** as an explicit counterexample to the original
Bernhard-Jablan conjecture. Montesinos obstructions handle `12n_288` and
`12n_501`; Greene's model supplies the correction terms needed for `12n_491`
and `13n_3370`. The published upper bounds are recorded separately.
The same method determines `u(13n_1587) = 2`, without deciding its status
under the strong Bernhard-Jablan equality.

For 959 alternating knots with range `[2,3]`, the recorded invariant tests and
verified knot identifications exclude a neighbour of unknotting number one
after every crossing change in a minimal diagram. Any member with `u=2` would
be an alternating counterexample; these conditional cases are separate from
the exact-value count.

## Repository structure

```text
.
├── consolidate_results.py        # Validate and combine bounds; JSON output
├── compare_knotinfo.py           # Compare intervals with the dated official ranges
├── audit_crossing_formulas.py    # Exact checks of recorded crossing-change formulas
├── requirements.txt
├── lower_bounds/
│   ├── linking_pairing.py        # Finite-group and linking-pairing calculations
│   ├── cyclic_cover_bound.py    # Generator bounds from cyclic branched covers
│   ├── owens/                   # Definite surgery forms in ranks 2, 3, and 4
│   ├── montesinos/              # Plumbing correction terms and reduced Floer homology
│   └── greene/                  # Kauffman states, correction terms, surgery tests
├── upper_bounds/                # Searches and crossing-change certificate verification
│   ├── xtait/                   # Diagram moves and reductions
│   ├── data/                    # Search inputs and reference diagrams
│   └── cluster/                 # Distributed searches and certificate replay
├── crossing_changes/            # Minimal-diagram neighbours and McCoy's criterion
├── tait_graphs/tait/             # Embedded graphs and topological invariants
├── results/                     # Frozen inputs, obstruction records, and certificates
└── tests/                       # Mathematical kernels and bound validation
```

The [results guide](results/README.md) identifies the datasets and their roles.
The subdirectory READMEs explain the algorithms, assumptions, and command-line
options. LaTeX sources, table generators, and manuscript-editing tools are not
part of this distribution.

## Reconstruct the computed intervals

This step requires only **Python 3.10 or later**. It uses deposited results,
without rerunning the topological computations or accessing the network.

```bash
git clone https://github.com/Pinocchio315/unknotting-tait-graphs.git
cd unknotting-tait-graphs
python3 consolidate_results.py
python3 compare_knotinfo.py --out /tmp/knotinfo_comparison.json
```

The first command writes three JSON files in the ignored `generated/` directory:

| Output | Contents |
| --- | --- |
| `u_table.json` | `[lower, upper]` for all 12,965 study knots |
| `consolidated.json` | Changed intervals, source records, and method tags |
| `counts.json` | Exact values and improved ranges, counted once per knot |

Input checksums are verified before bounds are combined. Only completed
obstructions raise lower bounds. Passing a test or reaching a resource limit
supplies no bound; upper bounds come from separate constructions.
Contradictory intervals stop the run.

Use `--outdir PATH` to keep separate outputs. The `pre-sweep` profile
reconstructs the ranges used to select the enlarged Greene cohort:

```bash
python3 consolidate_results.py --profile pre-sweep --outdir /tmp/unknotting-before-greene
```

The default `full` profile includes that sweep. The input manifests retain
their historical filenames (`paper_v1_1_manifest.json` and
`paper_v1_3_manifest.json`): they authenticate computational data, and no
manuscript is needed to use them.

## Environment for mathematical computations

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The main packages are SnapPy, Spherogram, Regina, knot_floer_homology, NumPy,
SymPy, and NetworkX. `database_knotinfo==2026.8.1` supplies the frozen diagram
and invariant inputs; it is distinct from the dated official range comparison.
Sage is not required. Slurm is optional for the distributed searches.

The recorded local verification environment used Python 3.13.7, SnapPy 3.3.2,
Spherogram 2.4.1, Regina 7.4.1, SymPy 1.14.0, NumPy 2.5.2, NetworkX 3.6.1,
and knot_floer_homology 1.2.2. This is not a complete environment manifest for
every historical cluster run. Only the KnotInfo input package is pinned.

## Run or verify individual calculations

Run these commands from the repository root with the environment active:

```bash
# Small rank-three and rank-four examples.
python lower_bounds/owens/owens_u3.py 7_1
python lower_bounds/owens/owens_u4.py 9_1

# Correction terms and surgery obstructions for the Greene/BJ knots.
python lower_bounds/greene/pin_dinv.py 12n_491 13n_3370

# Independently replay the quasi-alternating certificate.
python lower_bounds/greene/qa_search.py --verify results/bernhard_jablan/qa_certificate_12n_491_2026-09-08.json

# Replay the eight displayed crossing-change diagrams.
python upper_bounds/verify_presentation_diagrams.py

# Check all 69,632 recorded crossing changes with exact rational arithmetic.
python audit_crossing_formulas.py

# Apply McCoy's criterion to one reduced alternating diagram.
python crossing_changes/mccoy_u1_check.py --names 13a_1069 --out /tmp/mccoy_example.json
```

See the [lower-bound guide](lower_bounds/README.md), [Greene guide](lower_bounds/greene/README.md),
[upper-bound guide](upper_bounds/README.md), and [crossing-change guide](crossing_changes/README.md)
for full sweeps and their inputs. New runs should use separate output paths.
Large enumerations can take hours or longer.

The upper-bound search uses randomized diagram moves, permits up to 50
crossings, and defaults to at most 10,000 additional distinct diagrams per
target. Crossing subsets are filtered by determinant before simplification;
when needed, a random subset is selected within the reduction budget.
An unsuccessful search supplies no lower bound. Certificate verifiers retain
exact move replay and knot Floer genus detection; named hyperbolic knots also
use numerical SnapPy isometries.

## Tests

```bash
python -m unittest discover -s tests -v
python upper_bounds/selftest.py
python lower_bounds/owens/validate.py
```

These commands exercise mathematical kernels, known examples, diagram moves,
certificate rejection, and interval assembly. They do not rerun every stored
obstruction or the full search campaigns.

## Paper and citation

The arXiv identifier will be inserted after announcement. The title above
uses the reference link below so that its destination can be updated once.

```bibtex
@misc{Lee2026,
  author = {Lee, Seong-Jin},
  title = {Computation of Unknotting Numbers: Which Knot Breaks the Bernhard--Jablan Conjecture?},
  year = {2026},
  note = {Accompanying code and computational data},
  url = {https://github.com/Pinocchio315/unknotting-tait-graphs}
}
```

[paper]: #paper-and-citation
