# Upper bounds: diagrams and crossing changes

This directory contains the upper-bound constructions in Section 4 and the published bounds replayed in Section 5 of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*, manuscript v1.7. An upper bound is established by an authenticated diagram of the named knot, an explicit sequence of isotopies, and distinct marked crossing changes. The endpoint must be the unknot or a knot with an independently justified upper bound.

The search programs propose such constructions. `run_search.py` uses random walks in embedded Tait graphs, while `cluster/enum_search.py` enumerates diagrams reachable within a chosen crossing budget. `flipdet.py` screens crossing subsets with the matrix determinant lemma modulo a prime. A matching determinant is only a filter: it neither identifies a knot nor detects the unknot. Search exhaustion within a finite move budget gives no lower bound on the unknotting number.

## Results and deposited constructions

The eight new constructions establish `u = 2` for `13n_30`, `13n_45`,
`13n_80`, `13n_2379`, `13n_2809`, `13n_2907`, `13n_3033`, and `13n_3589`.
Each has lower bound two independently of this search. Five presentation
diagrams have 13 crossings and three have 14 crossings; their PD codes and
marked crossings appear in Appendix G. These eight exact values are separate
from the 2,711 obtained by lower bounds in the paper.

The manuscript v1.7 retains the frozen computational profile `v1.3`.
Verification of its deposited constructions can be run independently of the
search that found them. The [repository README](../README.md) describes the
input manifests and the current KnotInfo comparison.

## Verification

Run from the repository root in the topology environment described in the top-level README:

```sh
python upper_bounds/selftest.py
python upper_bounds/verify_presentation_diagrams.py
python upper_bounds/verify_braid_certificate.py
python upper_bounds/verify_bj_upper_bounds.py
python upper_bounds/verify_certificates.py results/witness_chains --paper --out generated/verified_witness_chains
python -m unittest discover -s tests -p test_upper_bounds.py -v
```

`verify_presentation_diagrams.py` requires all eight marked diagrams used in the paper and verifies that changing the two marked crossings gives the unknot. `verify_braid_certificate.py` checks the braid construction for `13n_1587`, including an independent unknotting crossing of its partner `10_113`. `verify_certificates.py --paper` replays the eight deposited campaign certificates; the flag excludes separate archived constructions. Omitting required diagrams, repeating marks, recording an uncounted crossing change as an isotopy, or claiming an unjustified partner bound causes verification to fail.

`verify_bj_upper_bounds.py` supports Section 5.1: it replays the published
Brittenham–Hermiller DT crossing changes for `12n_491` and `13n_3370`, and
checks an unknotting crossing of each partner (`6_3` and `11n_21`). These
existing upper bounds are combined with the new Greene lower bounds; they
are not new upper-bound constructions of this paper.

For `verify_certificates.py`, partner bounds default to the frozen manuscript comparison in `results/paper_v1_1_snapshot.json`; its `--known-bounds FILE` option supplies a different explicit source. The verifier does not infer exact unknotting numbers from certificate assertions. The scripts use exact combinatorial replay and knot Floer genus detection, together with numerical SnapPy isometries for named hyperbolic knots. These are reproducibility checks, not formal proof certificates for the numerical isometry step.

## Running a search

From the repository root, a small exploratory run uses:

```sh
python upper_bounds/run_search.py --names 13n1587 --diagrams 100 --n-target 30 --time-per-knot 30 --workers 1 --seed 1 --out-dir /tmp/unknotting_search_demo
```

Names are selected from `upper_bounds/data/targets.json`; `--names` does not
create a new target. The companion `known_u.json` supplies canonical diagram
codes and partner bounds for the exploratory search. Final certificate
verification checks those bounds against an explicit independent source.

The randomized search first tests the starting diagram, then generates up to
10,000 additional canonically distinct diagrams by Reidemeister II and III
moves and flypes. Its default crossing cutoff is 50. These moves preserve
crossing-number parity, so the effective target size for a 13-crossing start
is 49. The determinant screen visits all crossing subsets of the prescribed
sizes in each visited diagram. Reduction is restricted to a random sample if
too many subsets pass that screen.

| Option | Default | Meaning |
| --- | ---: | --- |
| `--diagrams` | 10000 | Maximum number of additional diagrams emitted by the random walk |
| `--n-target` | 50 | Crossing cutoff, adjusted to the starting parity |
| `--time-per-knot` | 1800 | Time budget in seconds, checked between search operations |
| `--full-per-diagram` | 15 | Maximum determinant-one subsets passed to reduction for a direct unknotting attempt |
| `--partial-per-diagram` | 6 | Maximum determinant-matching subsets passed to reduction for each smaller number of changes |
| `--seed` | 0 | Seed for diagram generation |
| `--workers` | 1 | Number of target knots processed concurrently |

Setting either subset cap to zero removes that cap. The time budget is not a
hard process timeout: a single diagram-generation or subset-enumeration step
can run beyond it. Other reduction limits are listed by `--help`. A run may
stop before the diagram limit because it succeeds, exhausts a move budget,
or reaches the time budget; its failure supplies no lower bound.

Each output directory receives `results_0_1.jsonl` for the default shard and
`SUCCESS_<name>.json` for each proposed construction. Result rows record the
numbers of diagrams, subsets, determinant survivors and reductions actually
processed. A rerun in the same directory skips names with completed result
rows. Use a new directory to compare seeds or budgets, and replay successful
constructions with `verify_certificates.py` before using them as bounds.

## Main files

| File | Purpose |
|---|---|
| `expand.py`, `run_search.py` | Search equivalent diagrams and screen marked crossings |
| `flipdet.py` | Modular determinant updates and an exact determinant cross-check |
| `certificate_checks.py` | Validate marked crossings and independent partner bounds |
| `verify_certificates.py` | Authenticate the source and replay the recorded construction |
| `verify_presentation_diagrams.py` | Check the eight diagrams displayed in the manuscript |
| `verify_braid_certificate.py` | Check the braid construction and the partner's unknotting crossing |
| `verify_bj_upper_bounds.py` | Replay the two published upper bounds used in Section 5.1 and verify the partners' unknotting crossings |
| `identify_candidates.py`, `resolve_unidentified.py`, `rehunt_certificates.py` | Identify proposed endpoints and recover replayable routes |
| `minimize_witness.py`, `reduce_witnesses.py` | Simplify marked diagrams while preserving their certificates |
| `xtait/` | Embedded graph, isotopy, reduction, and flype kernels |
| `bj_scan.py` | Auxiliary crossing-neighbour scan using historical selected `u=2` and `u=3` input lists; it does not regenerate the v1.6 consolidated analysis |
| `bj_check_13a647.py` | Historical diagram diagnostic; the current value is `u(13a647)=4`, and this script supplies no alternating counterexample |
| `make_targets.py`, `build_all_codes.py` | Build exploratory search data from an installed KnotInfo package |

The search-data builders read the installed database; they do not regenerate the frozen comparison used in the paper. Named endpoints are compared with KnotInfo PD codes; a table label or a Jones polynomial alone is insufficient for identification. Historical scripts and data remain archived alongside the paper constructions.

The historical BJ scanner uses lower bound two for an identified sum of two
nontrivial knots. Scharlemann's theorem does not raise this bound to three
when one summand has unknotting number greater than one. Thus factor values
one and two give `[2,3]` unless a separate obstruction applies. The older
`13a647` diagnostic now reports only its diagram checks: its former
counterexample announcement depended on the incorrect premise `u(13a647)=3`.
Neither auxiliary script supplies an input to the v1.6 result tables.
