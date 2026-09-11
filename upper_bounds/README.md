# Upper bounds: diagrams and crossing changes

This directory contains the upper-bound constructions in Section 4 and the published bounds replayed in Section 5 of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*, manuscript v1.8. An upper bound is established by an authenticated diagram of the named knot, an explicit sequence of isotopies, and distinct marked crossing changes. The endpoint must be the unknot or a knot with an independently justified upper bound.

The search programs propose such constructions. `run_search.py` uses random walks in embedded Tait graphs, while `cluster/enum_search.py` enumerates diagrams reachable within a chosen crossing budget. `flipdet.py` screens crossing subsets with the matrix determinant lemma modulo a prime. A matching determinant is only a filter: it neither identifies a knot nor detects the unknot. Search exhaustion within a finite move budget gives no lower bound on the unknotting number.

## Results and deposited constructions

The nineteen constructions in Section 4.3 establish `u=2`. Eight earlier
presentation PDs remain in `results/best/`; eleven additions are in
`results/extensions_2026-09-10/presentation.json`. Ten of the latter require
the new Greene lower bound two; 13n3733 already had that lower bound.
The nineteen PDs comprise fifteen with 13 crossings, three with 14 crossings,
and one with 18 crossings (13n4237).

Eighteen PDs have two marked crossings giving the unknot. The PD for 13n447
has one mark giving 10_91, whose unknotting number is one. The figures and
Appendix G use this distinction explicitly. For 13n3108 and 13n3733 the
original search certificates used a named partner; the presentation PDs
now use independently verified direct pairs. Both certificates are retained.

The manuscript's search description applies to all nineteen constructions.
Actual execution settings and reduction provenance are recorded below and in
the [v1.8 additions guide](../results/extensions_2026-09-10/README.md).

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

`verify_presentation_diagrams.py` defaults to v1.8 and requires all nineteen diagrams. It checks eighteen unknot endpoints and the identified partner 10_91 with its independent upper bound one. Use `--manuscript v1.3` for the original eight. `verify_braid_certificate.py` checks the braid construction for `13n_1587`, including an independent unknotting crossing of its partner `10_113`. `verify_certificates.py --paper` replays the eight deposited campaign certificates; the flag excludes separate archived constructions. Omitting required diagrams, repeating marks, recording an uncounted crossing change as an isotopy, or claiming an unjustified partner bound causes verification to fail.

`verify_bj_upper_bounds.py` supports Section 5.1: it replays the published
Brittenham–Hermiller DT crossing changes for `12n_491` and `13n_3370`, and
checks an unknotting crossing of each partner (`6_3` and `11n_21`). These
existing upper bounds are combined with the new Greene lower bounds; they
are not new upper-bound constructions of this paper.

For `verify_certificates.py`, partner bounds default to the frozen manuscript comparison in `results/paper_v1_1_snapshot.json`; its `--known-bounds FILE` option supplies a different explicit source. The verifier does not infer exact unknotting numbers from certificate assertions. The presentation verifier requires actual returned SnapPy isometries that extend over the meridians; it also uses exact knot Floer genus detection for the unknot endpoints. SnapPy documents positive isometry answers as rigorous. A failure to find an isometry is not evidence that two knots differ. The [SnapPy API documentation](https://snappy.computop.org/manifold.html#snappy.Manifold.is_isometric_to) describes that distinction.

## Running a search

From the repository root, a small exploratory run uses:

```sh
python upper_bounds/run_search.py --names 13n1587 --diagrams 100 --n-target 30 --time-per-knot 30 --workers 1 --seed 1 --out-dir /tmp/unknotting_search_demo
```

Names are selected from `upper_bounds/data/targets.json` by default;
`--targets FILE` selects an explicit frozen list. `--names` filters that list
and does not create a new target. The companion `known_u.json` supplies canonical diagram
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

## Recorded v1.8 run

The additional 190-target run used the same search implementation with
`--diagrams 10000 --n-target 25 --time-per-knot 600 --workers 3 --seed 0`.
The original execution replaced the target file in a working copy; the new
`--targets` option reads the same deposited list without changing old inputs:

```sh
python upper_bounds/run_search.py --targets results/extensions_2026-09-10/targets190.json --diagrams 10000 --n-target 25 --time-per-knot 600 --workers 3 --seed 0 --out-dir /tmp/unknotting_190
python upper_bounds/verify_certificates.py results/extensions_2026-09-10/runs190 --out /tmp/verified_v18
```

This is a full campaign command, not a short test. The 190 recorded rows
contain 11 successes and 1,799,825 diagram visits. A row can contain 10,001
visits because the initial diagram is tested before up to 10,000 additional
diagrams. The sum of recorded per-knot elapsed times is 37,684.2 seconds;
this is neither measured CPU time nor the wall time of the parallel run.

For compact presentation, marked diagrams are reduced by the existing move
engine and safe pass moves. `children_pass.py` only uses an over/under run
containing no marked crossing, so the move commutes with the crossing changes.
It tracks the surviving PD rows explicitly and rejects a changed determinant.
That determinant guard is only a consistency check; the endpoint is verified
separately. `minimize_witness_pass.py` stores the path it produces. The original
`min11/` records retain only path lengths, while the later 13n436 pass reduction
has a complete five-move path and a deposited replay.

```sh
python upper_bounds/minimize_witness_pass.py results/extensions_2026-09-10/min11/13n_436.json --endpoint unknot --out /tmp/13n_436_reduced.json
python upper_bounds/draw_witness.py --outdir /tmp/unknotting_figures
```

All nineteen drawings use the same orthogonal layout renderer, black strands,
underpass gaps, and red circles. PD row numbers are preserved as crossing
labels so the renderer can require every mark to be drawn.

## Main files

| File | Purpose |
|---|---|
| `expand.py`, `run_search.py` | Search equivalent diagrams and screen marked crossings |
| `flipdet.py` | Modular determinant updates and an exact determinant cross-check |
| `certificate_checks.py` | Validate marked crossings and independent partner bounds |
| `verify_certificates.py` | Authenticate the source and replay the recorded construction |
| `verify_presentation_diagrams.py` | Check all nineteen presentation diagrams and their endpoints |
| `verify_braid_certificate.py` | Check the braid construction and the partner's unknotting crossing |
| `verify_bj_upper_bounds.py` | Replay the two published upper bounds used in Section 5.1 and verify the partners' unknotting crossings |
| `identify_candidates.py`, `resolve_unidentified.py`, `rehunt_certificates.py` | Identify proposed endpoints and recover replayable routes |
| `minimize_witness.py`, `minimize_witness_pass.py`, `children_pass.py`, `reduce_witnesses.py`, `draw_witness.py` | Simplify and draw marked diagrams |
| `xtait/` | Embedded graph, isotopy, reduction, and flype kernels |
| `bj_scan.py` | Auxiliary crossing-neighbour scan using historical selected `u=2` and `u=3` input lists; it does not regenerate the v1.8 consolidated analysis |
| `make_targets.py`, `build_all_codes.py` | Build exploratory search data from an installed KnotInfo package |

The search-data builders read the installed database; they do not regenerate the frozen comparison used in the paper. Named endpoints are compared with KnotInfo PD codes; a table label or a Jones polynomial alone is insufficient for identification. Historical scripts and data remain archived alongside the paper constructions.

The historical BJ scanner uses lower bound two for an identified sum of two
nontrivial knots. Scharlemann's theorem does not raise this bound to three
when one summand has unknotting number greater than one. Thus factor values
one and two give `[2,3]` unless a separate obstruction applies.
This auxiliary scan supplies no input to the v1.8 result tables.
