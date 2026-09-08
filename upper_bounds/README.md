# Upper bounds: diagrams and crossing changes

This directory supports Section 4 of *How to Compute the Unknotting Number: Theory and Computation*. An upper bound is established by an authenticated diagram of the named knot, an explicit sequence of isotopies, and distinct marked crossing changes. The endpoint must be the unknot or a knot with an independently justified upper bound.

The search programs propose such constructions. `run_search.py` uses random walks in embedded Tait graphs, while `cluster/enum_search.py` enumerates diagrams reachable within a chosen crossing budget. `flipdet.py` screens crossing subsets with the matrix determinant lemma modulo a prime. A matching determinant is only a filter: it neither identifies a knot nor detects the unknot. Search exhaustion within a finite move budget gives no lower bound on the unknotting number.

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
| `minimize_witness.py`, `reduce_witnesses.py`, `draw_witness.py` | Simplify and draw marked diagrams |
| `xtait/` | Embedded graph, isotopy, reduction, and flype kernels |
| `bj_scan.py` | Examine crossing neighbors of prime alternating knots for the conditional Bernhard–Jablan analysis in Section 5.2 |
| `make_targets.py`, `build_all_codes.py` | Build exploratory search data from an installed KnotInfo package |

The search-data builders read the installed database; they do not regenerate the frozen comparison used in the paper. Named endpoints are compared with KnotInfo PD codes rather than accepted solely from a census label or a Jones polynomial. Historical scripts and data remain archived alongside the paper constructions.
