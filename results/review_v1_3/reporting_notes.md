# v1.3 reporting and provenance audit

The audit compares the unchanged v1.2 manuscript inputs with both the first
v1.3 submission (commit `6468494`) and the subsequent candidate-scan update
that the author explicitly included during review. The obstruction counts are
the same in both v1.3 states; the conditional alternating-knot list changes.

## Reproduced totals

| Quantity | v1.2 | Latest v1.3 | Change |
|---|---:|---:|---:|
| Exact unknotting numbers | 1,227 | 2,719 | +1,492 |
| Exact values equal to 2 | 840 | 2,149 | +1,309 |
| Exact values equal to 3 | 237 | 373 | +136 |
| Exact values equal to 4 | 129 | 172 | +43 |
| Exact values equal to 5 | 21 | 25 | +4 |
| Nonexact improved ranges | 176 | 390 | +214 |
| Conditional alternating-knot list | 806 | 959 | +153 |
| Incompatible September intervals | 7 | 10 | +3 |

The sweep supplies 1,707 new lower bounds: 1,492 make an interval exact and
215 improve a still-open interval. The net interval increase is 214 because
`13n_1587`, previously counted as an upper-bound improvement, is now exact.
The Greene method row contains 1,494 exact values: the new 1,492 plus the two
earlier computations for `12n_491` and `13n_3370`. These categories are not
interchangeable in the abstract or the sweep description.

Every frozen sweep entry uses the **improved v1.2 range**, not the original
comparison range. All 3,419 entry ranges agree with an independent reconstruction
of that v1.2 table. All 3,418 deposited job records have distinct names and valid
target/control roles. The new bounds coincide exactly with the 1,707 obstructed
target records; no control contributes a bound.

Of 1,967 targets, 1,964 have completed records (1,707 OBSTRUCTED, 256 PASS,
one UNDECIDED), two have ERROR records with return code -15, and one has no
record. The latter three are `13n_4876`, `13n_4973`, and `13n_5102`.
Of 1,452 controls, 1,421 PASS, nine are UNDECIDED, and 22 TIMEOUT. Six control
values were themselves obtained by the earlier v1.2 computations; the controls
should not all be described as independent published values.

The new conflicts are `13n_111`, `13n_142`, and `13n_196`: each has paper value
three and September value two. The full set has nine incompatible exact values
and the interval `[3,4]` for `13a_650`. Among the 2,719 exact values, 194 agree
with September exact entries, nine disagree, and 2,516 remain unresolved there.

## Candidate analysis

The latest population partitions as **959 certified + 46 heuristic + 22 with
identified candidate neighbours = 1,027**. The 46 heuristic cases consist of
39 in the zero-candidate list and seven moved out of the candidate list after
Jones-only composite identifications. Those seven are `13a_1154`, `13a_1158`,
`13a_1184`, `13a_120`, `13a_2098`, `13a_320`, and `13a_473`. They are not added
to the theorem list.

Of the historical 112 candidate child knots, the new bounds exclude unknotting
number one for 101 (82 now have exact value two; 19 have range `[2,3]`). Eleven
remain `[1,2]`. The detailed validity of the new crossing identifications is
covered by the separate updated-dichotomy audit, not by this counting script.

## Corrected reporting defects

- The verifier authenticated only `paper_v1_2` inputs, omitted uppercase `G`
  entries from the interval parser, and still expected exactly two Greene-marked
  knots. It now authenticates the selected version and verifies every `G` mark
  in Appendices A–E.
- Sweep ingestion checked only summary fields. It now validates the actual
  tabulated mod-2 vector, complete conjugation-symmetric candidate class sets,
  the number and outcomes of all recorded candidate tests, duplicate jobs,
  roles, signatures, tested lengths, and the absence of obstructed controls.
  This authenticates record consistency; it does not recalculate every complex.
- The expected `dichotomy` value in the manifest was silently ignored. It is
  now checked, and unknown expected-count keys are rejected.
- The release-audit interpretation text still hard-coded the v1.2 counts even
  when its numerical fields contained v1.3 counts. It now derives that prose
  from the same reconstruction.
- The stored sweep summary was stale about ERROR/TIMEOUT records. The Greene
  code audit regenerated the derived summary without changing the raw archive;
  the v1.3 manifest now authenticates that derived file too.
- Version metadata now says 1.3. Explicit `--manuscript` selection retains
  reproducibility of v1.2. All ten original v1.2 generated files remain
  byte-for-byte identical to their independently regenerated counterparts.

Concurrent work replaced the live candidate scan during review. To preserve the
historical paper, `open23/priority_u23_paper_v1_3.json` contains the exact bytes
of `open23/priority_u23_final.json` from commit `6468494`, with the same SHA-256.
The historical manifest points to that frozen alias. The optional
`v1.3-review` profile preserves the initial v1.3 audit state; the final paper
uses the current `v1.3` profile and its updated scan.

## Reproduction

From `code/`, with the project Python environment:

```text
python results/review_v1_3/reporting_audit.py
python -m unittest discover -s tests -p test_manuscript_inputs.py -v
python audit_knotinfo_releases.py --verify
python verify_manuscript.py --tex ../main_v1.2.tex
python verify_manuscript.py --tex ../main_v1.3.tex
```

The focused reporting suite passed 13 tests, including corruption of a candidate
set, an omitted comparison, a hidden surviving fit, a repeated job, mixed
manuscript versions, and attempted promotion of a Jones-only case into the
certified category. `reporting_audit.json` records all 1,707 transitions and
the source hashes. Independent mathematical computation checks are documented
by the other audit files in this directory.
