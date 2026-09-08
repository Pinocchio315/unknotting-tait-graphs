# Greene computation audit for v1.3

The reviewed additions are the cyclic rank-two through rank-four comparison
driven by Greene candidate correction terms, the boundary linking-form
pre-filter, and the enlarged sweep and its summary. The mathematical theorem
review is recorded separately in `theory_bj_notes.md`.

## Findings

No defect invalidating a deposited lower bound was found in this review. The
1,707 lower-bound records are unchanged. Input validation and failure handling
had weaknesses which were corrected before further use:

- `greene_sweep.validate_entry` now checks the actual reduced mod-2 Khovanov
  vector, rank/determinant equality, signature, range, selected test and control
  value before using the L-space argument. The worker previously trusted the
  builder to have checked these premises.
- `summarize_greene` now refuses to write bounds after conflicting completed
  verdicts or an obstructed control. It identifies controls from the frozen
  cohort, so omitting `control_ok` cannot hide a failed control. Malformed
  nonblank records are rejected rather than silently skipped.
- Resuming a sweep cannot bypass a previously obstructed control. A subprocess
  with nonzero exit status cannot publish a mathematical verdict merely
  because its stdout contains a JSON line. A lock serializes complete records
  written by workers in the same process.
- The rank test validates a complete exact rational correction-term vector and
  a generator pairing of full order. Its validation command now supplies the
  computed pairing, exercises the actual sweep pre-filter, and fails nonzero
  if a numerical comparison disagrees.
- The comments no longer assert `d(spin) = -sigma/4` for every knot whose
  branched cover is an L-space. Orientation comes from the signature and the
  same-sign crossing-change hypothesis. The spin-value equality is a further
  conservative implementation restriction: failure gives no bound. Likewise,
  the bounded search for a cyclic generator is no longer described as a
  generally exhaustive theorem; failure stops the computation.

These changes close paths through which future malformed inputs or failed
computations could be reported incorrectly. They do not alter any of the
archived successful results checked here.

## Full-cohort checks

`greene_code_audit.py` validated all 1,967 target and 1,452 control inputs
(3,419 total). Their raw mod-2 Khovanov vectors, determinants, signatures and PD
codes agree with the installed KnotInfo rows, and the tabulated double-cover
homology is cyclic. Every actual vector has mod-2 rank equal to the determinant.
This is a consistency check of the cited inputs, not an independent
recomputation of all Khovanov homology groups.

The deposited JSONL contains 3,418 distinct knot records. For every completed
record, the audit checks the correction-term domain, conjugation symmetry,
the number of candidate vectors, and agreement between the final verdict and
the candidate verdicts. Every recorded Ni--Wu candidate test has exactly
`2*phi(D)*D` tested affine identifications, with failure counters partitioning
that complete search.

The original derived summary had a stale description of unfinished targets.
The raw data contains `ERROR` for `13n4876` and `13n4973` and no record for
`13n5102`; the summary had instead counted one error and one timeout. The
summary was regenerated from the unchanged raw archive. All 1,707 bounds are
identical; only the unsuccessful-status bookkeeping changed. The control
counts remain 1,421 `PASS`, 9 `UNDECIDED` and 22 `TIMEOUT`. Thus “no obstructed
control among 1,430 completed records” must not be paraphrased as 1,430 passing
numerical tests.

## Bounded recomputations

All computations used one process at a time.

| Knot | Test | Fresh verdict |
|---|---|---|
| 13n111, 13n142, 13n196 | rank 2 | OBSTRUCTED for all three |
| 11n120 | unknotting number one | OBSTRUCTED |
| 11n180 | rank 3 | OBSTRUCTED |
| 8_20 | unknotting number one control | PASS |
| 9_43 | rank 2 control | PASS |
| 10_142 | rank 3 control | PASS |
| 12n474 | rank 4 control | PASS |

The fresh complete candidate sets and correction-term tables agree exactly
with the archived ones for all nine knots. Their combined replay and the
other audit checks took about 27 seconds in this environment.

The three new conflicts also received a separate rank-two comparison using
full rectangular characteristic-vector enumeration, adjugate class labels
and an independent enumeration of all cyclic isomorphisms, without using the
production `m_Q` or `form_admits` routines:

| Knot | Candidate forms | Excluded isomorphisms | Other rejection |
|---|---:|---:|---|
| 13n111 | 4 | 300 | one noncyclic discriminant group |
| 13n142 | 4 | 576 | none |
| 13n196 | 4 | 192 | none |

The full-box search contains a minimum in every characteristic class: changing
a vector by `2Qe_i` decreases its norm if its ith coordinate lies outside
`[-Qii,Qii]`; at the positive endpoint one may subtract `2Qe_i`, preserve the
norm and decrease `1^T Q^-1 xi`. A norm minimizer minimizing the latter linear
function therefore lies in the half-open box used by the audit. Each failed
isomorphism has an explicit violating class saved in `greene_code_audit.json`.

The new linking pre-filter agrees with exact matrix-derived primary linking
invariants on ten alternating Goeritz forms. The twelve published rank-two
and rank-three alternating controls were also rerun with the actual Greene
pairing supplied: all correction-term comparisons and verdicts agree. Their
output is in `alternating_validation.txt`.

The fresh replay of the new rank-four target `13n4361` was stopped by an
explicit 60-second subprocess timeout. Its archived run took 588.6 seconds;
the four newly obstructed rank-four targets historically took 588.6 to
3,709.4 seconds each. This review therefore does not claim to have repeated
those four full enumerations. The timeout supplies no mathematical result.
The rank-four positive control above completed in 10.9 seconds. See
`rank4_target_replay.json` for the bounded replay record.

## Reproduction

From the repository root, with its Python environment:

```sh
python results/review_v1_3/greene_code_audit.py
python -m unittest discover -s results/review_v1_3 -p test_greene_audit.py -v
python lower_bounds/greene/greene_ranks.py --validate 5_1 7_3 7_5 8_2 9_10 9_13 7_1 9_3 9_6 9_9 10_2 10_46
```

All 15 focused failure-handling regressions passed. The full server sweep was
not restarted, and no background process from this audit remains running.
