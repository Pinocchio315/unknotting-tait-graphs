# Comparison of the symmetric matching obstruction and full surgery conditions

This directory reproduces the comparison reported in Section 3.9 of manuscript v1.8.

**No knot in the tested collection is excluded by the full correction-term
conditions while passing the positive even symmetric matching conditions.**
After removing 59 overlaps, the comparison covers **2,776 distinct knots**.
This is a comparison of obstruction strength, not a new count of unknotting
numbers determined by the paper.

## Conditions compared

The comparison follows Owens--Strle, *Immersed disks, slicing numbers and
concordance unknotting numbers*, [arXiv:1311.6702, Theorem 3](https://arxiv.org/abs/1311.6702):

- **P**: nonnegative normalized correction terms;
- **E**: even integral normalized correction terms;
- **S**: the paired symmetry in Theorem 3(iii);
- **M**: monotonicity in Theorem 3(iv);
- **B**: the successive-increase bound in Theorem 3(v).

The primary comparison is **PES versus PESMB**. The independent comparison
**PES versus PEM** addresses the remark immediately following Theorem 3
about symmetry and monotonicity. Both comparisons give identical knot-level
outcomes in this collection.

The term "older symmetric matching obstruction" means the first three
conditions as identified by Owens--Strle. It does not mean that
Ozsvath--Szabo's paper contains no stronger result. Their
[*Knots with unknotting number one and Heegaard Floer homology*,
Theorem 8.4](https://arxiv.org/abs/math/0401426) already contains stronger
restrictions arising from L-space surgery. Nor do we directly apply their
alternating-knot theorem to every non-alternating knot: the comparison uses
the general half-integral-surgery necessary conditions in Owens--Strle.
No claim of an improvement over Theorem 8.4 is made.

## Results

Counts in the obstruction columns mean that **no complete candidate vector
admits any matching** satisfying the indicated conditions. A passing test is
inconclusive; it does not imply unknotting number one.

| Collection | Knots tested | PE excludes | PES excludes | PEM excludes | PESMB excludes | Extra exclusions from M/B beyond PES |
|---|---:|---:|---:|---:|---:|---:|
| Greene targets | 1,523 | 104 | 1,523 | 1,523 | 1,523 | 0 |
| Greene known-u=1 controls | 574 | 0 | 0 | 0 | 0 | 0 |
| Montesinos targets | 50 | 2 | 48 | 48 | 48 | 0 |
| Montesinos cyclic controls, known u=1 or u=2 | 688 | 314 | 485 | 485 | 485 | 0 |
| **Distinct knots, removing overlaps** | **2,776** | **420** | **2,056** | **2,056** | **2,056** | **0** |

All 59 overlaps are controls, and both constructions give the same outcome
for every compared profile. All known-u=1 controls pass the full test:
574 in the Greene cohort and 201 in the Montesinos cohort, or 716 distinct
knots after removing the overlap.

The Montesinos source contains 779 records. Of these, 41 have noncyclic
first homology, so they do not enter this comparison of cyclic correction-term
matchings. Such covers already fail the cyclicity required by half-integral
surgery. The remaining 738 include all 50 target records. Determinant-one
examples are handled as the trivial group, with one identification and no
successive-step conditions. The stronger reduced-Floer mapping-cone test in
the original Montesinos calculation is **not** part of this comparison.

The Greene input comprises the 2,095 completed u=1 records in the archived
sweep, with the v1.8 exact correction terms for 13n1619 replacing its older
6,561-vector candidate set. The separately deposited computations for 12n491
and 13n3370 add two distinct target knots. There are 1,984 exact-vector
records and 113 records with ambiguity, comprising 4,595 complete candidate
vectors in total. No expensive spanning-tree calculation was repeated here.
The Montesinos correction terms were freshly calculated from the archived
notation using the exact definite-plumbing enumeration and subsequently cached.

## Five candidate-level differences are not knot examples

For the following known-u=1 controls, an additional candidate correction-term
vector fails M or B after passing PES:

| Knot | Candidate vectors | Vectors passing PES | Vectors passing PESMB |
|---|---:|---:|---:|
| 12n275 | 9 | 2 | 1 |
| 12n340 | 9 | 2 | 1 |
| 13n1465 | 9 | 2 | 1 |
| 13n1598 | 729 | 2 | 1 |
| 13n2822 | 9 | 2 | 1 |

Each still has a vector satisfying all conditions. These records show that
the later conditions can remove ambiguous candidates; they do not exhibit a
knot whose unknotting number is obstructed only by the stronger test. The
discarded candidate vectors have not been identified as the actual
correction-term vectors. Their coordinates and failing steps are retained in
the comparison output for inspection.

## Arithmetic and verification

`compare.py` implements two separately indexed calculations:

1. The literal spin-centered formula of Owens--Strle Theorem 3, testing every
   unit modulo D and both orientations. It separates P/E, S, M, and B.
2. A replay using the recursive lens-space correction terms in Ni--Wu surgery
   coordinates, testing every affine bijection, including all translations.
   It separates nonnegative even gaps, equality of entries assigned the same
   V-index, and successive differences of the resulting V sequence.

Neither routine uses a linking-form filter. The second routine enlarges the
identification space further by not fixing the spin origin. For all compared
profiles, the two routines have identical knot-level outcomes. A negative
result therefore survives this conservative enlargement.

All computations use exact integers after multiplying d-invariants by 4D;
nonnegative even gaps are detected by nonnegativity and divisibility by 8D.
Conjugation symmetry, complete class coverage, and integral scaling are
asserted. Every Greene candidate's full-test verdict agrees with its deposited
Fraction-based result, including the current 13n1619 certificate. All cyclic
Montesinos full-test outcomes agree with the archived correction-term-only
verdicts, which also used a linking-form filter.

`verify.py` additionally compares the recursive and production lens formulas
at 2,499 entries, tests 49 lens-space vectors and the trivial-group case,
and replays 24 selected real candidate vectors with the original rational
surgery routine, including all five candidate-level differences. Synthetic
vectors separately exercise symmetry, monotonicity, and boundedness failures;
these are arithmetic tests, not proposed knot examples. These checks verify
the comparison code, not a fresh proof of every underlying Greene d-invariant.

The computations run sequentially. The final Greene comparison takes about
13 seconds on this laptop; the newly evaluated Montesinos lattice calculations
take about 24 seconds in aggregate. Preparation, source comparison, and
verification take additional time. Reusing saved Greene candidates and exact
integer scaling makes this substantially faster than recomputing correction
terms or using the earlier rational-only runtime estimate.

## Files and reproduction

- `compare.py`: the two obstruction comparisons and cohort readers.
- `verify.py`: arithmetic tests and production-code cross-checks.
- `report.py`: input-hash verification, deduplication, and CSV export.
- `greene_comparison.jsonl`, `montesinos_comparison.jsonl`: per-knot and
  per-candidate outcomes, including the candidate-level extra rejections.
- `montesinos_correction_terms.jsonl`: the newly computed exact vectors,
  plumbing matrices, notation, generator pairings, and calculation times.
- `greene_summary.json`, `montesinos_summary.json`, `summary.json`: counts
  and source hashes. Cohort timings reflect the most recent execution, which
  may reuse the Montesinos cache.
- `verification.json`: completed arithmetic and regression checks.
- `knots.csv`: one row per distinct knot with a completed cyclic comparison.

Python 3 suffices for the Greene comparison, verification, and reporting.
The Montesinos calculation additionally requires NumPy and SymPy,
along with the existing local repository modules. Run from the repository
root, using a Python environment with those dependencies:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python lower_bounds/obstruction_comparison/compare.py greene
python lower_bounds/obstruction_comparison/compare.py montesinos
python lower_bounds/obstruction_comparison/verify.py
python lower_bounds/obstruction_comparison/report.py
```

All generated comparison files and the optional Montesinos cache are written
to `generated/obstruction_comparison/`. The frozen per-knot CSV and summary
are in `results/obstruction_comparison/`. They are numerical evidence for the
paper; candidate-level logs are reproduced locally. Correction terms are
computed from the authenticated compact KnotInfo input, not a live download.
