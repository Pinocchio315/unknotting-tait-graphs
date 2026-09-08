# Correction-term obstructions for alternating knots

This directory accompanies **Computation of Unknotting Numbers: Which Knot Breaks the
Bernhard–Jablan Conjecture?**. It implements the definite surgery obstruction reviewed in
Section 2.4 and applied in Section 3. The signed Jones-polynomial test
also uses Section 2.6. These programs compute necessary conditions for an
unknotting sequence; finding an admissible surgery form does not construct one.

The principal source is Brendan Owens, [*Unknotting information from Heegaard
Floer homology*](https://arxiv.org/abs/math/0506485), *Advances in Mathematics*
217 (2008), 2353–2376, Theorems 1 and 5. The cases of three and four crossing
changes apply the existing general theorem, rather than a new extension of it.

## Mathematical conventions

Choose the representative of the knot with positive signature. A **negative
crossing change** changes a negative crossing to a positive crossing. Owens's
definite obstruction assumes that the number of negative changes is
`n = sigma(K)/2`.

- For a sequence of two changes with `|sigma| = 4`, the signature forces the
  pattern `(positive, negative) = (0, 2)`.
- For `|sigma| = 2`, the patterns `(1, 1)` and `(0, 2)` are possible. The definite
  obstruction tests only the first. `traczyk.py` must exclude the second before
  their combination can prove `u >= 3`. The program mirrors the Jones polynomial
  when the signature sign is reversed.
- For `|sigma| = 6` and a sequence of three changes, or `|sigma| = 8` and four
  changes, every change is negative. Excluding all admissible surgery forms
  proves `u >= 4` or `u >= 5`, respectively. An independently supplied upper
  bound is still needed to obtain an exact unknotting number.

For a proposed sequence of `r` changes, the integral surgery form has `r` linked
chains `(m_i, 2)`. Its Schur complement is half the rank-`r` form `Q' = 2M-I`,
with off-diagonal entries `2a_ij`.
Consequently `det(Qt) = det(Q')`, and positivity of `Q'` implies positivity of
`Qt`. The required negative changes correspond to even `m_i`, equivalently
diagonal entries of `Q'` congruent to 3 modulo 4.

For a positive definite integral matrix `A` of odd determinant, the core computes

```
m_A([xi]) = min (xi^T A^(-1) xi - rank(A))/4,
           xi_i = A_ii mod 2, [xi] in Z^r / A Z^r.
```

The class label is **`[xi]`**, without translating by the diagonal of `A`.
Because the discriminant group has odd order, this labels characteristic classes
bijectively, and class zero is the unique spin structure. Both forms in a
comparison must use this same convention. The two positive checkerboard
Laplacians describe opposite cover orientations. The spin value
`d(spin) = -|sigma|/4` selects the appropriate one. For alternating knots the
chosen Goeritz form is sharp, so `m_G` gives the cover's correction terms.

A candidate survives only if some group isomorphism `phi` satisfies, in every
class,

```
m_Qt(g) - m_G(phi(g)) is a nonnegative even integer.
```

The values are exact rational numbers; equality modulo two is not a comparison
of rounded or floating-point values.

## Algorithms and the meaning of completeness

`owens_obstruction.py` provides integral determinants, discriminant-group
coordinates, characteristic-vector minima, and exhaustive finite-group
isomorphism comparison. Integral diagonalization is followed by decomposition
into prime-power factors; the intermediate diagonal need not already obey the
usual Smith divisibility convention.

Characteristic covectors are enumerated in an ellipsoid using an exact rational
`LDL^T` decomposition. Integer square roots give rigorous coordinate bounds.
Once an entire ellipsoid contains a representative of every class, the recorded
minima are global: every omitted covector has larger norm. Exhausting the growth
limit raises `IncompleteCorrectionTerms`, and incomplete input tables cannot
produce a completed obstruction. This replaces the previous floating-point
Cholesky pruning, whose tolerances did not justify an exhaustive mathematical
comparison.

`owens_u3.py` and `owens_u4.py` enumerate a finite superset of Minkowski-reduced
forms, using the diagonal-product bounds `2 det(K)` in rank three and `4 det(K)`
in rank four. Reduction under the full integral general linear group is not
enough for the required meridians. The programs therefore lift every element of
`GL(3,F_2)` or `GL(4,F_2)`—168 and 20,160 elements—to recover the admissible bases
modulo `Gamma(2)`. This subgroup preserves diagonal residues modulo four.
The diagonal-product bounds follow from the classical Minkowski-reduction
bounds; see Barnes and Trenerry, [*The minimum determinant of Minkowski-reduced
quinary quadratic forms*](https://doi.org/10.1017/S1446788700024964),
Eqs. (1.2)–(1.4), where the rank-three and rank-four constants are 2 and 4.
Signed permutations identify equivalent surgery components. Positivity checks
use exact principal minors, rather than numerical eigenvalue thresholds.

`linkform.py` supplies a necessary preliminary filter: the primary cyclic factors
of the discriminant group must agree, as must the square class of each cyclic
odd-primary linking pairing. A deterministic scan of the coordinate generators
finds a generator for each such primary part. Noncyclic primary pairings are not
classified by this filter; omitting them only weakens it.

The enumeration routines are designed for the small determinants in the knot
tables. They retain NumPy integer arrays for the vectorized form enumeration;
this is not a general-purpose arbitrary-size lattice package. The rational
correction-term kernel itself uses Python integers and `Fraction` arithmetic.

## Verdicts

| Verdict | Interpretation |
| --- | --- |
| `OBSTRUCTED` | Every admissible candidate has failed a completed comparison; the proposed sequence is excluded. |
| `PASS` | At least one candidate or an untreated sign pattern remains possible; this does not prove the proposed upper bound. |
| `NOT_APPLICABLE` | The signature or supplied data do not meet this entry point's hypotheses. |
| `UNDECIDED` | An isomorphism-search limit or incomplete comparison prevents an obstruction. |
| `ERROR` | Input or computation failed; this supplies no lower bound. |

An external timeout or interrupted process supplies no verdict. A partial list
of processed knots is different from a partially processed knot: a stored
`OBSTRUCTED` row requires the latter's candidate search to be complete.

## Running and verification

Install the dependencies from the repository root, then run the small controls:

```bash
python -m pip install -r requirements.txt
python lower_bounds/owens/validate.py
python lower_bounds/owens/crosscheck.py 9_10 8_5 11a_63 13a_16
python -m unittest discover -s tests -p test_owens.py -v
```

`validate.py` checks six obstructed examples of Owens and eight known knots with
unknotting number two. It exits unsuccessfully if a verdict differs.
`crosscheck.py` independently builds checkerboard forms from diagram faces,
enumerates characteristic covectors in a full finite box, labels classes using
adjugates, and searches generator images. It does not share the primary core's
linear algebra or ellipsoid algorithm.

The unit tests compare the two correction-term algorithms on small forms, test
rational ellipsoid boundaries and incomplete computations, verify integral
basis invariance and discriminant-group structure, enumerate all mod-two lifts,
and compare the small reduced-form lists with a separate brute-force search.
They also test the number of even surgery coefficients and reject a zero
Jones-polynomial evaluation without hanging. The six added rank-four targets
have independent checks of both Goeritz correction-term tables; these input
checks are distinct from rerunning their complete surgery-form exclusions.

Individual higher-rank computations can be run as follows:

```bash
python lower_bounds/owens/owens_u3.py 7_1 12a_107
python lower_bounds/owens/owens_u4.py 9_1
```

For sweeps, always choose a fresh output location when verifying the archived
paper. The deposited JSON/JSONL files are retained unchanged:

```bash
python lower_bounds/owens/run_sweep.py --out /tmp/owens_rank2.jsonl --prefix 13a --workers 6
python lower_bounds/owens/run_sweep.py --out /tmp/owens_signed.jsonl --sigma 2 --workers 6
python lower_bounds/owens/run_u3_sweep.py apply /tmp/apply_u3.jsonl
python lower_bounds/owens/run_u4_local.py --targets-only --out /tmp/owens_rank4.jsonl --workers 3
```

The sweep drivers select targets using the installed `database_knotinfo` version;
a changed database can change that selection. The paper's counts instead refer
to its archived input ranges and deposited results. Rank-four targets and known
comparison knots may remain unfinished; report only completed rows when updating
the manuscript. A sweep may take hours, especially in rank four.

## Files

- `owens_obstruction.py`: shared exact correction terms and obstruction for two changes.
- `owens_u3.py`, `owens_u4.py`: surgery-form enumeration and obstruction for three/four changes.
- `linkform.py`, `traczyk.py`: linking-pairing filter and signed Jones criterion.
- `kinfo.py`, `_compat.py`, `taitgraph.py`: input access and compatibility with the shared diagram code.
- `validate.py`, `crosscheck.py`: known examples and an independent mathematical implementation.
- `run_sweep.py`, `run_u3_sweep.py`, `run_u4_local.py`: resumable batch drivers.
- `make_targets.py`, `owens_label_scan.py`, `slurm_*.sh`: auxiliary target and cluster commands.

Deposited results also appear in `results/owens_rank2`, `results/owens_rank3`, and
`results/owens_rank4`, relative to the repository root. Historical filenames
identify those records; running these programs does not update the manuscript.
