# Montesinos knots

These programs implement Section 3.7, **Montesinos knots**, of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*, manuscript v1.7. The frozen `v1.3` computational profile retains the 50 targets and deposited results used in this calculation.

A Montesinos knot is the numerator closure of a sum of rational tangles. Its double branched cover is a Seifert fibred space. The code normalizes its rational Seifert invariants and reverses orientation when necessary to obtain a negative-definite star-shaped plumbing. Each noncentral vertex has weight at most minus its valence, so at most the central vertex is bad. The graph hypotheses are checked before applying the plumbing formula.

The computation excludes unknotting number one for 49 targets. With the
independently known upper bounds, 44 have exact value two and five retain the
improved range `[2,3]`. Correction terms alone exclude 48 targets; reduced
Heegaard Floer homology supplies the remaining obstruction for `12n_457`.
The values for `12n_288` and `12n_501` are also used in the Bernhard–Jablan
argument in Section 5.1.

## Files and usage

| File | Purpose |
| --- | --- |
| `montesinos_u1.py` | Parse rational tangles, construct and check a plumbing, compute correction terms, and test the half-integral surgery constraints |
| `hf_red.py` | Compute reduced Floer ranks from the graded root and test the reduced-homology mapping cone |

From the repository root, after installing the dependencies described in the
[main README](../../README.md), run named examples without replacing the deposit:

```sh
python lower_bounds/montesinos/montesinos_u1.py 3_1 4_1 12n_288 12n_501 --out /tmp/montesinos_examples.json
python lower_bounds/montesinos/montesinos_u1.py 11n_102 12n_457 --out /tmp/montesinos_reduced_examples.json
```

The JSON report separates named results, target results (`apply`), and controls
(`validate`). Successful correction-term computations record the plumbing rank,
the number of classes, and a verdict; reduced ranks are included when computed.
`PASS` means the obstruction leaves unknotting number one possible.
`OBSTRUCTED` excludes it by correction terms or the preliminary cyclicity test,
and `OBSTRUCTED_HF` excludes it using reduced homology.

To rerun the target list or the larger control calculation, use:

```sh
python lower_bounds/montesinos/montesinos_u1.py --apply --workers 3 --out /tmp/montesinos_targets.json
python lower_bounds/montesinos/montesinos_u1.py --validate --workers 3 --out /tmp/montesinos_validation.json
```

Both batch modes require `--out`. Target ranges come from the frozen
`results/reference_table_2026-09-07.json`; Montesinos notation and diagrams
come from the installed `database_knotinfo` package. Keep the pinned package
when reproducing the paper. New runs do not update the aggregate manuscript
tables automatically.

## Correction terms and labels

For the negative plumbing matrix `Q_neg`, the code computes `d(Y) = -m_{-Q_neg}`, where `Y` is the oriented boundary of that plumbing. This sign is tied to the chosen plumbing orientation. Both orientations are examined in the surgery test. All quadratic values and correction terms use integers or rational numbers.

`ClassMap` labels characteristic vectors by `[k]` in the discriminant group, as in the Background's correction-term definition; it does not label them by `(k − diag(Q))/2`. For the odd determinants of knots, these labels give a bijection with the boundary Spinᶜ structures and put the spin structure at zero. A primary decomposition with two factors involving the same prime is not cyclic. The generator routine tests its actual element order and enumerates the whole group before the affine matching begins.

If `u(K) = 1`, the Montesinos trick realizes the cover as surgery with coefficient `±D/2`. The Ni–Wu formula imposes a nonnegative integral sequence `V_j` whose successive differences are zero or one. The code tests affine identifications and both orientations. Failure rules out unknotting number one; a match supplies no unknotting diagram. If an observed sequence has a nonzero final term, its unknown tail is not set to zero for the reduced-homology test.

## Reduced homology and termination

For a characteristic `k`, put `A = −Q_neg` and

```
chi_k(x) = (xᵀ A x − k·x)/2,   x ∈ Z^n.
```

The zero-dimensional lattice sublevel sets determine the graded root. For these plumbings, the sum over levels of `(number of connected components − 1)` is the dimension of reduced Heegaard Floer homology over `F₂`. Reversing orientation preserves that dimension in each conjugate class.

The sublevel enumeration uses a rational `LDLᵀ` decomposition. Integer bounds may enlarge the candidate interval, but the final rational inequality is exact. This avoids losing a boundary lattice point through floating-point rounding.

A weak local minimum under unit steps satisfies `|(2Ax − k)_i| ≤ A_ii`. The code enumerates this finite set exactly. Equal-weight components of weak minima that have an exit to a non-minimum can descend and create no new component. The remaining plateaux contain every possible birth of a sublevel component. Once a connected sublevel set contains all these births, all higher sublevel sets stay connected: every higher point has a finite path of nonincreasing unit steps to a lower level. Positive definiteness guarantees termination of descent. The rank is returned only with this certificate. A point, depth or representative limit returns an uncomputed result and cannot cause a Floer obstruction.

The half-integral mapping cone requires `rank HF_red(i) − T_i = c_i + c_(i−1)` with nonnegative integer `c_i`. Here `T_i` is the finite contribution of the tower kernel after its infinite `U`-tower is removed. For odd `D`, the cyclic system has a unique rational solution, so integrality and nonnegativity can be checked exactly.

## Independent controls and references

The fresh regression suite checks lens spaces against their recursive rational correction-term formula; twist knots must pass and their lens-space covers have zero reduced homology. Lidman's `11n_102` is a published obstruction control: its reduced ranks are `0, 0, 2`, while its correction terms alone allow the half-integral surgery pattern. The revised computation reproduces this distinction.

For `12n_457`, a separate exact calculation gives correction terms with numerators `11, 7, −5, 19, −9, −1, −1, −9, 19, −5, 7` over 22, up to permutation. Its reduced ranks are nine zeros and two ones. The revised algorithm reproduces these values and the reduced-homology obstruction. Replaying the 50 target knots gives the same 49 obstructions, of which only `12n_457` requires reduced homology.

Primary mathematical sources:

- [Manolescu–Owens, Section 4.1](https://www.maths.gla.ac.uk/~bowens/papers/conc.pdf): rational tangles and Montesinos plumbing conventions.
- [Ozsváth–Szabó, Theorem 1.2](https://arxiv.org/pdf/math/0203265): Floer homology of a negative-definite plumbing with at most one bad vertex.
- [Némethi, Section 11](https://arxiv.org/pdf/math/0310083): graded roots and the lattice description.
- [Ni–Wu, Proposition 1.6 and Remark 2.10](https://arxiv.org/pdf/1009.4720): the rational surgery correction terms.
- [Gainullin, Corollary 14 and Proposition 15](https://arxiv.org/pdf/1411.1275): tower and reduced parts of the mapping cone.
- [Lidman](https://arxiv.org/abs/2606.12431): the independent `11n_102` control.
