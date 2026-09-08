# Montesinos knots

These programs implement Section 3, **Montesinos knots**, of *Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?*.

A Montesinos knot is the numerator closure of a sum of rational tangles. Its double branched cover is a Seifert fibred space. The code normalizes its rational Seifert invariants and reverses orientation when necessary to obtain a negative-definite star-shaped plumbing. Each noncentral vertex has weight at most minus its valence, so at most the central vertex is bad. The graph hypotheses are checked before applying the plumbing formula.

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
