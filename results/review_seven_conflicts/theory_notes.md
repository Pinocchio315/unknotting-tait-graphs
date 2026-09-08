# Seven KnotInfo disagreements: theoretical audit

Scope: read-only examination of `main_v1.2.tex`, `lower_bounds/owens/owens_obstruction.py`, `owens_u3.py`, and the primary mathematical sources. This is an audit of the logical implications and completeness mechanisms; the separate numerical replays supply the per-knot outcomes.

## Conclusion

No mathematical defect was found in the signature-sharp Owens implications for the five rank-two knots or the two rank-three knots. In particular, square factors in determinants 245 and 225 do not require additional candidate determinants. The actual surgery trace, rather than an arbitrary bounding four-manifold, is being used.

## Primary-source anchors

- Brendan Owens, *Unknotting information from Heegaard Floer homology*, https://arxiv.org/pdf/math/0506485 : Theorem 1 (PDF pp.2–3), Theorems 3 and 5 (pp.3–4), Lemma 2.2 and following paragraph (pp.8–9), Theorem 4.1 (p.17), proof of Theorem 5 (p.19).
- E. S. Barnes and D. W. Trenerry, *The minimum determinant of Minkowski-reduced quinary quadratic forms*, J. Austral. Math. Soc. (Series A) 32 (1982), 405–411, https://doi.org/10.1017/S1446788700024964 : pp.405–406, Eqs.(1.2)–(1.4), lambda_3=2 and lambda_4=4.
- Ciprian Manolescu and Brendan Owens, *A concordance invariant from the Floer homology of double branched covers*, https://arxiv.org/pdf/math/0508065 : definition delta=2d and the alternating-knot formula delta=-sigma/2.

## Checks against the implementation

1. **Signs.** After mirroring to positive signature, signature 4 forces both changes to be negative if u=2; signature 6 forces all three to be negative if u=3. Thus the crucial hypothesis n=sigma/2 holds for all seven knots. The untreated mixed-sign issue for signature 2 is irrelevant here.

2. **Determinant and square factors.** Half-integral surgery is replaced by integral surgery on a link in S^3. Attaching only 2-handles to B^4 gives a simply connected trace. Its intersection form Qt presents the complete boundary H1, so det(Qt)=det(K), including when det(K) has square factors. Considering D/s^2 would be required in some more general filling obstructions, but not in this one. The full-D enumeration in `candidates` and `candidate_plumbings` is therefore appropriate.

3. **Rank-two list.** `candidates` implements the precise reduced domain of Owens Theorem 1: 0<=a<m1<=m2, (2m1-1)(2m2-1)-4a^2=D, with both m_i even. Its while bound is necessary because m2>=m1 and a<=m1-1, so no admissible pair is lost at termination.

4. **Rank-three reduced forms.** Every class has a Minkowski-reduced representative satisfying ordered diagonals a<=d<=f, pairwise bounds |2b|<=a, |2c|<=a, |2e|<=d, and adf<=2D. `reduced_ternary_forms` enumerates a superset of these representatives. Omitting other reduction inequalities adds forms and cannot create a false obstruction. The bound is documented in Barnes–Trenerry, p.406. Positivity follows from a>0, ad-b^2>0, and det=D>0.

5. **Mod-two classes.** Given an admissible Q=T^T F T and the enumerated lift P with P=T modulo 2, set R=P^{-1}T. Then R=I modulo 2 and Q=R^T(P^T F P)R. Thus all 168 residue lifts reach every Gamma(2) class. For Q=I modulo 2 and R=I+2B, each diagonal of R^TQR agrees with the original diagonal modulo 4. The m_i-even filter cannot discard a class that could become admissible through another lift. Signed permutations only relabel/orient the surgery components and preserve the obstruction.

6. **Correction-term orientation and class labels.** The two positive checkerboard forms give opposite cover orientations; for these nonzero signatures their spin correction terms have distinct signs. Selecting d(spin)=-|sigma|/4 fixes the orientation required by Owens. Because determinant is odd, labelling characteristic covectors by [xi] in coker(A) is a bijection on Spin^c classes, with spin at zero. Both source and target use the same convention. Testing all group isomorphisms already includes all possible identifications fixing spin; an extra translation is not necessary in these labels.

7. **Exact minima.** The LDL enumeration uses Fractions and integer square roots. Its coordinate endpoints implement the exact rational ellipsoid inequality. Once an exhaustively enumerated ellipsoid contains a representative of every class, no vector outside it can improve a recorded minimum. Missing classes raise `IncompleteCorrectionTerms`; they are not read as obstructions.

8. **All group isomorphisms.** The product of prime-power cyclic factors is canonical. Possible images of every generator satisfying its order relation are enumerated, followed by an injectivity check on the entire group. This reaches every isomorphism. The inequality and congruence are tested as a nonnegative even integer difference of exact rationals. Search limits return UNDECIDED. No arbitrary permutation of d-values is used as a substitute for group matching.

## Documentation corrections identified and subsequently applied

- At the start of the review, `owens_u3.py` lines 4–5 called the bounding trace a positive-definite *plumbing*. In general the original surgery components need not be unknots, so the theorem guarantees a positive-definite **2-handlebody**, not that the filling is a plumbing. This terminology does not affect the matrix obstruction. The function name `plumbing` can remain a legacy name, but its documentation should say integral surgery form.
- At the start of the review, `lower_bounds/owens/README.md` had the title typo `Conjecture?s` and stale section pointers. The reduction bound should receive the Barnes–Trenerry citation above.

The theoretical subreview was read-only. Its terminology, title typo, section pointer, and reduction-bound citation corrections were subsequently applied to the code documentation. The mathematical algorithms and manuscript values were not changed.
