# lower_bounds/greene — correction terms of the double branched cover from a knot diagram

These scripts decide the question "is `u(K) = 1`?" for a Khovanov-thin knot `K` whose double branched
cover is hyperbolic, where none of the lower-bound methods of the paper applies.  They were written for
`12n491`, the last undecided knot in the Bernhard–Jablan counterexample of Brittenham and Hermiller
(arXiv:1705.05985): `u(13n3370) ≤ 2`, and by their Theorem 1.3 the conjecture fails for `13n3370` if
`u(12n288) = u(12n491) = u(12n501) = 2`, and for one of these three knots otherwise.  The Montesinos
obstruction of §3.7 gives `u(12n288) = u(12n501) = 2`; the scripts here give `u(12n491) = 2`.

| script | what it does |
|---|---|
| `qa_search.py` | quasi-alternating certificate by resolving crossings (`det(L) = det(L_0) + det(L_1)`), with spherogram simplification of every resolved link; a certificate proves that `Σ₂(K)` is an L-space |
| `greene_dinv.py` | Greene's spanning-tree model (arXiv:0805.1381): Kauffman states of a marked diagram, their absolute gradings (Theorem 4.1), Spin^c structures (§4.5) and the solitary states that generate the `E_1` page (Theorem 6.8); the surgery test of §3.7 (`u1_admissible`) on the resulting correction terms |
| `pin_dinv.py` | intersects the solitary-state pages of all markings and colourings of the diagram: the correction term `d(t)` is one of the gradings of the solitary states in the class `t` for every marking, the classes of two markings are matched by the unique automorphism of `H² = Z/D` that preserves the self-linking of `c₁`, and `d(t) = d(-t)` |

Validation: on alternating knots the gradings reproduce the Ozsváth–Szabó formula exactly (every
marking, both colourings); on the Montesinos knots 8_20, 9_43, 9_44 the pinned correction terms
coincide with those of the star plumbing (`lower_bounds/montesinos`), and the surgery test passes for
the knots with `u = 1` (3_1, 4_1, 5_2, 8_20, 9_44) and obstructs 5_1, 7_4, 9_43 (`u = 2`).

Result for `12n491` (`results/bernhard_jablan/`): reduced Khovanov homology thin of rank 69 = det,
quasi-alternating (certificate deposited), all 69 correction terms pinned from 47 of the 48 markings
(the same multiset from the independent census diagram `K12n491`), `d(spin) = 0 = -σ/4`, and the
surgery test is obstructed for both orientations and both signs of the linking-form pairing.  Hence
`u(12n491) = 2`, and `13n3370` is a counterexample to the Bernhard–Jablan conjecture.

    python qa_search.py 12n_491
    python greene_dinv.py 3_1 5_1 8_20
    python pin_dinv.py 8_20 9_43 9_44 12n_491
