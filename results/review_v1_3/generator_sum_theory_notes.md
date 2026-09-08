# Review of the late connected-sum addition

The incoming paragraph correctly states that the double branched cover of a
connected sum is the connected sum of the covers. It incorrectly infers
additivity of the minimum number of generators of their first homology.
This breaks the displayed proof and invalidates the advertised pair counts
as consequences of that argument. It does not by itself establish a
counterexample to the separate unknotting-number statement.

## Exact finite-group calculation

For a finite abelian group H, let r_p(H)=dim_Fp(H/pH). The minimum number of
generators is

    g(H) = max_p r_p(H).

Every generating set maps to a spanning set of H/pH, giving the lower bound.
Conversely, primary decomposition and the Chinese remainder theorem let
generators of distinct prime-primary parts be combined, giving the maximum
rather than the sum of their numbers. Consequently,

    g(H_A direct_sum H_B) = max_p (r_p(H_A) + r_p(H_B)).

There need not be one prime that realizes the maxima of both summands.

Two examples use the actual sharp knots in the proposed dataset:

- 3_1 and 4_1 each have u=g=1. Their double covers have first homology Z/3
  and Z/5; the direct sum is cyclic Z/15 and has g=1, not 2. Scharlemann's
  theorem separately implies that this connected sum has u=2.
- 8_18 and 9_40 each have u=g=2. Their first homology factors are [3,15]
  and [5,15]. The 3-primary ranks are 2 and 1, and the 5-primary ranks are
  1 and 2. Their direct sum therefore has g=3, not 4. This generator
  argument alone does not determine the connected sum's unknotting number.

## A valid replacement

Fix one prime p. If

    u(A) = dim_Fp H_1(Sigma_2(A); F_p),
    u(B) = dim_Fp H_1(Sigma_2(B); F_p),

then u(A#B)=u(A)+u(B). Indeed, an unknotting sequence of r changes presents
the first homology of the double cover with at most r generators, so
u(K)>=r_p(K). The direct-sum identity makes r_p additive. Applying this
lower bound to A#B gives u(A)+u(B), while concatenating the two unknotting
sequences supplies the reverse inequality.

This proof also works for mirrors, since mirroring changes orientation but
not these mod-p dimensions. Counts should explicitly refer to unordered
pairs of tabulated knot types with repetition; they should not be called
counts of all oriented or chiral connected-sum types without an additional
counting convention.

## Independent bounded enumeration

`audit_generator_sum_claim.py` checked the first-homology factors from the
installed database and the reviewed u-table, saving its actual inputs in
`generator_sum_claim_audit.json`. The incoming count of 1,972 individual
sharp knots is reproduced; 1,504 have u=1. Counting nontrivial invariant
factors agrees with max_p r_p on these individual inputs, so the error is
specifically the passage to pairs.

| Quantity | Count |
|---|---:|
| All unordered pairs with repetition from the 1,972 knots | 1,945,378 |
| Pairs with a prime simultaneously attaining both unknotting numbers | 483,817 |
| Such common-prime pairs with at least one u>=2 summand | 269,032 |
| Pairs with both u=1, covered separately by Scharlemann | 1,131,760 |
| Pairs with some u>=2 not covered by the common-prime argument | 544,586 |

The last row is **not** a count of knots with unknown unknotting number:
other invariants or constructions might determine some of them. It says
only that this particular generator argument does not certify additivity.
These independent audit counts are kept outside production reporting here;
the reporting subtask authenticates and freezes the public source inputs
before incorporating a corrected common-prime result.

The further incoming paragraph about 297 small pairs, 26 settled pairs,
271 open pairs and the claim that no obstruction in the paper decides
5_1#mirror(5_1) is not established by this computation. It requires its own
input list and exhaustive evaluations. It should not be retained as an
unsupported extension of the corrected proposition.

No manuscript or production code was edited in this review subtask.
