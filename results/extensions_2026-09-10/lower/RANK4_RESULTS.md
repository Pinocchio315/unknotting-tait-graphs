# Three additional unknotting numbers

The completed calculations give

    u(13n4876) = u(13n4973) = u(13n5102) = 5.

All three knots were unresolved in the manuscript's results and have range
[4,5] in the KnotInfo responses retrieved on 10 September 2026. The responses,
timestamps, URLs, and checksums are retained in `rank4_live_retrieval.json`.

| Knot | Determinant | Signature | Enumerated surgery forms | New value |
|---|---:|---:|---:|---:|
| 13n4876 | 33 | -8 | 964 | 5 |
| 13n4973 | 105 | -8 | 3,849 | 5 |
| 13n5102 | 117 | -8 | 4,540 | 5 |

The form counts refer to the complete covering lists used in the obstruction,
not to pairwise non-isometric lattices. Every entry is accounted for in the
saved certificates, including any noncyclic discriminant groups.

Together with u(13n1619)=2, these are the four additional lower-bound
results integrated into manuscript v1.8. The full update also adds eleven
upper-bound constructions; see [the aggregate guide](../README.md).

## Correction terms and the lower bound

For each knot the archived reduced mod-2 Khovanov rank equals the determinant,
so the branched cover is an L-space. Greene candidates were freshly computed
from all 52 marked diagrams. Jones polynomials were independently computed
with Regina and checked against the archived polynomial.

For 13n4876 the Greene intersection leaves 27 conjugation-compatible vectors;
the spin identity leaves nine, and the Casson--Walker sum leaves one. Its
correction terms have sum 58/3. The other two knots already have unique
Greene vectors; both have sum 122/3, agreeing with the independent
Casson--Walker calculation. These identities use the normalization and
sources described in `README.md`.

After mirroring to signature +8, a sequence of four crossing changes must
have four negative changes. [Owens, Theorem 5](https://arxiv.org/abs/math/0506485)
then provides the positive definite rank-eight surgery forms tested here.
The theorem does not require the original knot to be alternating; correction
terms for these covers are supplied by the L-space/Greene calculation.

For each possible form Q and every group isomorphism phi, the required
conditions are

    m_Q(g) >= d(phi(g)),    m_Q(g) - d(phi(g)) in 2 Z.

The signature determines the orientation; the unique spin structure is the
origin of the odd-order discriminant group. The comparison includes every
unit of Z/D, with no linking-form prefilter.

The first implementation computes m_Q after an exact unimodular basis
change R=P^T Q P. The calculations for 13n4876 and 13n4973 finish by this
method and reject every form. The basis changes are verified by integer
matrix multiplication and det(P)=+/-1.

There is also a shorter certificate for all three knots. A characteristic
covector xi in class g supplies the upper bound

    m_Q(g) <= (xi^T Q^(-1) xi - 8)/4.

If this upper bound is below d(phi(g)), the isomorphism is impossible. A
wrong congruence modulo two also excludes it: changing characteristic
representatives in the same class changes this value by an even integer.
The saved proofs provide a specific xi for each isomorphism. In the reduced
bases, the vectors found here have coordinates in {-2,0,2}.

This direct verification uses exact matrix arithmetic and primitive linear
functionals from adj(Q) to identify the cyclic discriminant group. It does
not call the shortest-vector routine or its class-labelling implementation.
Every form is excluded. The longer m-function run for 13n5102 was therefore
stopped after the full covector proof was complete; its partial file is
explicitly marked `SUPERSEDED_INCOMPLETE`. The completed proof for this knot
is `13n_5102_rank4_direct.json.gz`.

## Matching upper bounds and validation

Regina independently simplifies the diagrams obtained by switching these
five zero-based crossings of the frozen PD to a one-component diagram with
zero crossings:

| Knot | Crossing indices |
|---|---|
| 13n4876 | 0, 1, 2, 5, 10 |
| 13n4973 | 0, 1, 2, 6, 9 |
| 13n5102 | 0, 1, 2, 3, 4 |

The verification uses actual Reidemeister simplification, not Jones-polynomial
matching as an unknot recognition procedure. Together with the lower bounds,
these witnesses establish the exact values five.

The final replay regenerates all three complete surgery-form lists, checks
every saved basis transformation, evaluates every covector proof, and repeats
the upper-bound simplifications. The known u=4 controls 12n474, 12n576, and
12n338 admit both the original comparison and the direct covector test. Exact
m-functions before and after basis changes also agree in additional test
forms of ranks one, two, four, and eight.

## Reproduction

To regenerate proofs, use a separate repository copy because the following
command writes certificates alongside the script:

```sh
python3 results/extensions_2026-09-10/lower/direct_rank4_extension.py 13n_4876 13n_4973 13n_5102
```

To verify the deposited bundle (including the independently generated
13n4876 certificate):

```sh
python3 results/extensions_2026-09-10/lower/replay_rank4_bundle.py
```

`rank4_summary.json` records the verified results, coverage counts, controls,
and hashes. `rank4_upper_bound_witnesses.json` records the original and changed
PDs. Full direct proofs for 13n4973 and 13n5102 are in their
`*_rank4_direct.json.gz` files. For 13n4876 the form list and reduction
matrices are in `13n_4876_rank4.json`, and its covector proofs are in
`13n_4876_rank4_covectors.json.gz`. No jobs from these calculations remain
running after the final replay.
