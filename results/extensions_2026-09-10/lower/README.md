# Further lower bounds in manuscript v1.8

The completed calculations give **u(13n1619)=2** and
**u(13n4876)=u(13n4973)=u(13n5102)=5**. These are four additions to the
frozen v1.3 profile, now integrated into Section 3.8. The full v1.8 totals,
including the eleven upper constructions, are explained [one level up](../README.md).
The rank-four proof is described in [RANK4_RESULTS.md](RANK4_RESULTS.md).

## Mathematical argument

Let Y be the double branched cover of K, and put D = det(K). The following
argument requires Y to be an L-space. In these computations this is certified
by the archived reduced mod-2 Khovanov rank equality, using the
[Ozsvath--Szabo spectral sequence](https://arxiv.org/abs/math/0309170).

Write lambda_CW for the Casson--Walker invariant in the normalization
lambda_CW(boundary of the negative E8 plumbing) = -1. By
[Rustamov, Theorem 3.3 and the definition in Section 2](https://arxiv.org/abs/math/0409294),

    sum_s d(Y,s) = -2 D lambda_CW(Y).

Indeed the renormalized Euler characteristic is -d(Y,s)/2 when the reduced
Floer group vanishes. Mullins's Jones-polynomial formula gives

    lambda_CW(Sigma(K)) = sigma(K)/8 - V'_K(-1)/(12 V_K(-1)).

The formula appears in [Greene--Watson, Theorem 13](https://arxiv.org/abs/1106.5559),
where their invariant is twice lambda_CW. Thus the directly implemented
identity is

    sum_s d(Sigma(K),s) = -D sigma(K)/4 + D V'_K(-1)/(6 V_K(-1)).

All fractions and polynomial coefficients are evaluated exactly. This
identity is orientation sensitive; the PD convention is checked against the
archived Jones polynomial, and reversing orientation negates all correction
terms and lambda_CW together. For non-L-spaces the reduced Floer Euler
characteristic cannot be omitted, so this filtering procedure is not used.

For 13n1619, D = 33, sigma = 0, and the reduced mod-2 Khovanov rank is 33.
Its Jones polynomial, independently computed from its PD with Regina, is

    t^(-2) - 2t^(-1) + 4 - 5t + 5t^2 - 5t^3 + 5t^4 - 3t^5 + 2t^6 - t^7.

Consequently V(-1) = 33, V'(-1) = -80, lambda_CW = 20/99, and

    sum_s d(Y,s) = -40/3.

The freshly recomputed [Greene model](https://arxiv.org/abs/0805.1381) agrees
with the frozen calculation on all 52 marked diagrams. In cyclic coordinates
on Spin^c(Y), the generator has self-pairing 1/33 and the spin structure has
coordinate zero. After intersecting marked-diagram candidate sets, 15 classes
remain ambiguous. Conjugation leaves eight free choices, each with three
possible values: 3^8 = 6,561 vectors. Six of these passed the old surgery test.

There is a particularly simple way to eliminate the ambiguity. If C_k is the
candidate set for class k, then d(Y,k) >= min(C_k). **The sum of these 33
minima is already -40/3.** Since the actual correction terms have that sum,
each individual minimum must be attained. This proves uniqueness without
having to assume that any particular differential vanishes. It also avoids
needing the separate spin identity for this knot. The implementation checks
this argument as well as enumerating all candidates.

The resulting values for k = 0,...,16 are below; the remaining values follow
from d(k) = d(33-k).

| k | d(Y,k) | k | d(Y,k) |
|---:|---:|---:|---:|
| 0 | 0 | 9 | 4/11 |
| 1 | -8/33 | 10 | -8/33 |
| 2 | -32/33 | 11 | 2/3 |
| 3 | -2/11 | 12 | -10/11 |
| 4 | 4/33 | 13 | -32/33 |
| 5 | -2/33 | 14 | -50/33 |
| 6 | -8/11 | 15 | -6/11 |
| 7 | 4/33 | 16 | -2/33 |
| 8 | -50/33 | | |

If u(K)=1, one orientation of Y must be 33/2 surgery on a knot in S^3.
The [Ni--Wu correction-term formula, Proposition 1.6](https://arxiv.org/abs/1009.4720)
then requires all differences between the lens-space terms and the terms for
that orientation of Y to be nonnegative even integers, under the induced
affine identification of Spin^c structures. We test all units a, all
translations b, and both orientations, without a linking-form prefilter:

    i -> a i + b mod 33; 2 phi(33) 33 = 2 * 20 * 33 = 1,320 maps.

No map works. The production implementation records 1,316 failures of the
gap condition and four earlier failures of the repeated-V condition. The
independent replay checks all gaps before the repeated-V condition and finds
that **all 1,320 maps already fail the gap condition**. These results agree:
the four early failures also have an invalid gap at a later index. Thus no
monotonicity condition on V_j is needed for this particular obstruction.

Finally, switching zero-based crossings [0,8] of the frozen PD gives a diagram
that Regina simplifies by Reidemeister moves to a zero-crossing knot. The pair
[0,12] gives another verified witness. Actual simplification is required;
Jones polynomial 1 alone is not used to recognize the unknot. Together these
calculations give u >= 2 and u <= 2.

## Verification and reproduction

Run from the repository root with NumPy, SymPy and Regina installed:

```sh
python3 results/extensions_2026-09-10/lower/replay_half_integral.py
python3 results/extensions_2026-09-10/lower/replay_rank4_bundle.py
python3 results/extensions_2026-09-10/lower/verify_casson_walker.py --outdir /tmp/casson_walker_recalculation
```

The first script independently excludes all 1,320 affine maps for 13n1619
using a lens-space recursion. It needs only the standard library and writes
`independent_half_integral_replay.json` beside the script. The rank-four replay
regenerates all 9,353 forms, verifies the exact basis changes and saved
covectors, checks the five-crossing upper bounds, and repeats three known
u=4 controls. It updates only `rank4_summary.json`, an output excluded from
the input manifest.

`verify_casson_walker.py` recomputes Greene candidates and Jones polynomials
for 13n1619 and a comparison knot, repeats the five known u=1 controls, and
checks the sum identity on all 3,204 complete vectors in the frozen sweep.
Its default output directory is `generated/v1_8_extension_recalculation/`;
`--outdir` keeps reruns separate from the deposited certificates. The same
18-knot screen of ambiguous candidates is included; only 13n1619 acquires a
new obstruction in that screen. Khovanov ranks are checked from archived
mod-2 homology vectors, not recomputed from chain complexes.

| File | Role |
| --- | --- |
| `13n_1619_certificate.json` | PD, L-space evidence, candidate sets, exact correction terms, surgery test, upper bounds |
| `verify_casson_walker.py` | Exact sum constraint, candidate filtering, fresh Greene and Jones calculations |
| `replay_half_integral.py` | Independent affine-surgery arithmetic |
| `rank4_extension.py` | Exact unimodular basis reduction and m-function computations |
| `certify_rank4.py` | Characteristic-covector proof and independent arithmetic checks |
| `direct_rank4_extension.py` | Full form enumeration with direct covector certificates |
| `replay_rank4_bundle.py` | Replay all three deposited rank-four proofs |
| `13n_4876_rank4.json`, `13n_4876_rank4_covectors.json.gz` | Complete first proof, transformations, covectors |
| `13n_4973_rank4_direct.json.gz`, `13n_5102_rank4_direct.json.gz` | Complete direct proofs for the other two knots |
| `rank4_upper_bound_witnesses.json` | Five-crossing upper constructions |
| `rank4_positive_controls.json` | Three admitted known-u=4 examples |

The generation scripts `rank4_extension.py` and `direct_rank4_extension.py`
write new certificates beside themselves. For fresh computations, run them
in a separate copy of the repository; the replay commands above preserve the
manifested proof files. The input manifest authenticates the preserved proof
records. Recorded code hashes inside an original certificate describe its
creation environment, before the scripts were relocated here.

Raw cyclic PD rows used by the Tait-graph routines need not use Regina's
oriented labelling. For an auxiliary Regina check, convert with
`Link(pd).PD_code(min_strand_index=1)` first. For the frozen lower-bound input
PDs, the oriented convention is already checked against the archive.

The replay reads the authenticated compact cells in `results/knotinfo_calculation_inputs_2026-09-09.json.gz`; no full spreadsheet or network access is needed.
