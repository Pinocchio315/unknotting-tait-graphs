# Lower-bound computations

This directory contains the lower-bound computations in Sections 3 and 5 of **Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?**, manuscript v1.8. The original `v1.3` input profile is preserved; four later lower-bound determinations are supplied by the v1.8 extension certificates. See the [repository README](../README.md) for installation and the comparison with the official KnotInfo snapshot retrieved on 9 September 2026.

Each method tests a necessary condition for an unknotting sequence. An obstruction raises the lower bound; an independently established upper bound is needed to obtain an exact value. A passing test or a resource limit supplies no new bound.

| Script | Location in the paper | Mathematical role |
|---|---|---|
| `scan_lower_bounds.py` | Linking pairings | Lickorish's cyclicity and self-linking conditions from a Goeritz presentation. Casson–Walker quantities remain auxiliary diagnostic output; signs from the two conventions are not intersected without an orientation identification. |
| `seifert_linking_check.py` | Linking pairings | The same pairing from `V + Vᵀ`, where `V` is a Seifert matrix. This gives a different presentation of the cover's homology and pairing. |
| `regina_linking_check.py` | Linking pairings; Summary and verification | Regina computes the pairing from a triangulated branched cover. The cover is constructed by lifting the order-two meridional orbifold filling in SnapPy, so the manifold is not selected by its homology order. |
| `linking_pairing.py`, `tait_tools.py` | Linking pairings | Exact finite-group algebra and Goeritz helpers. A generator is constructed by Chinese remaindering the primary components; a search among a few short vectors is not assumed exhaustive. |
| `scan_torsion_lower_bounds.py` | Torsion bounds | The maximal `U`-torsion order in knot Floer homology over `F₂[U]`, bounded above by the unknotting number. The reported torsion results are HFK results. Rational Lee `X`-torsion and characteristic-two Bar–Natan `h`-torsion are reviewed in the Background and are distinct invariants. |
| `owens/` | Two crossing changes; Three and four crossing changes | Correction-term obstructions for definite fillings and their lattice tests; see [owens/README.md](owens/README.md). |
| `generator_bound_sweep.py` | Two crossing changes; Cyclic-cover bounds | The minimum number of generators of `H₁(Σ₂)`, with the three-rank checked against the Lickorish–Millett Jones evaluation. |
| `cyclic_cover_bound.py` | Cyclic-cover bounds | Nakanishi's bound `u ≥ ceil(g_n/(n−1))`, using the tabulated homology of covers of degrees two through nine. Zero invariant factors are free summands and count toward `g_n`. |
| `montesinos/montesinos_u1.py`, `montesinos/hf_red.py` | Montesinos knots | Negative-definite star plumbings, correction terms, the Ni–Wu half-integral surgery pattern and the reduced Floer mapping cone. See [montesinos/README.md](montesinos/README.md) for hypotheses and exact termination. |
| `greene/qa_search.py`, `greene/greene_dinv.py`, `greene/pin_dinv.py` | Sections 2, 3 and 5 | Greene's Kauffman-state gradings and Spin^c structures give finite candidate sets for correction terms when the double branched cover is an L-space. Reduced Khovanov rank over `F_2` equal to the determinant certifies that premise. The half-integral surgery test excludes unknotting number one for `12n_491`, `13n_3370` and `13n_1587`; published upper bounds give value two for all three. See [greene/README.md](greene/README.md). |
| `greene/greene_sweep.py`, `greene/greene_ranks.py` | Section 3 | Apply the same half-integral and Owens rank-two, rank-three and rank-four comparisons to the frozen nonalternating target list. Ambiguous correction terms must be handled conservatively; no timeout, incomplete form enumeration or surviving candidate produces an obstruction. |

The deposited pairing list contains 815 obstructions. The HFK scan records order two for `13n_689`, `13n_1166`, `13n_2504` and `13n_2807`. The Montesinos computation excludes unknotting number one for 49 of its 50 targets: 44 have an independently recorded upper bound two and five have upper bound three. These raw records retain their original filenames in `../results/`; the aggregate table is assembled by `../consolidate_results.py`.

With the v1.8 extensions, the Greene sweep supplies 1,711 lower-bound
improvements: 1,496 additional exact values and 215 narrower ranges. The
exact values comprise 1,310 values two, 136 values three, 43 values four,
and seven values five. Including the two original BJ calculations, Greene's
model contributes 1,498 exact values before the new upper bounds are applied.

Across all methods, lower bounds determine 2,715 exact values and narrow 390
ranges. Nineteen explicit constructions then give 2,734 exact values and 380
remaining ranges; ten need both the lower-bound improvement and the new
construction. Against the official 9 September snapshot, 2,540 exact values
and 323 range improvements are new to that snapshot; 194 exact values and
57 ranges already agree. Run `python compare_knotinfo.py` for this comparison.

See the [four new lower-bound calculations](../results/extensions_2026-09-10/lower/README.md)
for the Casson--Walker constraint, rank-four covectors, and replay commands.
The [obstruction comparison](obstruction_comparison/README.md) reproduces
Section 3.9's comparison of symmetry and full surgery conditions on 2,776 knots.

Owens's general theorem uses the signed unknotting-sequence hypothesis and
the labelled correction terms of the branched cover; it does not require the
knot to be alternating or Montesinos. The L-space condition is used here to
extract correction-term candidates from solitary states. It does not guarantee
that every candidate set becomes a singleton. Spectral-sequence differentials
respect the Spin^c decomposition and lower the Maslov grading by one.

The reduced mod-2 Khovanov rank equality and cyclic-homology filters select
4,264 of the 6,236 nonalternating
prime knots with at most 13 crossings: 31.6% lie outside this implementation
on these grounds alone. The input ranges and signatures then select 1,967
undetermined targets for the half-integral and rank-two through rank-four tests.
For a rank-`n` obstruction, the spin correction term must be determined to be
`-n/2`, and every remaining candidate vector must fail every admissible form.

For the 959 alternating knots in Appendix F, the recorded neighbour bounds
exclude unknotting number one after every crossing change in a minimal diagram.
Their values remain in `[2,3]`; if any has value two, it fails the strong
Bernhard–Jablan equality. These knots are not added to the exact-value total.
Forty-six further cases require unverified identifications and remain outside
that list; 22 other parents have 11 named possible children.

## Running fresh checks

Run these commands from the repository root with the project's Python environment.
Start with the bounded regression suite and a few named Montesinos examples:

```sh
python -m unittest discover -s tests -p test_lower_bounds.py -v
python lower_bounds/montesinos/montesinos_u1.py 3_1 4_1 11n_102 12n_457 --out /tmp/montesinos_controls.json
```

The following commands perform larger computations. Their explicit destinations
keep fresh output separate from the deposited results:

```sh
python lower_bounds/seifert_linking_check.py --workers 4 --out /tmp/seifert_check.json
python lower_bounds/regina_linking_check.py --workers 4 --out /tmp/regina_check.json
python lower_bounds/cyclic_cover_bound.py --out /tmp/cyclic_cover_bound.json
python lower_bounds/montesinos/montesinos_u1.py --apply --workers 3 --out /tmp/montesinos_apply.json
```

The 19-test regression suite checks explicit abelian groups, exact ellipsoid boundaries, lens-space correction terms, known unknotting-number-one knots, known HFK torsion orders and Lidman's reduced-Floer obstruction for `11n_102`. It also checks `12n_457` against rational correction terms computed by a separate exact calculation. Agreement with a deposited list is a reproducibility check, not the independent mathematical control.

The September 2026 audit recomputed the 3,002 Seifert presentations: exactly 815 obstructions, with no undecided cases or determinant mismatches. The separate 1,516-knot Seifert control set with known unknotting number one gives no false obstruction or uncomputed case. All 12,965 cyclic-cover records were reproduced, and the four reported HFK torsion orders were recomputed as two. It also reran all 50 Montesinos targets with the exact lattice algorithm and reproduced all 49 obstructions, with no errors. The only target requiring reduced homology is `12n_457`; `12n_309` remains unexcluded. The fresh runs did not replace deposited raw files. Historical full-table logs and the bounded regression suite cover different computations; neither should be read as a rerun of every deposited result.

Reconstruct the full result table with `python consolidate_results.py`; use
`--profile pre-sweep` for the ranges used to select the enlarged Greene cohort.
