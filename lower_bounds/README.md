# Lower-bound computations

This directory accompanies **How to Compute the Unknotting Number: Theory and Computation**, Sections 3 and 5.1. The subsections of Section 3 are unnumbered; the names below identify the corresponding arguments. Every obstruction is a necessary condition for an unknotting sequence. A passing test does not prove an upper bound, and a resource limit supplies no obstruction.

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
| `greene/qa_search.py`, `greene/greene_dinv.py`, `greene/pin_dinv.py` | Section 5.1 | Greene's Kauffman-state gradings and Spin^c structures give candidate correction terms when the double branched cover is an L-space. For `12n_491` and `13n_3370`, the reduced Khovanov homology ranks over `F_2` equal the determinants and establish this premise. The half-integral surgery test excludes unknotting number one, and the published upper bounds give `u(12n_491) = u(13n_3370) = 2`. See [greene/README.md](greene/README.md). |

The deposited pairing list contains 815 obstructions. The HFK scan records order two for `13n_689`, `13n_1166`, `13n_2504` and `13n_2807`. The Montesinos computation excludes unknotting number one for 49 of its 50 targets: 44 have an independently recorded upper bound two and five have upper bound three. These raw records retain their original filenames in `../results/`; the aggregate table is assembled by `../consolidate_results.py`.

## Running fresh checks

Run these commands from the repository root with the project's Python environment. Explicit destinations keep fresh computations separate from deposited results.

```sh
python -m unittest discover -s tests -p test_lower_bounds.py -v
python lower_bounds/seifert_linking_check.py --workers 4 --out /tmp/seifert_check.json
python lower_bounds/regina_linking_check.py --workers 4 --out /tmp/regina_check.json
python lower_bounds/cyclic_cover_bound.py --out /tmp/cyclic_cover_bound.json
python lower_bounds/montesinos/montesinos_u1.py --apply --workers 3 --out /tmp/montesinos_apply.json
python lower_bounds/montesinos/montesinos_u1.py 3_1 4_1 11n_102 12n_457 --out /tmp/montesinos_controls.json
```

The 19-test regression suite checks explicit abelian groups, exact ellipsoid boundaries, lens-space correction terms, known unknotting-number-one knots, known HFK torsion orders and Lidman's reduced-Floer obstruction for `11n_102`. It also checks `12n_457` against rational correction terms computed by a separate exact calculation. Agreement with a deposited list is a reproducibility check, not the independent mathematical control.

The September 2026 audit recomputed the 3,002 Seifert presentations: exactly 815 obstructions, with no undecided cases or determinant mismatches. The separate 1,516-knot Seifert control set with known unknotting number one gives no false obstruction or uncomputed case. All 12,965 cyclic-cover records were reproduced, and the four reported HFK torsion orders were recomputed as two. It also reran all 50 Montesinos targets with the exact lattice algorithm and reproduced all 49 obstructions, with no errors. The only target requiring reduced homology is `12n_457`; `12n_309` remains unexcluded. The fresh runs did not replace deposited raw files. Historical full-table logs and the bounded regression suite cover different computations; neither should be read as a rerun of every deposited result.
