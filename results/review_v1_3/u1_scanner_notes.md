# Audit of the added unknotting-crossing scan

## Correct scope and mathematical conclusion

The program scans one KnotInfo PD per knot. It does not enumerate other
minimal diagrams, either by flypes or by any non-alternating diagram search.
The phrase “all minimal diagrams were scanned” would therefore overstate the
computation. The original module docstring already stated this limitation.

The distinction has different consequences for positive and negative results:

- An exhibited unknotting crossing in the selected minimal diagram is enough
  to prove `u_BJ^s(K) = 1` when `u(K) = 1`. No enumeration of other minimal
  diagrams is needed for this existential statement.
- Failure to find an unknotting crossing after deciding every crossing of
  the selected diagram proves only that this diagram has none. For a
  non-alternating knot it does not prove that all minimal diagrams have none,
  and cannot by itself produce a Bernhard--Jablan counterexample or prove
  `u(K) >= 2`.

The scan consequently supports a useful location for further searches, but
does not settle the unknotting numbers of the 546 open-range knots.

## Independent full-cohort verification

The two archived name sets agree exactly with the current consolidated
comparison table: 1,516 knots with range `[1,1]` and 546 knots with lower
endpoint one and larger upper endpoint. They are disjoint, have no duplicate
crossing labels and contain no undecided record. The known group consists of
710 alternating and 806 non-alternating knots; the open group is entirely
non-alternating. Every supplied PD has the tabulated crossing number.

`u1_scanner_audit.py` independently reproduced the actual positive and
negative statements, avoiding the scanner's Tait-graph reducer:

1. For all 1,516 positive records, one archived crossing was changed directly
   in the original PD. Spherogram reduced every resulting one-component
   diagram to zero crossings. Thus all 1,516 existential witnesses were
   independently verified; no polynomial-only knot identification or HFK
   call was needed for these positive checks.
2. For all 546 negative records, a separate exact Goeritz construction and
   rank-one determinant formula evaluated all 7,052 crossing changes. Of
   these, 6,542 have determinant different from one and hence are not unknots.
3. The remaining 510 direct PD crossing changes were independently simplified
   and passed to the knot Floer homology program. Their Seifert genera were
   all positive: 233 of genus two, 275 of genus three, and two of genus four.
   No determinant-one case was left untested by the bounded replay.

The combined verification took about 19 seconds in the repository environment.
The source hashes, all open-diagram crossing determinants, positive witness
indices and all 510 computed HFK rank tables are saved in
`u1_scanner_audit.json`. This evidence proves the stated results for the
chosen diagrams; it provides no coverage of other minimal diagrams.

The determinant-one filter is only a necessary condition. The negative
decision in the last step uses the genus-detection theorem, not `tau`, the
signature or the Alexander polynomial. The relevant primary sources are
[Ozsvath--Szabo, Holomorphic disks and genus bounds](https://arxiv.org/abs/math/0311496)
and the [HFK program's official repository](https://github.com/3-manifolds/knot_floer_homology).
For a knot in `S^3`, Seifert genus zero is equivalent to the unknot.

## Implementation corrections

The original top-level script ran immediately on import, and even a named
one-knot invocation overwrote the entire published 1,516-entry JSON file.
It also lacked a reproducible CLI selection for the separate 546-knot open
cohort. These workflow defects were corrected:

- Computation is isolated in `scan_row` and an import-safe `main`.
- `--cohort known` and `--cohort open` explicitly select the two tables;
  positional names remain available and duplicate arguments are removed.
- Without `--out`, the program prints its result and does not change an
  archived file. Writing requires an explicit output path.
- Each new result records its input PD, scope, source version/table hash and
  per-crossing evidence: exact determinant, replayable Tait reduction path,
  or HFK input and rank table. Existing published archives were preserved.
- The input is required to be a one-component PD with the tabulated crossing
  number and determinant. HFK output must have a nonnegative integer genus,
  the correct coefficient field, a complete rank count and support matching
  that genus. Genus-zero output must have unknot rank one.
- Recognition or reduction failures produce an undecided crossing; they
  cannot silently become a negative result. Determinants are computed before
  simplification, so the inexpensive exclusion does not depend on reduction.

All eight focused regressions in `test_u1_scanner.py` pass. They include a
determinant-one nontrivial crossing of `11n116`, a forced HFK failure,
malformed HFK data, input inconsistency, explicit open-cohort selection and
protection against accidental archive replacement.

## Reproduction

From the repository root, with its Python environment:

```sh
python results/review_v1_3/u1_scanner_audit.py
python -m unittest discover -s results/review_v1_3 -p test_u1_scanner.py -v
python crossing_changes/u1_minimal_diagram_scan.py 3_1 11n_116 --out /tmp/chosen_diagrams.json
```

No background computation from this audit remains running.
