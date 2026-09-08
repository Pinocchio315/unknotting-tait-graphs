# Recheck of the seven incompatible KnotInfo entries

These audit files supplement the manuscript's historical comparison. They do
not change its baseline or treat a database value as a mathematical certificate.

The official KnotInfo download page was consulted on 8 September 2026. Its
database archive was fetched at 10:13:55 UTC and compared byte for byte with the
archive pinned to official asset commit
`08265335f9f5463796c36fe7e5457af1f489179e` (2 September 2026, 11:22:11 UTC).
Both files have SHA-256
`6b97fae54e8e2e3f1e134a569f80176946f7b0b49920f2161b6049335f8ac275`.

The seven individual knots were also queried through the official website's
current single-knot GET search form at 10:24:42–10:25:07 UTC on the same day.
All seven web responses match the download, so the comparison does not merely
assume the downloadable spreadsheet matches the live search backend.

- `official_today_comparison.json` records the comparison and full provenance.
- `knotinfo_official_2026-09-08.json.gz` freezes all 12,965 study-population
  unknotting-number entries extracted from that official download.
- The header metadata and asset commit JSON record the HTTP response and source revision. Header line endings are normalized and transient session cookies are omitted.
- `live_web_comparison.json`, `live_web_*.html.gz`, and their HTTP headers retain
  the directly queried web results, URLs, retrieval dates, and response hashes.
- `verify_official_snapshot.py` reproduces the comparisons from the frozen
  extract and freshly reconstructed manuscript results. With `--xls-zip FILE`,
  it also verifies the archive hash and repeats the original Excel extraction;
  this optional operation needs `xlrd`.

All 12,965 unknotting-number ranges match the September package snapshot. There
remain **seven incompatible knots**, comprising **six exact manuscript values**
and **one manuscript interval**, namely `13a_650` with `[3,4]`. Thus “six exact
results” and “seven incompatible knots” count different, explicitly stated sets.
For all 1,227 exact manuscript results, the current official download has 194
matching exact entries, six different exact entries, and 1,027 unresolved entries.

## Independent upper-bound witnesses

`upper_bound_recheck.json` contains explicit crossing subsets in the official
minimal PD diagrams that establish all seven historical upper bounds:

| Knot | Verified upper bound | Changed PD rows (zero-based) |
| --- | ---: | --- |
| `12a_107` | 4 | 0, 2, 4, 7 |
| `13a_15` | 3 | 0, 6, 8 |
| `13a_55` | 3 | 3, 6, 11 |
| `13a_422` | 3 | 6, 7, 9 |
| `13a_568` | 3 | 6, 7, 8 |
| `13a_650` | 4 | 0, 3, 6, 7 |
| `13a_660` | 4 | 0, 3, 6, 7 |

Every source PD was compared exactly with the freshly downloaded official Excel
entry. Cyclically rotate each selected PD quadruple by one position to change
that crossing. Spherogram reduces every changed diagram to zero crossings using
Reidemeister moves. Independently, Regina reduces each to a one-component,
zero-crossing link. No knot identification based only on Jones polynomials is
used. `recheck_upper_bounds.py` repeats this bounded, single-process search and
both verifications; `--live-xls FILE` repeats the official PD comparison.

These witnesses establish upper bounds only. Their existence does not prove
minimality; the lower-bound obstruction audits must supply that separate step.

## Independent lower-bound audits

The rank-two audit uses independent face-based checkerboard graphs, complete
characteristic-vector boxes, adjugate labels, and direct generator-image
isomorphisms. All 22 candidate forms are excluded: two have a different
discriminant group, and the other 20 fail all 5,984 possible group isomorphisms.
The determinant-225 case includes its noncyclic group Z/3 + Z/75.

For rank three, `rank3_enumeration.py` solves the determinant equation for the
last diagonal entry and uses all 168 binary unimodular lifts. It independently
recovers 22 candidates at determinant 239 and 23 at determinant 155.
`rank3_audit.py` computes the Goeritz correction terms by a complete finite box
and excludes all 7,996 cyclic group isomorphisms with explicit characteristic
vectors. It does not reuse the production correction-term minimizer or matcher
for these exclusions. A vector's quadratic value is an upper bound for m_Q;
if it is smaller than the target d-value, the required inequality is impossible.
Values in one characteristic class differ by an even integer, so an incorrect
congruence is also decisive without finding a shortest vector.

All seven knot identities agree exactly with independently tabulated named
diagrams using Regina's canonical diagram signature (mirrors allowed).
Signatures and determinants also agree with exact recomputation from Seifert
matrices. The theory review, including determinant square factors, orientation,
and completeness of form enumeration, is in `theory_notes.md`.

Run with the repository's topology Python from the repository root:

```bash
python results/review_seven_conflicts/verify_official_snapshot.py
python results/review_seven_conflicts/rank2_audit.py
python results/review_seven_conflicts/rank3_enumeration.py
python results/review_seven_conflicts/rank3_audit.py
python results/review_seven_conflicts/recheck_upper_bounds.py
```

The summary JSON files are accompanied by compressed full correction-term
tables and failing-map witnesses. These single-process audits use the frozen
inputs and produce no changes to the manuscript or original result deposits.
The new audit outputs are kept separate from the paper's counting inputs.

The known example `7_1`, with unknotting number three, is not excluded by the independent rank-three witness scan; the surviving units are 2 and 5. Its bounded comparison is retained in `rank3_positive_control.json`.
