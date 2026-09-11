# Additions integrated into manuscript v1.8

These inputs add **15 exact values** to the earlier profile: four from the
lower-bound calculations and eleven from upper-bound constructions. The
aggregate is 2,734 exact values and 380 improved non-exact ranges. All fifteen
were unresolved in the official KnotInfo snapshot of 9 September 2026.

## Lower-bound additions: Section 3.8

| Knot | Previous range | New value | Obstruction |
| --- | --- | ---: | --- |
| 13n1619 | [1,2] | 2 | Casson–Walker sum, Greene candidates, half-integral surgery |
| 13n4876 | [4,5] | 5 | Sum and spin constraints, rank-four surgery forms |
| 13n4973 | [4,5] | 5 | Unique Greene vector, rank-four surgery forms |
| 13n5102 | [4,5] | 5 | Unique Greene vector, rank-four surgery forms |

The [lower-bound guide](lower/README.md) gives the hypotheses, normalization,
citations, exact certificates, and replay commands. The rank-four replay
regenerates all 9,353 forms; noncyclic discriminant groups are rejected before
testing cyclic identifications. Every remaining form and group identification
has a characteristic-covector obstruction. Independent upper constructions
give the matching values.

## Upper-bound additions: Section 4.3 and Appendix G

`presentation.json` is the final figure/PD source. Marks are **zero-based rows
of the specified PD**, not crossing labels of a separately relabelled diagram.

| Knot | Crossings | Marked rows | Result of changing the marks |
| --- | ---: | --- | --- |
| 13n221 | 13 | 3, 5 | unknot |
| 13n436 | 13 | 3, 4 | unknot |
| 13n447 | 13 | 8 | 10_91, with independently known u=1 |
| 13n636 | 13 | 3, 4 | unknot |
| 13n1439 | 13 | 6, 8 | unknot |
| 13n2251 | 13 | 6, 8 | unknot |
| 13n2639 | 13 | 8, 10 | unknot |
| 13n3108 | 13 | 2, 6 | unknot |
| 13n3733 | 13 | 9, 11 | unknot |
| 13n4025 | 13 | 2, 6 | unknot |
| 13n4237 | 18 | 1, 12 | unknot |

All eleven give u=2. For ten, the lower bound two comes from the Greene
sweep; these entries carry `G,w` in Appendix D. For 13n3733 it is already in
the archived comparison ranges and the official 9 September snapshot, both
of which record [2,3]. An earlier package's different value is not used.

The original certificates in `runs190/` have eight direct unknot endpoints
and three named partners. For 13n3108 and 13n3733, subsequent verification
found the direct pairs displayed above in the same reduced PDs. The original
partner proofs remain available. For the supplied PD of 13n447, all 78 pairs
of crossing changes have determinant different from one; the one-change
construction through 10_91 is the valid proof of u(K)≤2. This does not claim
anything about all minimal diagrams of that knot.

`reverification.json` records source identifications by actual
meridian-preserving isometries, endpoint verifications, and the 13 neighbours
of each of the reference diagrams of 13n2251, 13n4025, and 13n4237. These
neighbours all have u=2 in the consolidated data. The new minimal diagrams of
13n2251 and 13n4025 admit a neighbour with u=1, so the neighbour set can depend
on the minimal diagram. No exhaustive minimal-diagram enumeration is claimed
for 13n4237; it supplies an upper construction, not a new BJ counterexample.

## Files and reproduction

| File/directory | Contents |
| --- | --- |
| `bounds.json` | Four new lower bounds and eleven upper bounds, with source attribution |
| `lower/` | Lower-bound computation modules and exact certificates |
| `targets190.json` | Frozen list for the additional randomized search |
| `runs190/` | 190 execution rows and 11 successful search certificates |
| `min11/` | Reduced PDs, marks and reduction summaries |
| `min11p/` | Further pass-move reduction records for 13n436 and 13n4237 |
| `presentation.json` | The eleven final diagrams used in the manuscript |
| `partner_direct_pair_check.json` | Verified direct pairs for 13n3108 and 13n3733 |
| `13n_447_supplied_pd_pairs.json` | Determinants of all pairs in that one PD |
| `pass_replay_13n_436.json` | Exact five-move replay from 14 to 13 crossings |
| `reverification.json` | Independent final-PD, endpoint and neighbour identifications |

From the repository root:

```sh
python3 consolidate_results.py --manuscript v1.8 --outdir generated/v1_8
python upper_bounds/verify_presentation_diagrams.py --manuscript v1.8
python upper_bounds/verify_certificates.py results/extensions_2026-09-10/runs190 --out /tmp/verified_v18
python upper_bounds/draw_witness.py --outdir /tmp/unknotting_figures
```

The same random diagram exploration is used for the earlier eight and these
eleven constructions. Actual settings for the 190-target execution were
10,000 additional diagrams, a 25-crossing cutoff, a 600-second budget per
knot, three workers, and seed zero. The 190 rows sum to 1,799,825 diagram
visits and 37,684.2 seconds of per-knot elapsed time, not CPU time. The initial
diagram is tested separately, so a record can contain 10,001 visits. The
[upper-bound guide](../../upper_bounds/README.md#recorded-v18-run) gives the
full rerun command and safe marked-pass minimization procedure.

The input manifest is `../paper_v1_8_manifest.json`; the base inputs remain
under the v1.3 manifest. Aggregation keeps the lower-bound stage separate
from the final stage: 2,715 exact values plus 390 narrowed ranges first,
then nineteen constructions, ten of which promote those narrowed ranges.
This yields 2,734 exact values and 380 narrowed ranges, with no double count.
