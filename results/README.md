# Deposited records for manuscript version 1.1

These files accompany **How to Compute the Unknotting Number: Theory and
Computations**. They distinguish input ranges, completed lower-bound verdicts,
explicit upper-bound certificates, and exploratory data.

The current aggregation entry point is `../consolidate_results.py`. It reads
deposited inputs and writes only to `../generated/`; it does not replace the
historical aggregates in this directory.
The rebuilt counts are 1,225 exact values and 176 further improved ranges.

## Frozen comparison ranges

- `paper_v1_1_snapshot.json`: the 12,965 comparison intervals used in the paper.
- `paper_v1_1_manifest.json`: file checksums and per-record input provenance.
- `reference_table_2026-09-07.json`: the older archived aggregate on which the
  comparison is based. The 27 original input ranges in the deposited McCoy
  records are retained explicitly in the v1.1 snapshot. This counts their
  lower-bound contribution in Appendix E, rather than absorbing it into the
  initial range.
- `u_table_2026-09-07.json` and `consolidated_2026-09-07.json`: historical outputs.
  Use newly generated outputs for v1.1; the historical aggregate omits the
  27 McCoy contributions from its improvement count.

The frozen comparison is not a claim about what a live KnotInfo page records
today. The manifest is also not a proof of the mathematics: its checksums detect
changes to the input files used to assemble the paper's tables.

## Lower bounds — Section 3

| Files | Interpretation |
| --- | --- |
| `lickorish_cw_obstructed_815_knotinfo2026.8.1.json` | 815 recorded linking-pairing obstructions to unknotting number one; overlap with other methods is counted once |
| `torsion_lower_bounds_2026-08-24.json`, `torsion/` | Knot Floer torsion computations and supplementary historical torsion experiments; only the HFK contributions enter the paper's totals |
| `owens_u3_2026-08-24.json`, `owens_sweep/`, `generator_bound_2026-08-24.json` | Initial correction-term, Traczyk, and homology calculations excluding sequences of two crossing changes |
| `owens_rank2/` | Completed rank-two comparisons for 1,438 alternating knots with absolute signature four |
| `owens_rank3/` | Rank-three controls and applications; after consolidation, 129 exact values with `u=4` and 12 further improvements |
| `owens_rank4/` | Explicitly frozen partial output and two JSONL snapshots; completed exclusions give all 21 values with `u=5`, including the six added results |
| `cyclic_cover/` | Nakanishi's generator bounds from cyclic branched-cover homology |
| `montesinos/` | Correction terms and reduced Floer homology of Montesinos plumbings; 44 exact values and five further improved ranges as the primary method |
| `crossing_changes/mccoy_alternating_u1_2026-09-08.json` | 27 alternating knots, their original input ranges, and the completed McCoy obstruction |
| `comparison/gebel_prangley_2026.json` | Theorem lists used in the comparison with concurrent work in the conclusion |

`PASS` means that an obstruction is inconclusive. Only completed exclusions
supply lower bounds. Duplicate records from resumed jobs are deduplicated;
conflicting completed results are an error. The uncompleted rank-four targets
`13a_4770`, `13a_4787`, and `13a_4801` retain their input ranges.

## Upper bounds and minimal diagrams — Section 4

| Files | Interpretation |
| --- | --- |
| `summary.json`, `best/` | Eight presentation unknotting diagrams, with PD codes and the two marked rows used in Appendix G |
| `witness_chains/` | Replayable search certificates; the presentation verifier checks the eight final paper diagrams separately |
| `13n_1587_u_le_2_braid_certificate.json` | Braid crossing-change certificate following Brittenham–Hermiller |
| `cluster_runs/` | Archived search runs; an unsuccessful run establishes no lower bound |
| `crossing_changes/dataset_v2.json.gz` and accompanying analyses | Single crossing changes in the reference alternating diagrams, identified children, signatures, determinants, and crossing features |
| `open23/priority_u23_final.json` | Source of Appendix F: 806 knots without candidate crossings after excluding records that rely on heuristic identifications |
| Other files in `open23/` | Flype checks, child-knot data, and supplementary exploratory classifications |

The Appendix F conclusion is conditional on the Bernhard–Jablan property; its
806 knots are not part of the unconditional exact-value count. Historical
comparison and exploratory files, including those in `dkt/`, are preserved as
archival records and are not additional results of v1.1. Per-shard working files
are not deposited.

The verification commands in the directory READMEs exercise the mathematical
kernels, known examples, and deposited unknotting certificates. Rebuilding these
aggregate files does not rerun every obstruction or reidentify every crossing
neighbor used in the conditional Appendix F analysis.
