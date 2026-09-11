# Computational inputs and results

These records accompany [Computation of Unknotting Numbers: Which Knot Breaks
the Bernhard-Jablan Conjecture?](../README.md#paper-and-citation). They contain
input ranges, completed obstruction verdicts, crossing-change data, and
upper-bound certificates. `consolidate_results.py` combines them into a table
of intervals without accessing a manuscript or a live database.

## Input ranges and provenance

- `paper_v1_1_snapshot.json` contains the historical starting ranges for
  12,965 knots. It is a reconstructed research table: the August 2026 package
  ranges are retained except for 104 entries restored to their April ranges,
  including 27 lower bounds restored from two to one. It is not an unmodified
  KnotInfo release.
- `paper_v1_1_manifest.json` authenticates the inputs before the enlarged Greene
  sweep. `paper_v1_3_manifest.json` also includes that sweep and the updated
  alternating-diagram analysis. Their checksums and per-record provenance are
  computational inputs, retained with their original filenames.
- `reference_table_2026-09-07.json` is the earlier aggregate used when assembling
  those starting ranges. The manifests record the overrides.
- `knotinfo_2026-09-09.json` contains the official ranges for the same study
  population, retrieved on 9 September 2026. Its adjacent manifest records the
  source URL, source commit, retrieval time, and checksums. These ranges are
  used for numerical comparison, not to retroactively change the search inputs.
- `greene/greene_targets_2026-09-08.json` is the input to the enlarged Greene
  sweep: diagrams, determinants, signatures, reduced mod-2 Khovanov vectors,
  and pre-sweep intervals for 1,967 targets and 1,452 controls.

Reconstruct intervals from the repository root:

```bash
python consolidate_results.py
python compare_knotinfo.py --out /tmp/knotinfo_comparison.json
```

The default v1.8 profile gives 2,734 exact values and 380 narrower ranges
relative to the historical starting table. Of these, 2,540 exact values and
323 range improvements remain improvements over the dated official snapshot;
194 exact values and 57 ranges are already recorded there. Computed intervals
and per-knot sources are written to `generated/` as JSON.

To reconstruct the ranges used when selecting the enlarged Greene cohort:

```bash
python consolidate_results.py --profile pre-sweep --outdir /tmp/unknotting-before-greene
python lower_bounds/greene/build_greene_targets.py --table /tmp/unknotting-before-greene/u_table.json --out /tmp/unknotting-before-greene/greene_targets.json
```

The frozen worker input can be used directly for replays. A table containing
the sweep's final bounds must not be used to reconstruct its input selection.

## Lower-bound records

| Location | Calculation |
| --- | --- |
| `lickorish_cw_obstructed_815_knotinfo2026.8.1.json` | Double-cover linking-pairing obstructions |
| `torsion_lower_bounds_2026-08-24.json`, `torsion/` | Knot Floer torsion bounds and underlying homological output |
| `owens_u3_2026-08-24.json`, `owens_sweep/` | Initial correction-term comparisons and Traczyk's criterion |
| `owens_rank2/` | Rank-two comparisons for alternating knots |
| `owens_rank3/` | Rank-three obstructions and controls |
| `owens_rank4/` | Completed rank-four verdicts in frozen partial output and JSONL snapshots |
| `generator_bound_2026-08-24.json`, `cyclic_cover/` | Homology generator bounds for branched covers |
| `montesinos/` | Plumbing correction terms, reduced Floer homology, and surgery comparisons |
| `greene/` | Frozen inputs, worker records, and summary of the enlarged correction-term sweep |

`PASS` is an inconclusive obstruction test. Missing, timed-out, or incomplete
calculations do not supply bounds. Resumed jobs are deduplicated, and
conflicting completed records stop consolidation. For ambiguous correction
terms, every candidate vector must be excluded before raising the lower bound.
The Greene worker and summary code also check that controls are not obstructed.

## Upper-bound certificates and Bernhard-Jablan examples

| Location | Contents |
| --- | --- |
| `summary.json`, `best/` | The eight earlier displayed diagrams, with marked crossing rows and PD codes |
| `witness_chains/` | Replayable diagram-move and crossing-change certificates |
| `13n_1587_u_le_2_braid_certificate.json` | The published braid construction with upper bound two |
| `bernhard_jablan/brittenham_hermiller_upper_bounds.json` | Published upper bounds for `12n_491` and `13n_3370`, with diagram data and source locations |
| `bernhard_jablan/khovanov_inputs_2026-09-08.json` | Reduced mod-2 Khovanov ranks and diagrams certifying the L-space premise |
| `bernhard_jablan/greene_d_*.json` | Correction terms or finite candidate sets and the resulting surgery tests |
| `bernhard_jablan/qa_certificate_12n_491_2026-09-08.json` | A replayable quasi-alternating resolution tree |
| `cluster_runs/` | Deposited upper-bound search output |

For `13n_3370`, the archived starting range is `[1,3]`; the sharper published
upper bound two is a separate input. For `13n_1587`, the historical braid file
retains its original descriptive fields and checksum. Its construction is
from Brittenham-Hermiller Section 3; the Applebaum attribution in that archival
file is not the original source. The later Greene lower bound gives value two.

Use `upper_bounds/verify_presentation_diagrams.py` for all nineteen displayed
diagrams and the certificate-verification commands in the
[upper-bound guide](../upper_bounds/README.md) for move sequences. A stored
verification flag is not a substitute for replaying the certificate.

## Minimal diagrams and crossing changes

- `crossing_changes/dataset_v2.json.gz` records 69,632 crossing changes in
  5,546 alternating parent diagrams. Signatures, determinants, identifications,
  and crossing features are kept with the corresponding analysis records.
- `crossing_changes/mccoy_alternating_u1_2026-09-08.json` records the 27
  starting-table improvements. `mccoy_alternating_all_2026-09-08.json` records
  the full April lower-bound-one cohort: 710 diagrams with an unknotting
  crossing and 2,917 without one.
- `open23/priority_u23_final.json` supplies the current alternating analysis:
  959 certified exclusions of a neighbour with unknotting number one. The
  46 cases requiring unverified identifications remain outside that count.
- `open23/priority_u23_paper_v1_3.json` preserves the pre-sweep cohort required
  by its input manifest. Other records in `open23/` are intermediate neighbour
  classifications, flype-orbit data, or candidate-resolution results.
- `comparison/crossing_formula_verification.json` records the exact rational
  formula comparison reproduced by `audit_crossing_formulas.py`.
- `comparison/gebel_prangley_2026.json` is the published theorem list preserved
  by the computational source manifests.

A Jones-polynomial match alone does not identify a knot. The candidate code
keeps such hypotheses separate from bounds justified by identification or
invariants. The 959 conditional cases are not counted as exact determinations.

## v1.8 additions

`paper_v1_8_manifest.json` authenticates the new proof records and replay
code. [extensions_2026-09-10/README.md](extensions_2026-09-10/README.md)
lists the four lower-bound and eleven upper-bound additions, with full PDs,
marks, search certificates, covectors, and reproduction commands. These
increase the exact count by fifteen and resolve ten previously narrowed ranges.

`knotinfo_calculation_inputs_2026-09-09.json.gz` contains seven unaltered
spreadsheet columns for the 4,110 knots needed by the Greene/Montesinos
calculations, including the small controls. Its manifest gives the selected
fields, upstream commit, retrieval date, spreadsheet checksum, and extract
checksum. `knotinfo_inputs.py` authenticates the extract before loading it.
The complete spreadsheet and the old database-review web captures are omitted.

`obstruction_comparison/` contains the per-knot comparison and summary for
Section 3.9: 2,776 distinct knots, with no additional exclusion by the full
conditions over positive even symmetric matching. The code in
`../lower_bounds/obstruction_comparison/` regenerates all candidate-level
outputs from the original correction-term records and selected KnotInfo cells.
The generated cache and execution logs are not part of the distribution.
