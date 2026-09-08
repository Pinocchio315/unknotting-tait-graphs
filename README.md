# Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?

Companion code and deposited computational records for the manuscript by
Seong-Jin Lee, version 1.3 (September 2026).

The paper explains how topological obstructions and explicit crossing changes
work together. An obstruction raises a lower bound for `u(K)`; an unknotting
diagram gives an upper bound. An exact value follows when the two bounds agree.
The code follows this distinction: an unsuccessful search supplies no lower
bound, and a correction-term verdict `PASS` supplies no upper bound.

Relative to the archived comparison ranges used in the manuscript, the deposited
computations give **2,719 exact values** and **390 further improved ranges** for
prime knots with at most 13 crossings:

| Unknotting number | Exact values |
| --- | ---: |
| 2 | 2,149 |
| 3 | 373 |
| 4 | 172 |
| 5 | 25 |

Lower bounds account for 2,711 exact values and all 390 further improvements.
Eight unknotting diagrams account for the remaining exact values. The expanded
Greene computation produces 1,707 lower-bound improvements: 1,492 additional
exact values and 215 narrower ranges. Together with the two earlier Greene
cases, method `G` accounts for 1,494 exact values.
These counts include independently recovered results, as explained in the paper;
they are not a claim that every listed value first appears here.

The new correction-term computations give **`u(12n_491) = u(13n_3370) = 2`**.
Together with the Montesinos obstructions for `12n_288` and `12n_501`, they
identify **`13n_3370` as a counterexample to the Bernhard–Jablan conjecture**:
its strong Bernhard–Jablan unknotting number is three. The passage from the
three 12-crossing knots to this conclusion is
[Brittenham–Hermiller, Theorem 1.3](https://arxiv.org/abs/1705.05985v2).
Their published upper bound for `13n_3370` is combined with our lower bound;
it is not counted as a new upper-bound construction.

The same method now proves **`u(13n_1587) = 2`**. This confines the strong
counterexample in Brittenham–Hermiller's second further example to that knot
or one of its minimal-diagram neighbours with unknotting number one; it does
not identify which of those possibilities occurs. The complete BJ set of
`14a2539` has five prime members and `8_6 # 4_1`. Its unknotting number remains
in `[2,3]`; the composite possibility is retained in the conditional discussion.

## Relation to the paper

| Part of the paper | Mathematical method | Implementation |
| --- | --- | --- |
| Sections 2–3 | Linking pairings of double branched covers | [`lower_bounds/`](lower_bounds/README.md): Goeritz, Seifert, and Regina calculations |
| Sections 2–3 | Knot Floer torsion order | `lower_bounds/scan_torsion_lower_bounds.py` |
| Sections 2–3 | Owens's correction-term obstruction, Traczyk's criterion | [`lower_bounds/owens/`](lower_bounds/owens/README.md), for sequences of two, three, or four crossing changes |
| Sections 2–3 | Homology of cyclic branched covers | `lower_bounds/cyclic_cover_bound.py` |
| Sections 2–3 | Montesinos plumbings and reduced Heegaard Floer homology | [`lower_bounds/montesinos/`](lower_bounds/montesinos/README.md) |
| Sections 2, 3 and 5 | Greene's spanning-tree model, L-space correction terms, half-integral surgery and higher-rank obstructions | [`lower_bounds/greene/`](lower_bounds/greene/README.md) |
| Sections 3 and 5 | McCoy's theorem and minimal alternating diagrams | [`crossing_changes/`](crossing_changes/README.md) |
| Section 4 | Unknotting diagrams and crossing-change certificates | [`upper_bounds/`](upper_bounds/README.md) |
| Section 4 | Tait graphs, diagram moves, determinant and signature formulas | [`tait_graphs/`](tait_graphs/README.md) |
| Appendices A–H | Exact values, improved ranges, conditional results, PD codes, and correction terms | `consolidate_results.py`, deposited records in [`results/`](results/README.md) |

The torsion contributions in the result table use **knot Floer homology**. Other
torsion theories reviewed in Section 2 are not counted as additional computational
results. The rank-three and rank-four routines apply Owens's existing general
theorem. Their comments explain its signature and sign-pattern hypotheses.

For the **959 knots in Appendix F**, the conclusion is conditional: each has
`u = 3` unless it violates the Bernhard–Jablan property. These are not added to
the 2,719 exact values. The 1,027 alternating `[2,3]` parents split into these
959 cases, 22 parents with 11 named candidate children, and 46 further parents
whose exclusions would require unverified Jones-only identifications. The last
group consists of 39 entries in the heuristic list and seven further proposed
composites; none enters the proved dichotomy count.

## Inputs and provenance

[`results/paper_v1_1_snapshot.json`](results/paper_v1_1_snapshot.json) freezes a
**historical research baseline**, assembled on September 8, 2026, for the
12,965 knots in the study. It starts with the `database_knotinfo` 2026.8.1
package ranges and restores the April 2026 ranges for the 104 knots listed
in the deposited override file. Of these, 27 also have their lower bound
restored from two to one. It is **not an unmodified KnotInfo release**.
The [v1.3 input manifest](results/paper_v1_3_manifest.json) authenticates the
September 8, 2026 records used in the current consolidation. The
[historical manifest](results/paper_v1_1_manifest.json) preserves the earlier
profile; the
[release audit](results/comparison/baseline_comparison.json) explains and
reconstructs the comparison baseline.

The historical baseline differs from the public package snapshots 2026.6.1
and 2026.9.1 in 150 and 328 entries, respectively. Frozen extracts of the
April, June, August, and September releases include public wheel URLs,
SHA-256 hashes, package upload times, and original range fields. Reproduce
the comparison without installing or querying any live database:

```bash
python audit_knotinfo_releases.py --verify
```

Among the 2,719 historical-baseline exact determinations, the September
package has **194 agreeing exact values, nine conflicting exact values,
and 2,516 unresolved entries**. Thus 203 entries are already exact in that
package, but only 194 agree. These are database comparisons, not a claim
that all remaining values are previously unpublished. The 194 agreeing
exact values are precisely the values in Gebel–Prangley; their 30 improved
lower bounds also agree with the deposited computations. Including the
non-exact range for `13a_650`, there are ten disjoint intervals:

| Knot | Computed interval | Package 2026.9.1 |
| --- | ---: | ---: |
| `12a_107` | 4 | 3 |
| `13a_15` | 3 | 2 |
| `13a_55` | 3 | 2 |
| `13a_422` | 3 | 2 |
| `13a_568` | 3 | 2 |
| `13a_650` | [3,4] | 2 |
| `13a_660` | 4 | 3 |
| `13n_111` | 3 | 2 |
| `13n_142` | 3 | 2 |
| `13n_196` | 3 | 2 |

These disagreements are explicitly recorded in the paper and audit; the
pipeline never intersects incompatible intervals and calls the result an
improvement. Outside these ten conflicts, 94 retained upper bounds are
larger than those in the September package: 85 intervals are strictly
wider, and nine have a stronger lower bound but a weaker upper bound.
The historical comparison is retained for reproducibility, while the
release comparison makes these limitations visible.

The 27 McCoy improvements are improvements against this historical
baseline only: all 27 lower bounds were already two in the June package.
The larger deposited McCoy computation checks all 3,627 alternating knots
with April lower bound one. It finds an unknotting crossing for exactly the
710 knots recorded with `u = 1`, and none for the other 2,917 knots, whose
June lower bounds are two. This is an independent verification of those
recorded bounds, not 2,917 newly determined values. The separate crossing
dataset covers 3,105 of this cohort and agrees with every deposited
unknotting-crossing list in that overlap.

The comparison snapshot retains `[1,3]` for `13n_3370`. Its independently
published upper bound two is recorded in
[`brittenham_hermiller_upper_bounds.json`](results/bernhard_jablan/brittenham_hermiller_upper_bounds.json),
with the original DT code and the crossing change to `11n_21`. Consolidation
keeps this source distinct from the Greene lower bounds. Greene-derived
values use method tag `G`; the appendices mark them with superscript `G`.
Lower-case `g` continues to denote the homology generator bound.

The unchanged historical file
[`13n_1587_u_le_2_braid_certificate.json`](results/13n_1587_u_le_2_braid_certificate.json)
also supplies an upper bound two. **Erratum to its descriptive fields:** the
construction already appears in Brittenham–Hermiller, §3, rather than
originating in the later Applebaum paper. Its `[1,2]` interval and unresolved
status describe the upper-bound stage of the calculation. The new Greene
lower bound now makes the value exactly two. The original file is retained
unchanged to preserve historical checksums; its braid verification is valid.

The computational scripts use the `database_knotinfo==2026.8.1` snapshot for
knot diagrams and invariants. Reporting reads only the frozen deposited inputs;
it does not query a live database or silently merge unfinished worker logs.
See [`results/README.md`](results/README.md) for the individual files.

## Environment

Use Python 3.10 or later in a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The September 2026 audit ran on Python 3.13.7 with SnapPy 3.3.2,
Spherogram 2.4.1, Regina 7.4.1, SymPy 1.14.0, NumPy 2.5.2,
NetworkX 3.6.1, knot_floer_homology 1.2.2, and database_knotinfo 2026.8.1.
No Sage installation is needed for the verification commands below.
Matplotlib is needed only for drawing figures. Optional historical experiments
may require KnotJob (Java), scikit-learn, or a Slurm cluster; they are not needed
to reproduce the paper's aggregate data or verify the deposited unknotting diagrams.

## Regenerate the paper's data

Run from the repository root:

```bash
python consolidate_results.py --manuscript v1.3
```

This writes `generated/consolidated.json`, `generated/u_table.json`,
and `generated/counts.json`, including the enlarged Greene sweep and all 25
values with `u = 5`. It also regenerates the result and release-comparison macros,
and the rows of the lower-bound summary, the release-discrepancy table, and Appendices A–F
in `generated/paper_v1_3/`. The v1.3 manuscript uses these
files directly through `\input`; counts and table entries are not copied
into the manuscript. Layout, table legends, and mathematical explanations
remain in the manuscript. Neither a manuscript nor a deposited result file
is overwritten. Old files under `generated/appendix_tables/` are historical
outputs of a removed generator and are not used by this pipeline.

Version 1.3 is the default profile. To reproduce the preserved v1.2 totals of
1,227 exact values and 176 improved ranges in a separate output directory:

```bash
python consolidate_results.py --manuscript v1.2 --outdir /tmp/unknotting-v1.2
```

If the LaTeX manuscript is available separately, its actual entries and marked
PD codes and correction-term tables can also be compared with the deposited data:

```bash
python verify_manuscript.py --tex /path/to/main_v1.3.tex
```

The verifier resolves local `\input` files, rejects stale generated inputs,
and compares all appendix entries, marked PD codes, and correction-term
values with the deposited records. It also accepts the standalone v1.1 source.
To compile v1.3, keep the repository as the adjacent `code/` directory or
adjust the manuscript input paths after regenerating the files.

The consolidation pipeline needs only the Python standard library. It rejects a
changed source checksum, an unknown knot name, conflicting completed verdicts,
or a contradictory interval. This checks the assembly of the aggregate data; it
does **not** rerun the mathematical obstruction behind each deposited verdict.

## Verify the implementation and certificates

```bash
python -m unittest discover -s tests -v
python upper_bounds/selftest.py
python upper_bounds/verify_presentation_diagrams.py
python upper_bounds/verify_braid_certificate.py
python upper_bounds/verify_bj_upper_bounds.py
python audit_crossing_formulas.py
python audit_11a14_figure.py
python lower_bounds/owens/validate.py
python lower_bounds/owens/crosscheck.py 11a_63 13a_16 9_10 8_5
```

The presentation verifier requires all eight paper diagrams, authenticates their
knot identifications, checks two distinct marked crossings in each, and verifies
the unknot endpoints. The braid verifier checks the intermediate knot and its
unknotting crossing. The Bernhard–Jablan upper-bound verifier replays the
published constructions for `12n_491` and `13n_3370` and independently checks
the partners' unknotting crossings. The new lower-bound computations can be
repeated with `python lower_bounds/greene/pin_dinv.py 12n_491 13n_3370`; they
read the frozen diagram and mod-two Khovanov inputs and print their results.
The same command accepts `13n_1587`, reading its installed database row;
the result file and `results/review_v1_3/13n1587_replay.json` preserve the
actual diagram and Khovanov inputs used in its deposited and replayed calculations.
For a search campaign with replayable certificates:

```bash
python upper_bounds/verify_certificates.py results/witness_chains --paper --out generated/verified
```

Campaign verification and verification of the final presentation PD codes are
separate tasks. Some identifications use SnapPy's numerical isometry routines;
these should not be described as interval-certified proofs of homeomorphism.
`audit_crossing_formulas.py` compares all 69,632 stored crossing changes
from 5,546 parent knots with the determinant, signature, dual-resistance,
and self-linking formulas in exact rational arithmetic. The deposited
[audit summary](results/comparison/crossing_formula_verification.json)
records zero mismatches and hashes both input files. This is a comparison
of stored records; it does not recompute their invariants or certify their
knot identifications. The command prints its result and writes a file only
when `--out PATH` is supplied.

`audit_11a14_figure.py` freshly replays three exterior comparisons for the
published `11a_14` figure using the PD codes in
[`11a14_figure_verification.json`](results/comparison/11a14_figure_verification.json).
It compares the left diagram with `11a_14`, its boxed crossing change with
`8_8`, and the right diagram with `10_129`. These are SnapPy floating-point
canonical isometries, allowing mirrors, with at most 12 attempts per pair;
they are not interval-certified identifications. Exact integer Jones state
sums on the small 8- and 10-crossing reference diagrams separately verify
`V_8_8(t^-1) = V_10_129(t)`. The unknotting numbers two and one for those
references are explicitly dated tabulated inputs. The default command
reads the deposited inputs and writes nothing; `--write` deliberately
refreshes the record. No live KnotInfo access is needed. This diagnostic
concerns the displayed crossing change, not every upper bound in the
104-entry baseline override list.

The regression tests exercise mathematical kernels and known examples; they do
not rerun every deposited obstruction or identify every Appendix F crossing neighbor.

The [v1.3 review records](results/review_v1_3/) distinguish full arithmetic
replays, bounded diagram controls, and stored-result checks. The new four
rank-four successes were inspected in the deposited data, but their complete
form enumerations were **not freshly completed during this review**: a
bounded replay reached its 60-second limit. That limit supplies no verdict
and must not be described as a successful independent rerun. See
[the code audit](results/review_v1_3/greene_code_notes.md) and
[the updated dichotomy audit](results/review_v1_3/updated_dichotomy_notes.md).

Full searches and sweeps are described in the directory READMEs. They can take
hours or longer. Do not restart a long search merely to reproduce the aggregate counts.
Save new runs separately and review their completed records before deliberately
updating the frozen inputs and manuscript counts.

## Citation

```bibtex
@misc{Lee2026,
  author = {Lee, Seong-Jin},
  title = {Computation of Unknotting Numbers: Which Knot Breaks the Bernhard–Jablan Conjecture?},
  year = {2026},
  note = {Manuscript, version 1.3; accompanying code and computational data},
  url = {https://github.com/Pinocchio315/unknotting-tait-graphs}
}
```
