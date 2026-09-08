# How to Compute the Unknotting Number: Theory and Computations

Companion code and deposited computational records for the manuscript by
Seong-Jin Lee and Pavel Putrov, version 1.1 (September 2026).

The paper explains how topological obstructions and explicit crossing changes
work together. An obstruction raises a lower bound for `u(K)`; an unknotting
diagram gives an upper bound. An exact value follows when the two bounds agree.
The code follows this distinction: an unsuccessful search supplies no lower
bound, and a correction-term verdict `PASS` supplies no upper bound.

Relative to the archived comparison ranges used in the manuscript, the deposited
computations give **1,225 exact values** and **176 further improved ranges** for
prime knots with at most 13 crossings:

| Unknotting number | Exact values |
| --- | ---: |
| 2 | 838 |
| 3 | 237 |
| 4 | 129 |
| 5 | 21 |

Lower bounds account for 1,217 exact values and 175 further improvements. Eight
unknotting diagrams account for the remaining exact values. A braid certificate
for `13n_1587`, following Brittenham–Hermiller, improves one further upper bound.
These counts include independently recovered results, as explained in the paper;
they are not a claim that every listed value first appears here.

## Relation to the paper

| Part of the paper | Mathematical method | Implementation |
| --- | --- | --- |
| Sections 2–3 | Linking pairings of double branched covers | [`lower_bounds/`](lower_bounds/README.md): Goeritz, Seifert, and Regina calculations |
| Sections 2–3 | Knot Floer torsion order | `lower_bounds/scan_torsion_lower_bounds.py` |
| Sections 2–3 | Owens's correction-term obstruction, Traczyk's criterion | [`lower_bounds/owens/`](lower_bounds/owens/README.md), for sequences of two, three, or four crossing changes |
| Sections 2–3 | Homology of cyclic branched covers | `lower_bounds/cyclic_cover_bound.py` |
| Sections 2–3 | Montesinos plumbings and reduced Heegaard Floer homology | [`lower_bounds/montesinos/`](lower_bounds/montesinos/README.md) |
| Sections 3–4 | McCoy's theorem and minimal alternating diagrams | [`crossing_changes/`](crossing_changes/README.md) |
| Section 4 | Unknotting diagrams and crossing-change certificates | [`upper_bounds/`](upper_bounds/README.md) |
| Section 4 | Tait graphs, diagram moves, determinant and signature formulas | [`tait_graphs/`](tait_graphs/README.md) |
| Appendices A–G | Exact values, improved ranges, conditional results, PD codes | `consolidate_results.py`, deposited records in [`results/`](results/README.md) |

The torsion contributions in the result table use **knot Floer homology**. Other
torsion theories reviewed in Section 2 are not counted as additional computational
results. The rank-three and rank-four routines apply Owens's existing general
theorem. Their comments explain its signature and sign-pattern hypotheses.

For the **806 knots in Appendix F**, the conclusion is conditional: each has
`u = 3` unless it violates the Bernhard–Jablan property. These are not added to
the 1,225 exact values. The deposited rank-four computation contains completed
records for 21 of 24 targets; an unfinished target is not treated as obstructed.

## Inputs and provenance

[`results/paper_v1_1_snapshot.json`](results/paper_v1_1_snapshot.json) freezes the
12,965 comparison ranges used by the manuscript. It is an archived collection
of input ranges, **not an unmodified export of a currently installed database**.
The [manifest](results/paper_v1_1_manifest.json) records the source and SHA-256
checksum of every input used in consolidation.

The 27 McCoy records explicitly store their original `[1,b]` input ranges.
The frozen comparison retains those ranges, so the corresponding improvements
to `[2,b]` are counted once, as in Appendix E. An older deposited aggregate had
already incorporated these lower bounds and therefore counted only 149 further
improvements. Historical aggregates remain intact; use the v1.1 pipeline below
for the current paper.

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
python consolidate_results.py
```

This writes `generated/consolidated.json`, `generated/u_table.json`,
and `generated/counts.json`, including all 21 values with `u = 5`.
Neither manuscript version nor a deposited result file is overwritten.

If the LaTeX manuscript is available separately, its actual entries and marked
PD codes can also be compared with the deposited data:

```bash
python verify_manuscript.py --tex /path/to/main_v1.1.tex
```

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
python lower_bounds/owens/validate.py
python lower_bounds/owens/crosscheck.py 11a_63 13a_16 9_10 8_5
```

The presentation verifier requires all eight paper diagrams, authenticates their
knot identifications, checks two distinct marked crossings in each, and verifies
the unknot endpoints. The braid verifier checks the intermediate knot and its
unknotting crossing. For a search campaign with replayable certificates:

```bash
python upper_bounds/verify_certificates.py results/witness_chains --paper --out generated/verified
```

Campaign verification and verification of the final presentation PD codes are
separate tasks. Some identifications use SnapPy's numerical isometry routines;
these should not be described as interval-certified proofs of homeomorphism.
The regression tests exercise mathematical kernels and known examples; they do
not rerun every deposited obstruction or identify every Appendix F crossing neighbor.

Full searches and sweeps are described in the directory READMEs. They can take
hours or longer. Do not restart a long search merely to reproduce the aggregate counts.
Save new runs separately and review their completed records before deliberately
updating the frozen inputs and manuscript counts.

## Citation

```bibtex
@misc{LeePutrov2026,
  author = {Lee, Seong-Jin and Putrov, Pavel},
  title = {How to Compute the Unknotting Number: Theory and Computations},
  year = {2026},
  note = {Manuscript, version 1.1; accompanying code and computational data},
  url = {https://github.com/Pinocchio315/unknotting-tait-graphs}
}
```
