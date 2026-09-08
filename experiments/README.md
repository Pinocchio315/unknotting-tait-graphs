# experiments — exploratory computations that produced no new value

These scripts are kept for completeness; none of their outcomes is used in the paper.

| directory | what was tried | outcome |
|---|---|---|
| `ml/` | data mining of the Jones polynomial and other invariants for new unknotting-number rules (`jones_u_mining.py`, `ext_features.py`, `stage_ab.py`), and a probability ranking of the open alternating `[2,3]` knots (`u23_ranking_2026-09-07.json`, used only to order computations) | every exact rule found is a known theorem (signature, Lickorish, Traczyk); no new obstruction |
| `donaldson/` | refinements of the correction-term obstruction to `u = 2` by embedding the definite Goeritz lattice of the other checkerboard surface into the standard lattice (`donaldson_u2.py`, with the gluing identification `donaldson_u2_glue.py`) | validated on 747 knots with known `u`, no obstruction on any open knot (`results_donaldson_glue_2026-09-07.json`) |
| `sign_scan/` | combining the Lickorish, Casson–Walker and Traczyk sign conditions for `u = 1` (`u1_sign_scan.py`) | the combined rule is an identity; 0 new obstructions (`u1_crossing_sign_obstructions_2026-09-07.json`) |

Two further experiments were run interactively and are described in the project notes only: a
statistical model of the outcome of the rank-2 obstruction from the determinant and the linking
form (the verdict is not a function of these data), and the crossing-selection model of
`../crossing_changes/analyze_v2.py`, whose result is quoted in the paper as a remark.
