# Review of the additional crossing-neighbour update

The user explicitly asked to include the concurrent update committed as
`a511ba1`, after the first v1.3 review had started. This note extends the
earlier frozen review to the latest 959-knot dichotomy.

## Outcome

The 959-knot conditional theorem is supported by the updated classification.
An independent classifier, which never reads machine-learning scores,
reconstructs exactly the claimed set. Every source knot still has comparison
range [2,3]. Every crossing is present exactly once, with no missing row.
The theorem remains conditional: each knot has u=3, or has u=2 and fails the
strong Bernhard–Jablan equality. None of these 959 is counted as an unconditional
exact value solely because no candidate crossing survives.

The new set contains all 806 previous knots and 153 additional knots. Across
all 13,105 crossing rows for 1,027 parents, only two fields changed relative to
the starting commit `6468494`: 238 neighbour intervals, and 20 harmless
reorderings of connected-sum factors. All 238 intervals agree exactly with
the new Greene summary's `new_range` and `lower_bound` fields.

For the 153 newly included parents, 214 previously possible crossings are
now excluded by independently calculated neighbour lower bounds. There are
23 exact canonical-diagram identifications and 191 existing SnapPy isometry
identifications filtered by determinant/Jones polynomial. The latter are
geometric identifications, not Jones-polynomial-only identifications. Their
numerical rather than interval-certified nature is inherited from the
existing crossing-neighbour dataset. This update does not replace those
identifications or strengthen them silently.

Eight deterministic controls with distinct newly settled children were
replayed directly from KnotInfo PD codes. All eight crossing changes again
matched the stated child by SnapPy exterior isometry; their determinants
and absolute signatures were recomputed with exact arithmetic and agreed
with the stored rows. This is a bounded sample of the inherited geometric
identifications, not a fresh exhaustive geometric verification of all rows.

## The seven moved parents remain heuristic

The latest table has 29 parents in its candidate list. Eleven named child
knots cover 22 of them. For the other seven parents, the separate resolution
file records only `composite(jones)` guesses:

- 13a1154
- 13a1158
- 13a1184
- 13a120
- 13a2098
- 13a320
- 13a473

`resolve_unres.py` describes them as moved after resolution, but a product
Jones polynomial does not prove a connected-sum decomposition. They must
remain outside the 959 proved dichotomies. The reporting code currently does
exclude them from that count. There are also 39 Jones-dependent zero-candidate
parents in the distinct heuristic list. Thus the complete partition is

    1027 = 959 certified dichotomies + 39 heuristic zero-candidate parents
           + 22 parents with named candidates + 7 unresolved Jones-only parents.

It would be misleading to say that the 11 named children describe all
remaining possibilities, or that those seven parents have been proved to
have no candidate crossing. They are unresolved in the rigorous sense.

## Code and mathematical logic

`candidate_status` separates determinant/Jones filtering followed by actual
isometry from a Jones-only guess. It gives no proof status to an unidentified
knot or to an unknown identification method. A signature of absolute value
greater than two is an independent sufficient exclusion of unknotting number
one. Machine-learning scores affect ordering only.

For identified connected sums, the lower bound two requires at least two
nontrivial factors, as in Scharlemann's theorem. The helper does not assume
unknotting-number additivity. The rigorous set contains no crossing whose
exclusion depends solely on a `sum-jones` match.

For a prime alternating parent, the crossing-neighbour set is independent
of the chosen minimal diagram by the flyping lemma already stated in the
paper. Accordingly, once every crossing in one minimal diagram has a sound
exclusion, the conditional theorem follows. No heuristic search failure or
machine-learning classification is being substituted for this implication.

The change in `crossing_changes/paths.py` switches the default lower/upper
table from the older fixed file to `generated/u_table.json`, with an explicit
`CC_U_TABLE` override. This is appropriate for propagating the new Greene
bounds. Reproduction must generate that table first, or point `CC_U_TABLE`
to the intended immutable input; otherwise a missing or stale generated
file can fail or reproduce the wrong historical classification. The audit
records the exact crossing-row hash and checks the current intervals against
the separately audited Greene summary.

## Reproduction artifacts

- `audit_updated_dichotomy.py`: independent classification, coverage checks,
  old/new interval provenance and the complete list of newly decisive
  parent–child dependencies.
- `updated_dichotomy_audit.json`: successful audit, with 959 reconstructed
  rigorous parents, 153 additions and all changed intervals.
- `replay_new_dichotomy_controls.py` and `new_dichotomy_controls.json`:
  eight bounded diagram/invariant replays, all successful.
- `dichotomy_row_delta.json`: detailed raw-row comparison.

No production code or manuscript was modified in this subtask.
