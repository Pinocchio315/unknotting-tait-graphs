"""Proof-sensitive crossing-neighbor decisions for the Bernhard–Jablan discussion.

These are necessary tests for a neighbor with u=1. A parent with u=2 must
have such a *minimal-diagram* neighbor only if it also satisfies the
Bernhard–Jablan property. Scores and polynomial identifications are not proofs.
"""

VERIFIED_IDENTIFICATIONS = {'code', 'trivial', 'isometry', 'isometry-jones',
                            'prime-isometry', 'sum-blocks'}


def candidate_status(name, how, interval, signature):
    if signature is not None and abs(signature) > 2:
        return 'excluded:sigma'
    if name is None:
        return 'candidate:unidentified'
    if how == 'sum-ambiguous':
        return 'candidate:composite?(jones-ambiguous)'
    if how == 'sum-jones':
        # Preserve a distinct heuristic tier for historical comparisons, never
        # a certified witness or an unqualified exclusion.
        return 'excluded:u>=2(jones)' if interval and interval[0] >= 2 else 'candidate:unverified'
    if how not in VERIFIED_IDENTIFICATIONS:
        return 'candidate:unverified'
    if interval is None:
        return 'candidate:u-unknown-factor'
    if interval[1] <= 1:
        return 'WITNESS'
    if interval[0] >= 2:
        return 'excluded:u>=2'
    return 'candidate'


def connected_sum_interval(parts, known_ranges, signature=None):
    """Bounds for an identified connected sum, without assuming additivity.

    Scharlemann gives lower bound two only when at least TWO nontrivial
    summands remain. One nontrivial block and any unknot blocks give that
    block itself. Unknown partner bounds remain unknown.
    """
    parts = [name for name in parts if name != '0_1']
    if not parts:
        return [0, 0]
    if not all(name in known_ranges for name in parts):
        return None
    if len(parts) == 1:
        lower, upper = known_ranges[parts[0]]
    else:
        lower, upper = 2, sum(known_ranges[name][1] for name in parts)
    if signature is not None:
        lower = max(lower, abs(signature) // 2)
    if lower > upper:
        raise ValueError('identified-sum bounds are inconsistent')
    return [lower, upper]
