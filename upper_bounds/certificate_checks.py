"""Structural checks shared by upper-bound certificate verifiers.

An upper bound counts distinct marked crossings of a specified diagram. Neither
the claimed count nor a partner's unknotting number is authenticated by replaying
isotopy moves, so these data must be checked separately.
"""

def validate_marks(rows, n_crossings, claimed_count=None):
    """Validate zero-based crossing rows; reject duplicate, negative and noninteger rows."""
    if not isinstance(rows, (list, tuple)) or any(type(r) is not int for r in rows):
        raise ValueError('marked rows must be a list of integer indices')
    if len(set(rows)) != len(rows) or any(r < 0 or r >= n_crossings for r in rows):
        raise ValueError('marked rows must be distinct and within the diagram')
    if claimed_count is not None and (type(claimed_count) is not int or claimed_count != len(rows)):
        raise ValueError('claimed crossing count does not equal the number of marked rows')
    return tuple(rows)


def trusted_upper_bound(value):
    """Read an integer upper bound or an ordered [lower, upper] interval."""
    if isinstance(value, (list, tuple)):
        if len(value) != 2 or any(type(x) is not int for x in value) or not 0 <= value[0] <= value[1]:
            raise ValueError('expected a nonnegative ordered interval')
        return value[1]
    if type(value) is not int or value < 0:
        raise ValueError('expected a nonnegative integer upper bound')
    return value


def check_partner_bound(claimed, independently_known):
    """A certificate's partner_u may weaken, but cannot improve, its trusted bound."""
    known = trusted_upper_bound(independently_known)
    if claimed is not None and trusted_upper_bound(claimed) < known:
        raise ValueError('partner_u is smaller than the independently supplied upper bound')
    return known
