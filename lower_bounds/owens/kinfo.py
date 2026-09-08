"""Minimal KnotInfo access (database_knotinfo) for the owens folder — no external project code."""
import ast
import functools

@functools.lru_cache(maxsize=1)
def rows():
    from database_knotinfo import link_list
    out = {}
    for r in link_list():
        nm = r.get('name')
        if not nm or not r.get('crossing_number', '').strip().isdigit():
            continue
        out[nm] = r
    return out

def row(name):
    return rows()[name]

def pd_code(name):
    return [list(q) for q in ast.literal_eval(row(name)['pd_notation'].replace(';', ','))]

def parse_jones(name):
    """KnotInfo jones_polynomial string -> {exponent: int coefficient} (handles t^(-3), fractions a/t^k)."""
    import re
    s = row(name)['jones_polynomial'].replace(' ', '').replace('(', '').replace(')', '')
    s = s.replace('^-', '^~').replace('/-', '/~')
    s = s.replace('-', '+-').replace('~', '-')
    out = {}
    for term in s.split('+'):
        if not term:
            continue
        m = re.fullmatch(r'(-?\d*)\*?(?:t(?:\^(-?\d+))?)?', term)
        if m:
            c, e = m.group(1), m.group(2)
            coef = int(c) if c not in ('', '-') else (-1 if c == '-' else 1)
            exp = int(e) if e is not None else (1 if 't' in term else 0)
            out[exp] = out.get(exp, 0) + coef
            continue
        m = re.fullmatch(r'(-?\d*)/t\^?(-?\d+)?', term)
        if m:
            c, e = m.group(1), m.group(2)
            coef = int(c) if c not in ('', '-') else (-1 if c == '-' else 1)
            exp = -(int(e) if e is not None else 1)
            out[exp] = out.get(exp, 0) + coef
            continue
        raise ValueError(f'cannot parse jones term {term!r} of {name}')
    return {e: c for e, c in out.items() if c}
