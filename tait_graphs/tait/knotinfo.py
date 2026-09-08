"""KnotInfo access helpers (database_knotinfo package, '|'-delimited CSV) and Jones-string parsing.

KnotInfo names look like '3_1', '10_132', '11a_14', '12n_491', '13n_3370'.
"""
from __future__ import annotations

import ast
import csv
import os
import re
from functools import lru_cache


@lru_cache(maxsize=1)
def rows() -> dict[str, dict]:
    """name -> row dict (all 244 KnotInfo columns as strings)."""
    import database_knotinfo
    root = os.path.dirname(database_knotinfo.__file__)
    path = next(os.path.join(dp, f) for dp, dn, fn in os.walk(root) for f in fn
                if f == 'knotinfo_data_complete.csv')
    with open(path, newline='', encoding='utf-8') as h:
        out = {}
        for r in csv.DictReader(h, delimiter='|'):
            if not r.get('crossing_number', '').strip().isdigit():
                continue  # the file contains a second header-like row ("Crossing Number", "UnknottingNumber", ...)
            out[r['name']] = r
        return out


def row(name: str) -> dict:
    return rows()[normalize_name(name)]


def normalize_name(name: str) -> str:
    """Accept '12n491', 'K12n491', '12n_491', '10_132', 'K10a85'(HT names are NOT converted) -> KnotInfo name."""
    s = str(name).strip()
    m = re.fullmatch(r'K?(\d+)([an])_?(\d+)', s)
    if m:
        return f'{m.group(1)}{m.group(2)}_{m.group(3)}'
    return s


def pd_code(name: str) -> list[list[int]]:
    return ast.literal_eval(row(name)['pd_notation'])


def unknotting_interval(name: str) -> tuple[int, int]:
    """'2' -> (2,2); '[1,3]' -> (1,3)."""
    u = row(name)['unknotting_number'].replace(' ', '')
    m = re.fullmatch(r'\[(\d+),(\d+)\]', u)
    if m:
        return int(m.group(1)), int(m.group(2))
    if u.isdigit():
        return int(u), int(u)
    return (None, None)


def gapped_names() -> list[str]:
    """Knots whose unknotting number is not yet determined (interval with lower < upper)."""
    out = []
    for n, r in rows().items():
        u = r['unknotting_number'].replace(' ', '')
        m = re.fullmatch(r'\[(\d+),(\d+)\]', u)
        if m and int(m.group(1)) < int(m.group(2)):
            out.append(n)
    return out


def _parse_signed_int(value: str, default: int = 1) -> int:
    if value in ('', '+'):
        return default
    if value == '-':
        return -default
    return int(value)


def parse_jones_string(poly: str) -> dict[int, int]:
    """KnotInfo Jones polynomial string -> {exponent: coefficient}.
    Handles both 't^(-3)+ 2*t^(-2)' (<=12 crossings) and '2 - t^(-9) + 2/t^8 - 2/t' (13 crossings) formats."""
    s = str(poly).replace(' ', '')
    if not s:
        return {}
    terms, current, depth = [], [], 0
    for i, ch in enumerate(s):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
        if i > 0 and depth == 0 and ch in '+-':
            terms.append(''.join(current))
            current = [ch]
        else:
            current.append(ch)
    if current:
        terms.append(''.join(current))
    out: dict[int, int] = {}
    for term in terms:
        if not term:
            continue
        if 't' not in term:
            e, c = 0, int(term)
        elif '/t' in term:
            cp, tp = term.split('/t', 1)
            c = _parse_signed_int(cp, 1)
            if tp == '':
                e = -1
            else:
                ep = tp[1:] if tp.startswith('^') else tp
                if ep.startswith('(') and ep.endswith(')'):
                    ep = ep[1:-1]
                e = -int(ep)
        else:
            ti = term.index('t')
            cp = term[:ti].rstrip('*')
            c = _parse_signed_int(cp, 1)
            suffix = term[ti + 1:]
            if suffix == '':
                e = 1
            else:
                ep = suffix[1:] if suffix.startswith('^') else suffix
                if ep.startswith('(') and ep.endswith(')'):
                    ep = ep[1:-1]
                e = int(ep)
        out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}


def jones_dict(name: str) -> dict[int, int]:
    return parse_jones_string(row(name)['jones_polynomial'])
