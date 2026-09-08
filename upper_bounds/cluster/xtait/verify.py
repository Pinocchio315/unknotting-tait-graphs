"""Independent certificate replay and verification (stdlib only).

A counterexample certificate for the crossing-additivity conjecture is a JSON object:

    {
      "pair":            {"k1": ..., "k2": ..., "mirror2": bool, "arc1": int, "arc2": int, "variant": int},
      "c1": int, "c2": int,                # table crossing numbers of the summands
      "start_pd":        [[a,b,c,d], ...], # the spliced connected-sum diagram (c1+c2 crossings)
      "path":            [[kind, params], ...],   # isotopy moves (R1±/R2±/R3/flype/pass)
      "final_pd":        [[...], ...],
      "final_crossings": int               # == len(final_pd) < c1 + c2   <-- the contradiction
    }

`verify_certificate` replays the path from start_pd with the deterministic move engine and checks that the
result equals final_pd as a diagram (canonical code) and has the claimed crossing number.  Optionally the
Jones polynomial of final_pd is recomputed with the Kauffman bracket and compared with the product of the
summands' Jones polynomials (multiplicativity under connected sum) — a strong independent sanity check.
"""
from __future__ import annotations

import json

from . import moves as mv
from .canonical import canonical_code
from .graph import EmbeddedTaitGraph, from_pd, to_pd
from .jones import jones, jones_mirror, jones_mul


def _params_tuple(kind, params):
    if kind == 'pass':
        p = params[0] if isinstance(params, (list, tuple)) else params
        return ({'run_index': int(p['run_index']),
                 'route': [[int(a), int(b)] for a, b in p['route']]},)
    return tuple(params)


def replay(g0: EmbeddedTaitGraph, path):
    cur = g0
    for kind, params in path:
        if kind == 'cc':
            raise ValueError('crossing changes are not isotopy moves')
        cur = mv.apply_move(cur, kind, _params_tuple(kind, params))
    return cur


def verify_certificate(cert: dict, check_jones: bool = True, jones_of=None, log=print) -> bool:
    """jones_of: optional {name: {exp: coeff}} for the summands (else taken from cert['jones1/2'])."""
    start = [[int(x) for x in q] for q in cert['start_pd']]
    final = [[int(x) for x in q] for q in cert['final_pd']]
    c1, c2 = int(cert['c1']), int(cert['c2'])
    if len(start) != c1 + c2:
        log('verify: start diagram does not have c1+c2 crossings')
        return False
    if len(final) != int(cert['final_crossings']):
        log('verify: final_crossings mismatch')
        return False
    g0 = from_pd(start)
    g = replay(g0, cert['path'])
    if g.n_crossings() != len(final):
        log(f'verify: replay gives {g.n_crossings()} crossings, certificate claims {len(final)}')
        return False
    if canonical_code(g) != canonical_code(from_pd(final)):
        log('verify: replayed diagram differs from final_pd')
        return False
    if check_jones and len(final) <= 24:
        j1 = jones_of[cert['pair']['k1']] if jones_of else {int(e): c for e, c in cert['jones1'].items()}
        j2 = jones_of[cert['pair']['k2']] if jones_of else {int(e): c for e, c in cert['jones2'].items()}
        if cert['pair'].get('mirror2'):
            j2 = jones_mirror(j2)
        expected = jones_mul(j1, j2)
        got = jones(final)
        if got not in (expected, jones_mirror(expected)):
            log('verify: Jones polynomial of the final diagram is not V(K1)V(K2)')
            return False
    return True


def load_certificate(path: str) -> dict:
    with open(path) as h:
        return json.load(h)


def save_certificate(cert: dict, path: str) -> None:
    def enc(o):
        if isinstance(o, tuple):
            return list(o)
        raise TypeError(str(type(o)))
    with open(path, 'w') as h:
        json.dump(cert, h, indent=1, default=enc)
