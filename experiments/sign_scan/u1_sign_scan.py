"""u = 1 obstruction from the SIGN of the unknotting crossing (non-alternating rows with range [1,2]).

For a knot K with u(K) = 1 the sign of the unknotting crossing is constrained three times:
  * signature (Cochran–Lickorish): sigma(K) = -2 forces a positive crossing, sigma(K) = +2 a negative one
    (KnotInfo convention; validated below on every knot with known u = 1);
  * Lickorish's linking form <s*2/D> and the Casson–Walker integrality test each allow a set of surgery
    signs s (tait.invariants.lickorish_allowed / casson_walker_allowed, as in scan_lower_bounds.py);
  * Traczyk (skein relation at e^{i pi/3}, see owens/traczyk.py): when d = dim_3 H_1(Sigma_2) = 1 the
    unknotting crossing is positive iff eps(K) = +1, where V_K(e^{i pi/3}) = eps (i sqrt3)^d.
The dictionary between crossing signs and surgery signs is fixed empirically on the u = 1 knots and the
whole test must produce no obstruction on them; then a target knot whose allowed sign sets have empty
intersection has u >= 2.

    python u1_sign_scan.py validate out.jsonl     # all KnotInfo knots with u = 1
    python u1_sign_scan.py targets  out.jsonl     # non-alternating rows [1,2]
"""
from __future__ import annotations
import json, os, re, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))
sys.path.insert(0, os.path.join(HERE, 'owens'))


def worker(name):
    import tait  # noqa
    from spherogram import Link
    from tait import knotinfo, graph as tg, invariants as inv
    from traczyk import eps_and_dim3
    r = knotinfo.row(name)
    rec = {'name': name, 'u': str(r['unknotting_number']).strip(), 'sigma': int(r['signature']),
           'det': int(r['determinant']), 'alternating': r['alternating']}
    try:
        g = tg.from_link(Link(knotinfo.pd_code(name)), name)
        SL = inv.lickorish_allowed(g)
        J = inv.jones_regina(g)
        SC = inv.casson_walker_allowed(J, int(r['signature']))
        eps, d3 = eps_and_dim3(J)
        rec.update({'lick': sorted(SL) if SL is not None else None, 'cw': sorted(SC), 'eps': eps, 'd3': d3})
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {e}'
    return rec


def main():
    mode, out = sys.argv[1], sys.argv[2]
    from tait import knotinfo
    import database_knotinfo as dk
    rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()]
    def rng(u):
        m = re.findall(r'\d+', str(u)); return (int(m[0]), int(m[-1])) if m else None
    if mode == 'validate':
        names = [r['name'] for r in rows if rng(r['unknotting_number']) == (1, 1)]
    else:
        names = [r['name'] for r in rows if rng(r['unknotting_number']) and rng(r['unknotting_number'])[0] == 1
                 and rng(r['unknotting_number'])[1] > 1 and not str(r.get('alternating', '')).upper().startswith('Y')]
    done = set()
    try:
        for line in open(out): done.add(json.loads(line)['name'])
    except FileNotFoundError:
        pass
    names = [n for n in names if n not in done]
    print(f'{mode}: {len(names)} knots to do ({len(done)} done)', flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=8) as ex, open(out, 'a') as h:
        futs = [ex.submit(worker, n) for n in names]
        for i, f in enumerate(as_completed(futs), 1):
            h.write(json.dumps(f.result()) + '\n'); h.flush()
            if i % 250 == 0:
                print(f'  {i}/{len(names)} {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
