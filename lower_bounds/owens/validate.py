#!/usr/bin/env python
"""Validation of the Owens u=2 obstruction.
Positive controls (must OBSTRUCT; Owens Cor. 2, all with |sigma| = 4, u = 3): 9_10, 9_13, 9_38, 10_53, 10_101, 10_120.
Negative controls (must PASS): alternating knots with KnotInfo exact u = 2 and |sigma| = 4."""
import sys, time
sys.path.insert(0, '.')
import _compat  # noqa
import kinfo
from owens_obstruction import obstruct_u2

def run(nm, expect):
    r = kinfo.row(nm)
    sigma = int(r['signature']); det = int(r['determinant'])
    t0 = time.time()
    try:
        res = obstruct_u2(kinfo.pd_code(nm), sigma, det)
    except Exception as e:
        res = {'verdict': 'ERROR', 'reason': repr(e)}
    ok = res['verdict'] == expect
    print(f"{nm:8s} sigma {sigma:+d} det {det:4d} u={r['unknotting_number']:6s} -> {res['verdict']:10s} "
          f"(expected {expect}) {'OK' if ok else '** FAIL **'}  {time.time()-t0:.1f}s  {res.get('reason','')[:80]}")
    return ok

def main():
    pos = ['9_10', '9_13', '9_38', '10_53', '10_101', '10_120']
    neg = []
    for nm, r in kinfo.rows().items():
        if (r['alternating'] == 'Y' and r['unknotting_number'].strip() == '2' and abs(int(r['signature'])) == 4
                and int(r['crossing_number']) <= 10):
            neg.append(nm)
    neg = neg[:8]
    print('negative controls:', neg)
    good = sum(run(nm, 'OBSTRUCTED') for nm in pos) + sum(run(nm, 'PASS') for nm in neg)
    print(f'{good}/{len(pos)+len(neg)} controls correct')
    return 0 if good == len(pos) + len(neg) else 1


if __name__ == '__main__':
    sys.exit(main())
