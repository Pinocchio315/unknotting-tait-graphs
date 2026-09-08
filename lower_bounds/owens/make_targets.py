#!/usr/bin/env python
"""List the applicable targets (alternating, [2,3], |sigma| = 4) by family, for information."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import kinfo
from collections import Counter
c = Counter()
for nm, r in kinfo.rows().items():
    if (r['alternating'] == 'Y' and r['unknotting_number'].replace(' ', '') == '[2,3]'
            and abs(int(r['signature'])) == 4):
        c[nm.split('_')[0]] += 1
print(dict(sorted(c.items())))
