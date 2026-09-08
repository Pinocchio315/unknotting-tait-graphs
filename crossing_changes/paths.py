"""Repository-relative locations used by the crossing-change scripts (paper Section 4.4)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, 'results')
CC = os.path.join(RESULTS, 'crossing_changes')          # dataset of crossing changes (5,542 alternating knots)
OPEN23 = os.path.join(RESULTS, 'open23')                # the 833 open alternating [2,3] knots
WORK = os.environ.get('CC_WORK', os.path.join(CC, 'work'))   # per-shard intermediate files (not deposited)
U_TABLE = os.path.join(RESULTS, 'u_table_2026-09-07.json')   # name -> [lo, hi] (reference table + our bounds)
DKT_104 = os.path.join(RESULTS, 'dkt', 'dkt_104.json')
ALL_CODES = os.path.join(ROOT, 'upper_bounds', 'data', 'all_knot_codes.json')
RANKING = os.path.join(ROOT, 'experiments', 'ml', 'u23_ranking_2026-09-07.json')
for p in (os.path.join(ROOT, 'upper_bounds', 'cluster'), os.path.join(ROOT, 'upper_bounds'), HERE):
    if p not in sys.path: sys.path.insert(0, p)
os.makedirs(WORK, exist_ok=True)
