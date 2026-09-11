"""Shared paths for the crossing-change analysis (v1.8 Sections 4.4 and 5.2).

The deposited directory contains successive historical cohorts. Current paper
membership comes from the v1.3 computational profile in consolidate_results.py;
a path to this directory alone does not select the v1.8 analysis.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, 'results')
CC = os.path.join(RESULTS, 'crossing_changes')          # current dataset: 5,546 alternating parent knots
OPEN23 = os.path.join(RESULTS, 'open23')                # current 1,027-parent analysis and historical cohorts
WORK = os.environ.get('CC_WORK', os.path.join(CC, 'work'))   # per-shard intermediate files (not deposited)
U_TABLE = os.environ.get('CC_U_TABLE',                       # name -> [lo, hi] (comparison table + our bounds);
                        os.path.join(ROOT, 'generated', 'u_table.json'))   # the consolidation writes this one
DKT_104 = os.path.join(RESULTS, 'dkt', 'dkt_104.json')
ALL_CODES = os.path.join(ROOT, 'upper_bounds', 'data', 'all_knot_codes.json')
for p in (os.path.join(ROOT, 'upper_bounds', 'cluster'), os.path.join(ROOT, 'upper_bounds'), HERE):
    if p not in sys.path: sys.path.insert(0, p)
os.makedirs(WORK, exist_ok=True)
