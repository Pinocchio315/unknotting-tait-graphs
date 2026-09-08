#!/usr/bin/env python3
"""Replay upper bounds for the three non-alternating conflicts added in v1.3.

This is an upper-bound audit, not a search-based lower bound.  It tests subsets
of the tabulated minimal PD diagram.  A witness is accepted only if elementary
Reidemeister simplification removes every crossing; determinant one alone
is never accepted as recognition of the unknot.  The complete input PD and
changed crossing rows are saved so the upper bound can be replayed independently.

Run with the project's unknot-venv Python.  --live-xls optionally authenticates
each input PD against the official Excel download used in the companion audit.
xlrd must then be importable (a wheel on PYTHONPATH is sufficient).
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tait_graphs"))
sys.path.insert(0, str(ROOT / "upper_bounds"))

import tait  # noqa: E402,F401; package initializes topology dependencies
from tait import knotinfo  # noqa: E402
from spherogram import Link  # noqa: E402
import regina  # noqa: E402
from xtait.graph import from_pd  # noqa: E402
from flipdet import exact_det  # noqa: E402

BOUNDS = {"13n_111": 3, "13n_142": 3, "13n_196": 3}


def flipped(pd, rows):
    """A cyclic rotation exchanges over/under at that crossing only."""
    out = [list(q) for q in pd]
    for row in rows:
        out[row] = out[row][1:] + out[row][:1]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live-xls", type=Path)
    parser.add_argument("--seconds", type=float, default=120)
    args = parser.parse_args()
    live_pd = {}
    if args.live_xls:
        import xlrd
        sheet = xlrd.open_workbook(str(args.live_xls), on_demand=True).sheet_by_index(0)
        headers = sheet.row_values(0)
        j = headers.index("pd_notation")
        for i in range(2, sheet.nrows):
            name = sheet.cell_value(i, 0)
            if name in BOUNDS:
                live_pd[name] = ast.literal_eval(sheet.cell_value(i, j))
    started = time.monotonic()
    report = {"schema_version": 1, "scope": "independent upper-bound witnesses only",
              "recognition": "Spherogram Reidemeister I/II and, if needed, III simplification to zero crossings",
              "type_III_limit": 100, "random_seed": 0,
              "warning": "Failure of this bounded diagram search is not a lower-bound proof.",
              "regina_version": regina.versionString(),
              "knots": {}}
    if args.live_xls:
        report["official_xls_sha256"] = hashlib.sha256(args.live_xls.read_bytes()).hexdigest()
    for name, bound in BOUNDS.items():
        pd = knotinfo.pd_code(name)
        if args.live_xls:
            assert pd == live_pd[name], (name, "PD differs from live official source")
        assert len(Link(pd).link_components) == 1
        record = {"claimed_upper_bound": bound, "source_pd": pd,
                  "official_live_pd_exactly_equal": bool(args.live_xls),
                  "subsets_tested": 0, "determinant_one_candidates": 0,
                  "verified": False}
        for rows in itertools.combinations(range(len(pd)), bound):
            if time.monotonic() - started > args.seconds:
                record["stopped_at_time_limit"] = True
                break
            record["subsets_tested"] += 1
            changed = flipped(pd, rows)
            if exact_det(from_pd(changed)) != 1:
                continue
            record["determinant_one_candidates"] += 1
            link = Link(changed)
            assert len(link.link_components) == 1
            link.simplify("basic")
            mode = "basic"
            if link.crossings:
                random.seed(0)
                link.simplify("level", type_III_limit=100)
                mode = "level"
            if not link.crossings:
                # Replay in a second topology implementation.  Regina retains
                # the crossing-free component, so also require exactly one.
                independent = regina.Link.fromPD(
                    Link(changed).PD_code(min_strand_index=1))
                independent.simplify()
                assert independent.size() == 0
                assert independent.countComponents() == 1
                record.update(verified=True, changed_rows_zero_based=list(rows),
                              changed_pd=changed, remaining_crossings=0,
                              simplifying_mode=mode,
                              independent_regina_remaining_crossings=0,
                              independent_regina_components=1)
                break
        report["knots"][name] = record
        print(name, record["verified"], record["subsets_tested"],
              record.get("changed_rows_zero_based"), flush=True)
    report["elapsed_seconds"] = round(time.monotonic() - started, 3)
    report["all_three_verified"] = all(r["verified"] for r in report["knots"].values())
    report["independent_regina_all_three_verified"] = report["all_three_verified"]
    path = Path(__file__).with_name("upper_bound_recheck.json")
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
