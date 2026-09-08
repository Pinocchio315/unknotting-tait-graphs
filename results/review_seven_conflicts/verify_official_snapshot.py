#!/usr/bin/env python3
"""Reproduce the comparison with the official KnotInfo download of 2026-09-08.

The default check uses the deposited SHA-256-authenticated range extract.
Optionally pass --xls-zip to verify the full downloaded official archive and
repeat the extraction, with xlrd 2.0.2 importable.  Obtain the ZIP from the
immutable pinned_download_url in official_today_comparison.json.  No installed
database package is changed, and no network requests or knot searches are made.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xls-zip", type=Path)
    args = parser.parse_args()
    report = json.loads((HERE / "official_today_comparison.json").read_text())
    frozen_bytes = (HERE / report["frozen_file"]).read_bytes()
    assert sha256(frozen_bytes) == report["frozen_sha256"]
    frozen = json.loads(gzip.decompress(frozen_bytes))
    assert frozen["source"] == report["source"]
    current = frozen["rows"]
    baseline = json.loads((ROOT / "results/paper_v1_1_snapshot.json").read_text())
    assert set(current) == set(baseline)
    if args.xls_zip:
        import xlrd
        source = report["source"]
        archive = args.xls_zip.read_bytes()
        assert sha256(archive) == source["zip_sha256"]
        with zipfile.ZipFile(args.xls_zip) as zipped:
            raw = zipped.read(source["xls_member"])
        assert sha256(raw) == source["xls_sha256"]
        sheet = xlrd.open_workbook(file_contents=raw, on_demand=True).sheet_by_name(source["xls_sheet"])
        header = sheet.row_values(0)
        name_col, value_col, ref_col = (header.index(c) for c in source["extracted_fields"])
        extracted = {}
        for row in range(2, sheet.nrows):
            name = sheet.cell_value(row, name_col)
            if name not in baseline:
                continue
            value = sheet.cell_value(row, value_col)
            nums = [int(value)] if isinstance(value, float) else [int(n) for n in re.findall(r"\d+", value)]
            assert len(nums) in (1, 2)
            extracted[name] = {"range": [nums[0], nums[-1]],
                               "raw_unknotting_number": value,
                               "reference_html": sheet.cell_value(row, ref_col)}
        assert extracted == current
    with gzip.open(ROOT / "results/comparison/releases/knotinfo_2026.9.1.json.gz", "rt") as stream:
        september = json.load(stream)["rows"]
    assert set(september) == set(current)
    changed_ranges = {n for n in current if current[n]["range"] != september[n]["range"]}
    assert len(changed_ranges) == report["different_ranges_from_package_2026_9_1"] == 0
    from consolidate_results import reconstruct
    consolidated, _ = reconstruct()
    changed = consolidated["changed"]
    conflicts = {}
    partition = {"same_exact_value": 0, "different_exact_value": 0, "official_today_unresolved": 0}
    for name, row in changed.items():
        paper = row["new"]
        official = current[name]["range"]
        if paper[0] > official[1] or paper[1] < official[0]:
            conflicts[name] = {"paper": paper, "official_today": official,
                               "baseline": baseline[name], "exact_paper_value": paper[0] == paper[1]}
        if paper[0] == paper[1]:
            key = ("official_today_unresolved" if official[0] != official[1] else
                   "same_exact_value" if official == paper else "different_exact_value")
            partition[key] += 1
    assert conflicts == report["conflicts"]
    assert partition == report["exact_partition"]
    assert len(changed) == report["paper_changed_results"]
    assert sum(partition.values()) == report["exact_results_total"]
    assert len(conflicts) == report["conflicting_knots"] == 7
    exact_conflicts = sum(r["exact_paper_value"] for r in conflicts.values())
    assert exact_conflicts == report["conflicting_exact_paper_values"] == 6
    web = json.loads((HERE / "live_web_comparison.json").read_text())
    assert set(web["knots"]) == set(conflicts)
    for name, entry in web["knots"].items():
        raw = gzip.decompress((HERE / entry["retrieved_html_file"]).read_bytes())
        assert sha256(raw) == entry["html_sha256"]
        match = re.search(r"<td>unknotting_number:</td><td>(.*?)</td>", raw.decode(), re.S)
        assert match is not None
        value = html.unescape(re.sub(r"<[^>]*>", "", match[1])).strip()
        nums = [int(n) for n in re.findall(r"\d+", value)]
        assert [nums[0], nums[-1]] == entry["range"] == current[name]["range"]
    print(f"PASS: {len(current)} official ranges; zero differences from September package; "
          f"{len(conflicts)} incompatible knots = {exact_conflicts} exact results + 1 interval.")
    print("PASS: all seven directly queried live web pages match the official download.")
    print("Exact-result comparison:", partition)


if __name__ == "__main__":
    main()
