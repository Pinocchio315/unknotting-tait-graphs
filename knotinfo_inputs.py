"""Load the compact, authenticated KnotInfo cells used by the v1.8 calculations.

Only seven columns for the union of the Greene and Montesinos cohorts are
deposited. Cell values are unchanged from the dated official spreadsheet.
The manifest records the upstream commit and spreadsheet checksum. This
loader needs no spreadsheet package and does not contact the live database.
"""
import gzip
import hashlib
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent / 'results'
NAME = 'knotinfo_calculation_inputs_2026-09-09'


def load_metadata(results=RESULTS):
    results = Path(results)
    manifest = json.loads((results / (NAME + '_manifest.json')).read_text())
    data = (results / (NAME + '.json.gz')).read_bytes()
    if hashlib.sha256(data).hexdigest() != manifest['sha256']:
        raise ValueError('The frozen KnotInfo calculation cells have changed')
    rows = json.loads(gzip.decompress(data))
    if len(rows) != manifest['rows']:
        raise ValueError('The KnotInfo calculation cohort has changed')
    for name, row in rows.items():
        if row['name'] != name or set(row) != set(manifest['fields']):
            raise ValueError(f'Malformed KnotInfo calculation input: {name}')
    return rows
