#!/usr/bin/env python3
"""Compare the paper's crossing formulas with the deposited crossing data.

This is a short, single-process audit of stored numerical records. It does not
rebuild knot diagrams, compute their invariants afresh, or certify the recorded
knot identifications and unknotting numbers. All formula comparisons below use
exact rational arithmetic. The two input hashes identify precisely which
records were compared.

Run from the repository root:
    python audit_crossing_formulas.py
    python audit_crossing_formulas.py --out results/comparison/crossing_formula_verification.json
"""

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_records(path):
    """Hash the file as distributed, then read its possibly compressed JSON."""
    raw = path.read_bytes()
    payload = gzip.decompress(raw) if path.suffix == '.gz' else raw
    records = json.loads(payload)
    if not isinstance(records, list):
        raise ValueError(f'{path.name}: expected a JSON array')
    try:
        label = str(path.resolve().relative_to(ROOT))
    except ValueError:
        label = str(path.resolve())
    return records, {'path': label, 'sha256': hashlib.sha256(raw).hexdigest(),
                     'bytes': len(raw), 'records': len(records)}


def index_records(records, label):
    """Never silently discard repeated crossing identifiers during the join."""
    indexed = {}
    for row in records:
        key = (row['knot'], row['crossing'])
        if key in indexed:
            raise ValueError(f'{label}: duplicate crossing identifier {key}')
        indexed[key] = row
    return indexed


def rational(value):
    """Accept exact JSON integers or rational strings, never binary floats."""
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise ValueError(f'Expected an exact integer or fraction, got {value!r}')
    return Fraction(value)


def audit(dataset_path, electrical_path):
    dataset, dataset_info = load_records(dataset_path)
    electrical, electrical_info = load_records(electrical_path)
    if not dataset or not electrical:
        raise ValueError('Both inputs must contain crossing records')
    data_index = index_records(dataset, 'crossing dataset')
    electrical_index = index_records(electrical, 'electrical records')
    if data_index.keys() != electrical_index.keys():
        raise ValueError('The inputs do not contain the same crossing identifiers')

    rule_names = (
        'positive_tait_sign', 'oriented_crossing_sign',
        'resistance_between_zero_and_one', 'resistance_not_one_half',
        'dual_resistances_sum_to_one', 'determinant_from_resistance',
        'determinant_matches_electrical_record',
        'signature_change_from_resistance', 'self_linking_from_resistance',
        'recorded_u_one_passes_lickorish',
    )
    mismatches = {name: 0 for name in rule_names}
    examples = {name: [] for name in rule_names}
    signature_patterns = Counter()
    u_one_rows = 0

    def check(rule, condition, key):
        if not condition:
            mismatches[rule] += 1
            if len(examples[rule]) < 10:
                examples[rule].append({'knot': key[0], 'crossing': key[1]})

    for key, row in data_index.items():
        companion = electrical_index[key]
        resistance = rational(row['R_black'])
        dual_resistance = rational(row['R_white'])
        sign = row['sign']
        check('positive_tait_sign', row['tait_sign'] == 1, key)
        check('oriented_crossing_sign', sign in (-1, 1), key)
        check('resistance_between_zero_and_one', 0 < resistance < 1, key)
        check('resistance_not_one_half', resistance != Fraction(1, 2), key)
        check('dual_resistances_sum_to_one', resistance + dual_resistance == 1, key)
        check('determinant_from_resistance',
              rational(row['det_changed']) == rational(row['det']) * abs(1 - 2 * resistance), key)
        check('determinant_matches_electrical_record',
              row['det_changed'] == companion['det_changed_check'], key)

        # All deposited parent diagrams use the positive-definite checkerboard
        # form (Tait sign +1). In their oriented-crossing convention, sign +1
        # is type II and sign -1 is type I. The rank-one update changes the
        # form's signature by -2 exactly when R > 1/2. The type-II correction
        # contributes +2 to the knot signature change. Thus the prediction is
        # (1 + sign) - 2*[R > 1/2], the four cases of Proposition 4.1.
        actual_change = row['sigma_c'] - row['sigma']
        predicted_change = (1 + sign) - 2 * int(resistance > Fraction(1, 2))
        check('signature_change_from_resistance', actual_change == predicted_change, key)
        signature_patterns[(sign, resistance > Fraction(1, 2), actual_change)] += 1

        # This comparison includes non-generating classes [x]. No conclusion
        # about cyclicity follows solely from this self-linking calculation.
        predicted_linking = (resistance / (1 - 2 * resistance)) % 1 if resistance != Fraction(1, 2) else None
        check('self_linking_from_resistance',
              predicted_linking is not None and predicted_linking == rational(companion['self_linking_x']) % 1, key)
        if row.get('u_result') == 1:
            u_one_rows += 1
            # This is a consistency check against an archived label, not a
            # new proof of u=1. Passing Lickorish alone is inconclusive.
            check('recorded_u_one_passes_lickorish', companion['lickorish'] == 'pass', key)

    passed = not any(mismatches.values())
    return {
        'schema_version': 1,
        'status': 'PASSED' if passed else 'FAILED',
        'scope': 'Exact rational comparison of stored crossing-change data with the determinant, signature, and self-linking formulas.',
        'limitations': [
            'This audit does not recompute knot invariants from the underlying diagrams.',
            'Agreement between deposited records is a consistency test, not an independent certification of their original computation.',
            'The u_result=1 labels are archived inputs; passing the Lickorish obstruction does not establish unknotting number one.',
        ],
        'inputs': {'crossing_dataset': dataset_info, 'electrical_records': electrical_info},
        'implementation_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic': 'fractions.Fraction; exact integers and rational strings only',
        'execution': 'single process; no workers or external topology engines',
        'counts': {'parent_knots': len({key[0] for key in data_index}),
                   'crossing_changes': len(data_index), 'recorded_u_one_rows': u_one_rows},
        'mismatches': mismatches,
        'mismatch_examples': {name: values for name, values in examples.items() if values},
        'signature_convention': 'All parent Tait signs are +1; oriented crossing sign +1 is type II and -1 is type I.',
        'signature_patterns': [
            {'oriented_crossing_sign': sign, 'resistance_greater_than_one_half': above,
             'signature_change': change, 'count': count}
            for (sign, above, change), count in sorted(signature_patterns.items())
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', type=Path,
                        default=ROOT / 'results/crossing_changes/dataset_v2.json.gz')
    parser.add_argument('--electrical', type=Path,
                        default=ROOT / 'results/crossing_changes/lickorish_electrical_2026-09-07.json')
    parser.add_argument('--out', type=Path, help='Write the JSON summary to this explicit path')
    args = parser.parse_args()
    result = audit(args.dataset, args.electrical)
    output = json.dumps(result, indent=2) + '\n'
    print(output, end='')
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(output)
    return 0 if result['status'] == 'PASSED' else 1


if __name__ == '__main__':
    raise SystemExit(main())
