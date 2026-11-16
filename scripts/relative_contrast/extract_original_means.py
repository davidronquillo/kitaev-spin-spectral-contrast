#!/usr/bin/env python3
"""Reduce the saved scalar collection to the four mean-contrast input boxes.

Run from repository root: python3 scripts/relative_contrast/extract_original_means.py
Only saved scalar enclosures are read. No production module is imported.
CSV float endpoints are recovered as their exact binary64 rational values;
N24 orbit weights are exact integers over a common denominator. Output decimal
endpoints are rounded outwards, in 10^-3 absolute-spin intensity units.
"""
import csv
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "relative_contrast"
EXPECTED_INPUT = ROOT / "data" / "relative_contrast" / "original_mean_inputs.csv"


def number(s):
    return Q.from_float(float(s))


def box(row, key):
    result = [number(row[key + '_lo']), number(row[key + '_hi'])]
    assert result[0] <= result[1]
    return result


def outward(x, places, lower):
    x *= 10 ** places
    k = x.numerator // x.denominator if lower else -((-x.numerator) // x.denominator)
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k // 10**places}.{k % 10**places:0{places}d}'


def weighted(rows, key):
    den = sum(w for w, _ in rows)
    return [sum(w * r[key][e] for w, r in rows) / den for e in (0, 1)]


def main():
    root = ROOT / 'data' / 'scalars'
    OUT.mkdir(parents=True, exist_ok=True)
    used = []

    def read(relative, is_json=False):
        p = root / relative
        used.append({'path': 'data/scalars/' + relative,
                     'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
        return json.loads(p.read_text()) if is_json else list(csv.DictReader(p.open()))

    law = read('n24/n24_exact_law.json', True)
    classes = {int(r['probability_class']): int(r['exact_numerator'])
               for r in read('n24/n24_probability_classes.csv')}
    catalog = {int(r['id']): r for r in read('n24/n24_orbit_catalog.csv')}
    weights = {i: classes[int(r['probability_class'])] for i, r in catalog.items()}
    assert len(weights) == 8192 and sum(weights.values()) == int(law['common_denominator'])
    assert hashlib.sha256(b''.join(weights[i].to_bytes(law['bytes_per_numerator'], 'big')
                                  for i in range(8192))).hexdigest() == law['exact_Walsh_numerator_sha256']
    grouped = {a: [] for a in ('g005', 'g010')}
    ids = {a: set() for a in grouped}
    for row in read('n24/n24_odd_scalar_boxes.csv'):
        i, a = int(row['id']), row['axis']
        assert catalog[i]['classification'] == 'odd_modified' and i not in ids[a]
        ids[a].add(i)
        s1, s2, s3 = (box(row, k) for k in ('S1', 'S2', 'S3'))
        f = [(s1[0] + s3[0]) / 2 - s2[1], (s1[1] + s3[1]) / 2 - s2[0]]
        grouped[a].append((weights[i], {'F': f, 'M': box(row, 'D')}))
    assert all(len(v) == 3892 for v in ids.values()) and ids['g005'] == ids['g010']
    odd_weight = sum(weights[i] for i in ids['g005'])

    fresh = read('fresh_endpoint_boxes.csv')
    fresh = [r for r in fresh if r['geometry'] == 'C3n16']
    all_ids = {int(r['id']) for r in fresh}
    even_ids = {int(r['id']) for r in fresh if r['eligibility'] == 'excluded'}
    assert all_ids == set(range(320)) and len(even_ids) == 150
    output = []
    for a, width in [('g005', '0.05'), ('g010', '0.10')]:
        rows = grouped[a]
        output.append((24, width, len(rows), 'exact conditional odd law', weighted(rows, 'M'), weighted(rows, 'F')))
        selected = [r for r in fresh if r['kernel'] == a + '_D']
        assert len(selected) == 170 and len({int(r['id']) for r in selected}) == 170
        assert {int(r['id']) for r in selected} == all_ids - even_ids
        large = []
        for r in selected:
            assert r['eligibility'] == 'known' and r['parity'] == 'odd'
            s = box(r, 'S')
            # The stored signed kernel is S(2)-[S(1)+S(3)]/2 = -d(S).
            large.append((1, {'M': box(r, 'X'), 'F': [-s[1], -s[0]]}))
        output.append((6144, width, len(large), 'realized eligible odd cohort', weighted(large, 'M'), weighted(large, 'F')))

    output.sort(key=lambda r: (r[0], r[1]))
    with (OUT / 'original_mean_inputs.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['N', 'gamma_over_Jsigma', 'eligible', 'population', 'M_lo', 'M_hi', 'F_lo', 'F_hi', 'units'])
        for n, g, count, pop, m, full in output:
            # 1000 / 4 converts normalized intensity to 10^-3 absolute-spin units.
            values = [outward(b[e] * 250, 16, e == 0) for b in (m, full) for e in (0, 1)]
            w.writerow([n, g, count, pop, *values, '10^-3 absolute spin intensity'])
    provenance = {
        'n24_odd_probability_exact': str(Q(odd_weight, int(law['common_denominator']))),
        'n24_weight_digest_verified': True,
        'input_files': used,
        'scope': 'Exact reductions of saved enclosures. Inherited numerical guarantees remain conditional; no spectral solver, resampling or physics execution.',
        'F_construction': {'24': 'Exact conditional weighting of three full-ordinate boxes.',
                           '6144': 'Negated saved signed full-response D kernel, averaged over all 170 eligible draws.'},
        'M_construction': 'Exact conditional N24 or equal-weight N6144 average of saved paired D boxes; each-source matching.',
        'output_rounding': 'Outward to 16 decimals in 10^-3 absolute spin units; precision serves reproducibility, not an additional scientific accuracy claim.'}
    (OUT / 'input_bindings.json').write_text(json.dumps(provenance, indent=2) + '\n')
    assert (OUT / 'original_mean_inputs.csv').read_bytes() == EXPECTED_INPUT.read_bytes(), 'Reconstructed input snapshot differs'
    print('PASS: four original-record mean boxes; 3892 exact-weight N24 orbits and all 170 eligible N6144 draws at each width.')


if __name__ == '__main__':
    main()
