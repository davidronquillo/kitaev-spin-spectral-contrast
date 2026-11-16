#!/usr/bin/env python3
"""Offline ratio of odd-source mean contrasts; Python standard library only.

Run from repository root: python3 scripts/relative_contrast/reduce_relative_contrast.py
Inputs are decimal enclosure endpoints in 10^-3 absolute-spin intensity units.
All calculations use exact rational arithmetic. No independence is assumed
between numerator and denominator. The strictly positive full contrast is
checked before all four endpoint quotients are evaluated. Percentages are
rounded outwards to 0.1 percentage point. This is not a mean of source ratios,
a relative paired-median error or a large-cohort population confidence interval.
"""
import csv
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "relative_contrast"
OUT = ROOT / "outputs" / "relative_contrast"
EXPECTED = ROOT / "results" / "relative_contrast"


def outward(x, places, lower):
    x *= 10 ** places
    k = x.numerator // x.denominator if lower else -((-x.numerator) // x.denominator)
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k // 10**places}.{k % 10**places:0{places}d}'


def display(box, places):
    return [outward(v, places, e == 0) for e, v in enumerate(box)]


def interval(row, key):
    b = [Q(row[key + '_lo']), Q(row[key + '_hi'])]
    assert b[0] <= b[1]
    return b


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    output = []
    raw_results = {}
    for route, name in [('printed', 'printed_input_rows.csv'), ('original', 'original_mean_inputs.csv')]:
        for row in csv.DictReader((DATA / name).open()):
            n, width = int(row['N']), row['gamma_over_Jsigma']
            m = interval(row, 'M')
            if route == 'printed':
                a, b, c = (interval(row, k) for k in ('S1', 'S2', 'S3'))
                f = [(a[0] + c[0]) / 2 - b[1], (a[1] + c[1]) / 2 - b[0]]
            else:
                f = interval(row, 'F')
            assert f[0] > 0, (route, n, width, 'Denominator not strictly positive')
            q = [x / y for x in m for y in f]
            r = [min(q), max(q)]
            assert (m[0] > 0 and r[0] > 0) if n == 24 else (m[1] < 0 and r[1] < 0)
            # The spin-normalization factor must cancel in both numerator and denominator.
            assert [min(4*x/(4*y) for x in m for y in f), max(4*x/(4*y) for x in m for y in f)] == r
            raw_results[route, n, width] = (m, f, r)
            output.append({'input_route': route, 'N': n, 'gamma_over_Jsigma': width,
                           'population': row['population'], 'eligible': int(row['eligible']),
                           'M_absolute_x1000': display(m, 16), 'F_absolute_x1000': display(f, 16),
                           'denominator_strictly_positive': True,
                           'ratio_outward_16dp': display(r, 16),
                           'percent_outward_0p1pp': display([100*x for x in r], 1)})
    assert len(output) == 8
    for n, g in [(24, '0.05'), (24, '0.10'), (6144, '0.05'), (6144, '0.10')]:
        for p, o in zip(raw_results['printed', n, g], raw_results['original', n, g]):
            assert p[0] <= o[0] <= o[1] <= p[1], (n, g, 'Original interval outside conservative printed interval')
    result = {'status': 'PASS', 'definition': 'R = mean_odd[d(P)-d(S)] / d(mean_odd[S]); P is matched on each source.',
              'units': 'M,F: 10^-3 absolute spin intensity; R dimensionless.',
              'rounding': 'Exact rational endpoints; directed outward decimal rounding.',
              'coverage': 'Conditional on inherited saved-enclosure guarantees. No independence assumption or new population coverage.',
              'original_intervals_inside_printed_bounds': True, 'rows': output}
    (OUT / 'relative_contrast_results.json').write_text(json.dumps(result, indent=2) + '\n')
    with (OUT / 'relative_contrast_results.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['input_route', 'N', 'gamma_over_Jsigma', 'M_lo', 'M_hi', 'F_lo', 'F_hi', 'R_percent_lo', 'R_percent_hi'])
        for r in output:
            w.writerow([r['input_route'], r['N'], r['gamma_over_Jsigma'], *r['M_absolute_x1000'], *r['F_absolute_x1000'], *r['percent_outward_0p1pp']])
    for name in ('relative_contrast_results.json', 'relative_contrast_results.csv'):
        assert (OUT / name).read_bytes() == (EXPECTED / name).read_bytes(), ('Expected result differs', name)
    print('PASS: positive denominators, signs, units, endpoint division and containment verified.')
    for r in output:
        print(r['input_route'], r['N'], r['gamma_over_Jsigma'], r['percent_outward_0p1pp'], '%')


if __name__ == '__main__':
    main()
