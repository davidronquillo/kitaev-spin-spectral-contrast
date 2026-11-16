#!/usr/bin/env python3
"""Second arithmetic check using directed Decimal division, independent of reducer.

Checks the actual CSV inputs, all four corner quotients, expected printed-input
targets and outward 0.1-percentage-point output. Does not import the reducer.
"""
import csv
import json
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "relative_contrast"
OUT = ROOT / "outputs" / "relative_contrast"
EXPECTED = ROOT / "results" / "relative_contrast"


def main():
    published = json.loads((OUT/'relative_contrast_results.json').read_text())
    index = {(r['input_route'], r['N'], r['gamma_over_Jsigma']): r for r in published['rows']}
    targets = {(24, '0.05'): ['15.7', '15.8'], (24, '0.10'): ['17.4', '17.5'],
               (6144, '0.05'): ['-55.8', '-45.5'], (6144, '0.10'): ['-38.3', '-33.4']}
    checks = []
    for route, file in [('printed', 'printed_input_rows.csv'), ('original', 'original_mean_inputs.csv')]:
        for row in csv.DictReader((DATA/file).open()):
            key = (route, int(row['N']), row['gamma_over_Jsigma'])
            m = [D(row['M_lo']), D(row['M_hi'])]
            if route == 'printed':
                f = [(D(row['S1_lo'])+D(row['S3_lo']))/2-D(row['S2_hi']),
                     (D(row['S1_hi'])+D(row['S3_hi']))/2-D(row['S2_lo'])]
            else:
                f = [D(row['F_lo']), D(row['F_hi'])]
            assert f[0] > 0
            expected = []
            for rounding, take in [(ROUND_FLOOR, min), (ROUND_CEILING, max)]:
                with localcontext() as ctx:
                    ctx.prec = 70; ctx.rounding = rounding
                    q = take(100*a/b for a in m for b in f)
                    expected.append(str(q.quantize(D('0.1'), rounding=rounding)))
            assert expected == index[key]['percent_outward_0p1pp'], key
            if route == 'printed':
                assert expected == targets[key[1:]], (key, 'Printed Table S7 target differs')
            checks.append({'route': route, 'N': key[1], 'width': key[2], 'percent': expected, 'status': 'PASS'})
    result = {'status': 'PASS', 'method': 'Independent directed 70-digit Decimal arithmetic; no reducer import.',
              'checks': checks, 'new_physics': False,
              'limits': 'Arithmetic check on supplied boxes, not independent validation of their spectral construction.'}
    (OUT/'independent_arithmetic_check.json').write_text(json.dumps(result, indent=2)+'\n')
    assert (OUT/'independent_arithmetic_check.json').read_bytes() == (EXPECTED/'independent_arithmetic_check.json').read_bytes(), 'Expected arithmetic-check result differs'
    print('PASS: independent Decimal check agrees for all eight ratios and all four published printed-input targets.')


if __name__ == '__main__':
    main()
