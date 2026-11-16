#!/usr/bin/env python3
"""Reduce saved finite-cohort scalar enclosures. Standard library only.

No production module is imported. No spectra, physical states, random samples,
bootstrap samples or fitted models are generated. Exact rational arithmetic
is used for registered binomial ranks, means, quotients and printed rounding.
"""
import csv
import json
import math
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'data/scalars'
OUT = ROOT/'outputs/scalars'
OUT.mkdir(parents=True,exist_ok=True)

def interval(row, key):
    return [Q.from_float(float(row[key+'_lo'])), Q.from_float(float(row[key+'_hi']))]

def outward(value, places, lower):
    n = value * 10**places
    k = n.numerator // n.denominator if lower else -((-n.numerator)//n.denominator)
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k//10**places}.{k%10**places:0{places}d}'

def rounded(box, factor, places):
    return [outward(box[0]*factor, places, True), outward(box[1]*factor, places, False)]

def mean(boxes):
    return [sum(b[e] for b in boxes)/len(boxes) for e in (0,1)]

def ranks(n, denominator):
    tail = 0
    best = None
    for k in range(1, n//2+1):
        tail += math.comb(n, k-1)
        if 2*denominator*tail > 2**n:
            break
        best = (k, n-k+1, Q(2**n-2*tail, 2**n))
    if best is None:
        raise ValueError('No registered finite rank exists')
    return best

def main():
    rows = list(csv.DictReader((DATA/'fresh_endpoint_boxes.csv').open()))
    expected = {(x['geometry'],x['axis'],x['endpoint']):x for x in
                csv.DictReader((DATA/'protected_medians.csv').open())}
    groups = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        groups[row['geometry']][int(row['id'])][row['kernel']] = row
    results = []
    counts = []
    coverage_failure = Q(1,100)
    for geometry, sources in sorted(groups.items()):
        assert set(sources) == set(range(320))
        eligible = [s for s in sources.values() if '' not in s]
        excluded = [s for s in sources.values() if '' in s]
        assert all(s['']['eligibility']=='excluded' and s['']['parity']=='even' for s in excluded)
        required = {a+'_'+k for a in ('g005','g010') for k in ('c180','c200','c220','D')}
        assert all(set(s)==required for s in eligible)
        assert all(v['eligibility']=='known' and v['parity']=='odd' for s in eligible for v in s.values())
        counts.append({'geometry':geometry,'N':int(next(iter(next(iter(sources.values())).values()))['N']),
                       'total':320,'odd':len(eligible),'even':len(excluded),'ambiguous':0})
        for axis in ('g005','g010'):
            for metric in ('c180','c200','c220','D'):
                boxes = [interval(s[axis+'_'+metric],'X') for s in eligible]
                denominator = 800 if metric in ('D','c200') else 2400
                lo_rank,hi_rank,coverage = ranks(len(boxes),denominator)
                lower = sorted(b[0] for b in boxes)[lo_rank-1]
                upper = sorted(b[1] for b in boxes)[hi_rank-1]
                endpoint = 'D' if metric=='D' else 'X_'+metric
                e = expected[geometry,axis,endpoint]
                assert (lower,upper)==(Q.from_float(float(e['lo'])),Q.from_float(float(e['hi']))), (geometry,axis,metric)
                assert int(e['known'])==len(boxes) and int(e['ambiguous'])==0
                coverage_failure += 1-coverage
                results.append({'geometry':geometry,'N':counts[-1]['N'],'axis':axis,'endpoint':endpoint,
                                'eligible':len(boxes),'ranks':[lo_rank,hi_rank],
                                'alpha':f'1/{denominator}','saved_bounds_normalized':[e['lo'],e['hi']],
                                'outward_bounds_1e_minus3_normalized':rounded([lower,upper],1000,3),
                                'matches_original_protected_table':True})
    normalization = []
    sources = groups['C3n16']
    for axis in ('g005','g010'):
        data = [s[axis+'_D'] for s in sources.values() if '' not in s]
        ms,mc,mw,md = [mean([interval(r,k) for r in data]) for k in ('S','C1','w1','X')]
        assert mw[0] > 0
        ratios = [c/w for c in mc for w in mw]
        pooled = [ms[0]-max(ratios),ms[1]-min(ratios)]
        # Signed D kernel is F(2)-[F(1)+F(3)]/2 = minus the mean-shoulder depth.
        # Hence D_each = X_D and D_pool = mean(S_D)-mean(C1_D)/mean(w1).
        each_printed = rounded(md,250,4)
        all_from_printed = [Q(x)*Q(17,32) for x in each_printed]
        normalization.append({'N':6144,'axis':axis,'odd':170,'total':320,
                              'units':'10^-3 absolute spin intensity',
                              'each_source':each_printed,'pooled':rounded(pooled,250,4),
                              'complete_cohort_from_original_rows':rounded([x*Q(17,32) for x in md],250,4),
                              'complete_cohort_from_printed_S6_bounds':rounded(all_from_printed,1,4),
                              'interpretation':'Numerical bounds for this realized cohort, not population confidence intervals.'})
    expect = {'g005':(['-0.9767','-0.8567'],['-1.2502','-1.1292']),
              'g010':(['-1.3953','-1.2656'],['-1.7007','-1.5704'])}
    for row in normalization:
        assert (row['each_source'],row['pooled']) == expect[row['axis']]
    report = {'status':'PASS','source_count':1920,'kept_all_sources':True,
              'retained_wider_numerical_boxes':True,'groups':counts,
              'fresh_Table_S2_intervals_reproduced':len(results),'fresh_Table_S2_scalar_endpoints_reproduced':2*len(results),
              'Table_S2':results,'N6144_Table_S6':normalization,
              'joint_coverage_floor_fraction':str(1-coverage_failure),
              'joint_coverage_floor_decimal':str(float(1-coverage_failure)),
              'scope':'Reductions of supplied finite-cohort numerical enclosures using exact ranks, means, and directed rounding.',
              'new_physics_calculations':0}
    (OUT/'fresh_scalar_reduction.json').write_text(json.dumps(report,indent=2)+'\n')
    with (OUT/'Table_S2_fresh_reproduced.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['geometry','N','width','endpoint','eligible','lower_rank','upper_rank','lo_1e_minus3','hi_1e_minus3'])
        for x in results:writer.writerow([x['geometry'],x['N'],'.05' if x['axis']=='g005' else '.10',x['endpoint'],x['eligible'],*x['ranks'],*x['outward_bounds_1e_minus3_normalized']])
    print('PASS: all 48 fresh Table S2 intervals (96 scalar endpoints) and both N6144 Table S6 rows reproduced.')
    print(json.dumps(normalization,indent=2))

if __name__=='__main__':
    main()
