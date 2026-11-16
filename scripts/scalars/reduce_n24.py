#!/usr/bin/env python3
"""Independent reductions of saved N24 scalar boxes and exact orbit weights.

No production module imports, network access, eigensolves, state construction or
spectral calculation. Numeric CSV cells round-trip their original binary64
endpoints; Decimal.from_float preserves those endpoints exactly. Directed
90-digit decimal arithmetic encloses weighted sums and normalization operations.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import csv, hashlib, json

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/scalars/n24'
OUT=ROOT/'outputs/scalars'
OUT.mkdir(parents=True,exist_ok=True)
def number(s): return D.from_float(float(s))
def down(fn):
    with localcontext() as c:
        c.prec=90;c.rounding=ROUND_FLOOR;return fn()
def up(fn):
    with localcontext() as c:
        c.prec=90;c.rounding=ROUND_CEILING;return fn()
def add(a,b):return [down(lambda:a[0]+b[0]),up(lambda:a[1]+b[1])]
def neg(a):return [a[1].copy_negate(),a[0].copy_negate()]
def sub(a,b):return add(a,neg(b))
def mul(a,b):return [down(lambda:min(x*y for x in a for y in b)),up(lambda:max(x*y for x in a for y in b))]
def div(a,b):
    assert b[0]>0
    return [down(lambda:min(x/y for x in a for y in b)),up(lambda:max(x/y for x in a for y in b))]
def depth(row,k):return sub(mul(add(row[k+'1'],row[k+'3']),[D('.5')]*2),row[k+'2'])
def weighted_mean(rows,key,den=None):
    den=sum(r['n'] for r in rows) if den is None else den
    return [down(lambda:sum(D(r['n'])*r[key][0] for r in rows)/D(den)),
            up(lambda:sum(D(r['n'])*r[key][1] for r in rows)/D(den))]
def weighted_median(rows,key):
    total=sum(r['n'] for r in rows);result=[]
    for end in (0,1):
        cumulative=0
        for r in sorted(rows,key=lambda r:(r[key][end],r['id'])):
            cumulative+=r['n']
            if 2*cumulative>total or (end==0 and 2*cumulative==total):
                result.append(r[key][end]);break
    return result
def scaled(a,factor):return mul(a,[D(factor)]*2)
def exported(a,places):
    q=D(1).scaleb(-places)
    return [str(a[0].quantize(q,rounding=ROUND_FLOOR)),str(a[1].quantize(q,rounding=ROUND_CEILING))]
def main():
    law=json.loads((DATA/'n24_exact_law.json').read_text())
    classes={int(r['probability_class']):int(r['exact_numerator']) for r in csv.DictReader((DATA/'n24_probability_classes.csv').open())}
    catalog={int(r['id']):r for r in csv.DictReader((DATA/'n24_orbit_catalog.csv').open())}
    weights={i:classes[int(r['probability_class'])] for i,r in catalog.items()}
    denominator=int(law['common_denominator'])
    assert len(catalog)==8192 and sum(weights.values())==denominator and min(weights.values())>0
    raw=b''.join(weights[i].to_bytes(law['bytes_per_numerator'],'big') for i in range(8192))
    assert hashlib.sha256(raw).hexdigest()==law['exact_Walsh_numerator_sha256']
    grouped={a:[] for a in ['g005','g010']}
    for rawrow in csv.DictReader((DATA/'n24_odd_scalar_boxes.csv').open()):
        row={'id':int(rawrow['id']),'n':weights[int(rawrow['id'])]}
        assert catalog[row['id']]['classification']=='odd_modified'
        for k in ['D','X18','X2','X22','w1','P1','P2','P3','S1','S2','S3']:
            row[k]=[number(rawrow[k+'_lo']),number(rawrow[k+'_hi'])]
            assert row[k][0]<=row[k][1]
        row['dP']=depth(row,'P');row['dS']=depth(row,'S');row['w1dP']=mul(row['w1'],row['dP'])
        grouped[rawrow['axis']].append(row)
    expected=json.loads((DATA/'n24_published_reference_summary.json').read_text())
    result={'status':'PASS','scope':'Reductions of supplied N24 row enclosures and exact integer orbit probabilities.',
            'all_orbits':8192,'eligible_odd_orbits':3892,'exact_law_numerator_digest_verified':True,'axes':{}}
    for axis,rows in grouped.items():
        assert len(rows)==3892 and len({r['id'] for r in rows})==3892
        med={k:weighted_median(rows,k) for k in ['D','X18','X2','X22']}
        for k,b in med.items():
            assert b==[number(x) for x in expected['axes'][axis]['conditional_odd']['median'][k]],(axis,k)
        mean=weighted_mean(rows,'D')
        oldmean=[number(x) for x in expected['axes'][axis]['conditional_odd']['mean']['D']]
        assert oldmean[0]<=mean[0]<=mean[1]<=oldmean[1]
        pooled=sub(div(weighted_mean(rows,'w1dP'),weighted_mean(rows,'w1')),weighted_mean(rows,'dS'))
        allmean=weighted_mean(rows,'D',denominator)
        oddnumerator=sum(r['n'] for r in rows)
        result['axes'][axis]={
            'odd_probability_exact':str(Fraction(oddnumerator,denominator)),
            'odd_probability_decimal':str(down(lambda:D(oddnumerator)/D(denominator))),
            'table_s2_normalized_x1000':{k:exported(scaled(b,1000),3) for k,b in med.items()},
            'table_s2_original_endpoint_identity':True,
            'table_s6_absolute_x1000':{'each_source':exported(scaled(mean,250),4),'after_averaging':exported(scaled(pooled,250),4)},
            'table_s6_unrounded_normalized':{'each_source':mean,'after_averaging':pooled},
            'all_source_D_absolute_x1000':{'outward_4dp':exported(scaled(allmean,250),4),'unrounded':scaled(allmean,250)},
            'all_source_median_exact':[0,0],
            'all_source_median_reason':'The exactly weighted even class has probability greater than one half, and its correction is identically zero.'}
        assert 2*oddnumerator < denominator
    (OUT/'n24_reduction_results.json').write_text(json.dumps(result,default=str,indent=2)+'\n')
    with (OUT/'n24_reproduced_tables.csv').open('w',newline='') as f:
        out=csv.writer(f);out.writerow(['table_or_figure','axis','quantity','lo','hi','units'])
        for a,v in result['axes'].items():
            for k,b in v['table_s2_normalized_x1000'].items():out.writerow(['Table S2',a,k,*b,'10^-3 normalized intensity'])
            for k,b in v['table_s6_absolute_x1000'].items():out.writerow(['Table S6',a,k,*b,'10^-3 absolute spin intensity'])
            out.writerow(['Figure 3 complete N24',a,'all-source D',*v['all_source_D_absolute_x1000']['outward_4dp'],'10^-3 absolute spin intensity'])
    print('PASS: exact N24 weighted medians match original reference endpoints; both normalization orders and complete-ensemble annotations reduced.')
    for a,v in result['axes'].items():
        print(a,'Table S2',v['table_s2_normalized_x1000'],'Table S6',v['table_s6_absolute_x1000'],'Fig3',v['all_source_D_absolute_x1000']['outward_4dp'])

if __name__=='__main__':main()
