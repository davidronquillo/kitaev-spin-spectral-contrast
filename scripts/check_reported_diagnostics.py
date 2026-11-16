#!/usr/bin/env python3
"""Check stored diagnostic/plotting scalars; no physics or plot regeneration."""
import csv
from decimal import Decimal as D
from pathlib import Path
import json
import re

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/diagnostics'
OUT=ROOT/'outputs/diagnostics'
OUT.mkdir(parents=True,exist_ok=True)
read=lambda name:list(csv.DictReader((DATA/name).open()))
geom={'24, all':('C3n01','exact_all_mean'),'24, odd':('C3n01','exact_odd_mean'),
      '96':('C3n02','fixed_odd_mean'),'384':('C3n04','fixed_odd_mean'),
      '1536':('C3n08','fixed_odd_mean'),'3456':('C3n12','fixed_odd_mean'),
      '6144':('C3n16','fixed_odd_mean'),'1536 rect.':('R24x32','fixed_odd_mean')}
raw=read('absolute_model_ordinates.csv')
lookup={(x['geometry'],x['estimand'],x['width'],x['c_Jsigma']):x for x in raw}
checked=0
for line in (DATA/'paper_tables/absolute_rows.txt').read_text().splitlines():
    if '&' not in line: continue
    label,width,*_=map(str.strip,line.split('&'))
    g,law=geom[label]; axis='g005' if width=='0.05' else 'g010'
    boxes=[(D(a),D(b)) for a,b in re.findall(r'\[([0-9.]+),([0-9.]+)\]',line)]
    for j,c in enumerate(('1','2','3')):
        row=lookup[g,law,axis,c]
        for k,key in enumerate(('P','S')):
            lo,hi=boxes[2*j+k]
            assert lo<=D(row[key+'_abs_lo'])*1000<=D(row[key+'_abs_hi'])*1000<=hi,(label,width,c,key)
            checked+=1
assert checked==96
window=read('window_fraction_bounds_N6144.csv')
assert len(window)==170 and len({x['id'] for x in window})==170
counts={'higher':0,'one_particle':0,'unresolved':0}; forced=[]
for x in window:
    lo,hi=D(x['f_h_minus_f_P_lo']),D(x['f_h_minus_f_P_hi'])
    assert lo<=hi
    label='higher' if lo>0 else ('one_particle' if hi<0 else 'unresolved')
    counts[label]+=1
    if label=='higher':forced.append(int(x['id']))
    if label=='unresolved':assert x['forced_order']=='unresolved'
assert counts=={'higher':37,'one_particle':0,'unresolved':133}
printed_ids=[int(x) for x in re.findall(r'\b\d{3}\b',(DATA/'paper_tables/forced_ids.txt').read_text())]
assert sorted(forced)==sorted(printed_ids)
candidates=read('physical_minimum_candidates.csv')
expected={'raw':set(map(D,['.080','.085','.090','.095','.105','.110','.115','.120'])),
          'inelastic':set(map(D,['.075','.080','.085','.090','.095','.110','.115','.120']))}
displacements={}
for layer in expected:
    low=[D(x['center_Jsigma']) for x in candidates if x['field']=='050' and x['layer']==layer]
    high=[D(x['center_Jsigma']) for x in candidates if x['field']=='055' and x['layer']==layer]
    shifts={b-a for a in low for b in high};assert shifts==expected[layer]
    assert D('.100') not in shifts
    displacements[layer]=[str(x) for x in sorted(shifts)]
curves=read('physical_field_curves.csv')
for field in ['050','055']:
    rows=sorted((x for x in curves if x['field']==field),key=lambda x:D(x['center_Jsigma']))
    assert len(rows)==401
    assert [D(x['center_Jsigma']) for x in rows]==[D('1.5')+D('.005')*i for i in range(401)]
    assert all(D(x['Y_lo'])<=D(x['Y_hi']) for x in rows)
report={'status':'PASS','table_S1_numeric_boxes_enclose_original_ordinate_rows':checked,
        'stored_window_bounds_classification':counts,'forced_ids_match_printed_list':True,
        'field_displacement_sets_Jsigma':displacements,'field_curve_points':len(curves),
        'scope':'Saved-value compatibility checks, window classifications, field grids, and discrete displacement sets.'}
(OUT/'reported_diagnostic_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
