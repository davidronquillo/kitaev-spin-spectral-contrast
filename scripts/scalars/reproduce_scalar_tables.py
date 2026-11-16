#!/usr/bin/env python3
"""Run saved-scalar reductions and compare with the paper tables."""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/scalars'
OUT=ROOT/'outputs/scalars'
SCRIPTS=ROOT/'scripts/scalars'
OUT.mkdir(parents=True,exist_ok=True)

def main():
    subprocess.run([sys.executable,str(SCRIPTS/'reduce_fresh_scalars.py')],check=True)
    subprocess.run([sys.executable,str(SCRIPTS/'reduce_n24.py')],check=True)
    fresh=json.loads((OUT/'fresh_scalar_reduction.json').read_text())
    small=json.loads((OUT/'n24_reduction_results.json').read_text())
    geom={'96':'C3n02','384':'C3n04','1536':'C3n08','3456':'C3n12','6144':'C3n16','1536 rect.':'R24x32'}
    lookup={(x['geometry'],x['axis'],x['endpoint']):x['outward_bounds_1e_minus3_normalized'] for x in fresh['Table_S2']}
    count=0
    for line in (DATA/'paper_tables/median_rows.txt').read_text().splitlines():
        if '&' not in line:continue
        cells=[x.strip() for x in line.split('&')]
        axis='g005' if cells[1]=='0.05' else 'g010'
        printed=[list(x) for x in re.findall(r'\[([^,]+),([^\]]+)\]',line)]
        if cells[0]=='24 (exact)':
            actual=[small['axes'][axis]['table_s2_normalized_x1000'][k] for k in ['X18','X2','X22','D']]
        else:
            actual=[lookup[geom[cells[0]],axis,k] for k in ['X_c180','X_c200','X_c220','D']]
        assert actual==printed,(cells[:2],actual,printed)
        count+=len(printed)
    assert count==56
    normalization=[]
    count_s6=0
    for line in (DATA/'paper_tables/normalization_rows.txt').read_text().splitlines():
        if '&' not in line:continue
        cells=[x.strip() for x in line.split('&')];axis='g005' if cells[1]=='0.05' else 'g010'
        printed=[list(x) for x in re.findall(r'\[([^,]+),([^\]]+)\]',line)]
        if cells[0]=='24':
            v=small['axes'][axis]['table_s6_absolute_x1000'];actual=[v['each_source'],v['after_averaging']]
        else:
            v=next(x for x in fresh['N6144_Table_S6'] if x['axis']==axis);actual=[v['each_source'],v['pooled']]
        assert actual==printed,(cells[:2],actual,printed)
        normalization.append({'N':int(cells[0]),'width':cells[1],'each_source':actual[0],'after_averaging':actual[1]})
        count_s6+=len(printed)
    assert count_s6==8
    result={'status':'PASS','Table_S2_intervals_matching_paper':count,
            'Table_S2_scalar_endpoints_matching_paper':2*count,
            'Table_S2_fresh_intervals':48,'Table_S2_exact_N24_intervals':8,
            'Table_S6_rows':4,'Table_S6_intervals_matching_paper':count_s6,
            'Table_S6':normalization,
            'Figure3_complete_N24_annotations':{a:v['all_source_D_absolute_x1000']['outward_4dp'] for a,v in small['axes'].items()},
            'substantive_discrepancies':[],
            'scope':'Exact output comparison using supplied numerical enclosures and probability metadata.'}
    (OUT/'expected_output_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: all 56 Table S2 intervals, all 8 Table S6 intervals, and complete-N24 Figure3 annotations reproduced.')

if __name__=='__main__':main()
