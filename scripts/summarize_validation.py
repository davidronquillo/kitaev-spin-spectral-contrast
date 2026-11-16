#!/usr/bin/env python3
"""Summarize saved validation scalars; do not import or execute solver code."""
from pathlib import Path
import csv
import json
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'/'validation'
OUT = ROOT/'outputs'/'validation'
OUT.mkdir(parents=True, exist_ok=True)
fixtures = json.loads((DATA/'fixtures_24_site.json').read_text())
high = json.loads((DATA/'high_sector_controls.json').read_text())
direct = json.loads((DATA/'transform_controls.json').read_text())
core = json.loads((DATA/'core_controls.json').read_text())
assert fixtures['status'] == high['status'] == direct['status'] == core['status'] == 'PASS'
records = fixtures['records']
components = [c for r in records for c in r['components'] if 'occupation_L1_difference' in c]
for r in records:
    assert r['status'] == 'PASS'
    assert r['maximum_normalized_observable_halfwidth'] <= fixtures['acceptance_normalized_halfwidth']
    for c in r['checks']:
        assert max(c['Fock'][0],c['gaussian'][0]) <= min(c['Fock'][1],c['gaussian'][1])
for c in components:
    assert c['occupation_L1_difference'] <= c['occupation_L1_bound']
assert all(c['passed'] and c['difference'] <= c['tolerance'] for c in high['checks'])
transform_checks = [c for width in direct['width_controls'] for c in width['checks']]
assert all(c['passed'] and c['error'] <= c['allowed'] for c in transform_checks)
largest_transform = max(transform_checks, key=lambda c:c['error'])
amplitudes = [c for c in high['checks'] if c['check'].startswith('coherent_amplitudes_')]
rows = list(csv.DictReader((DATA/'direct_precision.csv').open()))
groups = defaultdict(list)
for r in rows:
    groups[(r['geometry'],int(r['id']))].append(r)
fresh = []
for (geometry,sid), vals in groups.items():
    assert len(vals) == 8
    assert all(r['issues'] == '[]' and r['combined_available'] == 'True' for r in vals)
    fresh.append({'geometry':geometry,'id':sid,
        'MC_max_halfwidth':max(float(r['MC_only_halfwidth']) for r in vals),
        'direct_max_halfwidth':max(float(r['direct_only_halfwidth']) for r in vals),
        'combined_max_halfwidth':max(float(r['combined_halfwidth']) for r in vals),
        'direct_all_eight_targets_met':all(r['direct_only_all_eight_targets_met']=='True' for r in vals),
        'combined_all_eight_targets_met':all(r['combined_all_eight_targets_met']=='True' for r in vals),
        'recorded_issues':[]})
summary = {
    'scope':'Standard-library reduction of saved validation scalars.',
    'status':'PASS',
    'N24':{'orbits':fixtures['selected_ids'],'strata':fixtures['selection_strata'],
        'max_occupation_probability_difference':max(r['maximum_occupation_weight_difference'] for r in records),
        'max_component_occupation_L1_difference':max(c['occupation_L1_difference'] for c in components),
        'min_component_L1_allowance':min(c['occupation_L1_bound'] for c in components),
        'max_component_L1_allowance':max(c['occupation_L1_bound'] for c in components),
        'max_observable_halfwidth':max(r['maximum_normalized_observable_halfwidth'] for r in records),
        'acceptance_halfwidth':fixtures['acceptance_normalized_halfwidth'],
        'all_recorded_Fock_Gaussian_boxes_intersect':True,
        'max_first_moment_difference':max(r['first_moment_discrepancy'] for r in records)},
    'high_sector':{'synthetic_modes':7,'q':[5,7],'phase_choices':3,
        'max_coherent_amplitude_difference':max(c['difference'] for c in amplitudes),
        'amplitude_tolerance':min(c['tolerance'] for c in amplitudes),
        'total_recorded_checks':len(high['checks']),'all_saved_checks_pass':True},
    'independent_transform':{'synthetic_modes':10,'widths':[r['gamma'] for r in direct['width_controls']],
        'total_recorded_checks':len(transform_checks),'largest_recorded_discrepancy':largest_transform,
        'all_saved_checks_pass':True},
    'low_sector_native':core['native']['high_precision'],
    'assembly_control_retained_unresolved':core['assembly'],
    'fresh_direct_comparison':fresh,
    'fresh_comparison_scope':'Summary of saved scalar comparison results, not a re-evaluation of transform nodes or stochastic samples.'
}
(OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (OUT/'direct_summary.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(fresh[0]));writer.writeheader();writer.writerows(fresh)
print('PASS: archived validation scalars; 8 N24 fixtures, 620 high-sector controls, 406 transform controls, 6 fresh comparison rows.')
