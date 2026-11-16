"""Ten prescribed rational Lorentzians and 24 paired kernels."""
import hashlib
import json
import math
import numpy as np

EPS=float(np.finfo(float).eps)
AXES=('g005','g010')
CENTERS=('c100','c180','c200','c220','c300')
CENTER_VALUES=(1.,1.8,2.,2.2,3.)
SCALED_CENTERS=(20,36,40,44,60)
GAMMAS={'g005':.05,'g010':.10}
BASE_KEYS=tuple(a+'_'+c for a in AXES for c in CENTERS)
_rows=[]
def _add(name,terms,bounds):
    coef=[0.]*10
    for index,value in terms:coef[index]+=value
    _rows.append((name,coef,bounds))
for j,axis in enumerate(AXES):
    b=5*j
    for i,c in enumerate(CENTERS):_add(axis+'_'+c,[(b+i,1.)],(0.,1.))
    _add(axis+'_D',[(b+2,1.),(b,-.5),(b+4,-.5)],(-1.,1.))
    _add(axis+'_L',[(b,1.),(b+2,-1.)],(-1.,1.))
    _add(axis+'_R',[(b+4,1.),(b+2,-1.)],(-1.,1.))
for i,c in ((1,'c180'),(2,'c200'),(3,'c220')):
    _add('width_'+c,[(5+i,1.),(i,-1.)],(-1.,1.))
_add('width_D',[(7,1.),(5,-.5),(9,-.5),(2,-1.),(0,.5),(4,.5)],(-2.,2.))
for j,axis in enumerate(AXES):
    for i,c in ((1,'c180'),(3,'c220')):
        _add(axis+'_neighbor_'+c,[(5*j+i,1.),(5*j+2,-1.)],(-1.,1.))
KEYS=tuple(r[0] for r in _rows)
KERNEL_KEYS=KEYS
COEFFICIENTS=np.asarray([r[1] for r in _rows],dtype=float)
REGISTRY={}
for name,coef,bounds in _rows:
    absolute=[abs(x) for x in coef]
    REGISTRY[name]=dict(coefficients=coef,range=list(bounds),sup=max(abs(x) for x in bounds),
        oscillation=max(0.,bounds[1])-min(0.,bounds[0]),l1=math.fsum(absolute),
        energy_lipschitz=math.fsum(x/(.05 if i<5 else .10) for i,x in enumerate(absolute)),
        scalar_lipschitz=math.fsum(x*3*math.sqrt(3)/(8*(.05 if i<5 else .10)) for i,x in enumerate(absolute)))
    REGISTRY[name]['derivative']=REGISTRY[name]['scalar_lipschitz']
ALPHA=.01/(1920*2*12*32)
TARGET=.00025
PRIMARY_KEYS=tuple(a+'_'+c for a in AXES for c in ('c180','c200','c220','D'))
def binding():
    body=dict(base_keys=BASE_KEYS,keys=KEYS,registry=REGISTRY,
        filters='d=20*omega-C; C=20,36,40,44,60; narrow1/(d*d+1), broad4/(d*d+4)',
        required_statistics=24,MC_alpha_per_statistic_look=ALPHA,
        optional_mass_slot=None,remaining_slots='eight unused; no mass or squared-filter observable endpoints')
    body['sha256']=hashlib.sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return body
def base_values(energy):
    e=np.asarray(energy,dtype=float)
    d=20*e[...,None]-np.asarray(SCALED_CENTERS)
    z=d*d
    return np.concatenate((1/(z+1),4/(z+4)),axis=-1)
def values(energy):
    f=base_values(energy)
    return f@COEFFICIENTS.T
evaluate=values
def registry_sha256():return binding()['sha256']
def second_upper(key,positive_upper,mass_upper):
    r=REGISTRY[key]
    return min(r['sup']**2*mass_upper,r['sup']*math.fsum(abs(a)*positive_upper[b] for a,b in zip(r['coefficients'],BASE_KEYS)))
def controls():
    if len(KEYS)!=24 or len(set(KEYS))!=24 or len(BASE_KEYS)!=10:raise ArithmeticError('kernel inventory')
    if REGISTRY['width_D']['range']!=[-2.,2.]:raise ArithmeticError('paired depth range')
    return dict(status='PASS',binding=binding(),physical_arrays_read=0,H_actions=0)
