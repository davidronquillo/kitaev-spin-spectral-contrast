"""Complete-frame source and noncommuting energy-operator certificates.

No eigensolver is called here. The owner supplies one positive eigensystem of
H=iA, and all diagnostics refer to that same represented frame. Particle-hole
conjugation supplies its complete negative partner. Public receipts are scalar.
"""
import math
import numpy as np
from scipy import linalg
from scipy.sparse import csr_matrix

EPS=float(np.finfo(float).eps)

def up(x): return math.nextafter(float(x),math.inf)
def down(x): return math.nextafter(float(x),-math.inf)
def gamma(k):
    x=float(k)*EPS
    if x>=.25: raise ArithmeticError('roundoff model outside admitted domain')
    return up(x/(1-x))
def fsumup(x): return up(math.fsum(float(z) for z in x)*(1+4*EPS)+1e-290)
def normF(a):
    return up(float(linalg.norm(a))*(1+gamma(4*np.asarray(a).size+16))+1e-145)
def public(cert):
    return {k:v for k,v in cert.items() if k not in ('O','Gram') and not k.startswith('_')}

def compact(value):
    """Scalar public receipt; full per-mode data remain hash-bound privately."""
    if isinstance(value,dict):return {k:compact(v) for k,v in value.items() if k not in ('O','Gram') and not k.startswith('_')}
    if isinstance(value,(list,tuple)) and len(value)>16:
        if all(isinstance(v,(int,float)) for v in value):
            return dict(count=len(value),minimum=min(value),maximum=max(value),representation='full values in bound private certificate')
        if all(isinstance(v,(list,tuple)) and len(v)==2 and all(isinstance(z,(int,float)) for z in v) for v in value):
            return dict(count=len(value),lower_extrema=[min(v[0] for v in value),max(v[0] for v in value)],
                        upper_extrema=[min(v[1] for v in value),max(v[1] for v in value)],representation='full intervals in bound private certificate')
    if isinstance(value,(list,tuple)):return [compact(v) for v in value]
    return value

def _parity(O,g):
    """LU determinant sign with an explicit backward-stability enclosure."""
    N=len(O);lu,piv=linalg.lu_factor(O,overwrite_a=False,check_finite=False)
    d=np.diag(lu)
    if np.any(d==0):return None,dict(status='UNRESOLVED',reason='zero LU pivot')
    L=np.tril(lu,-1);np.fill_diagonal(L,1.);U=np.triu(lu)
    # Standard partial-pivot LU backward bound |Delta O|<=gamma(3N)|L||U|.
    err=up(gamma(3*N+16)*normF(L)*normF(U))
    lower=down(math.sqrt(max(0.,1-g)))
    sign=int((-1 if np.count_nonzero(piv!=np.arange(N))%2 else 1)*np.prod(np.sign(d)))
    valid=err<lower
    return (sign if valid else None),dict(status='PASS' if valid else 'UNRESOLVED',
        determinant_sign=sign,LU_backward_F=err,singular_value_lower=lower,
        rule='homotopy preserves determinant sign when LU backward norm is below sigma_min(O)')

def analyze(A_sparse,E,Vplus,coefficient_error=0.,gram_receipt=None,residual_receipt=None,check=lambda:None):
    """Certify the complete PH frame with one sparse residual and real Gram.

    E contains sorted positive represented frequencies. coefficient_error is an
    operator bound from the exact prescribed A to its supplied binary matrix.
    Optional receipts must already include every arithmetic/formation allowance.
    """
    check();A=csr_matrix(A_sparse);E=np.asarray(E,float);V=np.asarray(Vplus,complex)
    N,n=V.shape
    if N!=2*n or E.shape!=(n,) or np.any(np.diff(E)<0) or np.any(E<0):raise ValueError('positive PH frame dimensions/order')
    if coefficient_error<0 or not np.all(np.isfinite(V)) or not np.all(np.isfinite(E)):raise ValueError('nonfinite frame')
    nnz=int(max(np.diff(A.indptr),default=0));arow=float(np.max(np.asarray(abs(A).sum(axis=1))))
    Hnorm=up(arow*(1+gamma(2*nnz+8))+coefficient_error)
    nv=np.linalg.norm(V,axis=0);pad=gamma(4*N+16);vlo=nv*(1-pad)-1e-145;vhi=nv*(1+pad)+1e-145
    if np.min(vlo)<=0:raise ArithmeticError('nonpositive vector norm bound')
    if residual_receipt is None:
        R=(1j*A)@V-V*E
        measured=np.linalg.norm(R,axis=0)
        err=(measured*(1+pad)+gamma(8*nnz+32)*(Hnorm+np.abs(E))*vhi+coefficient_error*vhi)*(1+32*EPS)+1e-145
        del R
    else:
        err=np.asarray(residual_receipt['residual_columns'],float)
        if err.shape!=(n,) or np.any(err<0):raise ValueError('invalid residual receipt')
    rad=(err/vlo)*(1+32*EPS)+1e-145
    Rplus=up(math.sqrt(fsumup(x*x for x in err))*(1+16*EPS))
    if gram_receipt is None:
        O=np.empty((N,N));O[:,0::2]=math.sqrt(2)*V.real;O[:,1::2]=-math.sqrt(2)*V.imag
        formation=up(gamma(8)*math.sqrt(2)*normF(V))
        check();G=O.T@O;G.flat[::N+1]-=1.
        vf=normF(O);g=up(normF(G)+gamma(N+8)*(vf*vf+math.sqrt(N)))
    else:
        O=np.asarray(gram_receipt['O']);G=gram_receipt.get('Gram')
        g=float(gram_receipt['gram_bound']);formation=float(gram_receipt['real_frame_formation_F'])
        if O.shape!=(N,N) or min(g,formation)<0:raise ValueError('invalid canonical Gram receipt')
    if not g<1:raise ArithmeticError('complete real frame Gram bound >=1')
    sourceg=up(g+2*math.sqrt(1+g)*formation+formation*formation)
    if not sourceg<1:raise ArithmeticError('complete complex frame Gram bound >=1')
    polar=up(g/(1+math.sqrt(1-g))*(1+32*EPS))
    sourcepolar=up(sourceg/(1+math.sqrt(1-sourceg))*(1+32*EPS))
    fullR=up(math.sqrt(2)*Rplus*(1+8*EPS))
    weyl=up(fullR+(Hnorm+float(np.max(E)))*sourcepolar)
    fullvalues=np.concatenate((-E[::-1],E));fullrad=np.concatenate((rad[::-1],rad))
    spacing=gamma(16)*(1+float(np.max(E)))
    disjoint=bool(np.all(np.diff(fullvalues)>fullrad[:-1]+fullrad[1:]+spacing))
    picked=np.minimum(rad,weyl) if disjoint else np.full(n,weyl)
    intervals=[[down(e-r),up(e+r)] for e,r in zip(E,picked)]
    sign,parity=_parity(O,g)
    check()
    return dict(schema='complete-PH-frame-v1',status='PASS',dimension=N,modes=n,
        E=E.tolist(),O=O,Gram=G,gram_bound=g,source_gram_bound=sourceg,
        real_frame_formation_F=formation,polar_displacement=polar,
        complex_positive_frame_displacement_F=up((polar+formation)/math.sqrt(2)),
        residual_columns=err.tolist(),residual_F_positive=Rplus,residual_F_complete=fullR,
        positive_radii=picked.tolist(),positive_intervals=intervals,
        full_spectral_bounds=[[-b,-a] for a,b in intervals[::-1]]+intervals,
        all_sorted_disjoint=disjoint,complete_Weyl_radius=weyl,
        H_norm_upper=Hnorm,coefficient_error=float(coefficient_error),max_NNZ_per_row=nnz,
        positive_vector_norm_upper=vhi.tolist(),parity=sign,parity_certificate=parity,
        source_positive_gap_lower=min(x[0] for x in intervals),
        lowest_occupation_gap_lower=(intervals[1][0]-intervals[0][1]) if n>1 else None,
        rounding_model='binary64 bounded sparse products, real Gram, Frobenius reductions and complete polar frame; per-mode residual disks only when all complete PH disks are disjoint')

def _leakage(cert,occupied_indices=()):
    E=np.asarray(cert['E'],float);r=np.asarray(cert['residual_columns'],float)
    occupied={int(i) for i in occupied_indices};n=len(E)
    if any(i<0 or i>=n for i in occupied):raise ValueError('occupation index')
    iv=np.asarray(cert['positive_intervals'],float)
    if float(np.min(iv[:,0]))<=0:
        return dict(valid=False,status='UNRESOLVED',reason='positive/negative spectral membership not separated')
    selected=np.asarray([-e if i in occupied else e for i,e in enumerate(E)])
    selected_intervals=np.asarray([[-iv[i,1],-iv[i,0]] if i in occupied else iv[i] for i in range(n)])
    complement=np.asarray([iv[i] if i in occupied else [-iv[i,1],-iv[i,0]] for i in range(n)])
    separations=[];isolation=[]
    for x,si in zip(selected,selected_intervals):
        isolation.append(down(float(np.min(np.maximum(np.maximum(complement[:,0]-si[1],si[0]-complement[:,1]),0.)))))
        distances=np.maximum(np.maximum(complement[:,0]-x,x-complement[:,1]),0.)
        separations.append(down(float(np.min(distances))))
    if min(isolation)<=0:
        return dict(valid=False,status='UNRESOLVED',reason='occupied target cluster not isolated from complementary spectral intervals',selected_interval_isolation=isolation)
    if min(separations)<=0:return dict(valid=False,status='UNRESOLVED',reason='selected/complement spectral separation not positive',separations=separations)
    g=float(cert['source_gram_bound'])
    raw=up(math.sqrt(fsumup((a/d)**2 for a,d in zip(r,separations))))
    leakage=up(raw/math.sqrt(1-g)*(1+16*EPS))
    return dict(valid=True,status='PASS',selected_complement_separations=separations,
                selected_interval_isolation=isolation,
                raw_leakage_F=raw,polar_subspace_leakage_F=leakage,selected_Gram_F=g,
                occupied_indices=sorted(occupied))

def assembly_allowance(cert,occupied_indices=()):
    """Direct S S* dense construction, real/skew projection and c0 reflection.

    Real/skew orthogonal projection is contractive. The operation bound covers
    forming i(2 S S*−I), real extraction, averaging and the exact sign reflection.
    """
    n=int(cert['modes']);N=2*n;vf2=fsumup(v*v for v in cert['positive_vector_norm_upper'])
    return up(2*gamma(8*n+32)*vf2+gamma(24)*(2*vf2+math.sqrt(N)))

def occupied(cert,occupied_indices=(),assembly_error_F=None):
    """Physical occupied covariance -> unchanged raw represented covariance.

    With S the selected raw columns and W=polar(S), complement leakage is at
    most sqrt(sum_j (r_j/d_j)^2)/sqrt(1-g). Projector distance is sqrt(2)
    times leakage, while ||SS*−WW*||F<=g. Multiply by2 for the Majorana
    covariance. The source reflection is an exact orthogonal conjugation.
    """
    out=_leakage(cert,occupied_indices)
    if not out['valid']:return dict(out,initial_covariance_F=None,prechart_trace_budget=1.)
    err=assembly_allowance(cert,occupied_indices) if assembly_error_F is None else float(assembly_error_F)
    if err<0:raise ValueError('negative source assembly allowance')
    projector=up(math.sqrt(2)*out['polar_subspace_leakage_F'])
    covariance=up(2*(projector+out['selected_Gram_F'])+err)
    return dict(out,occupied_projector_F=projector,raw_projector_polar_F=out['selected_Gram_F'],
                source_assembly_F=err,initial_covariance_F=covariance,
                prechart_trace_budget=min(1.,up(covariance/math.sqrt(8))),
                norm_units='projector_F and covariance_F are Frobenius; prechart_trace_budget is only an additive budget, not a trace-distance certificate for the nonpure raw covariance; include pure-chart reconstruction before division by sqrt(8)')

def final_alignment(cert,canonical_displacement_F=None):
    """Aligned exact positive subspace AND its full compressed Hamiltonian.

    If B is the exact canonical frame underlying the represented chart, let
    s>=||(I-P)B||op and U be the canonical direct rotation to exact P. With
    C=(B*PB)^1/2 and Y=(I-P)B, B*HB=C h C+Y*H_-Y. Hence
    ||h−B*HB||op <=3||H||s². The retained B*HB−diag(E) residual bounds
    within-positive rotations without dividing by interior positive gaps.
    """
    out=_leakage(cert,())
    if not out['valid']:return dict(out,final_alignment_covariance_F=None,one_body_operator_error=None)
    d=float(cert['complex_positive_frame_displacement_F'] if canonical_displacement_F is None else canonical_displacement_F)
    if d<0:raise ValueError('negative canonical frame displacement')
    # Raw V leakage is bounded before polarizing; B is already exactly canonical.
    leakage=up(out['raw_leakage_F']+d)
    if leakage>=1:return dict(out,valid=False,status='UNRESOLVED',reason='canonical direct rotation overlap not certified',final_alignment_covariance_F=None,one_body_operator_error=None)
    rotation=up(2*leakage/math.sqrt(1+math.sqrt(1-leakage*leakage))*(1+32*EPS))
    H=float(cert['H_norm_upper']);emax=max(cert['E'])
    RBF=up(cert['residual_F_positive']+(H+emax)*d)
    leakage_energy=up(3*H*leakage*leakage)
    eh=up(RBF+leakage_energy)
    return dict(out,canonical_displacement_F=d,canonical_leakage_F=leakage,
        direct_rotation_F=rotation,final_alignment_covariance_F=up(2*rotation),
        canonical_positive_residual_F=RBF,compression_leakage_operator_allowance=leakage_energy,
        one_body_operator_error=eh,
        q_sector_rule='||dGamma_q(h)-dGamma_q(diag E)||op <= q*one_body_operator_error',
        filter_rule='noncommuting Lorentzian resolvent identity gives ||F_gamma(A)-F_gamma(B)||op <= min(1,||A-B||op/gamma); signed coefficient sums and squared-filter factors are applied separately')

def transport(initial,final,reconstruction_F,old_t=None,old_e_h=None):
    if not initial.get('valid') or not final.get('valid'):
        return dict(status='UNRESOLVED',valid=False,t=1.,reason='source or final alignment not certified')
    total=up(initial['initial_covariance_F']+final['final_alignment_covariance_F']+float(reconstruction_F))
    t=min(1.,up(total/math.sqrt(8)))
    if old_t is not None:t=min(t,float(old_t))
    eh=final['one_body_operator_error']
    if old_e_h is not None:eh=max(eh,float(old_e_h))
    return dict(status='PASS',valid=True,t=t,initial_covariance_F=initial['initial_covariance_F'],
        final_alignment_covariance_F=final['final_alignment_covariance_F'],reconstruction_F=float(reconstruction_F),
        e_h=eh,derived_energy_operator_error=final['one_body_operator_error'],
        inherited_energy_bound_sufficient=None if old_e_h is None else float(old_e_h)>=final['one_body_operator_error'],
        identity='unchanged represented x/B/E and source assembly; new physical-to-represented bound only')

def from_saved(row,spectra):
    """Revalidate an unchanged source using authenticated scalar receipts.

    The saved complete sparse residual upper bound encloses every positive
    column collectively. Its maximum normalized eigenpair radius, multiplied
    by a complete-frame norm upper bound, encloses each raw residual. The
    measured real PH Gram is used, not just Vplus*Vplus. No smaller eigenvalue
    radius is substituted for a projector perturbation certificate.
    """
    certificates={}
    N=int(row['N']);n=N//2
    for side,key in (('initial','Ei'),('final','Ef')):
        E=np.asarray(spectra[key],float);radii=np.asarray(spectra[key+'_radii'],float)
        d=row[side+'_diagnostics'];p=row[side+'_spectrum_enclosure']
        if E.shape!=(n,) or radii.shape!=(n,) or not np.all(np.isfinite(E)):raise ValueError('saved spectrum dimensions')
        completeg=float(p['polar_Gram_bound'])
        if not 0<=completeg<1:raise ArithmeticError('saved complete Gram unqualified')
        # The PH real frame is formed from positive columns, so its own Gram
        # defect, including Vplus^T Vplus, is indispensable near zero modes.
        measured=float(d['canonical_O_orthogonality_F'])
        g=up(measured*(1+gamma(4*N*N+16))+gamma(N+8)*(2*N+math.sqrt(N)))
        formation=up(gamma(8)*math.sqrt(N*(1+completeg)))
        sourceg=up(g+2*math.sqrt(1+g)*formation+formation*formation)
        if sourceg>=1:raise ArithmeticError('saved complete PH frame unresolved')
        vhi=up(math.sqrt(1+completeg))
        residual=up(float(p['max_sparse_eigenpair_radius'])*vhi)
        iv=[[down(e-r),up(e+r)] for e,r in zip(E,radii)]
        c=dict(schema='saved-complete-PH-frame-v1',status='PASS',dimension=N,modes=n,E=E.tolist(),
            gram_bound=g,source_gram_bound=sourceg,real_frame_formation_F=formation,
            polar_displacement=up(g/(1+math.sqrt(1-g))),
            residual_columns=[residual]*n,residual_F_positive=float(p['sparse_residual_F_upper']),
            residual_F_complete=up(math.sqrt(2)*float(p['sparse_residual_F_upper'])),
            positive_radii=radii.tolist(),positive_intervals=iv,
            full_spectral_bounds=[[-b,-a] for a,b in iv[::-1]]+iv,
            all_sorted_disjoint=bool(p['all_sorted_residual_disks_disjoint']),
            complete_Weyl_radius=float(p['independent_Weyl_radius']),
            H_norm_upper=float(p['sparse_absolute_row_norm']),coefficient_error=0.,
            positive_vector_norm_upper=[vhi]*n,parity=d.get('parity'),
            parity_certificate=dict(status='inherited',scope='unchanged authenticated classification'),
            sparse_residual_provenance='complete cached sparse residual certificate, not the generic dense delta_A',
            real_Gram_provenance='same positive frame measured by pinned gaussian.analyze canonical_O_orthogonality_F')
        c['complex_positive_frame_displacement_F']=up((c['polar_displacement']+formation)/math.sqrt(2))
        certificates[side]=c
    occ=[] if not row['initial_occupancy'] else [0]
    # The historical source used vacuum SS* plus a rank-one occupation
    # update before reflection. Cover its separate arithmetic and cleanups.
    cleanup=float(row['initial_diagnostics']['covariance_defects']['cleanup_F'])
    cleanup+=float(row.get('source_covariance',{}).get('occupied_covariance_defects',{}).get('cleanup_F',0.))
    assembly=up(assembly_allowance(certificates['initial'],occ)+cleanup+gamma(64*N+64))
    initial=occupied(certificates['initial'],occ,assembly_error_F=assembly)
    dg=row['chart_diagnostics']
    displacement=up(float(dg['total_frame_displacement_allowance_F'])+gamma(32)*math.sqrt(N))
    final=final_alignment(certificates['final'],displacement)
    reconstruction=fsumup(dg[k] for k in ('covariance_reconstruction_F','reconstruction_block_calculation_F',
                                        'reconstruction_polar_tail_F','reconstruction_site_arithmetic_F'))
    tr=transport(initial,final,reconstruction,row['t'],row['e_h'])
    return dict(initial=initial,final=final,transport=tr,
                certificates={k:public(v) for k,v in certificates.items()},
                assumptions=['pinned source and scalar diagnostics are the same calculation as bound spectra and chart bytes',
                             'fixed dyadic sparse hopping representation is unchanged',
                             'no new source/chart/proposal/kernel/stream is substituted'])

def controls(charge=None):
    """Owner-only synthetic scalar/matrix controls; no Hamiltonian eigensolves."""
    rows=[]
    def test(name,ok):
        rows.append(dict(name=name,passed=bool(ok)))
        if not ok:raise ArithmeticError('certificate control: '+name)
    for energies in ([.125,.5,1.],[2**-20,.5,1.]):
        E=np.asarray(energies);N=2*len(E);A=np.zeros((N,N));V=np.zeros((N,len(E)),complex)
        for j,e in enumerate(E):A[2*j,2*j+1]=e;A[2*j+1,2*j]=-e;V[2*j,j]=1/math.sqrt(2);V[2*j+1,j]=-1j/math.sqrt(2)
        c=analyze(A,E,V);o=occupied(c,[0]);f=final_alignment(c)
        test('exact PH complete disks '+str(E[0]),c['all_sorted_disjoint'])
        test('exact PH parity '+str(E[0]),c['parity']==1)
        test('occupied-source certificate '+str(E[0]),o['valid'] and o['initial_covariance_F']>=0)
        test('full operator certificate '+str(E[0]),f['valid'] and f['one_body_operator_error']>=0)
        # Pure within-positive mixing leaves the vacuum projector invariant,
        # but the old diagonal energy representation changes noncommutatively.
        R=np.eye(len(E));theta=1e-4;R[-2:,-2:]=[[math.cos(theta),-math.sin(theta)],[math.sin(theta),math.cos(theta)]]
        W=V@R;cc=analyze(A,E,W);ff=final_alignment(cc)
        actual=float(linalg.norm(W.conj().T@(1j*A)@W-np.diag(E)))
        test('within-positive rotation has nonzero operator error '+str(E[0]),ff['one_body_operator_error']>1e-6 and actual<=ff['one_body_operator_error']+1e-13)
    # Degenerate lowest occupied-mode choice must remain unresolved; vacuum
    # positive-vs-negative separation is nevertheless valid via full Weyl.
    E=np.array([.5,.5,1.]);A=np.zeros((6,6));V=np.zeros((6,3),complex)
    for j,e in enumerate(E):A[2*j,2*j+1]=e;A[2*j+1,2*j]=-e;V[2*j,j]=1/math.sqrt(2);V[2*j+1,j]=-1j/math.sqrt(2)
    c=analyze(A,E,V)
    test('degenerate complete frame uses polar Weyl',not c['all_sorted_disjoint'])
    test('degenerate occupied singleton not silently chosen',not occupied(c,[0])['valid'])
    test('degenerate vacuum remains separated',occupied(c,[])['valid'])
    return dict(status='PASS',checks=rows,quadratic_factorization_attempts=0,physical_state_reads=0,physical_H_actions=0)
