"""Gaussian-chart construction with direct source contractions and occupation spectrum.

Phase-consistent two-sector Gaussian chart and ordinary numerical diagnostics.

Inputs are an admitted allowed-initial site covariance after the source
reflection and the already validated final positive-energy columns. All
factorizations here concern the frame, occupation, or pairing chart.
"""
import math

import numpy as np
from scipy import linalg


EPS = float(np.finfo(np.float64).eps)


class ChartUnresolved(ValueError):
    """Intrinsic numerical obstruction, carrying the attained diagnostics."""

    def __init__(self, reason, diagnostics):
        super().__init__(reason)
        self.diagnostics = diagnostics
        self.reason = reason


def frob(a):
    return float(np.linalg.norm(a, ord="fro"))


def gamma_n(k):
    ku = k * EPS
    return ku / (1.0 - ku)


def canonical_frame(V, diag, certificate=None):
    """One Newton polar step on the real canonical frame, with a posteriori error.

    An exact nearby orthogonal polar frame Q* is within
    ||Q^T Q-I||_F/(1+sqrt(1-||Q^T Q-I||_F)) of a near-orthogonal Q.
    The complex positive block built from its interleaved real pairs preserves
    particle-hole conjugation. We carry the remaining correction explicitly.
    """
    m,n=V.shape
    if certificate is None:
        frame=np.empty((m,m));frame[:,::2]=math.sqrt(2)*V.real;frame[:,1::2]=-math.sqrt(2)*V.imag
        defect=frame.T@frame-np.eye(m);measured=frob(defect)
        radius=measured+gamma_n(m+8)*(frob(frame)**2+math.sqrt(m))
        formation=gamma_n(16)*math.sqrt(m)
        if radius>=1:raise ChartUnresolved("canonical polar frame unresolved",diag)
        polar=radius/(1+math.sqrt(1-radius))
        displacement=(polar+formation)/math.sqrt(2)
    else:
        frame=certificate["O"];measured=frob(certificate["Gram"])
        radius=certificate["gram_bound"];polar=certificate["polar_displacement"]
        displacement=certificate["complex_positive_frame_displacement_F"]
    # The exact polar frame is counterfactual, not another matrix product.
    # Returned V is the represented contraction frame; every discrepancy from
    # its exact canonical polar is explicitly charged in reconstruction.
    diag.update(frame_dimension=m,positive_modes=n,original_frame_gram_defect_F=measured,
        cleaned_frame_gram_defect_F=measured,frame_newton_displacement_F=0.,
        exact_polar_tail_allowance_F=polar,total_frame_displacement_allowance_F=displacement,
        canonical_frame_rule="raw positive frame with certified counterfactual real polar; no Newton products")
    return V,displacement,math.sqrt(2)*displacement


def construct_chart(gamma, V, energies, odd, initial_covariance_F=0.0, final_certificate=None):
    """Return x (odd only), projected skew B, Z_b and the separate error budgets.

    t is the trace-distance allowance for the ideal normalized Gaussian
    represented by the returned chart versus the input ray in matching final
    coordinates. e_h is the aligned one-mode energy-operator allowance.
    chart_mass_roundoff covers representing the exact orthogonal blocked/skew
    amplitudes by stored entries; normalization_relative_error is the separate
    relative error in sqrt(Z_b). The streaming routine handles its own
    amplitude enumeration, line energy, filter, and reductions.
    """
    m,n=V.shape;energies=np.asarray(energies,dtype=float)
    if (m!=2*n or n==0 or gamma.shape!=(m,m) or energies.shape!=(n,) or
        not np.all(np.isfinite(energies)) or not np.all(energies>0) or
        not np.all(np.isfinite(gamma)) or not np.all(np.isfinite(V)) or
        not np.isfinite(initial_covariance_F) or initial_covariance_F<0):
        raise ValueError("invalid chart input")
    diag={"odd":bool(odd),"initial_covariance_F":float(initial_covariance_F)}
    positive,frame_shift,polar_tail=canonical_frame(V,diag,final_certificate)
    if final_certificate is not None:
        from numerical_certificates import final_alignment
        align=final_alignment(final_certificate,frame_shift)
        if not align['valid']:raise ChartUnresolved(align.get('reason','final alignment unresolved'),diag)
        final_cov=align['final_alignment_covariance_F'];energy_error=align['one_body_operator_error']
        diag['final_alignment_certificate']=align
    else:
        final_cov=0.;energy_error=0.
    diag['final_alignment_covariance_F']=final_cov

    source_matrix = 0.5 * (np.eye(m) + 1j*gamma)
    left = positive.conj().T @ source_matrix
    raw_R = left @ positive
    raw_K = left @ positive.conj()
    r_defect = frob(raw_R-raw_R.conj().T)
    k_defect = frob(raw_K+raw_K.T)
    R = 0.5*(raw_R+raw_R.conj().T)
    K = 0.5*(raw_K-raw_K.T)
    occ_values, occ_vectors = linalg.eigh(R, driver="evd", check_finite=True)
    occ_residual = frob(R @ occ_vectors - occ_vectors * occ_values)
    occ_orth = frob(occ_vectors.conj().T @ occ_vectors-np.eye(n))
    contract_roundoff = gamma_n(8*m+16)*(2*m+frob(gamma)+frob(R))
    tau = (0.5*(initial_covariance_F+final_cov)+polar_tail*2+0.5*(r_defect+k_defect)
           +occ_residual+frob(R)*occ_orth+contract_roundoff)
    diag.update(R_hermiticity_defect_F=r_defect, K_skew_defect_F=k_defect,
                R_eigen_residual_F=occ_residual, R_eigen_orthogonality_F=occ_orth,
                chart_contraction_arithmetic=contract_roundoff,
                R_spectral_resolution=tau, min_R_eigenvalue=float(occ_values[0]),
                max_R_eigenvalue=float(occ_values[-1]))
    if occ_values[0] < -tau or occ_values[-1] > 1+tau:
        raise ChartUnresolved("occupation spectrum outside probability range", diag)

    eye = np.eye(n, dtype=complex)
    x = None
    P = np.zeros_like(eye)
    if odd:
        x = occ_vectors[:, -1].copy()
        x /= np.linalg.norm(x)
        P = np.outer(x,x.conj())
        separation = float(occ_values[-1]-occ_values[-2]) if n>1 else 1.
        block_residual = float(np.linalg.norm(R@x-x))
        diag.update(blocked_eigenvalue_gap=separation, blocked_eigen_residual=block_residual)
        if separation <= 2*tau or block_residual > tau:
            raise ChartUnresolved("unique occupied blocked direction unresolved", diag)
        smin = min(2-float(occ_values[-1]), 1-float(occ_values[-2]) if n>1 else 1.)
    else:
        smin = 1-float(occ_values[-1])
    system = eye-R+P
    diag.update(deflated_smin_estimate=smin, deflated_smax_upper=2+tau,
                deflated_condition_estimate=(2+tau)/smin if smin>0 else None)
    if smin <= 2*tau:
        raise ChartUnresolved("unblocked pairing chart unresolved", diag)
    unclean_B = linalg.solve(system,K,assume_a="pos",check_finite=True)
    skew = (unclean_B-unclean_B.T)*0.5
    B = skew
    if odd:
        B = skew - np.outer(x,x.conj()@skew)
        B -= np.outer(B@x.conj(),x)
        B = (B-B.T)*0.5
    bnorm = frob(B)
    skew_error = frob(B+B.T)
    blocked_leak = float(np.linalg.norm(x.conj()@B)) if odd else 0.
    linear_residual = frob(system@B-K)
    diag.update(raw_B_skew_defect_F=frob(unclean_B+unclean_B.T), final_B_skew_defect_F=skew_error,
                B_projection_cleanup_F=frob(B-unclean_B), final_x_dagger_B_norm=blocked_leak,
                B_F=bnorm, deflated_linear_residual_F=linear_residual)
    if linear_residual > 10*math.sqrt(n)*tau*(1+bnorm):
        raise ChartUnresolved("inconsistent deflated solve after projection", diag)

    # Reconcile stored x,B with a mathematical unit x and exactly projected
    # skew B. Spectral norm of the omitted block is bounded by its Frobenius
    # norm. Exterior multiplication by a normalized x has operator norm 1.
    x_error = abs(float(np.linalg.norm(x))-1)+gamma_n(8*n+8) if odd else 0.
    b_error = 0.5*skew_error+2*blocked_leak+4*x_error*bnorm+gamma_n(8*n+16)*(1+bnorm)
    amplitude_error = x_error*(1+bnorm)+b_error
    mass_roundoff = (2+amplitude_error)*amplitude_error

    singular = linalg.svdvals(B,check_finite=True)
    log_z = -0.5*math.fsum(math.log1p(float(s)**2) for s in singular)
    Z_b = math.exp(log_z)
    if Z_b == 0:
        raise ChartUnresolved("pairing chart normalization underflow", dict(diag,log_Z_b=log_z))
    # d log(sqrt Z_b)/ds = -s/[2(1+s²)], absolute value <= 1/4.
    # Singular-value changes in l2 are <= ||dB||_F. Include SVD backward
    # error and fsum/log1p/exp arithmetic in the logarithmic radius.
    svd_error = 64*n*EPS*bnorm
    log_amp_error = 0.25*math.sqrt(n)*(svd_error+b_error)+gamma_n(8*n+16)*(1+abs(log_z))
    relative_norm_error = math.expm1(log_amp_error)+4*EPS
    diag.update(log_Z_b=log_z, normalization_svd_F_allowance=svd_error,
                normalization_log_sqrt_Z_error=log_amp_error, x_representation_error=x_error,
                projected_B_representation_error_F=b_error, chart_mass_roundoff=mass_roundoff)

    G = eye+B@B.conj().T
    chol = linalg.cho_factor(G,lower=True)
    invG = linalg.cho_solve(chol,eye)
    RB = P + eye - invG
    KB = linalg.cho_solve(chol,B)
    inv_res = frob(G@invG-eye)
    k_res = frob(G@KB-B)
    # G(true)>=I. Residual bounds for the two solves account also for the
    # correction to the exactly projected B and for assembling the Gram.
    gram_error = (2*bnorm+b_error)*b_error + gamma_n(8*n+16)*(1+bnorm**2)
    solve_arithmetic = gamma_n(8*n+16)*(1+frob(G))*(1+frob(invG)+frob(KB))
    r_reconstruction_error = inv_res+gram_error*frob(invG)+3*x_error+solve_arithmetic
    k_reconstruction_error = k_res+gram_error*frob(KB)+b_error+solve_arithmetic
    block_error = 4*(r_reconstruction_error+k_reconstruction_error)
    # For Wick arithmetic only: compare the actual floating R,K inputs with
    # the pure chart in exactly the same energy labels. Neither covariance
    # assembly nor the physical-model displacement needs charging here: R,K
    # are the stored input being compared. The distinct t below covers the
    # chart-to-aligned-physical-source displacement for a positive observable.
    normal_diff, anomalous_diff = frob(R-RB), frob(K-KB)
    measure_roundoff = 64*n*EPS*(1 + frob(R)+frob(RB)+frob(K)+frob(KB))
    wick_F = math.sqrt(2)*math.hypot(normal_diff+r_reconstruction_error+measure_roundoff,
                                     anomalous_diff+k_reconstruction_error+measure_roundoff)
    diag.update(direct_normal_minus_chart_F=normal_diff, direct_anomalous_minus_chart_F=anomalous_diff,
                normal_reconstruction_arithmetic_F=r_reconstruction_error,
                anomalous_reconstruction_arithmetic_F=k_reconstruction_error,
                direct_chart_comparison_roundoff_F=measure_roundoff,
                wick_input_vs_pure_chart_F=wick_F*(1+16*EPS))

    # In an exact canonical frame the full covariance Frobenius difference
    # is twice the full Nambu density difference, i.e.2sqrt2*hypot(dR,dK).
    # Reality/skew symmetry fixes the lower blocks. This is unitarily invariant
    # and avoids reconstructing a site matrix through two full dense products.
    # Both direct contractions were computed in represented V, so include their
    # arithmetic, symmetrization and displacement from its exact polar frame.
    # Completing only the upper raw Nambu blocks to the exact PH structure
    # can enlarge their full-frame Frobenius discrepancy by sqrt(2).
    polar_cov_error=2*math.sqrt(2)*(2+polar_tail)*polar_tail*(1+.5*initial_covariance_F)
    direct_arithmetic=2*math.sqrt(2)*(2*contract_roundoff+.5*(r_defect+k_defect))
    recons=2*math.sqrt(2)*math.hypot(normal_diff,anomalous_diff)
    recon_error=2*wick_F+polar_cov_error+direct_arithmetic
    model_cov_error=initial_covariance_F+final_cov
    trace_radius=min(1.,(model_cov_error+recon_error)*(1+64*EPS)/math.sqrt(8))
    diag.update(covariance_reconstruction_F=recons,gram_formation_allowance=gram_error,
        G_inverse_residual_F=inv_res,G_KB_residual_F=k_res,
        reconstruction_block_calculation_F=max(0.,2*wick_F-recons),
        reconstruction_polar_tail_F=polar_cov_error,reconstruction_site_arithmetic_F=direct_arithmetic,
        reconstruction_total_F=recon_error,aligned_model_covariance_F=model_cov_error,
        trace_distance_allowance=trace_radius,normalization_relative_error=relative_norm_error,
        energy_allowance_terms={"certified_noncommuting_one_body_operator":energy_error},
        error_scope="Frobenius occupied-covariance certificate plus canonical-coordinate reconstruction; one-body operator error includes within-positive mixing")
    return {"singular":singular,"R":R,"K":K,"occ_values":occ_values,"B":B,"x":x,"Z_b":Z_b,"t":trace_radius,"e_h":energy_error,
            "chart_mass_roundoff":mass_roundoff,
            "normalization_relative_error":relative_norm_error,"diagnostics":diag}
