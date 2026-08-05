"""
Local pressure loss correlations for pipe fittings.

Implements Idelchik curved bend correlations across all flow regimes
(laminar, transition, turbulent) with smooth regime blending.

Functions:
    - compute_zeta_LOC_curved: Local resistance coefficient
    - compute_cf_fri: Roughness correction factor
    - compute_zeta_tot_curved: Total loss coefficient
    - compute_pressure_drop: Pressure drop from ζ_tot
    - compute_re_lam_leave: Laminar–transition boundary
    - compute_reynolds_correction_curved: Reynolds correction for curvature
    - compute_coefficients_curved: Idelchik coefficients (A1, A2, B1, C1)
"""

from labothappy.toolbox.helper_function.step_smoother import stepsmoother
from correlations.properties.dimensionless import compute_reynolds
from correlations.pressure_drop.straight_pipe_DP import friction_factor_swamee_jain

import math

PI = math.pi
EPS = 1e-12

RE_MIN = 1.0
RE_LAM_MAX = 6.5e3
RE_TURB_MIN = 4e4
RE_TURB_MAX = 3e5
RE_TURB_CONST = 1e6

def compute_re_lam_leave(K, d_hyd):
    """Compute Reynolds number at which laminar regime ends.
    sPECIFICALLY FOR CURVED BENDS, TAKES INTO ACCOUNT ROUGHNESS EFFECTS."""
    k_rel = K / max(EPS, d_hyd) # Relative roughness
    exponent = 0.0065 / max(k_rel, 0.007)
    return min(RE_LAM_MAX, max(1e2, 754.0 * math.exp(exponent)))

def compute_coefficients_curved(d_hyd, R0, delta_rad):
    """Compute Idelchik coefficients for curved bends."""
    MIN = EPS
    delta_deg = abs(delta_rad) * 180 / PI
    
    R0_d_hyd = R0 / max(MIN, d_hyd) # Ratio of bend radius to hydraulic diameter
    
    sin_half_delta = math.sin(PI / 180 * delta_deg / 2)
    A1 = max(MIN, 0.95 * sin_half_delta**2 + 2.05 * sin_half_delta**4)
    
    if R0_d_hyd > 2.0:
        A2 = 600.0
    elif R0_d_hyd > 0.55:
        if R0_d_hyd > 1.0:
            A2 = 1000.0
        elif R0_d_hyd > 0.7:
            A2 = 3000.0
        else:
            A2 = 6000.0
    else:
        A2 = 4000.0
    
    if R0_d_hyd >= 1.0:
        B1 = 0.21 * (R0_d_hyd)**(-0.5)
    else:
        B1 = 0.21 * (R0_d_hyd)**(-2.5)
    B1 = max(MIN, B1)
    
    C1 = 1.0
    
    return {
        'A1': A1,
        'A2': A2,
        'B1': B1,
        'C1': C1,
        'R0_d_hyd': R0_d_hyd
    }

def compute_reynolds_correction_curved(Re, R0_d_hyd):
    """Compute Reynolds-number correction factor for curved bends in turbulent regime."""
    if R0_d_hyd > 0.7:
        return 11.5 / (max(EPS, Re)**0.19)
    elif R0_d_hyd >= 0.55:
        return 5.45 / (max(EPS, Re)**0.131)
    else:
        return 1.0 + 4400.0 / max(EPS, Re)


def compute_zeta_LOC_curved(d_hyd, R0, delta_rad, K, rho, mu, m_dot):
    """Compute local resistance coefficient ζ_LOC for curved bends across all flow regimes."""
    MIN = EPS
    
    A_cross = PI * d_hyd**2 / 4.0
    v = abs(m_dot) / max(MIN, rho * A_cross)
    Re = compute_reynolds(rho, v, d_hyd, mu)
    
    coeffs = compute_coefficients_curved(d_hyd, R0, delta_rad)
    A1 = coeffs['A1']
    A2 = coeffs['A2']
    B1 = coeffs['B1']
    C1 = coeffs['C1']
    R0_d_hyd = coeffs['R0_d_hyd']
    
    zeta_base = A1 * B1 * C1
    Re_lam_leave = compute_re_lam_leave(K, d_hyd)
    Re_turb_min = 4e4
    Re_turb_max = 3e5
    
    zeta_lam = A2 / max(MIN, Re) + zeta_base
    
    k_Re = compute_reynolds_correction_curved(Re, R0_d_hyd) # Reynolds correction factor for turbulent regime
    zeta_turb_inter = max(MIN, k_Re * zeta_base)
    
    zeta_turb_high = zeta_base
    
    s_lam = stepsmoother(Re, Re_lam_leave, Re_turb_min)
    zeta_mid = (1.0 - s_lam) * zeta_lam + s_lam * zeta_turb_inter
    
    s_turb = stepsmoother(Re, Re_turb_min, Re_turb_max)
    zeta_LOC = (1.0 - s_turb) * zeta_mid + s_turb * zeta_turb_high
    
    return zeta_LOC

def compute_cf_fri(zeta_LOC, d_hyd, R0, delta_rad, K, rho, eta, m_dot):
    """Compute surface roughness correction factor CF_fri."""
    MIN = EPS
    
    A_cross = PI * d_hyd**2 / 4.0
    v = abs(m_dot) / max(MIN, rho * A_cross)
    Re = compute_reynolds(rho, v, d_hyd, eta)
    
    L = abs(delta_rad) * R0
    f_fri = friction_factor_swamee_jain(K, d_hyd, Re) # Friction factor for straight pipe, used for roughness correction
    Re_lam_leave = compute_re_lam_leave(K, d_hyd)
    
    raw_effect = min(1.4, (f_fri * L / max(MIN, d_hyd)) / max(MIN, zeta_LOC))
    s = stepsmoother(Re, RE_LAM_MAX, Re_lam_leave)
    CF_fri = 1.0 + s * raw_effect
    
    return CF_fri

def compute_zeta_tot_curved(d_hyd, R0, delta_rad, K, rho, mu, m_dot):
    """Compute total pressure loss coefficient ζ_tot for curved bends."""
    zeta_LOC = compute_zeta_LOC_curved(d_hyd, R0, delta_rad, K, rho, mu, m_dot)
    CF_fri = compute_cf_fri(zeta_LOC, d_hyd, R0, delta_rad, K, rho, mu, m_dot)
    zeta_tot = max(EPS, CF_fri) * zeta_LOC
    return zeta_tot


def compute_pressure_drop(zeta_tot, rho, v):
    """Compute pressure drop from total loss coefficient."""
    return zeta_tot * 0.5 * rho * v**2
