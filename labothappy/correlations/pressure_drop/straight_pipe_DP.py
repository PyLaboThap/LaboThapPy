"""
Pressure drop correlations for single-phase and two-phase flow in straight pipes.

This module implements:
- Single-phase friction factor correlations (Churchill, Swamee-Jain)
- Two-phase friction multipliers (Friedel)
- Two-phase pressure drop components (friction, acceleration, gravity)
- Total pressure drop calculation

References
----------
Churchill, S.W. (1977). Friction-factor equation spans all fluid-flow regimes.
    Chemical Engineering, 84(24), 91-92.

Swamee, P.K., Jain, A.K. (1976). Explicit equations for pipe flow problems.
    Journal of the Hydraulics Division, ASCE, 102(HY5), 657-664.

Friedel, L. (1979). Improved friction pressure drop correlations for horizontal
    and vertical two-phase pipe flow. European Two-Phase Flow Group Meeting,
    Ispra, Italy.

Darcy, H. (1857). Recherches expérimentales sur le mouvement de l'eau dans
    les tuyaux. Paris: Mallet-Bachelier.

Moody, L.F. (1944). Friction factors for pipe flow. Transactions of the ASME,
    66(8), 671-684.
"""

import math
import CoolProp.CoolProp as CP

from labothappy.correlations.properties.dimensionless import (
    compute_reynolds,
    compute_froude,
    compute_weber,
)
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction
from labothappy.correlations.properties.two_phase import compute_two_phase_density

EPS = 1e-12
PI = math.pi
G_GRAVITY = 9.81  # m/s²


# ============================================================================
# FRICTION FACTOR CORRELATIONS
# ============================================================================


def friction_factor_churchill(K, d_hyd, Re):
    """
    Darcy friction factor using Churchill correlation.

    Single formula valid for all Reynolds numbers (laminar, transition, turbulent)
    and accounts for surface roughness. Smooth transition at Re ≈ 2300 without
    branching logic.

    Formula (SOURCE: Churchill 1977)
    --------------------------------
    f = 8 · [(8/Re)^12 + 1/(A+B)^1.5]^(1/12)

    where:
        A = [2.457 · ln(1/((7/Re)^0.9 + 0.27·K/d_hyd))]^16
        B = (37530/Re)^16

    Parameters
    ----------
    K : float
        Absolute roughness [m]
    d_hyd : float
        Hydraulic diameter [m]
    Re : float
        Reynolds number [-]

    Returns
    -------
    f : float
        Darcy friction factor [-]

    Notes
    -----
    - Valid for Re > 0 (all flow regimes)
    - Valid for any relative roughness K/d_hyd
    - Smooth transition from laminar to turbulent (no discontinuities)

    References
    ----------
    Churchill, S.W. (1977). Friction-factor equation spans all fluid-flow regimes.
        Chemical Engineering, 84(24), 91-92.
    """
    Re_safe = max(EPS, Re)

    # Laminar term
    laminar_term = (8.0 / Re_safe) ** 12

    # Turbulent term with roughness correction
    k_rel = K / max(EPS, d_hyd)
    ln_arg = 1.0 / ((7.0 / Re_safe) ** 0.9 + 0.27 * k_rel)
    A = (2.457 * math.log(max(EPS, ln_arg))) ** 16
    B = (37530.0 / Re_safe) ** 16

    turbulent_term = 1.0 / (A + B) ** 1.5

    return 8.0 * (laminar_term + turbulent_term) ** (1.0 / 12.0)


def friction_factor_swamee_jain(K, d_hyd, Re):
    """
    Darcy friction factor using Swamee-Jain approximation.

    Explicit approximation of Colebrook-White correlation. Valid for turbulent
    flow (Re ≥ 4000, recommended Re ≥ 5000) and any relative roughness.

    Formula (SOURCE: Swamee & Jain 1976)
    -----------------------------------
    f = 0.25 / [log₁₀(K/(3.7·d_hyd) + 5.74/Re^0.9)]²

    Parameters
    ----------
    K : float
        Absolute roughness [m]
    d_hyd : float
        Hydraulic diameter [m]
    Re : float
        Reynolds number [-]

    Returns
    -------
    f : float
        Darcy friction factor [-]

    Notes
    -----
    - Valid for Re ≥ 5000 (turbulent, recommended)
    - Valid for Re ≥ 4000 (transition, acceptable)
    - Accuracy: ±2% vs. Colebrook-White
    - For laminar flow (Re < 2300), use f = 64/Re
    - For general use across all regimes, prefer Churchill correlation

    References
    ----------
    Swamee, P.K., Jain, A.K. (1976). Explicit equations for pipe flow problems.
        Journal of the Hydraulics Division, ASCE, 102(HY5), 657-664.
    """
    Re_safe = max(EPS, Re)
    d_hyd_safe = max(EPS, d_hyd)

    k_rel = K / d_hyd_safe
    log_arg = k_rel / 3.7 + 5.74 / (Re_safe ** 0.9)

    return 0.25 / (math.log10(max(EPS, log_arg)) ** 2)


# ============================================================================
# TWO-PHASE FRICTION MULTIPLIER
# ============================================================================

def friedel_multiplier(m_dot, x, D, rho_l, rho_g, mu_l, mu_g, sigma, K=0.0):
    """
    Compute Friedel two-phase friction multiplier Φ_l².

    The Friedel correlation is a separated flow model that uses a friction
    multiplier to correct single-phase friction factors for two-phase flow.
    Accounts for quality, density ratio, viscosity ratio, and interfacial
    tension effects.

    Formula (SOURCE: Friedel 1979)
    ----------------------------
    Φ_l² = E + (3.24 · F · H) / (Fr^0.0454 · We^0.035)

    where:
        E = (1-x)² + x² · (ρ_l·f_g / ρ_g·f_l)
        F = x^0.78 · (1-x)^0.224
        H = (ρ_l/ρ_g)^0.91 · (μ_g/μ_l)^0.19 · (1 - μ_g/μ_l)^0.7
        Fr = G² / (g·D·ρ_h²)  [Froude number]
        We = G²·D / (σ·ρ_h)   [Weber number]
        G = m_dot / A_cross   [mass velocity, kg/(m²·s)]

    Parameters
    ----------
    m_dot : float
        Mass flow rate [kg/s]
    x : float
        Vapor quality (dryness fraction) [0, 1]
    D : float
        Pipe inner diameter [m]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    mu_l : float
        Liquid dynamic viscosity [Pa·s]
    mu_g : float
        Vapor dynamic viscosity [Pa·s]
    sigma : float
        Surface tension [N/m]
    K : float, optional
        Absolute roughness [m] (default: 0.0, smooth pipe)

    Returns
    -------
    Phi_l2 : float
        Two-phase friction multiplier Φ_l² [-]

    Notes
    -----
    - Uses homogeneous void fraction model (slip ratio = 1)
    - Valid for separated flow regimes
    - Froude and Weber numbers computed with homogeneous density
    - Multiplier is typically 2-10 for two-phase flow

    References
    ----------
    Friedel, L. (1979). Improved friction pressure drop correlations for
        horizontal and vertical two-phase pipe flow. European Two-Phase Flow
        Group Meeting, Ispra, Italy.
    """
    # Safeguard quality
    x = max(EPS, min(1.0 - EPS, x))

    # Cross-sectional area and mass velocity
    A_cross = PI * D ** 2 / 4.0
    G = m_dot / A_cross  # kg/(m²·s)

    # ====================================================================
    # Homogeneous two-phase properties
    # ====================================================================
    alpha_h = compute_void_fraction(x, rho_l, rho_g, slip_model=None)
    rho_h = compute_two_phase_density(x, rho_l, rho_g, alpha=alpha_h)

    # ====================================================================
    # Single-phase friction factors
    # ====================================================================
    # Liquid (at total mass flow)
    v_l = m_dot / (rho_l * A_cross)
    Re_l = compute_reynolds(rho_l, v_l, D, mu_l)
    f_l = friction_factor_churchill(K, D, Re_l)

    # Vapor (at total mass flow)
    v_g = m_dot / (rho_g * A_cross)
    Re_g = compute_reynolds(rho_g, v_g, D, mu_g)
    f_g = friction_factor_churchill(K, D, Re_g)

    # ====================================================================
    # Dimensionless numbers (using homogeneous density)
    # ====================================================================
    Fr = compute_froude(G, D, rho_h)
    We = compute_weber(G, D, sigma, rho_h)

    # ====================================================================
    # Friedel multiplier components
    # ====================================================================
    # F factor: F = x^0.78 · (1-x)^0.224
    F = (x ** 0.78) * ((1.0 - x) ** 0.224)

    # H factor: H = (ρ_l/ρ_g)^0.91 · (μ_g/μ_l)^0.19 · (1 - μ_g/μ_l)^0.7
    rho_ratio = rho_l / rho_g
    mu_ratio = mu_g / mu_l
    H = (rho_ratio ** 0.91) * (mu_ratio ** 0.19) * ((1.0 - mu_ratio) ** 0.7)

    # E factor: E = (1-x)² + x² · (ρ_l·f_g / ρ_g·f_l)
    E = (1.0 - x) ** 2 + x ** 2 * (rho_l * f_g) / (rho_g * f_l)

    # ====================================================================
    # Two-phase multiplier Φ_l²
    # ====================================================================
    denominator = (Fr ** 0.0454) * (We ** 0.035)
    Phi_l2 = E + (3.24 * F * H) / max(EPS, denominator)

    return Phi_l2


# ============================================================================
# SINGLE-PHASE PRESSURE DROP
# ============================================================================


def pressure_drop_single_phase(L, d_hyd, rho, v, K, mu, theta=0.0):
    """
    Compute total pressure drop in a straight pipe (single-phase flow).

    Combines frictional (Darcy-Weisbach) and gravitational pressure drops.
    Uses Churchill correlation for friction factor, valid across all Reynolds
    numbers and surface roughness conditions.

    Formula
    -------
    ΔP_friction = f · (L/D) · (ρ·v²/2)
    ΔP_gravity = ρ·g·L·sin(θ)
    ΔP_total = ΔP_friction + ΔP_gravity

    where f is computed using Churchill (1977) correlation.

    Parameters
    ----------
    L : float
        Pipe length [m]
    d_hyd : float
        Hydraulic diameter [m]
    rho : float
        Fluid density [kg/m³].
        Use average density if pressure drop > 5% of inlet pressure.
    v : float
        Mean flow velocity [m/s]
    K : float
        Absolute surface roughness [m]
    mu : float
        Dynamic viscosity [Pa·s]
    theta : float, optional
        Pipe inclination angle from horizontal [degrees] (default: 0)

    Returns
    -------
    deltaP_total : float
        Total pressure drop [Pa]
    f : float
        Darcy friction factor [-] (for diagnostics)
    Re : float
        Reynolds number [-] (for diagnostics)

    Notes
    -----
    - For large pressure drops (ΔP > 0.05·P_inlet), use mean density
    - Positive theta = downward flow (increases ΔP)
    - Negative theta = upward flow (decreases ΔP)

    References
    ----------
    Darcy, H. (1857). Recherches expérimentales sur le mouvement de l'eau dans
        les tuyaux. Paris: Mallet-Bachelier.

    Moody, L.F. (1944). Friction factors for pipe flow. Transactions of the
        ASME, 66(8), 671-684.
    """
    Re = compute_reynolds(rho, v, d_hyd, mu)
    f = friction_factor_churchill(K, d_hyd, Re)

    # Frictional pressure drop (Darcy-Weisbach)
    dP_friction = f * (L / d_hyd) * (rho * v ** 2 / 2.0)

    # Gravitational pressure drop (hydrostatic)
    dP_gravity = G_GRAVITY * rho * L * math.sin(math.radians(theta))

    # Total pressure drop
    dP_total = dP_friction + dP_gravity

    return dP_total, f, Re


# ============================================================================
# TWO-PHASE PRESSURE DROP COMPONENTS
# ============================================================================


def pressure_drop_frictional_two_phase(m_dot, D, L, rho_l, rho_g, x, 
                                        mu_l, mu_g, sigma, K=0.0):
    """
    Compute frictional pressure drop in two-phase flow.

    Uses Friedel correlation: ΔP_f = ΔP_lo · Φ_l²

    where ΔP_lo is single-phase pressure drop (liquid only, at total mass flow)
    and Φ_l² is the Friedel two-phase multiplier.

    Parameters
    ----------
    m_dot : float
        Mass flow rate [kg/s]
    D : float
        Pipe inner diameter [m]
    L : float
        Pipe length [m]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    x : float
        Vapor quality [0, 1]
    mu_l : float
        Liquid dynamic viscosity [Pa·s]
    mu_g : float
        Vapor dynamic viscosity [Pa·s]
    sigma : float
        Surface tension [N/m]
    K : float, optional
        Absolute roughness [m] (default: 0.0, smooth pipe)

    Returns
    -------
    dP_friction : float
        Frictional pressure drop [Pa]

    References
    ----------
    Friedel, L. (1979). Improved friction pressure drop correlations for
        horizontal and vertical two-phase pipe flow. European Two-Phase Flow
        Group Meeting, Ispra, Italy.
    """
    A_cross = PI * D ** 2 / 4.0

    # Single-phase pressure drop (liquid only, at total mass flow)
    v_l = m_dot / (rho_l * A_cross)
    Re_l = compute_reynolds(rho_l, v_l, D, mu_l)
    f_l = friction_factor_churchill(K, D, Re_l)
    dP_lo = f_l * (L / D) * (rho_l * v_l ** 2 / 2.0)

    # Two-phase multiplier
    Phi_l2 = friedel_multiplier(m_dot, x, D, rho_l, rho_g, mu_l, mu_g, sigma, K)

    # Frictional pressure drop
    dP_friction = dP_lo * Phi_l2

    return dP_friction


def pressure_drop_acceleration_two_phase(m_dot, D, rho_l, rho_g, x_inlet, 
                                          x_outlet, slip_model=None):
    """
    Compute acceleration pressure drop in two-phase flow.

    Uses separated flow model with selectable slip ratio.

    Formula (SOURCE: ACHP methodology)
    ---------------------------------
    ΔP_A = G² · [f(x_outlet) - f(x_inlet)]

    where:
        f(x) = x²·v_g/α + (1-x)²·v_l/(1-α)
        α = void fraction (depends on slip model)
        G = m_dot / A_cross (mass velocity)

    Parameters
    ----------
    m_dot : float
        Mass flow rate [kg/s]
    D : float
        Pipe inner diameter [m]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    x_inlet : float
        Vapor quality at inlet [0, 1]
    x_outlet : float
        Vapor quality at outlet [0, 1]
    slip_model : str, optional
        Slip ratio model: None/'Homogeneous' (default), 'Zivi'

    Returns
    -------
    dP_acceleration : float
        Acceleration pressure drop [Pa]

    Notes
    -----
    - Positive when quality increases (vapor accelerates, opposes flow)
    - Negative when quality decreases (flow decelerates, aids flow)
    - Typically 5-15% of total ΔP for evaporators

    References
    ----------
    ACHP Documentation: https://achp.readthedocs.io/en/latest/ACHPComponents/FluidCorrelations.html
    """

    def f_acceleration(x, rho_l, rho_g, slip_model):
        """Helper: f(x) = x²·v_g/α + (1-x)²·v_l/(1-α)"""
        if abs(x) < EPS:
            return 1.0 / rho_l
        elif abs(1.0 - x) < EPS:
            return 1.0 / rho_g
        else:
            alpha = compute_void_fraction(x, rho_l, rho_g, slip_model=slip_model)
            v_l = 1.0 / rho_l
            v_g = 1.0 / rho_g

            term1 = x ** 2 * v_g / max(EPS, alpha)
            term2 = (1.0 - x) ** 2 * v_l / max(EPS, 1.0 - alpha)

            return term1 + term2

    A_cross = PI * D ** 2 / 4.0
    G = m_dot / A_cross  # kg/(m²·s)

    f_inlet = f_acceleration(x_inlet, rho_l, rho_g, slip_model)
    f_outlet = f_acceleration(x_outlet, rho_l, rho_g, slip_model)

    dP_acceleration = G ** 2 * (f_inlet - f_outlet)

    return dP_acceleration


def pressure_drop_gravity_two_phase(m_dot, D, L, rho_l, rho_g, x_inlet, 
                                     x_outlet, theta):
    """
    Compute gravitational pressure drop in two-phase flow.

    Uses mean quality to estimate mean density.

    Formula
    -------
    ΔP_gravity = ρ_mean · g · L · sin(θ)

    Parameters
    ----------
    m_dot : float
        Mass flow rate [kg/s]
    D : float
        Pipe inner diameter [m]
    L : float
        Pipe length [m]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    x_inlet : float
        Vapor quality at inlet [0, 1]
    x_outlet : float
        Vapor quality at outlet [0, 1]
    theta : float
        Pipe inclination angle from horizontal [degrees]

    Returns
    -------
    dP_gravity : float
        Gravitational pressure drop [Pa]

    Notes
    -----
    - Positive theta = upward flow (increases ΔP)
    - Negative theta = downward flow (decreases ΔP)
    - For horizontal pipes (theta=0), ΔP_gravity = 0
    """
    x_mean = (x_inlet + x_outlet) / 2.0
    rho_mean = compute_two_phase_density(x_mean, rho_l, rho_g, alpha=None)

    dP_gravity = G_GRAVITY * rho_mean * L * math.sin(math.radians(theta))

    return dP_gravity


# ============================================================================
# TOTAL TWO-PHASE PRESSURE DROP
# ============================================================================


def pressure_drop_two_phase(m_dot, D, L, rho_l, rho_g, x_inlet, AS,
                            mu_l, mu_g, sigma, K=0.0, theta=0.0, 
                            slip_model=None, include_acceleration=True, 
                            include_gravity=True):
    """
    Compute total two-phase pressure drop.

    Combines frictional, acceleration, and gravitational components.

    Formula
    -------
    ΔP_total = ΔP_friction + ΔP_acceleration + ΔP_gravity

    Parameters
    ----------
    m_dot : float
        Mass flow rate [kg/s]
    D : float
        Pipe inner diameter [m]
    L : float
        Pipe length [m]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    x_inlet : float
        Vapor quality at inlet [0, 1]
    x_outlet : float
        Vapor quality at outlet [0, 1]
    mu_l : float
        Liquid dynamic viscosity [Pa·s]
    mu_g : float
        Vapor dynamic viscosity [Pa·s]
    sigma : float
        Surface tension [N/m]
    K : float, optional
        Absolute roughness [m] (default: 0.0, smooth pipe)
    theta : float, optional
        Pipe inclination angle from horizontal [degrees] (default: 0)
    slip_model : str, optional
        Slip ratio model for void fraction: None/'Homogeneous' (default), 'Zivi'
    include_acceleration : bool, optional
        Include acceleration pressure drop (default: True)
    include_gravity : bool, optional
        Include gravitational pressure drop (default: True)

    Returns
    -------
    dict with keys:
        'dP_total' : float
            Total pressure drop [Pa]
        'dP_friction' : float
            Frictional component [Pa]
        'dP_acceleration' : float
            Acceleration component [Pa]
        'dP_gravity' : float
            Gravitational component [Pa]

    Notes
    -----
    - Friction typically 85-95% of total ΔP
    - Acceleration typically 5-15% of total ΔP
    - Gravity important for vertical pipes only

    References
    ----------
    Friedel, L. (1979). Improved friction pressure drop correlations for
        horizontal and vertical two-phase pipe flow. European Two-Phase Flow
        Group Meeting, Ispra, Italy.

    ACHP Documentation: https://achp.readthedocs.io/en/latest/ACHPComponents/FluidCorrelations.html
    """

    # ========== STEP 1: Compute friction ΔP ==========

    # Frictional pressure drop
    dP_friction = pressure_drop_frictional_two_phase(
        m_dot, D, L, rho_l, rho_g, x_inlet,
        mu_l, mu_g, sigma, K
    )

    # ========== STEP 2: Estimate outlet pressure & compute x_outlet ==========
    p_inlet = AS.p()  # Inlet pressure
    h = AS.hmass()   # Specific enthalpy (assumed constant)
    p_outlet = p_inlet - dP_friction
    
    # Get outlet quality at (p_outlet, h)
    AS.update(CP.HmassP_INPUTS, p_outlet, h)
    x_outlet = AS.Q()  # Quality at outlet

    # Acceleration pressure drop
    if include_acceleration:
        dP_acceleration = pressure_drop_acceleration_two_phase(
            m_dot, D, rho_l, rho_g, x_inlet, x_outlet, slip_model=slip_model
        )
    else:
        dP_acceleration = 0.0

    # Gravitational pressure drop
    if include_gravity:
        dP_gravity = pressure_drop_gravity_two_phase(
            m_dot, D, L, rho_l, rho_g, x_inlet, x_outlet, theta
        )
    else:
        dP_gravity = 0.0
    # Total pressure drop
    dP_total = dP_friction + dP_acceleration + dP_gravity

    return dP_total