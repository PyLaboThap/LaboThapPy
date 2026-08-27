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

from labothappy.correlations.void_fraction.void_fraction import (
    compute_void_fraction,
    void_fraction_homogeneous,
    void_fraction_zivi,
)
from labothappy.correlations.properties.two_phase import compute_two_phase_density, get_saturated_phase_properties

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

    f = 8.0 * (laminar_term + turbulent_term) ** (1.0 / 12.0)

    return f


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

def friedel_multiplier(AS, m_dot, D, K=0.0):
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
    AS : CoolProp.AbstractState
        Two-phase fluid state (quality and pressure already set). Liquid
        and vapor properties (ρ_l, ρ_v, μ_l, μ_v, σ) are derived from it via
        `get_saturated_phase_properties`.
    m_dot : float
        Mass flow rate [kg/s]
    D : float
        Pipe inner diameter [m]
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
    props = get_saturated_phase_properties(AS)
    x = props["x"]
    rho_l = props["rho_l"]
    rho_v = props["rho_v"]
    mu_l = props["mu_l"]
    mu_v = props["mu_v"]
    sigma = props["sigma"]

    # Cross-sectional area and mass velocity
    A_cross = PI * D ** 2 / 4.0
    G = m_dot / A_cross  # kg/(m²·s)

    # ====================================================================
    # Homogeneous two-phase properties
    # ====================================================================
    rho_h = compute_two_phase_density(x, rho_l, rho_v)

    # ====================================================================
    # Single-phase friction factors
    # ====================================================================
    # Liquid (at total mass flow)
    v_l = m_dot / (rho_l * A_cross)
    Re_l = compute_reynolds(D, mu_l, rho_l, v_l)
    f_l = friction_factor_churchill(K, D, Re_l)

    # Vapor (at total mass flow)
    v_v = m_dot / (rho_v * A_cross)
    Re_v = compute_reynolds(D, mu_v, rho_v, v_v)
    f_v = friction_factor_churchill(K, D, Re_v)

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
    rho_ratio = rho_l / rho_v
    mu_ratio = mu_v / mu_l
    H = (rho_ratio ** 0.91) * (mu_ratio ** 0.19) * ((1.0 - mu_ratio) ** 0.7)

    # E factor: E = (1-x)² + x² · (ρ_l·f_g / ρ_g·f_l)
    E = (1.0 - x) ** 2 + x ** 2 * (rho_l * f_v) / (rho_v * f_l)

    # ====================================================================
    # Two-phase multiplier Φ_l²
    # ====================================================================
    denominator = (Fr ** 0.0454) * (We ** 0.035)
    Phi_l2 = E + (3.24 * F * H) / max(EPS, denominator)

    return Phi_l2


# ============================================================================
# SINGLE-PHASE PRESSURE DROP
# ============================================================================


def pressure_drop_pipe_single_phase(AS, params, m_dot, P_su=None, T_su=None, correlation='Churchill'):
    """
    Compute total pressure drop in a straight pipe (single-phase flow).

    Combines frictional (Darcy-Weisbach) and gravitational pressure drops.
    The friction factor is computed with the correlation selected via
    `correlation` ('Churchill' or 'Swamee-Jain').

    Formula
    -------
    ΔP_friction = f · (L/D) · (ρ·v²/2)
    ΔP_gravity = ρ·g·L·sin(θ)
    ΔP_total = ΔP_friction + ΔP_gravity

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Single-phase fluid state. If `P_su`/`T_su` are not given, AS must
        already have its state set (e.g. via HmassP_INPUTS) by the caller.
    params : dict
        Pipe geometry: 'D' (or 'd_hyd') [m], 'L' [m], and optionally
        'K' (absolute roughness [m], default 0.0) and 'theta' (inclination
        from horizontal [degrees], default 0.0).
    m_dot : float
        Mass flow rate [kg/s]
    P_su : float, optional
        Inlet pressure [Pa]. If given together with `T_su`, AS is updated
        to this (P, T) state before computing properties.
    T_su : float, optional
        Inlet temperature [K]. See `P_su`.
    correlation : str, optional
        Friction factor correlation: 'Churchill' (default, valid across all
        Reynolds numbers) or 'Swamee-Jain' (explicit, turbulent flow only).

    Returns
    -------
    dP_total : float
        Total pressure drop [Pa]

    Notes
    -----
    - For large pressure drops (ΔP > 0.05·P_inlet), consider using mean density
    - Positive theta = downward flow (increases ΔP)
    - Negative theta = upward flow (decreases ΔP)

    References
    ----------
    Darcy, H. (1857). Recherches expérimentales sur le mouvement de l'eau dans
        les tuyaux. Paris: Mallet-Bachelier.

    Moody, L.F. (1944). Friction factors for pipe flow. Transactions of the
        ASME, 66(8), 671-684.
    """
    if (P_su is None) != (T_su is None):
        raise ValueError("P_su and T_su must be given together, or not at all.")
    if P_su is not None and T_su is not None:
        AS.update(CP.PT_INPUTS, P_su, T_su)

    d_hyd = params.get('d_hyd', params.get('D'))
    A_cross = PI * d_hyd ** 2 / 4.0
    L = params['L']
    theta = params.get('theta', 0.0)
    K = params.get('K', 0.0)

    mu = AS.viscosity()
    rho = AS.rhomass()

    v = m_dot / (rho * A_cross)  # Mean velocity [m/s]

    Re = compute_reynolds(d_hyd, mu, rho, v)

    if correlation == 'Churchill':
        f = friction_factor_churchill(K, d_hyd, Re)
    elif correlation == 'Swamee-Jain':
        f = friction_factor_swamee_jain(K, d_hyd, Re)
    else:
        raise ValueError(f"Unknown single-phase friction correlation: {correlation!r}")

    # Frictional pressure drop (Darcy-Weisbach)
    dP_friction = f * (L / d_hyd) * (rho * v ** 2 / 2.0)

    # Gravitational pressure drop (hydrostatic)
    dP_gravity = G_GRAVITY * rho * L * math.sin(math.radians(theta))

    # Total pressure drop
    dP_total = dP_friction + dP_gravity

    return dP_total


# ============================================================================
# TWO-PHASE PRESSURE DROP COMPONENTS
# ============================================================================

def pressure_drop_pipe_frictional_two_phase(AS, params, m_dot, correlation='Friedel'):
    """
    Compute frictional pressure drop in two-phase flow.

    Uses ΔP_f = ΔP_lo · Φ_l², where ΔP_lo is the single-phase pressure drop
    (liquid only, at total mass flow) and Φ_l² is the two-phase friction
    multiplier from the selected `correlation`.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Two-phase fluid state (quality and pressure already set). Liquid
        properties are derived from it via `get_saturated_phase_properties`.
    params : dict
        Pipe geometry: 'D' (or 'd_hyd') [m], 'L' [m], and optionally
        'K' (absolute roughness [m], default 0.0). Also passed through to
        `friedel_multiplier` / `compute_void_fraction`.
    m_dot : float
        Mass flow rate [kg/s]
    correlation : str, optional
        Two-phase friction multiplier correlation. Only 'Friedel' (default)
        is currently implemented.

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
    props = get_saturated_phase_properties(AS)
    rho_l = props["rho_l"]
    mu_l = props["mu_l"]

    D = params.get('D', params.get('d_hyd'))
    A_cross = PI * D ** 2 / 4.0
    L = params['L']
    K = params.get('K', 0.0)

    # Single-phase pressure drop (liquid only, at total mass flow)
    v_l = m_dot / (rho_l * A_cross)
    Re_l = compute_reynolds(D, mu_l, rho_l, v_l)
    f_l = friction_factor_churchill(K, D, Re_l)
    dP_l = f_l * (L / D) * (rho_l * v_l ** 2 / 2.0)

    # Two-phase multiplier
    if correlation == 'Friedel':
        Phi_l2 = friedel_multiplier(AS, m_dot, D, K)
    else:
        raise ValueError(f"Unknown two-phase friction correlation: {correlation!r}")

    # Frictional pressure drop
    dP_friction = dP_l * Phi_l2

    return dP_friction


def pressure_drop_pipe_acceleration_two_phase(m_dot, D, rho_l, rho_g, x_inlet, 
                                          x_outlet, void_fraction_model=None):
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
    void_fraction_model : str, optional
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

    def f_acceleration(x, rho_l, rho_g, void_fraction_model):
        """Helper: f(x) = x²·v_g/α + (1-x)²·v_l/(1-α)"""
        if abs(x) < EPS:
            return 1.0 / rho_l
        elif abs(1.0 - x) < EPS:
            return 1.0 / rho_g
        else:
            if void_fraction_model in (None, 'Homogeneous'):
                alpha = void_fraction_homogeneous(rho_l, rho_g, x)
            elif void_fraction_model == 'Zivi':
                alpha = void_fraction_zivi(rho_l, rho_g, x)
            else:
                raise ValueError(
                    f"Unsupported void_fraction_model for acceleration pressure drop: {void_fraction_model!r} "
                    "(only 'Homogeneous'/None and 'Zivi' are supported here)"
                )
            v_l = 1.0 / rho_l
            v_g = 1.0 / rho_g

            term1 = x ** 2 * v_g / max(EPS, alpha)
            term2 = (1.0 - x) ** 2 * v_l / max(EPS, 1.0 - alpha)

            return term1 + term2

    A_cross = PI * D ** 2 / 4.0
    G = m_dot / A_cross  # kg/(m²·s)

    f_inlet = f_acceleration(x_inlet, rho_l, rho_g, void_fraction_model)
    f_outlet = f_acceleration(x_outlet, rho_l, rho_g, void_fraction_model)

    dP_acceleration = G ** 2 * (f_inlet - f_outlet)

    return dP_acceleration


def pressure_drop_pipe_gravity_two_phase(L, rho_l, rho_g, x_inlet, 
                                     x_outlet, theta):
    """
    Compute gravitational pressure drop in two-phase flow.

    Uses mean quality to estimate mean density.

    Formula
    -------
    ΔP_gravity = ρ_mean · g · L · sin(θ)

    Parameters
    ----------
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


def pressure_drop_pipe_two_phase(AS, params, m_dot, P_su=None, T_su=None, correlation='Friedel',
                            void_fraction_model=None, include_acceleration=True,
                            include_gravity=True):
    """
    Compute total two-phase pressure drop in a straight pipe.

    Combines frictional, acceleration, and gravitational components.

    Formula
    -------
    ΔP_total = ΔP_friction + ΔP_acceleration + ΔP_gravity

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Two-phase fluid state at the pipe inlet (quality and pressure already
        set). Liquid/vapor properties and inlet quality are derived from it
        via `get_saturated_phase_properties`. NOTE: this function updates
        `AS` in place to the outlet state (used to obtain x_outlet) — it is
        left at the outlet state when this function returns.
    params : dict
        Pipe geometry: 'D' (or 'd_hyd') [m], 'L' [m], and optionally
        'K' (absolute roughness [m], default 0.0) and 'theta' (inclination
        from horizontal [degrees], default 0.0). Also passed through to
        `pressure_drop_pipe_frictional_two_phase`.
    m_dot : float
        Mass flow rate [kg/s]
    P_su : float, optional
        Inlet pressure [Pa]. If given, AS is updated to this pressure at its
        current quality (via PQ_INPUTS) before computing properties.
    T_su : float, optional
        Unused for two-phase flow (state is already fixed by pressure and
        quality at saturation); accepted only so the single-phase and
        two-phase functions share the same call signature.
    correlation : str, optional
        Two-phase frictional pressure drop correlation, passed to
        `pressure_drop_pipe_frictional_two_phase`. Only 'Friedel' (default)
        is currently implemented.
    void_fraction_model : str, optional
        Slip ratio model for void fraction: None/'Homogeneous' (default), 'Zivi'
    include_acceleration : bool, optional
        Include acceleration pressure drop (default: True)
    include_gravity : bool, optional
        Include gravitational pressure drop (default: True)

    Returns
    -------
    dP_total : float
        Total pressure drop [Pa]

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
    if P_su is not None:
        AS.update(CP.PQ_INPUTS, P_su, AS.Q())

    D = params.get('D', params.get('d_hyd'))
    L = params['L']
    theta = params.get('theta', 0.0)

    # ========== STEP 1: Compute friction ΔP ==========

    props = get_saturated_phase_properties(AS)
    x_inlet = props["x"]
    rho_l = props["rho_l"]
    rho_g = props["rho_v"]

    dP_friction = pressure_drop_pipe_frictional_two_phase(
        AS, params, m_dot, correlation=correlation
    )

    # ========== STEP 2: Estimate outlet pressure & compute x_outlet ==========
    p_inlet = AS.p()  # Inlet pressure
    h = AS.hmass()   # Specific enthalpy (assumed constant)
    p_outlet = p_inlet - dP_friction

    # Get outlet quality at (h, p_outlet) — HmassP_INPUTS takes (Hmass, P)
    AS.update(CP.HmassP_INPUTS, h, p_outlet)
    x_outlet = AS.Q()  # Quality at outlet

    # Acceleration pressure drop
    if include_acceleration:
        dP_acceleration = pressure_drop_pipe_acceleration_two_phase(
            m_dot, D, rho_l, rho_g, x_inlet, x_outlet, void_fraction_model=void_fraction_model
        )
    else:
        dP_acceleration = 0.0

    # Gravitational pressure drop
    if include_gravity:
        dP_gravity = pressure_drop_pipe_gravity_two_phase(
            L, rho_l, rho_g, x_inlet, x_outlet, theta
        )
    else:
        dP_gravity = 0.0
    # Total pressure drop
    dP_total = dP_friction + dP_acceleration + dP_gravity

    return dP_total