"""
Pressure drop correlations for single-phase and two-phase flow in straight pipes.

author: Elise Neven (elise.neven@uliege.be)
"""

import math
import warnings
import numpy as np
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

from labothappy.correlations.properties.supercritical import (
    pseudo_boiling_bounds,
    pseudo_critical_temperature_CO2,
)

EPS = 1e-12
PI = math.pi
G_GRAVITY = 9.81  # m/s²

# ============================================================================
# SINGLE PHASE PRESSURE DROP
# ============================================================================

# 1) Friction factor correlations for single-phase flow in pipes (Darcy-Weisbach)

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


def friction_factor_konakov(Re):
    """
    Darcy friction factor using Konakov formula.

    Formula
    ------------------------
    f= (1.8 log10​(Re)−1.5)^−2

    Parameters
    ----------
    Re : float
        Reynolds number [-]

    Returns
    -------
    f : float
        Darcy friction factor [-]

    Notes
    -----
    - Valid for 4×10³ ≤ Re ≤ 3×10⁶
    - Valid for smooth pipes
    
    References
    ----------
    Khlapuk, M., Bezusyak, O., Volk, L., Zhang, Z. (2021). "Theoretical 
    research of friction factor in hydraulically smooth pipes." E3S Web 
    of Conferences, 280, 10009. 
    https://doi.org/10.1051/e3sconf/202128010009
    """

    f = (1.8*math.log10(Re) - 1.5)**(-2)

    return f


def friction_factor_petukhov(Re):
    """
    Darcy friction factor using Petukhov formula.

    Parameters
    ----------
    Re : float
        Reynolds number [-]

    Returns
    -------
    f : float
        Darcy friction factor [-]

    Notes
    -----
    - Valid for 3×10³–10⁴ ≤ Re ≤ 5×10⁶ (turbulent)
    - Valid for 0.5 ≤ Pr ≤ 200 
    - Valid for smooth pipes

    References
    ----------
    Petukhov, B.S. (1970). "Heat transfer and friction in turbulent pipe 
    flow with variable physical properties." Advances in Heat Transfer, 
    Vol. 6, pp. 503-564. Academic Press.
    (as reproduced in, e.g., Incropera et al., Fundamentals of Heat and 
    Mass Transfer, for smooth tubes: 3000 < Re < 5e6)
    """

    f = (0.79*np.log(Re) - 1.64)**(-2)

    return f

def friction_factor_haaland(K, d_hyd, Re):
    """
    Darcy friction factor using Haaland formula.

    Explicit approximation of Colebrook-White correlation.

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
    - Valid for 4×10³ ≤ Re ≤ 1×10^8 (turbulent)
    - Valid for 1×10−6 ≤ K/d_hyd ≤ 0.05

    References
    ----------
    Haaland, S.E. (1983). "Simple and Explicit Formulas for the Friction
    Factor in Turbulent Pipe Flow." Journal of Fluids Engineering, 105(1),
    89-90.
    """

    f = (1/(-1.8*np.log10((K/(3.75*d_hyd))**1.11 + 6.9/Re)))**2

    return f


def friction_factor_cheng_CO2(G, d_hyd, P, h, mu):
    """
    Friction factor for supercritical CO2 in a horizontal tube, based on
    a Froude number analysis distinguishing liquid-like, vapor-like, and
    two-phase-like (pseudo-boiling) flow regimes.

    Self-contained: creates its own CoolProp AbstractState for CO2
    internally to compute the properties it needs (density, pseudo-boiling
    bounds, pseudo-vapor quality), so no AbstractState is needed from the
    caller.

    Assumptions
    -----------
    d_hyd = 10 mm
    G = 496.7-1346.2 kg/(m^2*s)
    heat flux (qw): 97.4 ~ 400.3 kW/m^2
    Pressure (P): 7.53-23.51 MPa
    fluid : CO2

    Froude number (Fr): 7.58e-5 ~ 1834
    Reynolds number (Re): 6.18e4 ~ 5.35e5

    Inputs
    ------
    G     : Mass flux, flow rate per cross-section area [kg/(m^2*s)]
    d_hyd : Hydraulic diameter [m]
    P     : Pressure [Pa]
    h     : Specific enthalpy [J/kg]
    mu    : Dynamic viscosity [Pa*s]

    Outputs
    -------
    f : Darcy friction factor [-]

    Reference
    ---------
    Supercritical carbon dioxide heat transfer in horizontal tube based
    on the Froude number analysis (2024)

    Liangyuan Cheng, Jinliang Xu, Wenxuan Cao, Kaiping Zhou, Guanglin Liu
    """
    fluid = 'CO2'
    AS = CP.AbstractState('HEOS', fluid)

    AS.update(CP.HmassP_INPUTS, h, P)
    rho = AS.rhomass()

    # Pseudo-boiling transition bounds (T_minus = liquid-like,
    # T_plus = vapor-like) and the resulting pseudo-vapor quality at inlet
    T_pc, T_minus, T_plus = pseudo_boiling_bounds(
        P, fluid, T_pc_func=pseudo_critical_temperature_CO2
    )
    AS.update(CP.PT_INPUTS, P, T_minus)
    h_LL = AS.hmass() # Liquid-like enthalpy at T_minus
    AS.update(CP.PT_INPUTS, P, T_plus)
    h_LV = AS.hmass() # Vapor-like enthalpy at T_plus
    x_pb = (h - h_LL) / (h_LV - h_LL) # Pseudo-vapor quality at inlet (0 = liquid-like, 1 = vapor-like)

    Re = G * d_hyd / mu

    # 1) VL (Vapor-Like) and LL (Liquid-Like) regimes
    if x_pb < 0 or x_pb > 1:
        C = 1.5
        f = C * (1.82 * np.log10(Re) - 1.64)**(-2)

    # 2) TPL (Two-Phase-Like) regime
    else:
        Fr_VL = G**2 * x_pb / (rho**2 * G_GRAVITY * d_hyd)
        f = 4.379 * Re**(-0.388) * Fr_VL**(-0.0167)

    return f

# 2) Compute the pressure drop in a straight pipe (single-phase flow)

def pressure_drop_pipe_single_phase(AS, pipe_geom, m_dot, correlation='Churchill'):
    """
    Compute total pressure drop in a straight pipe (single-phase flow).

    Combines frictional (Darcy-Weisbach) and gravitational pressure drops.
    The friction factor is computed with the correlation selected via
    `correlation`.

    Formula
    -------
    ΔP_friction = f · (L/d_hyd) · (ρ·v²/2)
    ΔP_gravity = ρ·g·L·sin(θ)
    ΔP_total = ΔP_friction + ΔP_gravity

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Single-phase fluid state. AS must already have its state set (e.g.
        via HmassP_INPUTS) by the caller.
    pipe_geom : dict
        Pipe geometry: 'D' [m] (pipe inner diameter), 'L' [m], and
        optionally 'K' (absolute roughness [m], default 0.0) and 'theta'
        (inclination from horizontal [degrees], default 0.0). For a round
        pipe the hydraulic diameter equals D; this is where that identity
        is applied (`d_hyd = D`).
    m_dot : float
        Mass flow rate [kg/s]
    correlation : str, optional
        Friction factor correlation:
        'Churchill' (default, valid across all Reynolds numbers, accounts
        for roughness), 'Swamee-Jain' (explicit, turbulent flow only,
        accounts for roughness), 'Haaland' (explicit, turbulent flow only,
        accounts for roughness), 'Konakov' or 'Petukhov' (turbulent,
        smooth pipes only — roughness `K` is ignored), 'Cheng-CO2'
        (supercritical CO2 only — `AS` must be a CO2 state; roughness is
        not used).

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

    D = pipe_geom['D']
    d_hyd = D
    A_cross = PI * D ** 2 / 4.0
    L = pipe_geom['L']
    theta = pipe_geom.get('theta', 0.0)
    K = pipe_geom.get('K', 0.0)

    mu = AS.viscosity()
    rho = AS.rhomass()

    v = m_dot / (rho * A_cross)  # Mean velocity [m/s]

    Re = compute_reynolds(d_hyd, mu, rho, v)

    if correlation == 'Churchill':
        f = friction_factor_churchill(K, d_hyd, Re)
    elif correlation == 'Swamee-Jain':
        f = friction_factor_swamee_jain(K, d_hyd, Re)
    elif correlation == 'Haaland':
        f = friction_factor_haaland(K, d_hyd, Re)
    elif correlation == 'Konakov':
        if K > 0:
            warnings.warn(
                f"K={K} was given but the 'Konakov' correlation is "
                "smooth-pipe-only and ignores roughness.",
                stacklevel=2,
            )
        f = friction_factor_konakov(Re)
    elif correlation == 'Petukhov':
        if K > 0:
            warnings.warn(
                f"K={K} was given but the 'Petukhov' correlation is "
                "smooth-pipe-only and ignores roughness.",
                stacklevel=2,
            )
        f = friction_factor_petukhov(Re)
    elif correlation == 'Cheng-CO2':
        if AS.name() != 'CarbonDioxide':
            raise ValueError(
                "The 'Cheng-CO2' correlation is specific to supercritical "
                f"CO2; got a fluid state for {AS.name()!r}."
            )
        if K > 0:
            warnings.warn(
                f"K={K} was given but the 'Cheng-CO2' correlation does not "
                "use roughness.",
                stacklevel=2,
            )
        G = m_dot / A_cross  # Mass flux [kg/(m²·s)]
        P = AS.p()
        h = AS.hmass()

        f = friction_factor_cheng_CO2(G, d_hyd, P, h, mu)
    else:
        raise ValueError(
            f"Unknown single-phase friction correlation: {correlation!r}. "
            "Available: 'Churchill', 'Swamee-Jain', 'Haaland', 'Konakov', "
            "'Petukhov', 'Cheng-CO2'"
        )

    # Frictional pressure drop (Darcy-Weisbach)
    dP_friction = f * (L / d_hyd) * (rho * v ** 2 / 2.0)

    # Gravitational pressure drop (hydrostatic)
    dP_gravity = G_GRAVITY * rho * L * math.sin(math.radians(theta))

    # Total pressure drop
    dP_total = dP_friction + dP_gravity

    return dP_total


# ============================================================================
# TWO-PHASE PRESSURE DROP
# ============================================================================

# Two phase pressure drops: in ACHP, Ian bell used lockhart correlation
# In fluid library: used Muller Steinhaguen Heck

def pressure_drop_muller_steinhagen_heck(m_dot, x, rho_l, rho_v, mu_l, mu_v, d_hyd, L, K=45e-6):
    """
    Two-phase frictional pressure drop, Muller-Steinhagen and Heck (1986) 
    correlation.

    Unlike two-phase multiplier models, it does not treat the two phases as 
    physically interacting streams; it is a direct empirical curve fit.

    Parameters
    ----------
    m_dot : float
        Total (liquid + gas) mass flow rate [kg/s]
    x : float
        Vapor quality [-], 0 < x < 1
    rho_l : float
        Liquid density [kg/m^3]
    rho_v : float
        Vaport density [kg/m^3]
    mu_l : float
        Liquid dynamic viscosity [Pa*s]
    mu_v : float
        Vapor dynamic viscosity [Pa*s]
    d_hyd : float
        Hydraulic (or pipe) diameter [m]
    L : float
        Tube length [m]
    K : float
        Absolute roughness [m]

    Output
    ------
    dP_tp = Two-phase frictional pressure drop over length L [Pa]

    Notes
    -----
    - Applicable for 0 < x < 1.
    - Developed from a data bank of 9300 measurements: air-water,
    water-hydrocarbon, and refrigerant two-phase flows in channels with
    diameters from 4 to 392 mm.
    - Reported as one of the most accurate simple correlations, especially
    for annular flow and mini/micro-channel turbulent flow.
    - Used in the Fluid library

    Reference
    ---------
    Muller-Steinhagen, H., Heck, K. (1986). "A Simple Friction Pressure 
    Drop Correlation for Two-Phase Flow in Pipes." Chemical Engineering 
    and Processing: Process Intensification, 20(6), 297-308.
    """
    if not (0 < x < 1):
        raise ValueError(f"Muller-Steinhagen-Heck correlation requires 0 < x < 1, got x = {x}")

    A_cross = np.pi / 4 * d_hyd**2  # cross-sectional flow area, circular channel assumed

    # Mass flux for each single-phase limit: total mass flow rate m, as if
    # it were entirely liquid or entirely gas, through the same flow area.
    G_lo = m_dot / A_cross
    G_go = m_dot / A_cross

    v_l = m_dot / (rho_l * A_cross)  # Mean velocity [m/s]
    Re_l = compute_reynolds(d_hyd, mu_l, rho_l, v_l)
    f_l = friction_factor_swamee_jain(K, d_hyd, Re_l)
    dP_lo = f_l * (L / d_hyd) * (rho_l * v_l ** 2 / 2.0)

    v_v = m_dot / (rho_v * A_cross)  # Mean velocity [m/s]
    Re_v = compute_reynolds(d_hyd, mu_v, rho_v, v_v)
    f_v = friction_factor_swamee_jain(K, d_hyd, Re_v)
    dP_go = f_v * (L / d_hyd) * (rho_v * v_v ** 2 / 2.0)

    # Muller-Steinhagen and Heck interpolation
    G_MSH = dP_lo + 2 * (dP_go - dP_lo) * x
    dP_tp = G_MSH * (1 - x)**(1 / 3) + dP_go * x**3

    return dP_tp

def pressure_drop_choi(AS, G, rho_su, rho_ex, P_sat, x_su, x_ex, L, d_hyd):
    """
    Two-phase pressure drop for evaporation and condensation in smooth
    and micro-fin tubes.

    A modification of the Bo Pierre (1964) homogeneous-flow correlation:
    the smooth-tube diameter is replaced by the hydraulic diameter,
    the two-phase specific volume average now includes the liquid phase
    (Pierre neglected it), and the friction factor was re-regressed
    against NIST micro-fin tube pressure drop data.

    Inputs
    ------
    AS      : CoolProp AbstractState object for the working fluid
    G       : Mass flux, flow rate per cross-section area [kg/(m^2*s)]
    rho_ex : Two-phase (quality-weighted) density at tube outlet [kg/m^3]
    rho_su  : Two-phase (quality-weighted) density at tube inlet [kg/m^3]
    P_sat   : Saturation pressure, evaluated at the linearly-averaged
            refrigerant temperature between inlet and outlet [Pa]
    x_ex     : Vapor quality at tube outlet [-]
    x_su     : Vapor quality at tube inlet [-]
    L       : Tube length [m]
    d_hyd      : Hydraulic diameter [m] (use the actual tube ID for a
                smooth tube, or the micro-fin hydraulic diameter otherwise)

    Outputs
    -------
    dP_tp : Two-phase pressure drop over length L [Pa]

    Notes
    -----
    - Applicable for both evaporation and condensation.
    - Developed and validated for refrigerants R125, R134a, R32, R410A,
        R22, R407C, and R32/R134a (25/75 % mass) in smooth and micro-fin
        tubes.
    - Predicted the NIST micro-fin database with an average absolute
        residual of 10.8 %, and smooth-tube data with a mean deviation of
        15.0 %.
    - rho_ex and rho_su must be the quality-weighted two-phase
        densities (1/rho = x*v_v + (1-x)*v_l) at the outlet/inlet, not
        the saturated-liquid or -vapor density alone.

    Reference
    ---------
    Choi, J.Y., Kedzierski, M.A., Domanski, P.A. (2001). "Generalized
    Pressure Drop Correlation for Evaporation and Condensation in
    Smooth and Micro-Fin Tubes." IIR Commission B1, Paderborn, Germany.
    (Full derivation and validation in NISTIR 6333, Choi, Kedzierski,
    Domanski, 1999.)
    """

    # Saturated liquid and vapor enthalpies at P_sat, for latent heat
    AS.update(CP.PQ_INPUTS, P_sat, 0)
    mu_l = AS.viscosity()   # Liquid dynamic viscosity [Pa*s]
    h_l = AS.hmass()        # Saturated liquid enthalpy [J/kg]

    AS.update(CP.PQ_INPUTS, P_sat, 1)
    h_v = AS.hmass()        # Saturated vapor enthalpy [J/kg]

    h_lv = h_v - h_l        # Latent heat of vaporization [J/kg]

    # All-liquid Reynolds number: entire mass flux G as if it were
    # 100% liquid, using the hydraulic diameter and liquid viscosity
    Re_fo = G * d_hyd / mu_l

    # Pierre's two-phase (boiling) number
    K_f = abs(x_ex - x_su) * h_lv / (G_GRAVITY * L)

    # New two-phase friction factor (Choi, Kedzierski, Domanski, 2001)
    f_N = 0.00506 * Re_fo**(-0.0951) * K_f**(0.1554)

    v_ex = 1 / rho_ex
    v_su = 1 / rho_su

    dP_friction = f_N * L * (v_ex + v_su) / d_hyd * G**2
    dP_acceleration = (v_ex - v_su) * G**2

    dP_tp = dP_friction + dP_acceleration

    return dP_tp



# ============================================================================
# TWO-PHASE FRICTION MULTIPLIER
# ============================================================================

def friedel_multiplier(AS, m_dot, d_hyd, K=0.0):
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
        Fr = G² / (g·d_hyd·ρ_h²)  [Froude number]
        We = G²·d_hyd / (σ·ρ_h)   [Weber number]
        G = m_dot / A_cross   [mass velocity, kg/(m²·s)]

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Two-phase fluid state (quality and pressure already set). Liquid
        and vapor properties (ρ_l, ρ_v, μ_l, μ_v, σ) are derived from it via
        `get_saturated_phase_properties`.
    m_dot : float
        Mass flow rate [kg/s]
    d_hyd : float
        Hydraulic diameter [m]
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
    A_cross = PI * d_hyd ** 2 / 4.0
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
    Re_l = compute_reynolds(d_hyd, mu_l, rho_l, v_l)
    f_l = friction_factor_churchill(K, d_hyd, Re_l)

    # Vapor (at total mass flow)
    v_v = m_dot / (rho_v * A_cross)
    Re_v = compute_reynolds(d_hyd, mu_v, rho_v, v_v)
    f_v = friction_factor_churchill(K, d_hyd, Re_v)

    # ====================================================================
    # Dimensionless numbers (using homogeneous density)
    # ====================================================================
    Fr = compute_froude(G, d_hyd, rho_h)
    We = compute_weber(G, d_hyd, sigma, rho_h)

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
# TWO-PHASE PRESSURE DROP COMPONENTS
# ============================================================================

def pressure_drop_pipe_frictional_two_phase(AS, pipe_geom, m_dot, correlation='Friedel'):
    """
    Compute frictional pressure drop in two-phase flow.

    Dispatches on `correlation` between two families of two-phase
    frictional correlations, both of which only need the pipe's inlet
    state (no outlet/exit information required):
    - 'Friedel': a multiplier model. ΔP_f = ΔP_lo · Φ_l², where ΔP_lo is
      the single-phase pressure drop (liquid only, at total mass flow)
      and Φ_l² is the Friedel friction multiplier (`friedel_multiplier`).
    - 'MSH' (Muller-Steinhagen & Heck): a direct empirical curve fit that
      returns the two-phase frictional pressure drop directly, with no
      multiplier step (`pressure_drop_muller_steinhagen_heck`).

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Two-phase fluid state (quality and pressure already set). Liquid
        and vapor properties are derived from it via
        `get_saturated_phase_properties`.
    pipe_geom : dict
        Pipe geometry: 'D' [m] (pipe inner diameter), 'L' [m], and
        optionally 'K' (absolute roughness [m], default 0.0). For a round
        pipe the hydraulic diameter equals D (`d_hyd = D`).
    m_dot : float
        Mass flow rate [kg/s]
    correlation : str, optional
        'Friedel' (default) or 'MSH'.

    Returns
    -------
    dP_friction : float
        Frictional pressure drop [Pa]

    References
    ----------
    Friedel, L. (1979). Improved friction pressure drop correlations for
        horizontal and vertical two-phase pipe flow. European Two-Phase Flow
        Group Meeting, Ispra, Italy.

    Muller-Steinhagen, H., Heck, K. (1986). "A Simple Friction Pressure
        Drop Correlation for Two-Phase Flow in Pipes." Chemical Engineering
        and Processing: Process Intensification, 20(6), 297-308.
    """
    props = get_saturated_phase_properties(AS)
    x = props["x"]
    rho_l = props["rho_l"]
    rho_v = props["rho_v"]
    mu_l = props["mu_l"]
    mu_v = props["mu_v"]

    D = pipe_geom['D']
    d_hyd = D
    L = pipe_geom['L']
    K = pipe_geom.get('K', 0.0)

    if correlation == 'Friedel':
        # Single-phase pressure drop (liquid only, at total mass flow)
        A_cross = PI * d_hyd ** 2 / 4.0
        v_l = m_dot / (rho_l * A_cross)
        Re_l = compute_reynolds(d_hyd, mu_l, rho_l, v_l)
        f_l = friction_factor_churchill(K, d_hyd, Re_l)
        dP_l = f_l * (L / d_hyd) * (rho_l * v_l ** 2 / 2.0)

        Phi_l2 = friedel_multiplier(AS, m_dot, d_hyd, K)
        dP_friction = dP_l * Phi_l2

    elif correlation == 'MSH':
        dP_friction = pressure_drop_muller_steinhagen_heck(
            m_dot, x, rho_l, rho_v, mu_l, mu_v, d_hyd, L, K=K
        )

    else:
        raise ValueError(
            f"Unknown two-phase friction correlation: {correlation!r}. "
            "Available: 'Friedel', 'MSH'"
        )

    return dP_friction


def pressure_drop_pipe_acceleration_two_phase(m_dot, d_hyd, rho_l, rho_g, x_inlet,
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
    d_hyd : float
        Hydraulic diameter [m]
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

    A_cross = PI * d_hyd ** 2 / 4.0
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


def pressure_drop_pipe_two_phase(AS, pipe_geom, m_dot, correlation='Friedel', void_fraction_model=None):
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
    pipe_geom : dict
        Pipe geometry: 'D' [m] (pipe inner diameter), 'L' [m], and
        optionally 'K' (absolute roughness [m], default 0.0) and 'theta'
        (inclination from horizontal [degrees], default 0.0). For a round
        pipe the hydraulic diameter equals D (`d_hyd = D`). Also passed
        through to `pressure_drop_pipe_frictional_two_phase`.
    m_dot : float
        Mass flow rate [kg/s]
    correlation : str, optional
        Two-phase frictional pressure drop correlation, passed to
        `pressure_drop_pipe_frictional_two_phase`: 'Friedel' (default) or
        'MSH' (Muller-Steinhagen & Heck).
    void_fraction_model : str, optional
        Slip ratio model used in the acceleration term:
        None/'Homogeneous' (default), 'Zivi'

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

    Muller-Steinhagen, H., Heck, K. (1986). "A Simple Friction Pressure
        Drop Correlation for Two-Phase Flow in Pipes." Chemical Engineering
        and Processing: Process Intensification, 20(6), 297-308.

    ACHP Documentation: https://achp.readthedocs.io/en/latest/ACHPComponents/FluidCorrelations.html
    """
    D = pipe_geom['D']
    d_hyd = D
    L = pipe_geom['L']
    theta = pipe_geom.get('theta', 0.0)

    # ========== STEP 1: Compute friction ΔP ==========

    props = get_saturated_phase_properties(AS)
    x_inlet = props["x"]
    rho_l = props["rho_l"]
    rho_g = props["rho_v"]

    dP_friction = pressure_drop_pipe_frictional_two_phase(
        AS, pipe_geom, m_dot, correlation=correlation
    )

    # ========== STEP 2: Estimate outlet pressure & compute x_outlet ==========
    p_inlet = AS.p()  # Inlet pressure
    h = AS.hmass()   # Specific enthalpy (assumed constant)
    p_outlet = p_inlet - dP_friction

    # Get outlet quality at (h, p_outlet) — HmassP_INPUTS takes (Hmass, P)
    AS.update(CP.HmassP_INPUTS, h, p_outlet)
    x_outlet = AS.Q()  # Quality at outlet

    # Acceleration pressure drop
    dP_acceleration = pressure_drop_pipe_acceleration_two_phase(
        m_dot, d_hyd, rho_l, rho_g, x_inlet, x_outlet, void_fraction_model=void_fraction_model
    )

    # Gravitational pressure drop
    dP_gravity = pressure_drop_pipe_gravity_two_phase(
        L, rho_l, rho_g, x_inlet, x_outlet, theta
    )

    # Total pressure drop
    dP_total = dP_friction + dP_acceleration + dP_gravity

    return dP_total