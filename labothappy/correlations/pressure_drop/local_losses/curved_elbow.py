"""
Local pressure loss correlations for curved pipe bends.

Implements Idelchik curved bend correlations across all flow regimes
(laminar, transition, turbulent) with smooth regime blending.

The Idelchik method uses a local resistance coefficient ζ that depends on:
- Bend geometry (radius ratio R0/D, deflection angle)
- Flow regime (laminar, transition, turbulent)
- Surface roughness (via friction factor correction)

Functions
---------
- compute_zeta_LOC: Local resistance coefficient (geometry + flow regime)
- compute_cf_fri: Roughness correction factor
- compute_zeta_tot: Total loss coefficient (with roughness)
- pressure_drop_curved_elbow: Pressure drop from bend

References
----------
Idelchik, I. E. (1986). Handbook of Hydraulic Resistance: Coefficients of
    Local Resistance and of Friction. Hemisphere Publishing Corporation.

Modelica Standard Library: Fluid.Dissipative.PressureLoss.Bend.CurvedBend_LDP
"""

import math
from labothappy.correlations.properties.dimensionless import compute_reynolds
from labothappy.correlations.pressure_drop.straight_pipe_DP import (
    friction_factor_swamee_jain,
)
from labothappy.toolbox.helper_function.step_smoother import stepsmoother

PI = math.pi
EPS = 1e-12

# Reynolds number boundaries for flow regime transitions
RE_LAM_MAX = 6.5e3  # Upper limit of laminar regime
RE_TURB_MIN = 4e4  # Lower limit of turbulent regime
RE_TURB_MAX = 3e5  # Upper limit of intermediate turbulent regime


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def compute_re_lam_transition(K, d_hyd):
    """
    Compute Reynolds number at laminar–transition boundary.

    For curved bends, the transition depends on relative roughness.
    Rougher pipes transition earlier to turbulent flow.

    Parameters
    ----------
    K : float
        Absolute roughness [m]
    d_hyd : float
        Hydraulic diameter [m]

    Returns
    -------
    Re_transition : float
        Reynolds number at transition [-]

    Notes
    -----
    Based on Idelchik correlations for curved bends.
    """
    k_rel = K / max(EPS, d_hyd) # relative roughness
    exponent = 0.0065 / max(k_rel, 0.007)
    return min(RE_LAM_MAX, max(1e2, 754.0 * math.exp(exponent)))


def compute_idelchik_coefficients(d_hyd, R0, delta_rad):
    """
    Compute Idelchik coefficients for curved bends.

    These coefficients (A1, A2, B1, C1) characterize the bend geometry
    and define the pressure loss across different flow regimes.

    Parameters
    ----------
    d_hyd : float
        Hydraulic diameter [m]
    R0 : float
        Bend radius (centerline) [m]
    delta_rad : float
        Deflection angle [radians]

    Returns
    -------
    dict with keys:
        'A1' : float
            Laminar/turbulent base coefficient (geometry-dependent)
        'A2' : float
            Laminar regime coefficient (Re-dependent term)
        'B1' : float
            Reynolds correction factor (curvature-dependent)
        'C1' : float
            Constant factor (typically 1.0)
        'R0_d_hyd' : float
            Radius ratio R0/d_hyd

    Notes
    -----
    A1 depends on bend angle (δ).
    A2 depends on radius ratio (R0/D) – tighter bends have larger A2.
    B1 depends on radius ratio – affects turbulent regime.

    References
    ----------
    Idelchik, I. E. (1986). Handbook of Hydraulic Resistance.
    """
    delta_deg = abs(delta_rad) * 180 / PI
    R0_d_hyd = R0 / max(EPS, d_hyd)

    # ====================================================================
    # A1: Laminar/turbulent base coefficient (depends on angle)
    # ====================================================================
    sin_half_delta = math.sin(PI / 180 * delta_deg / 2)
    A1 = max(EPS, 0.95 * sin_half_delta ** 2 + 2.05 * sin_half_delta ** 4)

    # ====================================================================
    # A2: Laminar regime coefficient (depends on radius ratio)
    # ====================================================================
    # Tighter bends (smaller R0/D) have larger laminar losses
    if R0_d_hyd > 2.0:
        A2 = 600.0
    elif R0_d_hyd > 1.0:
        A2 = 1000.0
    elif R0_d_hyd > 0.7:
        A2 = 3000.0
    elif R0_d_hyd > 0.55:
        A2 = 6000.0
    else:
        A2 = 4000.0

    # ====================================================================
    # B1: Reynolds correction factor (depends on radius ratio)
    # ====================================================================
    if R0_d_hyd >= 1.0:
        B1 = 0.21 * (R0_d_hyd) ** (-0.5)
    else:
        B1 = 0.21 * (R0_d_hyd) ** (-2.5)
    B1 = max(EPS, B1)

    # ====================================================================
    # C1: Constant factor
    # ====================================================================
    C1 = 1.0

    return {
        "A1": A1,
        "A2": A2,
        "B1": B1,
        "C1": C1,
        "R0_d_hyd": R0_d_hyd,
    }


def compute_reynolds_correction_turbulent(Re, R0_d_hyd):
    """
    Compute Reynolds-number correction factor for turbulent regime.

    In turbulent flow through curved bends, the loss coefficient decreases
    with increasing Reynolds number due to reduced viscous effects relative
    to inertial forces.

    Formula (SOURCE: Idelchik)
    -------------------------
    k_Re = 11.5 / Re^0.19  (for R0/D > 0.7)
    k_Re = 5.45 / Re^0.131 (for 0.55 ≤ R0/D ≤ 0.7)
    k_Re = 1 + 4400/Re     (for R0/D < 0.55)

    Parameters
    ----------
    Re : float
        Reynolds number [-]
    R0_d_hyd : float
        Radius ratio R0/d_hyd [-]

    Returns
    -------
    k_Re : float
        Reynolds correction factor [-]

    Notes
    -----
    Tighter bends (smaller R0/D) have stronger Reynolds dependence.
    """
    if R0_d_hyd > 0.7:
        return 11.5 / (max(EPS, Re) ** 0.19)
    elif R0_d_hyd >= 0.55:
        return 5.45 / (max(EPS, Re) ** 0.131)
    else:
        return 1.0 + 4400.0 / max(EPS, Re)


# ============================================================================
# LOCAL RESISTANCE COEFFICIENT
# ============================================================================


def compute_zeta_LOC(d_hyd, R0, delta_rad, K, rho, mu, m_dot):
    """
    Compute local resistance coefficient ζ_LOC for curved bends.

    Uses Idelchik correlations with smooth blending across flow regimes:
    - Laminar (Re < Re_transition): ζ = A2/Re + A1·B1·C1
    - Transition: Smooth blending between laminar and turbulent
    - Turbulent (Re > Re_turb_max): ζ = A1·B1·C1·k_Re

    The loss coefficient accounts for:
    - Bend geometry (radius ratio, deflection angle)
    - Flow regime (viscous vs. inertial effects)
    - Surface roughness (via friction factor)

    Parameters
    ----------
    d_hyd : float
        Hydraulic diameter [m]
    R0 : float
        Bend radius (centerline) [m]
    delta_rad : float
        Deflection angle [radians]
    K : float
        Absolute roughness [m]
    rho : float
        Fluid density [kg/m³]
    mu : float
        Dynamic viscosity [Pa·s]
    m_dot : float
        Mass flow rate [kg/s]

    Returns
    -------
    zeta_LOC : float
        Local resistance coefficient (geometry + flow regime) [-]

    Notes
    -----
    - ζ_LOC does NOT include roughness correction (use compute_zeta_tot)
    - Smooth transition between regimes via stepsmoother
    - Valid for all Reynolds numbers

    References
    ----------
    Idelchik, I. E. (1986). Handbook of Hydraulic Resistance.
    """
    A_cross = PI * d_hyd ** 2 / 4.0
    v = abs(m_dot) / max(EPS, rho * A_cross)
    Re = compute_reynolds(rho, v, d_hyd, mu)

    # Get Idelchik coefficients
    coeffs = compute_idelchik_coefficients(d_hyd, R0, delta_rad)
    A1 = coeffs["A1"]
    A2 = coeffs["A2"]
    B1 = coeffs["B1"]
    C1 = coeffs["C1"]
    R0_d_hyd = coeffs["R0_d_hyd"]

    # ====================================================================
    # Base loss coefficient (turbulent, smooth transition)
    # ====================================================================
    zeta_base = A1 * B1 * C1

    # ====================================================================
    # Laminar regime (Re < Re_transition)
    # ====================================================================
    Re_transition = compute_re_lam_transition(K, d_hyd)
    zeta_lam = A2 / max(EPS, Re) + zeta_base

    # ====================================================================
    # Intermediate turbulent regime (Re_turb_min < Re < Re_turb_max)
    # ====================================================================
    k_Re = compute_reynolds_correction_turbulent(Re, R0_d_hyd)
    zeta_turb_inter = max(EPS, k_Re * zeta_base)

    # ====================================================================
    # High turbulent regime (Re > Re_turb_max)
    # ====================================================================
    zeta_turb_high = zeta_base

    # ====================================================================
    # Smooth blending between regimes
    # ====================================================================
    s_lam = stepsmoother(Re, Re_transition, RE_TURB_MIN)
    zeta_mid = (1.0 - s_lam) * zeta_lam + s_lam * zeta_turb_inter

    s_turb = stepsmoother(Re, RE_TURB_MIN, RE_TURB_MAX)
    zeta_LOC = (1.0 - s_turb) * zeta_mid + s_turb * zeta_turb_high

    return zeta_LOC


# ============================================================================
# ROUGHNESS CORRECTION
# ============================================================================


def compute_cf_fri(zeta_LOC, d_hyd, R0, delta_rad, K, rho, mu, m_dot):
    """
    Compute surface roughness correction factor CF_fri.

    Accounts for the effect of pipe surface roughness on pressure loss
    in curved bends. Roughness becomes important in turbulent flow.

    Formula
    -------
    raw_effect = (f_fri · L_arc / d_hyd) / ζ_LOC
    CF_fri = 1 + S · raw_effect

    where:
        f_fri: friction factor for equivalent straight pipe
        L_arc: arc length of bend = R0 · δ
        S: smoothing function (0 in laminar, 1 in turbulent)

    Parameters
    ----------
    zeta_LOC : float
        Local resistance coefficient (geometry only) [-]
    d_hyd : float
        Hydraulic diameter [m]
    R0 : float
        Bend radius [m]
    delta_rad : float
        Deflection angle [radians]
    K : float
        Absolute roughness [m]
    rho : float
        Fluid density [kg/m³]
    mu : float
        Dynamic viscosity [Pa·s]
    m_dot : float
        Mass flow rate [kg/s]

    Returns
    -------
    CF_fri : float
        Roughness correction factor (typically 1.0–1.4) [-]

    Notes
    -----
    - CF_fri ≈ 1.0 in laminar flow (roughness has little effect)
    - CF_fri > 1.0 in turbulent flow (roughness increases loss)
    - Capped at 1.4 to avoid unrealistic corrections

    References
    ----------
    Idelchik, I. E. (1986). Handbook of Hydraulic Resistance.
    """
    A_cross = PI * d_hyd ** 2 / 4.0
    v = abs(m_dot) / max(EPS, rho * A_cross)
    Re = compute_reynolds(rho, v, d_hyd, mu)

    # Arc length of the bend
    L_arc = abs(delta_rad) * R0

    # Friction factor for straight pipe (used for roughness correction)
    f_fri = friction_factor_swamee_jain(K, d_hyd, Re)

    # Raw roughness effect
    raw_effect = min(1.4, (f_fri * L_arc / max(EPS, d_hyd)) / max(EPS, zeta_LOC))

    # Smoothing: roughness effect only in turbulent regime
    Re_transition = compute_re_lam_transition(K, d_hyd)
    s = stepsmoother(Re, RE_LAM_MAX, Re_transition)
    CF_fri = 1.0 + s * raw_effect

    return CF_fri


# ============================================================================
# TOTAL LOSS COEFFICIENT & PRESSURE DROP
# ============================================================================


def compute_zeta_tot(d_hyd, R0, delta_rad, K, rho, mu, m_dot):
    """
    Compute total pressure loss coefficient ζ_tot for curved bends.

    Combines geometry-dependent loss (ζ_LOC) with roughness correction (CF_fri).

    Formula
    -------
    ζ_tot = CF_fri · ζ_LOC

    Parameters
    ----------
    d_hyd : float
        Hydraulic diameter [m]
    R0 : float
        Bend radius [m]
    delta_rad : float
        Deflection angle [radians]
    K : float
        Absolute roughness [m]
    rho : float
        Fluid density [kg/m³]
    mu : float
        Dynamic viscosity [Pa·s]
    m_dot : float
        Mass flow rate [kg/s]

    Returns
    -------
    zeta_tot : float
        Total loss coefficient (geometry + roughness) [-]

    References
    ----------
    Idelchik, I. E. (1986). Handbook of Hydraulic Resistance.
    """
    zeta_LOC = compute_zeta_LOC(d_hyd, R0, delta_rad, K, rho, mu, m_dot)
    CF_fri = compute_cf_fri(zeta_LOC, d_hyd, R0, delta_rad, K, rho, mu, m_dot)
    zeta_tot = max(EPS, CF_fri) * zeta_LOC
    return zeta_tot


def pressure_drop_curved_elbow(d_hyd, R0, delta_rad, K, rho, mu, m_dot):
    """
    Compute pressure drop across a curved pipe bend.

    Uses Idelchik correlations with smooth blending across all flow regimes.

    Formula
    -------
    ΔP = ζ_tot · (ρ · v²) / 2

    where:
        ζ_tot: total loss coefficient
        v: mean flow velocity = m_dot / (ρ · A_cross)

    Parameters
    ----------
    d_hyd : float
        Hydraulic diameter [m]
    R0 : float
        Bend radius (centerline) [m]
    delta_rad : float
        Deflection angle [radians]
    K : float
        Absolute roughness [m]
    rho : float
        Fluid density [kg/m³]
    mu : float
        Dynamic viscosity [Pa·s]
    m_dot : float
        Mass flow rate [kg/s]

    Returns
    -------
    dP : float
        Pressure drop [Pa]

    Examples
    --------
    >>> # 90° bend, R0/D = 1.0, steel pipe
    >>> dP = pressure_drop_curved_elbow(
    ...     d_hyd=0.01,           # 10 mm diameter
    ...     R0=0.01,              # 1:1 radius ratio
    ...     delta_rad=math.pi/2,  # 90 degrees
    ...     K=45e-6,              # Steel roughness
    ...     rho=1000,             # Water
    ...     mu=1e-3,
    ...     m_dot=0.1
    ... )
    >>> print(f"ΔP = {dP:.2f} Pa")

    References
    ----------
    Idelchik, I. E. (1986). Handbook of Hydraulic Resistance: Coefficients of
        Local Resistance and of Friction. Hemisphere Publishing Corporation.

    Modelica Standard Library: Fluid.Dissipative.PressureLoss.Bend.CurvedBend_LDP
    """
    A_cross = PI * d_hyd ** 2 / 4.0
    v = abs(m_dot) / max(EPS, rho * A_cross)

    zeta_tot = compute_zeta_tot(d_hyd, R0, delta_rad, K, rho, mu, m_dot)
    dP = zeta_tot * 0.5 * rho * v ** 2

    return dP
