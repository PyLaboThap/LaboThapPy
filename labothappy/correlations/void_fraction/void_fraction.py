"""
Void fraction correlations for two-phase flow.

This module implements a set of void fraction (alpha) correlations commonly
used for two-phase refrigerant/ORC working fluid flow, ranging from simple
slip-ratio models to iterative and drift-flux formulations. A single
dispatcher, `compute_void_fraction`, selects and evaluates the requested
model from a CoolProp AbstractState and a params dict of model-specific
inputs (mass flux, surface tension, hydraulic diameter, inclination, etc.).

Correlations and sources
-------------------------
Zivi, Premoli, Lockhart-Martinelli, Hughmark and Graham were implemented
following their presentation in:
    Dickes, R. (n.d.). Charge-sensitive methods for the off-design
    performance characterization of organic Rankine cycle (ORC) power
    systems. PhD thesis.

Fauske, Armand-Treschev, Bankoff, Rouhani-Axelsson, Dix,
Woldesemayat-Ghajar and Cioncolini-Thome were implemented following:
    Thome, J. R., Cioncolini, A. (2015). "Void Fraction." In: Woodhead
    Publishing Series in Energy, pp. 85-112.
    https://doi.org/10.1142/9789814623216_0021

Quick reference
----------------
    Model                   Type            Best suited for
    ----------------------- --------------- --------------------------------
    Homogeneous             no-slip         theoretical baseline, x > 0.7
    Zivi                    slip-ratio      annular, high void fraction
    Fauske                  slip-ratio      steam-water, high mass flux
    Premoli                 slip-ratio      vertical adiabatic, upward flow
    Lockhart-Martinelli     empirical       turbulent-turbulent, horizontal
    Hughmark                iterative       vertical & horizontal, low-mid x
    Graham                  empirical       condensation only
    Armand-Treschev         K*alpha_h       vertical, heated, high pressure
    Bankoff                 K*alpha_h       steam-water, broad but coarse
    Rouhani-Axelsson        drift-flux      vertical, quality boiling, x>0.1
    Dix                     drift-flux      vertical, wide range, robust
    Woldesemayat-Ghajar     drift-flux      horizontal & upward inclined
    Cioncolini-Thome        empirical       annular, macro & microscale

@author: elise.neven@uliege.be
"""

import numpy as np
import CoolProp.CoolProp as CP

from labothappy.correlations.properties.dimensionless import compute_reynolds, compute_weber
from labothappy.toolbox.solvers.zero_brent import zero_brent

EPS = 1e-12
G_GRAVITY = 9.81  # Gravitational acceleration [m/s^2]
PI = np.pi

#======================================================================
# Homogeneous void fraction model
#======================================================================

def void_fraction_homogeneous(rho_l, rho_v, x):
    """
    Homogeneous void fraction model (no-slip, S = 1).

    Formula
    -------
    alpha = 1 / (1 + S * (rho_v / rho_l) * (1-x) / x), with S = 1

    Parameters
    ----------
    rho_l : float
        Liquid density [kg/m3]
    rho_v : float
        Saturated vapor density [kg/m3]
    x : float
        Saturated vapor quality [-]

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity/Assumptions
    --------------------
    Theoretically applicable to all two-phase conditions (no-slip assumption).
    Tends to overpredict void fraction at low quality (x < 0.3) where 
    significant slip occurs. Best accuracy at high quality (x > 0.7).

    References
    ----------
    The homogeneous model is foundational theory, not a specific paper. 
    For more general information:
    Thome, John R., et Andrea Cioncolini. Void Fraction. 2015. https://doi.org/10.1142/9789814623216_0021.
    """
    # Compute slip ratio
    s = 1

    # Compute void fraction
    # alpha = 1 / (1 + S * (rho_v / rho_l) * (1-x) / x)
    denominator = 1.0 + s * (rho_v / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha

#======================================================================
# Zivi void fraction model
#======================================================================

def void_fraction_zivi(rho_l, rho_v, x):
    """
    Zivi void fraction correlation (slip-ratio method).
    
    Formula
    -------
    alpha = 1 / (1 + S * (rho_v / rho_l) * (1-x) / x), with S = (rho_l / rho_v)^(1/3)

    Parameters
    ----------
    rho_l : float
        Liquid density [kg/m3]
    rho_v : float
        Vapor density [kg/m3]
    x : float
        Vapor quality [-]

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity/Assumptions
    --------------------
    Best-suited to annular, high-quality/high-void-fraction, low-friction conditions. 
    Not validated for bubbly/slug flow.

    References
    ----------
    S. M. Zivi. "Estimation of steady-state steam void-fraction by means of
    the principle of minimum entropy production".
    In: Journal of Heat Transfer 86.2 (1964), pp. 247-251.
    """
    # Compute slip ratio
    s_zivi = (rho_l / max(EPS, rho_v)) ** (1.0 / 3.0)

    # Compute void fraction
    denominator = 1.0 + s_zivi * (rho_v / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha


#======================================================================
# Faukse void fraction model
#======================================================================

def void_fraction_fauske(rho_l, rho_v, x):
    """
    Fauske void fraction correlation (slip-ratio method).

    Formula
    -------
    alpha = 1 / (1 + S * (rho_v / rho_l) * (1-x) / x), with S = (rho_l / rho_v)^(1/2)
    
    Parameters
    ----------
    rho_l : float
        Saturated iquid density [kg/m3]
    rho_v : float
        Saturated vapor density [kg/m3]
    x : float
        Vapor quality [-]

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity/Assumptions
    --------------------
    Best suited for watersteam two-phase flows at high mass flux conditions.

    References
    ----------
    Fauske, H.K., (1962). Contribution to the theory of two-phase, one component critical flow.  
    Argonne National Laboratory Report ANL‒6633.
    """
    # Compute slip ratio
    s = (rho_l / max(EPS, rho_v)) ** (1.0 / 2.0)

    # Compute void fraction
    denominator = 1.0 + s * (rho_v / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha

#======================================================================
# Premoli void fraction model
#======================================================================

def void_fraction_premoli(rho_l, rho_v, x, mu_l, sigma, d_hyd, G):
    """
    Premoli (also called CISE) void fraction correlation (slip-ratio method).
    
    Formula
    -------
    alpha = 1 / (1 + S * (rho_v / rho_l) * (1-x) / x), with 

    S = 1 + F1 * (y / (1 + y*F2) - y*F2)^0.5
    y = alpha_h / (1 - alpha_h)
    F1 = 1.578 * Re_l^-0.19 * (rho_l/rho_v)^0.22
    F2 = 0.0273 * We_l * Re_l^-0.51 * (rho_l/rho_v)^-0.08
        
    Parameters
    ----------
    rho_l : float
        Saturated liquid density [kg/m3]
    rho_v : float
        Saturated vapor density [kg/m3]
    x : float
        Vapor quality [-]
    mu_l : float
        Saturated liquid dynamic viscosity [Pa.s]
    sigma : float
            Surface tension [N/m]
    d_hyd : float
            Hydraulic diameter [m]
    G : float
        Mass flux [kg/(m2.s)]

    Returns
    -------
    alpha : float
        Void fraction [-]
    
    Validity/Assumptions
    --------------------
    Empirical model. These correlations were developed with a large experimental 
    dataset of two-phase mixtures flowing upwardly in vertical adiabatic channels. 
    Account for the influence of the mass flux.

    References
    ----------
    A. Premoli, D.D. Francesco, and A. Prina. “A dimensional correlation 
    for evaluating two-phase mixture density”. In: La Termotecnica 25 (1971), pp. 17–26.
    """
    Re_l = compute_reynolds(d_hyd=d_hyd, mu=mu_l, G=G)
    We_l = compute_weber(G=G, sigma=sigma, L=d_hyd, rho=rho_l)

    alpha_h = void_fraction_homogeneous(rho_l, rho_v, x)
    y = alpha_h / (1 - alpha_h)
    F1 = 1.578 * (Re_l ** (-0.19)) * (rho_l / rho_v) ** 0.22
    F2 = 0.0273 * We_l * (Re_l ** (-0.51)) * (rho_l / rho_v) ** (-0.08)
    s = 1 + F1 * np.sqrt(max(0, (y / (1 + y * F2)) - y * F2))

    denominator = 1.0 + s * (rho_v / max(EPS, rho_l)) * (1.0 - x) / max(EPS, x)
    alpha = 1.0 / max(EPS, denominator)

    return alpha

#======================================================================
# Lockhart Martinelli void fraction model
#======================================================================

def void_fraction_lockhart_martinelli(rho_l, rho_v, x, mu_l, mu_v):
    """
    Lockhart-Martinelli void fraction correlation (turbulent-turbulent form).

    Computes the local void fraction from the Martinelli parameter X_tt,
    assuming both liquid and vapor phases are turbulent (the most common
    case for refrigerant flows in evaporators/condensers). The correlation
    was originally presented as an empirical curve (no closed-form equation)
    fitting isothermal, adiabatic two-component (air-liquid) flow data in
    pipes; the piecewise closed-form expression used here is a later
    curve-fit approximation to that original data (commonly attributed to
    Chisholm, 1967), not an equation from the original 1949 paper.

    Formula:
    --------
    X_tt = ((1-x)/x)^0.9 * (mu_l/mu_v)^0.1 * (rho_v/rho_l)^0.5
    alpha = (1 + X_tt^0.8)^-0.378              for X_tt <= 10
    alpha = 0.823 - 0.157*ln(X_tt)             for X_tt > 10

    Parameters
        ----------
    rho_l, rho_v : float
        Saturated liquid and vapor densities [kg/m^3].
    x : float
        Vapor quality [-].
    mu_l, mu_v : float
        Saturated liquid and vapor dynamic viscosities [Pa.s].

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity / Applicability
    -------------------------
    - Assumes turbulent-turbulent flow regime (both phases turbulent);
        not valid as-is for laminar-liquid or laminar-vapor regimes, which
        the original LM framework treats with separate parameter forms.
    - Original 1949 data: isothermal, adiabatic, two-component (mostly
        air-liquid) flow in pipes at near-atmospheric pressure, horizontal
        orientation, small pipe diameters.
    - Not derived from a single-component boiling/condensing two-phase
        system (no heat transfer, no phase change) - applying it to
        refrigerant evaporation/condensation is an extrapolation common
        in HVAC literature but not directly validated by the original data.
    - The X_tt <= 10 / > 10 split point (and the fitted constants) come
        from the later closed-form approximation, not the 1949 paper.

    References
    ----------
    Lockhart, R. W., Martinelli, R. C. (1949). "Proposed correlation of
    data for isothermal two-phase two-component flow in pipes."
    Chemical Engineering Progress, 45(1), pp. 39-48.
    Chisholm, D. (1967). Closed-form approximation to the Lockhart-Martinelli
    correlation curves (turbulent-turbulent case).
    """
    X_tt = (((1 - x) / x) ** 0.9) * ((mu_l / mu_v) ** 0.1) * ((rho_v / rho_l) ** 0.5)

    if X_tt <= 10:
        alpha = min(1, max(0, (1 + X_tt ** 0.8) ** (-0.378)))
    else:
        alpha = min(1, max(0, 0.823 - 0.157 * np.log(X_tt)))

    return alpha

#======================================================================
# Hughmark void fraction model
#======================================================================

def void_fraction_hughmark(rho_l, rho_v, x, mu_l, mu_v, d_hyd, G):
    """
    Hughmark void fraction correlation (iterative, slip-ratio method).
    
    Corrects the homogeneous void fraction alpha_h with an empirical factor
    Kh. The relationship used here is a polynomial fit to Hughmark's original
    correlation chart/tabulated data, not a closed-form expression from the
    1962 paper itself (which presented Kh as a graphical correlation).

    Formula:
    --------
    Z = (d_hyd*G / mu_mix)^(1/6) * ((G*x/rho_v)^2 / (g*D_h*alpha_h*(1-alpha_h)))^(1/8)
    ln(Kh) = polynomial fit in ln(Z)  [degree-4, fitted to Hughmark's chart]
    alpha = Kh * alpha_h              (solved iteratively/implicitly)

    Parameters
    ----------
    rho_l, rho_v : float
        Saturated liquid and vapor densities [kg/m^3].
    x : float
        Vapor quality [-], clipped to [0.001, 0.99] to avoid singularities.
    mu_l, mu_v : float
        Saturated liquid and vapor dynamic viscosities [Pa.s].
    d_hyd : float
        Hydraulic diameter [m].
    G : float
        Total mass flux [kg/(m^2 s)].

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity / Applicability
    -------------------------
    - Developed from air-water and other two-component gas-liquid holdup
      data covering both vertical and horizontal pipe orientations.
    - No official numeric quality/pressure/diameter bounds are reproduced
      here since the original chart-based correlation doesn't state hard
      cutoffs; the quality clipping to [0.001, 0.99] here is a numerical
      safeguard against the correlation's own asymptotes, not a stated
      physical limit from the paper.

    References
    ----------
    Hughmark, G. A. (1962). "Holdup in gas-liquid flow." Chemical
    Engineering Progress, 58(4), pp. 62-65.
    """

    def residual_void_fraction_hughmark(alpha_guess, x, alpha_h, rho_v, mu_v, mu_l, d_hyd, G):
        Z = (((d_hyd * G) / (mu_l + alpha_guess * (mu_v - mu_l))) ** (1/6)) * \
            (((1 / (G_GRAVITY * d_hyd)) * (G * x / (rho_v * alpha_h * (1 - alpha_h))) ** 2) ** (1/8))
        ln_Z = np.log(Z)
        p1 = -0.010060658854755
        p2 = 0.155594796014726
        p3 = -0.870912508715887
        p4 = 2.167004115373165
        p5 = -2.224608445535130
        ln_Kh = p1 * ln_Z ** 4 + p2 * ln_Z ** 3 + p3 * ln_Z ** 2 + p4 * ln_Z + p5
        Kh = np.exp(ln_Kh)
        alpha_new = Kh * alpha_h
        res = (alpha_guess - alpha_new)
        return alpha_new, res

    x_min = 0.001
    x_max = 0.99
    if x < x_min:
        x = x_min
    elif x > x_max:
        x = x_max

    alpha_h = void_fraction_homogeneous(rho_l=rho_l, rho_v=rho_v, x=x)

    def f(alpha_guess):
        _, res = residual_void_fraction_hughmark(
            alpha_guess, x, alpha_h, rho_v, mu_v, mu_l, d_hyd, G
        )
        return res

    machep = np.finfo(float).eps
    t = 1e-10
    tol_f = 1e-8

    alpha, res_alpha = zero_brent(1e-6, 1 - 1e-6, machep, t, f, tol_f)

    return alpha

#======================================================================
# Graham void fraction model (condensation only)
#======================================================================

def void_fraction_graham(rho_v, x, d_hyd, G):
    """
    Graham void fraction correlation (condensation only).

    Formula:
    --------
    Ft = sqrt( x^3 * G^2 / (g * rho_v^2 * D * (1-x)) )
    alpha = 1 - exp(-1 - 0.3*ln(Ft) - 0.0328*ln(Ft)^2)   for Ft > 0.01032
    alpha = 0                                             otherwise

    Parameters
    ----------
    rho_v : float
        Vapor density [kg/m^3].
    x : float
        Vapor quality [-].
    d_hyd : float
        Hydraulic diameter [m].
    G : float
        Total mass flux [kg/(m^2 s)].

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity / Applicability
    -------------------------
    - Explicitly a condensation-only correlation.
    - The cutoff Ft = 0.01032 acts as a regime switch rather than a
    smooth physical limit - worth flagging in your flow-regime map
    as a hard discontinuity in the correlation's output.

    References
    ----------
    Graham B Wallis. One-dimensional two-phase flow. 1969.

    """
    Ft = (((x ** 3) * (G ** 2)) / (G_GRAVITY * (rho_v ** 2) * d_hyd * (1 - x))) ** 0.5

    if Ft > 0.01032:
        alpha = min(1, max(0, 1 - np.exp(-1 - 0.3 * np.log(Ft) - 0.0328 * (np.log(Ft)) ** 2)))
    else:
        alpha = 0

    return alpha

#======================================================================
# Armand Treschev void fraction model (K-epsilon method)
#======================================================================

def void_fraction_armand_treschev(rho_l, rho_v, x):
    """
    Armand-Treschev void fraction correlation (K-epsilon method).
    
    Corrects the homogeneous void fraction alpha_h with a simple linear
    factor K that depends only on vapor quality x.

    Formula:
    --------
    K = 0.833 + 0.167*x
    alpha = K * alpha_h

    Parameters
    ----------
    rho_v, rho_l : float
        Saturated vapor and liquid densities [kg/m^3].
    x : float
        Vapor quality [-].
    
    Returns
    -------
    alpha : float
        Void fraction [-]
    
    Validity / Applicability
    -------------------------
    - Originates from steam-water mixture flow in heated boiler pipes at
        high pressure - a specifically vertical, high-pressure, boiling
        (heated) flow context, quite different from adiabatic refrigerant
        flow in evaporator/condenser piping.
    - Less used and accurate correlation
    """
    alpha_h = void_fraction_homogeneous(rho_l=rho_l, rho_v=rho_v, x=x)
    K = 0.833 + 0.167 * x

    alpha = K * alpha_h
    return alpha

#======================================================================
# Bankoff void fraction model (K-epsilon method)
#======================================================================

def void_fraction_bankoff(rho_l, rho_v, x, P):
    """
    Bankoff void fraction correlation (K-epsilon method).

    Corrects the homogeneous void fraction alpha_h with a pressure-dependent
    factor K, using Bankoff's "variable density single-fluid" model. 

    Formula:
    --------
        K = 0.71 + 0.0145 * (P / 1e6)     [P in Pa, converted to MPa]
        alpha = K * alpha_h

    Parameters
    ----------
    rho_l, rho_v : float
        Saturated liquid and vapor densities [kg/m^3].
    x : float
        Vapor quality [-].
    P : float
        System pressure [Pa].
    
    Returns
    -------
    alpha : float
        Void fraction [-]
    
    Validity / Applicability
    -------------------------
    - Developed for steam-water flow; the model's entrainment/distribution
    parameter was tuned empirically against three independent sets of
    void fraction data.
    - Simple and not the most accurate correlation.

    References
    ----------
    Bankoff, S. G. (1960). "A variable density single fluid model for
    two-phase flow with particular reference to steam-water flow." ASME
    Journal of Heat Transfer, 82, pp. 265-272.
    """
    alpha_h = void_fraction_homogeneous(rho_l=rho_l, rho_v=rho_v, x=x)
    K = 0.71 + 0.0145 * (P / 1e6)  # pressure term in MPa

    alpha = K * alpha_h

    return alpha

#======================================================================
# Rouhani-Axelsson void fraction model (Drift-flux correlations)
#======================================================================

def void_fraction_rouhani_axelsson(rho_l, rho_v, x, G, sigma):
    """
    Rouhani-Axelsson void fraction correlation (drift-flux method, vertical flow).

    Formula:
    --------
    J_l = (1-x)*G/rho_l          (liquid superficial velocity)
    J_g = x*G/rho_g              (vapor superficial velocity)
    C0 = 1 + 0.2*(1-x)           (distribution parameter)
    V_drift = 1.18 * (g*sigma*(rho_l-rho_g)/rho_l^2)^0.25
    alpha = J_g / (C0*(J_l+J_g) + V_drift)

    Parameters
    ----------
    rho_l, rho_v : float
        Saturated liquid and vapor densities [kg/m^3].
    x : float
        Vapor quality [-].
    G : float
        Total mass flux [kg/(m^2 s)].
    sigma : float
        Surface tension [N/m].

    Returns
    -------
    alpha : float
        Void fraction [-]
    
    Validity / Applicability
    -------------------------
    - Originally developed for subcooled and quality boiling of steam-water
    in vertical channels.
    - Reported valid for void fractions larger than 0.1; below that the
    correlation is not recommended.
    - Independent comparison studies found it has a tendency to slightly
    under-predict void fractions higher than 0.6.

    References
    ----------
    Rouhani, S. Z., Axelsson, E. (1970). "Calculation of void volume
    fraction in the subcooled and quality boiling regions." International
    Journal of Heat and Mass Transfer, 13, pp. 383-393.
    """
    J_l = (1 - x) * G / rho_l
    J_g = x * G / rho_v

    C0 = 1 + 0.2 * (1 - x)
    V_drift = 1.18 * ((G_GRAVITY * sigma * (rho_l - rho_v)) / (rho_l ** 2)) ** 0.25

    alpha = J_g / (C0 * (J_l + J_g) + V_drift)

    return alpha

#======================================================================
# Chexal/Dix void fraction model (Drift-flux correlations)
#======================================================================

def void_fraction_dix(x, rho_l, rho_g, G, sigma):
    """


    References
    ----------
    Dix, G. E. (1971). "Vapor void fractions for forced convection with
    subcooled boiling at low flow rates." PhD thesis (or ANS/Berkeley
    report), University of California, Berkeley.
    Chexal, B., Horowitz, J., Lellouche, G. (1986). "An assessment of
    eight void fraction models for vertical flows." EPRI Report
    NSAC-107, Palo Alto, USA.
    """


def void_fraction_dix(rho_l, rho_v, x, G, sigma):
    """
    Dix void fraction correlation (drift-flux method, vertical flow).

    Formula:
    --------
    J_l = (1-x)*G/rho_l          (liquid superficial velocity)
    J_g = x*G/rho_g              (vapor superficial velocity)
    n = (rho_g/rho_l)^0.1
    C0 = (J_g/(J_l+J_g)) * (1 + (J_l/J_g)^n)
    V_drift = 2.9 * (g*sigma*(rho_l-rho_g)/rho_l^2)^0.25
    alpha = J_g / (C0*(J_l+J_g) + V_drift)

    Parameters
    ----------
    rho_l, rho_v : float
        Liquid and vapor densities [kg/m^3].
    x : float
        Vapor quality [-].
    G : float
        Total mass flux [kg/(m^2 s)].
    sigma : float
        Surface tension [N/m].

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity / Applicability
    -------------------------
    - In broad multi-source comparisons, the Dix drift-flux model was
    one of the best-performing correlations for void fraction
    prediction across a large consolidated database.
    - Reported to have outstanding consistency across the entire void
    fraction range.
    - Vertical-flow-oriented; applying it to horizontal
    pipe segments in your library is an extrapolation.

    References
    ----------
    Dix, G. E. (1971). "Vapor void fractions for forced convection with
    subcooled boiling at low flow rates." PhD thesis (or ANS/Berkeley
    report), University of California, Berkeley.
    Chexal, B., Horowitz, J., Lellouche, G. (1986). "An assessment of
    eight void fraction models for vertical flows." EPRI Report
    NSAC-107, Palo Alto, USA.
    """
    J_l = (1 - x) * G / rho_l
    J_g = x * G / rho_v

    n = (rho_v / rho_l) ** 0.1
    C0 = J_g / (J_l + J_g) * (1 + (J_l / J_g) ** n)
    V_drift = 2.9 * ((G_GRAVITY * sigma * (rho_l - rho_v)) / (rho_l ** 2)) ** 0.25

    alpha = J_g / (C0 * (J_l + J_g) + V_drift)

    return alpha


#======================================================================
# Woldesemayat Ghajar void fraction model (Drift-flux correlations)
#======================================================================

def void_fraction_woldesemayat_ghajar(rho_l, rho_v, x, G, sigma, D, P, theta=np.pi/2):
    """
    Woldesemayat and Ghajar (2007) void fraction correlation (drift-flux method).
    
    Extends the Dix (1971) distribution parameter C0 with a drift velocity
    that additionally accounts for pipe diameter, inclination angle, and
    system pressure.

    Formula:
    -------
    J_l = (1-x)*G/rho_l          (liquid superficial velocity)
    J_g = x*G/rho_g              (vapor superficial velocity)
    n = (rho_g/rho_l)^0.1
    C0 = (J_g/(J_l+J_g)) * (1 + (J_l/J_g)^n)
    V_drift = 2.9 * (g*sigma*D*(1+cos(theta))*(rho_l-rho_g)/rho_l^2)^0.25
                * (1.22 + 1.22*sin(theta))^(101325/P)
    alpha = J_g / (C0*(J_l+J_g) + V_drift)

    Parameters
    ----------
    rho_l, rho_v : float
        Saturated liquid and vapor densities [kg/m^3].
    x : float
        Vapor quality [-].
    G : float
        Total mass flux [kg/(m^2 s)].
    sigma : float
        Surface tension [N/m].
    D : float
        Pipe internal diameter [m].
    P : float
        System pressure [Pa].
    theta : float, optional
        Pipe inclination angle from horizontal [rad], default pi/2
        (vertical). theta=0 is horizontal.

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity / Applicability
    -------------------------
    - Built as an improvement of the Dix correlation, explicitly adding
    diameter, inclination and pressure dependence.
    - Calibrated against a large unbiased database.
    - Tube diameters in the underlying database ranged from 12.7 to
    102.26 mm.
    - NOT validated for downward-inclined or vertical-downward flow -
    restricted to horizontal and upward orientations only.

    References
    ----------
    Woldesemayat, M. A., Ghajar, A. J. (2007). "Comparison of void
    fraction correlations for different flow patterns in horizontal and
    upward inclined pipes." International Journal of Multiphase Flow,
    33, pp. 347-370.
    """
    J_l = (1 - x) * G / rho_l
    J_g = x * G / rho_v

    n = (rho_v / rho_l) ** 0.1
    C0 = J_g / (J_l + J_g) * (1 + (J_l / J_g) ** n)

    V_drift = (2.9 * ((G_GRAVITY * sigma * D * (1 + np.cos(theta)) * (rho_l - rho_v)) / (rho_l ** 2)) ** 0.25
                   * (1.22 + 1.22 * np.sin(theta)) ** (101325 / P))

    alpha = J_g / (C0 * (J_l + J_g) + V_drift)

    return alpha

#======================================================================
# Cioncolini-Thome void fraction model (Drift-flux correlations)
#======================================================================

def void_fraction_cioncolini_thome(rho_l, rho_v, x):
    """
    Cioncolini and Thome (2012). Annular flow correlation.
    Best suited for: 0 < x < 1, 1e-3 < rho_g/rho_l < 1, 0.7 < eps < 1.

    Formula:
    -------
    rho_ratio = rho_v / rho_l
    h = -2.129 + 3.129 * rho_ratio^-0.2186
    n = 0.3487 + 0.6513 * rho_ratio^0.515
    alpha = (h * x^n) / (1 + (h - 1) * x^n)

    Parameters
    ----------
    rho_l, rho_v : float
        Liquid and vapor (gas) densities [kg/m^3].
    x : float
        Vapor quality [-].

    Returns
    -------
    alpha : float
        Void fraction [-]

    Validity / Applicability
    -------------------------
    - Developed and validated specifically for annular flow, in both
    macro- and microscale channels.
    - Studies for 8 different gas-liquid and vapor-liquid combinations 
    (water-steam, R410a, water-air, water-argon, water-nitrogen, water plus
    alcohol-air, alcohol-air and kerosene-air)
    - Tube diameters 1.05 mm to 45.5 mm and for both circular 
    and non-circular channels.

    References
    ----------
    Cioncolini, A., Thome, J. R. (2012a). "Void fraction prediction in
    annular two-phase flow." International Journal of Multiphase Flow,
    43, pp. 72-84.
    """
    rho_ratio = rho_v / rho_l

    h = -2.129 + 3.129 * rho_ratio ** (-0.2186)
    n = 0.3487 + 0.6513 * rho_ratio ** 0.515

    return (h * x ** n) / (1 + (h - 1) * x ** n)


#=======================================================================
# VOID FRACTION COMPUTATION
#=======================================================================

def compute_void_fraction(AS, params, m_dot = None, void_fraction_model='Homogeneous'):
    """
    Compute void fraction for two-phase flow using the chosen correlation.

    Parameters
    ----------
    AS : object
        Fluid state with attributes Q (quality), rhol, rhog, mu_l, mu_g, P
        (only the ones needed by the chosen model are accessed).
    params : dict
        May contain 'sigma', 'd_hyd', 'D', 'G', 'theta' depending on model.
    void_fraction_model : str
        One of: 'Homogeneous', 'Zivi', 'Fauske', 'Premoli',
        'Lockhart-Martinelli', 'Hughmark',
        'Graham', 'Armand-Treschev', 'Bankoff', 'Rouhani-Axelsson',
        'DIX', 'Woldesemayat-Ghajar', 'Cioncolini-Thome'.

    Returns
    -------
    float
        Void fraction [-]
    """
    x = AS.Q()
    x = max(EPS, min(1.0 - EPS, x))

    P = AS.p()

    AS_l = CP.AbstractState(AS.backend_name(), AS.fluid_names()[0])
    AS_l.update(CP.PQ_INPUTS, P, 0)
    rho_l = AS_l.rhomass()

    AS_v = CP.AbstractState(AS.backend_name(), AS.fluid_names()[0])
    AS_v.update(CP.PQ_INPUTS, P, 1)
    rho_v = AS_v.rhomass()

    if void_fraction_model == 'Homogeneous':
        alpha = void_fraction_homogeneous(rho_l, rho_v, x)

    elif void_fraction_model == 'Zivi':
        alpha = void_fraction_zivi(rho_l, rho_v, x)

    elif void_fraction_model == 'Fauske':
        alpha = void_fraction_fauske(rho_l, rho_v, x)

    elif void_fraction_model == 'Premoli':
        mu_l = AS_l.viscosity()
        sigma = AS.surface_tension()
        d_hyd = params.get('d_hyd')
        A_cross = PI * d_hyd ** 2 / 4.0
        G = m_dot / A_cross  # kg/(m²·s)
        alpha = void_fraction_premoli(rho_l, rho_v, x, mu_l, sigma, d_hyd, G)

    elif void_fraction_model == 'Lockhart-Martinelli':
        mu_l = AS_l.viscosity()
        mu_v = AS_v.viscosity()
        alpha = void_fraction_lockhart_martinelli(rho_l, rho_v, x, mu_l, mu_v)

    elif void_fraction_model == 'Hughmark':
        mu_l = AS_l.viscosity()
        mu_v = AS_v.viscosity()
        d_hyd = params.get('d_hyd')
        A_cross = PI * d_hyd ** 2 / 4.0
        G = m_dot / A_cross  # kg/(m²·s)
        alpha = void_fraction_hughmark(rho_l, rho_v, x, mu_l, mu_v, d_hyd, G)

    elif void_fraction_model == 'Graham':
        d_hyd = params.get('d_hyd')
        A_cross = PI * d_hyd ** 2 / 4.0
        G = m_dot / A_cross  # kg/(m²·s)
        alpha = void_fraction_graham(rho_v, x, d_hyd, G)

    elif void_fraction_model == 'Armand-Treschev':
        alpha = void_fraction_armand_treschev(rho_l, rho_v, x)

    elif void_fraction_model == 'Bankoff':
        alpha = void_fraction_bankoff(rho_l, rho_v, x, P)

    elif void_fraction_model == 'Rouhani-Axelsson':
        d_hyd = params.get('d_hyd')
        A_cross = PI * d_hyd ** 2 / 4.0
        G = m_dot / A_cross  # kg/(m²·s)
        sigma = AS.surface_tension()
        alpha = void_fraction_rouhani_axelsson(rho_l, rho_v, x, G, sigma)

    elif void_fraction_model == 'DIX':
        d_hyd = params.get('d_hyd')
        A_cross = PI * d_hyd ** 2 / 4.0
        G = m_dot / A_cross  # kg/(m²·s)
        sigma = AS.surface_tension()
        alpha = void_fraction_dix(rho_l, rho_v, x, G, sigma)

    elif void_fraction_model == 'Woldesemayat-Ghajar':
        d_hyd = params.get('d_hyd')
        A_cross = PI * d_hyd ** 2 / 4.0
        G = m_dot / A_cross  # kg/(m²·s)
        sigma = AS.surface_tension()
        theta = params.get('theta', np.pi / 2)
        alpha = void_fraction_woldesemayat_ghajar(rho_l, rho_v, x, G, sigma, d_hyd, P, theta)

    elif void_fraction_model == 'Cioncolini-Thome':
        alpha = void_fraction_cioncolini_thome(rho_l, rho_v, x)

    else:
        raise ValueError(f"Unknown void fraction model: {void_fraction_model}")

    # Safeguard output
    alpha = max(EPS, min(1.0 - EPS, alpha))

    return alpha
