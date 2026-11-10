# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:51:41 2026

@author: basil
"""

#%% IMPORTS
import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import warnings
from labothappy.correlations.properties.thermal_conductivity import conducticity_R1233zd

#%% HELPER fun

def chevron_deg(angle):
    """Chevron angle in degrees, whether given in radians or degrees.
    Chevron angles lie in [0, 90] deg = [0, pi/2] rad, so a value above
    pi/2 can only be in degrees."""
    if angle > np.pi / 2:
        if angle > 90:
            raise ValueError(f"Chevron angle {angle} is above 90 deg.")
        return float(angle)          # already in degrees
    return float(np.degrees(angle))  # radians -> degrees

def _zeta_martin(Re, alpha):
    """Martin (1996) / VDI Heat Atlas N6 Darcy-type friction factor.
    alpha: corrugation angle from the main flow direction [rad]."""
    Re = max(Re, 1.0)
    # Straight-channel (phi = 0) and cross-corrugated (phi = 90 deg) factors
    if Re < 2000:
        zeta0, zeta1_0 = 64 / Re, 597 / Re + 3.85
    else:
        zeta0 = (1.8 * np.log10(Re) - 1.5) ** -2     # Konakov; check vs. the source
        zeta1_0 = 39 / Re ** 0.289

    a, b, c = 3.8, 0.18, 0.36
    zeta1 = a * zeta1_0
    rhs = (np.cos(alpha) / np.sqrt(b * np.tan(alpha) + c * np.sin(alpha) + zeta0 / np.cos(alpha))
           + (1 - np.cos(alpha)) / np.sqrt(zeta1))
    return 1.0 / rhs**2

_SAT_CACHE = {}

def _sat_props(AS, p_sat):
    """Saturated liquid / vapour properties at p_sat [Pa].

    They depend only on (fluid, backend, p): results are cached with p rounded
    to 1 Pa, which avoids recomputing them at every residual evaluation.
    Falls back on (p, T) evaluations next to saturation when the backend fails
    on transport properties. AS is left in a saturated state on return: the
    caller must re-set it before reusing it.
    Keys: rho_l, rho_v, mu_l, mu_v (nan if unavailable), k_l, cp_l, Pr_l, h_lv, sigma."""
    fluid = AS.fluid_names()[0]
    key = (fluid, AS.backend_name(), round(p_sat))
    props = _SAT_CACHE.get(key)
    if props is None:
        if len(_SAT_CACHE) > 20000:
            _SAT_CACHE.clear()
        props = _SAT_CACHE[key] = _compute_sat_props(AS, fluid, float(round(p_sat)))
    return props

def _compute_sat_props(AS, fluid, p):
    # Saturated liquid
    AS.update(CP.PQ_INPUTS, p, 0)
    T_sat = AS.T()
    rho_l, h_l, cp_l = AS.rhomass(), AS.hmass(), AS.cpmass()
    try:
        sigma = AS.surface_tension()
    except ValueError:
        sigma = PropsSI('surface_tension', 'P', p, 'Q', 0, fluid)

    if fluid in ('R1233zd(E)', 'R1233ZD(E)'):
        k_l = conducticity_R1233zd(T_sat, p)
        k_known = True
    else:
        k_known = False
    try:
        mu_l = AS.viscosity()
        if not k_known:
            k_l = AS.conductivity()
    except ValueError:
        AS.update(CP.PT_INPUTS, p, T_sat - 0.001)
        mu_l = AS.viscosity()
        if not k_known:
            k_l = AS.conductivity()
    Pr_l = mu_l * cp_l / k_l

    # Saturated vapour
    AS.update(CP.PQ_INPUTS, p, 1)
    rho_v, h_v = AS.rhomass(), AS.hmass()
    try:
        mu_v = AS.viscosity()
    except ValueError:
        try:
            AS.update(CP.PT_INPUTS, p, T_sat + 0.001)
            mu_v = AS.viscosity()
        except ValueError:
            mu_v = np.nan

    return {"rho_l": rho_l, "rho_v": rho_v, "mu_l": mu_l, "mu_v": mu_v, "k_l": k_l,
            "cp_l": cp_l, "Pr_l": Pr_l, "h_lv": h_v - h_l, "sigma": sigma}

#%% SINGLE PHASE

def BPHE_sinnott_SP_HTC(AS, geom, cond):
    """
    Coulson & Richardson's Chemical Engineering Design Vol. 6
    
    R. K. Sinnott (1983)
    
    Eq. 12.77
    """
    # 1) Condition Extraction
    G = cond['G']
    mu = cond['mu']
    mu_w = cond['mu_w']
    Pr = cond['Pr']
    k = cond['k']

    # 2) Geometry Parameters
    Dh = geom['Dh']
    
    # 3) Reynolds number
    Re = G*Dh/mu
    
    # 4) Nusselt number
    Nu = 0.26*Re**(0.65)*Pr**(0.4)*(mu/mu_w)**0.14
    
    # 5) HTC
    h_conv = (k/Dh)*Nu
    
    return h_conv

def BPHE_kumar_SP_HTC(AS, geom, cond):
    """Kumar, H., The Plate Heat Exchanger: Construction and Design,
    IChemE Symp. Ser. 86, 1275-1288, 1984.
    Angle measured from the direction perpendicular to the flow (geom['beta'] is
    in the library convention, from the flow direction, and is converted here);
    values between tabulated angles are interpolated linearly, outside are clamped."""
    
    def _kumar_row(rows, Re):
        for Re_max, C, m in rows:
            if Re <= Re_max:
                return C * Re**m
            
    # Kumar (1984): (angle [deg], [(Re_max, C, m), ...]); last row has Re_max = inf
    _KUMAR = [
        (30, [(10, 0.718, 0.349), (np.inf, 0.348, 0.663)]),
        (45, [(10, 0.718, 0.349), (100, 0.400, 0.598), (np.inf, 0.300, 0.663)]),
        (50, [(20, 0.630, 0.333), (300, 0.291, 0.591), (np.inf, 0.130, 0.732)]),
        (60, [(20, 0.562, 0.326), (400, 0.306, 0.529), (np.inf, 0.108, 0.703)]),
        (65, [(20, 0.562, 0.326), (500, 0.331, 0.503), (np.inf, 0.087, 0.718)]),
    ]
    
    G, mu, mu_w, Pr, k = cond['G'], cond['mu'], cond['mu_w'], cond['Pr'], cond['k']
    Dh = geom['Dh']
    beta = 90.0 - chevron_deg(geom['beta'])     # from the perpendicular to the flow

    Re = G * Dh / mu
    if not 0.1 <= Re <= 10000 or not 30 <= beta <= 65:
        warnings.warn("Kumar HTC: outside the range of the data (0.1 <= Re <= 10000, "
                      "30-65 deg): extrapolated or clamped.")
    angles = [a for a, _ in _KUMAR]
    Nu_tab = [_kumar_row(rows, Re) for _, rows in _KUMAR]
    Nu = np.interp(beta, angles, Nu_tab) * Pr**(1/3) * (mu / mu_w)**0.17
    return k / Dh * Nu

def BPHE_thonon_SP_HTC(AS, geom, cond):
    """Thonon, B., Design Method for Plate Evaporators and Condensers,
    BHR Group Conf. Ser. 18, 37-47, 1995.  Valid for 50 <= Re <= 15000.
    Angle measured from the direction perpendicular to the flow (geom['beta'] is
    in the library convention, from the flow direction, and is converted here)."""
    # Thonon (1995): angle [deg] -> (C, m)
    _THONON = [(30, 0.2946, 0.700), (45, 0.2998, 0.645), (60, 0.2267, 0.631), (75, 0.1, 0.687)]
    
    G, mu, Pr, k = cond['G'], cond['mu'], cond['Pr'], cond['k']
    Dh = geom['Dh']
    beta = 90.0 - chevron_deg(geom['beta'])     # from the perpendicular to the flow

    Re = G * Dh / mu
    if not 50 <= Re <= 15000:
        warnings.warn(f"Re = {Re:.1f} outside the Thonon range [50, 15000]: extrapolated.")
    angles = [a for a, _, _ in _THONON]
    Nu_tab = [C * Re**m for _, C, m in _THONON]
    Nu = np.interp(beta, angles, Nu_tab) * Pr**(1/3)
    return k / Dh * Nu

def BPHE_muley_manglik_SP_HTC(AS, geom, cond):
    """
    Muley, A. & Manglik, R.M. (1999). Experimental Study of Turbulent Flow Heat
    Transfer and Pressure Drop in a Plate Heat Exchanger With Chevron Plates.
    J. Heat Transfer 121(1):110-117. DOI: 10.1115/1.2825923

    Validity: Re >= 1000, 30 <= beta <= 60 deg, 1 <= phi <= 1.5, 2 < Pr < 6.
    beta measured from the main flow direction (higher beta = 'harder' plate).
    geom['phi'] (surface enlargement factor) is optional; default 1.17.
    """
    G, mu, mu_w, Pr, k = cond['G'], cond['mu'], cond['mu_w'], cond['Pr'], cond['k']
    Dh = geom['Dh']
    beta = chevron_deg(geom['beta'])
    if 'phi' not in geom:
        warnings.warn("Muley-Manglik HTC: 'phi' not given in geom, using 1.17.")
    phi = geom.get('phi', 1.17)

    Re = G * Dh / mu
    if Re < 1000 or not 30 <= beta <= 60 or not 1.0 <= phi <= 1.5:
        warnings.warn("Muley-Manglik HTC: outside the validated range (Re >= 1000, "
                      "30-60 deg, 1 <= phi <= 1.5): extrapolated.")

    C = 0.2668 - 0.006967*beta + 7.244e-5*beta**2
    D = 20.7803 - 50.9372*phi + 41.1585*phi**2 - 10.1507*phi**3
    p = 0.728 + 0.0543*np.sin(np.pi*beta/45 + 3.7)

    Nu = C * D * Re**p * Pr**(1/3) * (mu/mu_w)**0.14
    return k / Dh * Nu

def BPHE_martin_SP_HTC(AS, geom, cond):
    """Martin (1996), VDI Heat Atlas N6. Angle from the main flow direction.
    cond['mu_w'] is optional: without it the viscosity correction is 1 (gases)."""
    
    G, mu, Pr, k = cond['G'], cond['mu'], cond['Pr'], cond['k']
    mu_w = cond.get('mu_w')
    phi = geom.get('phi', 1.17)
    Dh = geom['Dh']/phi
    beta = chevron_deg(geom['beta'])
    alpha = np.radians(beta)

    Re = G * Dh / mu
    
    if not 200 <= Re <= 10000 or not 0 <= beta <= 80:
        warnings.warn("Martin HTC: outside the range of the data (200 <= Re <= 10000, "
                      "0-80 deg): extrapolated.")
    zeta = _zeta_martin(Re, alpha)                       # log10 version
    visc = (mu / mu_w) ** (1/6) if mu_w else 1.0
    Nu = 0.122 * Pr**(1/3) * visc * (zeta * Re**2 * np.sin(2*alpha))**0.374

    return k / Dh * Nu

#%% TWO-PHASE : BOILING 

def BPHE_cooper_boiling_HTC(AS, geom, cond):
    """
    Nucleate pool boiling heat transfer coefficient, after M.G. Cooper (1984),
    "Heat flow rates in saturated nucleate boiling - A wide-ranging examination
    using reduced properties", Advances in Heat Transfer, Vol. 16, pp. 157-239.
    """
    p_sat = cond['p']
    q = cond['q']
    
    roughness = geom.get("roughness", 0)
    
    p_crit = AS.p_critical()
    M_molar = AS.molar_mass()  # kg/mol
    
    # Unit conversions to what Cooper's correlation expects
    p_red = p_sat / p_crit                       # reduced pressure [-]
    M_gmol = M_molar * 1e3                   # [kg/mol] -> [g/mol]

    if not 0.0 < p_red < 1.0:
        raise ValueError(f"Cooper: reduced pressure {p_red:.3f} outside ]0, 1[.")
    if not (np.isfinite(q) and q > 0.0):
        raise ValueError(f"Cooper: heat flux q = {q} must be finite and > 0.")
    
    if roughness is None or roughness <= 0:
        Rp_um = 1.0                          # Cooper's default when unknown [um]
    else:
        Rp_um = roughness * 1e6              # [m] -> [um]

    htc = 55 * p_red ** (0.12 - 0.2 * np.log10(Rp_um)) * (-np.log10(p_red)) ** (-0.55) * q ** 0.67 * M_gmol ** (-0.5)
    
    return htc

def BPHE_amalfi_boiling_HTC(AS, geom, cond):
    """
    Amalfi, Vakili-Farahani & Thome (2016), Int. J. Refrig. 61, 185-203 (Part 2).
    Flow boiling in plate heat exchangers.

    cond : G (mass flux in ONE channel) [kg/m2/s], q (local heat flux) [W/m2],
           x (mean quality of the cell) [-], p (saturation pressure) [Pa]
    geom : Dh (hydraulic diameter, 2b/Phi) [m],
           beta (chevron angle from the main flow direction, rad or deg)
    Returns h [W/m2/K].

    TODO check: beta_max = 70 deg here; the ht library uses 45 deg.
    """
    beta_max_deg=70.0
    g = 9.81
    
    p_sat = cond['p']
    x = cond['x']
    G = cond['G']
    q = cond['q']

    if not np.isfinite(x):
        raise ValueError("Amalfi boiling: cond['x'] must be finite (two-phase cell expected).")
    x = min(max(x, 1e-4), 1.0 - 1e-4)
    
    beta = chevron_deg(geom['beta'])
    Dh = geom['Dh']
    
    # Saturated properties
    s = _sat_props(AS, p_sat)
    k_l = s["k_l"]
    rho_l = s["rho_l"]
    mu_l = s["mu_l"]
    rho_v = s["rho_v"]
    mu_v = s["mu_v"]
    sigma = s["sigma"]
    h_lv = s["h_lv"]

    beta_star = beta / beta_max_deg
    rho_star = rho_l / rho_v
    rho_m = 1.0 / (x / rho_v + (1.0 - x) / rho_l)       # homogeneous density
    Bo = max(q / (G * h_lv), 1e-6)                       # guard: q = 0 would give h = 0
    Bd = (rho_l - rho_v) * g * Dh**2 / sigma

    if Bd < 4.0:
        We = G**2 * Dh / (rho_m * sigma)
        Nu = 982.0 * beta_star**1.101 * We**0.315 * Bo**0.320 * rho_star**-0.224
    else:
        if not np.isfinite(mu_v):
            raise ValueError("Amalfi boiling (Bd >= 4) needs the vapour viscosity, "
                             "not available for this fluid.")
        Re_g = G * x * Dh / mu_v
        Re_lo = G * Dh / mu_l
        Nu = (18.495 * beta_star**0.248 * Re_g**0.135 * Re_lo**0.351
              * Bd**0.235 * Bo**0.198 * rho_star**-0.223)
    return Nu * k_l / Dh

#%% TWO-PHASE : CONDENSING 

def BPHE_shah_condensing_HTC(AS, geom, cond):
    """
    Shah (2021) condensation in plate heat exchangers (to be checked against the paper),
    forced-convection branch from Longo et al. (2015).

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality of the cell) [-],
           p (saturation pressure) [Pa]
    geom : Dh (hydraulic diameter) [m],
           phi (surface enlargement factor, NOT the chevron angle) [-], default 1.17
    Returns h [W/m2/K].
    """
    G = cond['G']    
    x = cond['x']
    p_sat = cond['p']
    
    phi = geom.get('phi', 1.17)

    Dh = geom['Dh']/phi
    
    if not 2.3 <= G <= 165:
        warnings.warn("Shah condensation: mass flux outside the validated range "
                      "(2.3-165 kg/m2/s): extrapolated.")

    # Saturated properties
    s = _sat_props(AS, p_sat)
    Pr_l = s["Pr_l"]
    k_l = s["k_l"]
    rho_l = s["rho_l"]
    mu_l = s["mu_l"]
    rho_v = s["rho_v"]
    
    g = 9.81
    x = min(max(x, 0.0), 0.99)                    # Re_LS -> 0 at x = 1
    Re_l = G * (1 - x) * Dh / mu_l
    G_eq = G * ((1 - x) + x * (rho_l / rho_v) ** 0.5)
    Re_eq = G_eq * Dh / mu_l

    h_fc = 1.875 * (k_l / Dh) * phi * Re_eq**0.445 * Pr_l**(1/3)
    if Re_eq >= 1600:
        return h_fc
    h_grav = (1.32 * phi * Re_l**(-1/3)
              * (rho_l * (rho_l - rho_v) * g * k_l**3 / mu_l**2) ** (1/3))
    return max(h_grav, h_fc)

def BPHE_han_condensing_HTC(AS, geom, cond):
    """
    Han, Lee & Kim (2003), J. Korean Phys. Soc. 43(1), 66-73 (to be checked).
    Condensation in brazed plate HX, Nu = Ge1 Re_eq^Ge2 Pr^(1/3).

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-], p [Pa]
    geom : Dh (hydraulic diameter) [m], p_co (corrugation pitch) [m],
           beta (chevron angle from the main flow direction, rad or deg;
           Han's own angle is measured from the horizontal axis: pi/2 - beta)
    Returns h [W/m2/K].
    """
    G = cond['G']    
    x = cond['x']
    p_sat = cond['p']
    
    Dh = geom['Dh']
    
    theta = np.radians(chevron_deg(geom["beta"]))   # from the flow direction = Han's pi/2 - beta
    p_co = geom["p_co"]

    if not 0.0 < theta <= np.pi / 2:
        raise ValueError(f"Han condensation: angle {theta} rad from the flow direction "
                         "outside ]0, pi/2].")
    if not 13 <= G <= 34 or not 45 <= np.degrees(theta) <= 70:
        warnings.warn("Han condensation: outside the range of the data (G 13-34 kg/m2/s, "
                      "plates of 45/55/70 deg from the flow direction): extrapolated.")
    
    # Saturated properties
    s = _sat_props(AS, p_sat)
    Pr_l = s["Pr_l"]
    k_l = s["k_l"]
    rho_l = s["rho_l"]
    mu_l = s["mu_l"]
    rho_v = s["rho_v"]

    G_eq = G * ((1 - x) + x * (rho_l / rho_v) ** 0.5)
    Re_eq = G_eq * Dh / mu_l

    Ge1 = 11.22 * (p_co / Dh) ** (-2.83) * theta ** (-4.5)
    Ge2 = 0.35 * (p_co / Dh) ** 0.23 * theta ** 1.48
    Nu = Ge1 * Re_eq ** Ge2 * Pr_l ** (1 / 3)
    return Nu * k_l / Dh

def BPHE_longo_condensing_HTC(AS, geom, cond):
    """
    Condensation HTC in a brazed plate heat exchanger, in the simplified form used
    by ACHP (Bell, ACHP docs, Plate-Heat-Exchangers, eqs. 21-23), built on:
      Longo, Gasparella & Sartori (2004), IJHMT 47, 4125-4136;
      Longo (2010), IJHMT 53, 1079-1087; Longo (2010), IJR 33, 944-953.
    j = 60 for Re_eq < 1750, linear in Re_eq above (no upper limit in ACHP);
    Nu = j Pr_L^(1/3).  NOT the Longo et al. (2015) model.
    D_h: ACHP uses d_h = 4*amplitude/Phi (2b/Phi).
    """
    G = cond['G']    
    x = cond['x']
    p_sat = cond['p']
    
    Dh = geom['Dh']

    # Saturated properties
    s = _sat_props(AS, p_sat)
    Pr_l = s["Pr_l"]
    k_l = s["k_l"]
    rho_l = s["rho_l"]
    mu_l = s["mu_l"]
    rho_v = s["rho_v"]

    # Equivalent Reynolds number (Akers): two-phase flow converted to an all-liquid flow
    Re_eq = G * ((1 - x) + x * np.sqrt(rho_l / rho_v)) * Dh / mu_l

    if Re_eq < 1750:
        Nu = 60 * Pr_l ** (1 / 3)
    else:
        if Re_eq > 3000:
            warnings.warn("Longo condensation: Re_eq > 3000, linear law extrapolated.")
        Nu = ((75 - 60) / (3000 - 1750) * (Re_eq - 1750) + 60) * Pr_l ** (1 / 3)

    return Nu * k_l / Dh