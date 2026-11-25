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
from labothappy.correlations.convection.correlation_warnings import (
    AssumedParameterWarning, ValidityRangeWarning)

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
                      "30-65 deg): extrapolated or clamped.", ValidityRangeWarning)
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
        warnings.warn(f"Re = {Re:.1f} outside the Thonon range [50, 15000]: extrapolated.",
                      ValidityRangeWarning)
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
        warnings.warn("Muley-Manglik HTC: 'phi' not given in geom, using 1.17.",
                      AssumedParameterWarning)
    phi = geom.get('phi', 1.17)

    Re = G * Dh / mu
    if Re < 1000 or not 30 <= beta <= 60 or not 1.0 <= phi <= 1.5:
        warnings.warn("Muley-Manglik HTC: outside the validated range (Re >= 1000, "
                      "30-60 deg, 1 <= phi <= 1.5): extrapolated.", ValidityRangeWarning)

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
                      "0-80 deg): extrapolated.", ValidityRangeWarning)
    zeta = _zeta_martin(Re, alpha)                       # log10 version
    visc = (mu / mu_w) ** (1/6) if mu_w else 1.0
    Nu = 0.122 * Pr**(1/3) * visc * (zeta * Re**2 * np.sin(2*alpha))**0.374

    return k / Dh * Nu

def BPHE_khan_khan_SP_HTC(AS, geom, cond):
    """
    Single-phase heat transfer coefficient for chevron plate heat exchangers.
    Khan, T.S., Khan, M.S., Chyu, M.-C., Ayub, Z.H. (2010). Experimental
    investigation of single phase convective heat transfer coefficient in a
    corrugated plate heat exchanger for multiple plate configurations.
    Applied Thermal Engineering 30(8), 1058-1065.
    DOI: 10.1016/j.applthermaleng.2010.01.021

    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.conv_plate.Nu_plate_Khan_Khan, which cites the reference above.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Not used; kept for common signature.
    geom : dict
        'Dh' : hydraulic diameter of the channel [m]
        'beta' : chevron angle from the main flow direction [rad or deg]
    cond : dict
        'G' : mass flux in ONE channel [kg/m2/s]
        'mu' : dynamic viscosity (bulk) [Pa s]
        'Pr' : Prandtl number (bulk) [-]
        'k' : thermal conductivity (bulk) [W/m/K]

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K].

    Definitions
    -----------
    Reynolds number based on Dh with bulk properties.
    Nusselt number based on Dh.
    No wall-viscosity correction (the paper only recommends a blanket
    Sieder-Tate (mu/mu_w)^0.14; the ht implementation has none).

    Validity (warning if outside)
    -----------------------------
    500 <= Re <= 2500, 30 <= beta <= 60 deg, 3.5 <= Pr <= 6.
    Tested fluid: water.
    """
    # 1) Extract conditions
    G = cond['G']
    mu = cond['mu']
    Pr = cond['Pr']
    k = cond['k']

    # 2) Extract geometry
    Dh = geom['Dh']
    beta = chevron_deg(geom['beta'])

    # 3) Reynolds number (Khan et al. 2010, as implemented in ht)
    Re = G * Dh / mu

    # 4) Angle ratio (Khan et al. 2010, as implemented in ht)
    beta_max = 60.0
    beta_ratio = beta / beta_max

    # 5) Coefficient C (Khan et al. 2010, as implemented in ht)
    C = 0.0161 * beta_ratio + 0.1298

    # 6) Exponent m (Khan et al. 2010, as implemented in ht)
    m = 0.198 * beta_ratio + 0.6398

    # 7) Nusselt number (Khan et al. 2010, as implemented in ht)
    Nu = C * Re**m * Pr**0.35

    # 8) Heat transfer coefficient (Khan et al. 2010, as implemented in ht)
    h = Nu * k / Dh

    # 9) Validity warnings (one warning for all three conditions)
    if not (500 <= Re <= 2500 and 30 <= beta <= 60 and 3.5 <= Pr <= 6):
        warnings.warn(
            "Khan-Khan HTC: outside validated range "
            "(500 <= Re <= 2500, 30-60 deg, 3.5 <= Pr <= 6, water): extrapolated.",
            ValidityRangeWarning
        )

    return h

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
    Validity: ranges of the heat transfer database, Amalfi (2016) EPFL thesis 6856,
    Table 7.3 (microscale, Bd < 4): Bo 2.97e-5-4.05e-3, Bd 1.89-3.76, rho* 77.5-1350,
    We_m 0.0267-41.5, beta* 0.429-0.929;
    Table 7.2 (macroscale, Bd >= 4): Bo 1.15e-4-3.75e-3, Bd 4.33-78.9, rho* 19.1-128,
    Re_lo 83.8-5360, Re_v 7.94-34500, beta* 0.400-1.00. Warning outside (beta* compared at the 3
    digits of the tables: 0.429 = 30/70 deg).
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
        # Validity: microscale database (thesis Table 7.3)
        out = []
        if not 2.97e-5 <= Bo <= 4.05e-3:
            out.append(f"Bo = {Bo:.3g} (2.97e-5-4.05e-3)")
        if not 1.89 <= Bd <= 3.76:
            out.append(f"Bd = {Bd:.3g} (1.89-3.76)")
        if not 77.5 <= rho_star <= 1350:
            out.append(f"rho* = {rho_star:.3g} (77.5-1350)")
        if not 0.0267 <= We <= 41.5:
            out.append(f"We_m = {We:.3g} (0.0267-41.5)")
        if not 0.429 <= round(beta_star, 3) <= 0.929:
            out.append(f"beta* = {beta_star:.3g} (0.429-0.929)")
        if out:
            warnings.warn("Amalfi boiling (microscale): outside the database: " + ", ".join(out)
                          + ": extrapolated.", ValidityRangeWarning)
        Nu = 982.0 * beta_star**1.101 * We**0.315 * Bo**0.320 * rho_star**-0.224
    else:
        if not np.isfinite(mu_v):
            raise ValueError("Amalfi boiling (Bd >= 4) needs the vapour viscosity, "
                             "not available for this fluid.")
        Re_g = G * x * Dh / mu_v
        Re_lo = G * Dh / mu_l
        # Validity: macroscale database (thesis Table 7.2)
        out = []
        if not 1.15e-4 <= Bo <= 3.75e-3:
            out.append(f"Bo = {Bo:.3g} (1.15e-4-3.75e-3)")
        if not 4.33 <= Bd <= 78.9:
            out.append(f"Bd = {Bd:.3g} (4.33-78.9)")
        if not 19.1 <= rho_star <= 128:
            out.append(f"rho* = {rho_star:.3g} (19.1-128)")
        if not 83.8 <= Re_lo <= 5360:
            out.append(f"Re_lo = {Re_lo:.3g} (83.8-5360)")
        if not 7.94 <= Re_g <= 34500:
            out.append(f"Re_v = {Re_g:.3g} (7.94-34500)")
        if not 0.400 <= round(beta_star, 3) <= 1.00:
            out.append(f"beta* = {beta_star:.3g} (0.400-1.00)")
        if out:
            warnings.warn("Amalfi boiling (macroscale): outside the database: " + ", ".join(out)
                          + ": extrapolated.", ValidityRangeWarning)
        Nu = (18.495 * beta_star**0.248 * Re_g**0.135 * Re_lo**0.351
              * Bd**0.235 * Bo**0.198 * rho_star**-0.223)
    return Nu * k_l / Dh

def BPHE_han_lee_kim_boiling_HTC(AS, geom, cond):
    """
    Han, Lee & Kim (2003) flow boiling (evaporation) in brazed plate heat exchangers.
    Applied Thermal Engineering 23(10), 1209-1225. doi:10.1016/S1359-4311(03)00061-9

    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.boiling_plate.h_boiling_Han_Lee_Kim, which cites the reference above
    (reviewed in Amalfi et al. 2016 Part 1, IJR 61, 166-184, doi:10.1016/j.ijrefrig.2015.07.010).

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state (left in a saturated state on return; caller must re-set it).
    geom : dict
        Geometry dictionary with keys:
        - 'Dh' : hydraulic diameter [m]
        - 'Lambda' : corrugation pitch (wavelength) [m]
        - 'beta' : chevron angle from the main flow direction [rad or deg]
          (Han's paper measures its angle beta_H from the horizontal axis and writes
          (pi/2 - beta_H): that is exactly theta used here.)
    cond : dict
        Operating conditions with keys:
        - 'G' : mass flux in ONE channel [kg/m2/s]
        - 'x' : vapour quality of the cell [-]
        - 'q' : local heat flux [W/m2]
        - 'p' : saturation pressure [Pa]

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K].

    Validity (warning if outside)
    -----------------------------
    G 13-34 kg/m2/s; x 0.15-0.9; q 2500-8500 W/m2; theta 45-70 deg from the flow
    direction (plates of 45, 35 and 20 deg in Han's convention); fluids R410A and R22;
    evaporating temperature 5-15 degC (not checked).

    Definitions
    -----------
    G is the mass flux in one channel (flow area = 2 x amplitude x plate width).
    Re_eq and Bo_eq use the equivalent mass flux G_eq (Akers-type) and saturated-liquid
    properties. Pr is the Prandtl number of the saturated liquid.
    """
    # 1) Read conditions
    G = cond['G']           # [kg/m2/s] mass flux in one channel
    x = cond['x']           # [-] vapour quality
    q = cond['q']           # [W/m2] local heat flux
    p_sat = cond['p']       # [Pa] saturation pressure

    # 2) Read geometry
    Dh = geom['Dh']         # [m] hydraulic diameter
    p_co = geom['Lambda']     # [m] corrugation pitch (wavelength)
    # Chevron angle from the main flow direction [rad]
    # Han's beta_H is from horizontal; theta = pi/2 - beta_H = angle from flow
    theta = np.radians(chevron_deg(geom['beta']))

    # 3) Angle validity (same as Han condensation)
    if not 0.0 < theta <= np.pi / 2:
        raise ValueError(f"Han-Lee-Kim boiling: angle {theta:.3f} rad from the flow "
                         "direction outside ]0, pi/2].")

    # 4) Quality check and clamp (as in Amalfi)
    if not np.isfinite(x):
        raise ValueError("Han-Lee-Kim boiling: cond['x'] must be finite "
                         "(two-phase cell expected).")
    x = min(max(x, 1e-4), 1.0 - 1e-4)

    # 5) Saturated properties at p_sat
    s = _sat_props(AS, p_sat)
    k_l = s["k_l"]          # [W/m/K] liquid thermal conductivity
    rho_l = s["rho_l"]      # [kg/m3] liquid density
    rho_v = s["rho_v"]      # [kg/m3] vapour density
    mu_l = s["mu_l"]        # [Pa s] liquid dynamic viscosity
    Pr_l = s["Pr_l"]        # [-] liquid Prandtl number
    h_lv = s["h_lv"]        # [J/kg] latent heat

    # 6) Equivalent mass flux (Akers-type)
    G_eq = G * ((1.0 - x) + x * (rho_l / rho_v) ** 0.5)   # [kg/m2/s]

    # 7) Equivalent Reynolds number
    Re_eq = G_eq * Dh / mu_l                               # [-]
    
    # 8) Equivalent Boiling number with guard (q=0 would give h=0)
    Bo_eq = q / (G_eq * h_lv)                              # [-]
    Bo_eq = max(Bo_eq, 1e-6)

    # 9) Geometry coefficients (Eq. from paper)
    Ge1 = 2.81 * (p_co / Dh) ** (-0.041) * theta ** (-2.83)
    Ge2 = 0.746 * (p_co / Dh) ** (-0.082) * theta ** (0.61)

    # 10) Heat transfer coefficient (direct, not through Nu)
    h = Ge1 * (k_l / Dh) * (Re_eq ** Ge2) * (Pr_l ** 0.4) * (Bo_eq ** 0.3)

    # 11) Validity warnings
    if not (13.0 <= G <= 34.0):
        warnings.warn(f"Han-Lee-Kim boiling: G = {G:.1f} kg/m2/s outside "
                      "validated range [13, 34]: extrapolated.", ValidityRangeWarning)
    if not (0.15 <= x <= 0.9):
        warnings.warn(f"Han-Lee-Kim boiling: x = {x:.3f} outside "
                      "validated range [0.15, 0.9]: extrapolated.", ValidityRangeWarning)
    if not (2500.0 <= q <= 8500.0):
        warnings.warn(f"Han-Lee-Kim boiling: q = {q:.0f} W/m2 outside "
                      "validated range [2500, 8500]: extrapolated.", ValidityRangeWarning)
    theta_deg = np.degrees(theta)
    if not (45.0 <= theta_deg <= 70.0):
        warnings.warn(f"Han-Lee-Kim boiling: theta = {theta_deg:.1f} deg from flow "
                      "direction outside validated range [45, 70]: extrapolated.",
                      ValidityRangeWarning)

    return h

def BPHE_yan_lin_boiling_HTC(AS, geom, cond):
    """
    Yan & Lin (1999) flow boiling heat transfer coefficient for plate heat exchangers.
    Journal of Heat Transfer 121(1), 118-127. DOI: 10.1115/1.2825924
    
    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.boiling_plate.h_boiling_Yan_Lin, which cites the reference above.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state object (left in a saturated state on return).
    geom : dict
        Geometry dictionary with keys:
        - 'Dh' : hydraulic diameter [m]
        (chevron angle 'beta' is present in geom but NOT used by this correlation)
    cond : dict
        Operating conditions with keys:
        - 'G' : mass flux in ONE channel [kg/m2/s]
        - 'x' : mean vapour quality of the cell [-]
        - 'q' : local heat flux [W/m2]
        - 'p' : saturation pressure [Pa]

    Returns
    -------
    h : float
        Flow boiling heat transfer coefficient [W/m2/K].

    Validity range (claimed by authors)
    -----------------------------------
    2000 <= Re_eq <= 10000.
    Tested with R134a, one plate with chevron angle 60 deg, x 0.1-0.8,
    q 11-15 kW/m2, G 55-70 kg/m2/s.

    Definitions
    -----------
    - G is the mass flux in ONE channel.
    - Re_eq uses G_eq and saturated-liquid properties.
    - Bo_eq uses G_eq and saturated-liquid properties.
    - Re_lo uses the total mass flux G as if all liquid.
    - Pr is the Prandtl number of the saturated liquid.
    - geom['beta'] is not used in this correlation.
    """
    # 1) Read inputs from dictionaries
    G = cond['G']           # mass flux in one channel [kg/m2/s]
    x = cond['x']           # vapour quality [-]
    q = cond['q']           # heat flux [W/m2]
    p_sat = cond['p']       # saturation pressure [Pa]
    Dh = geom['Dh']         # hydraulic diameter [m]

    # 2) Quality check and clamp (same approach as Amalfi)
    if not np.isfinite(x):
        raise ValueError("Yan-Lin boiling: cond['x'] must be finite (two-phase cell expected).")
    x = min(max(x, 1e-4), 1.0 - 1e-4)

    # 3) Saturated properties at p_sat
    s = _sat_props(AS, p_sat)
    k_l = s["k_l"]          # liquid thermal conductivity [W/m/K]
    rho_l = s["rho_l"]      # liquid density [kg/m3]
    rho_v = s["rho_v"]      # vapour density [kg/m3]
    mu_l = s["mu_l"]        # liquid dynamic viscosity [Pa*s]
    Pr_l = s["Pr_l"]        # liquid Prandtl number [-]
    h_lv = s["h_lv"]        # latent heat [J/kg]

    # 4) Equivalent mass flux (Eq. in Yan & Lin 1999)
    G_eq = G * ((1.0 - x) + x * (rho_l / rho_v) ** 0.5)

    # 5) Equivalent Reynolds number
    Re_eq = G_eq * Dh / mu_l

    # 6) All-liquid Reynolds number (total mass flux as liquid)
    Re_lo = G * Dh / mu_l

    # 7) Equivalent Boiling number with guard against zero heat flux
    Bo_eq = q / (G_eq * h_lv)
    Bo_eq = max(Bo_eq, 1e-6)

    # 8) Heat transfer coefficient (Yan & Lin 1999 correlation)
    h = (1.926
         * (k_l / Dh)
         * Re_eq
         * Pr_l ** (1.0 / 3.0)
         * Bo_eq ** 0.3
         * Re_lo ** (-0.5))

    # 9) Validity warning (claimed range 2000 <= Re_eq <= 10000)
    if not (2000 <= Re_eq <= 10000):
        warnings.warn(
            f"Yan-Lin boiling: Re_eq = {Re_eq:.1f} outside the claimed range "
            "[2000, 10000]. Correlation extrapolated.", ValidityRangeWarning
        )

    return h

def BPHE_lee_kang_kim_boiling_HTC(AS, geom, cond):
    """
    Lee, Kang & Kim (2014) flow boiling heat transfer correlation for plate heat exchangers.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state object (must be set to the working fluid).
    geom : dict
        Geometry dictionary with key:
        - 'Dh' : hydraulic diameter [m]. In this library Dh = 2*b (plate gap).
          The correlation uses the equivalent diameter D_eq = 4*a = 2*b,
          which is exactly geom['Dh'].
    cond : dict
        Operating conditions with keys:
        - 'G' : mass flux in ONE channel [kg/m2/s]
        - 'x' : mean vapour quality of the cell [-]
        - 'q' : local heat flux [W/m2]
        - 'p' : saturation pressure [Pa]

    Returns
    -------
    h : float
        Flow boiling heat transfer coefficient [W/m2/K].

    Reference
    ---------
    Lee, E., Kang, H., Kim, Y. (2014). Flow boiling heat transfer and pressure drop
    of water in a plate heat exchanger with corrugated channels at low mass flux
    conditions. International Journal of Heat and Mass Transfer, 77, 37-45.
    doi:10.1016/j.ijheatmasstransfer.2014.05.019

    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.boiling_plate.h_boiling_Lee_Kang_Kim, which cites the reference
    above (reviewed in Amalfi et al. 2016 Part 1, IJR 61, 166-184,
    doi:10.1016/j.ijrefrig.2015.07.010).

    Validity range (warning if outside)
    -----------------------------------
    G: 14.5-33.6 kg/m2/s; q: 15000-30000 W/m2; x: 0.09-0.6;
    Re_ratio: 2.3-32.1; Bo: 0.00019-0.001; X_tt: 0.028-0.3.
    Fluid: water. Also stated: 200 < Re < 600, mean absolute deviation 4.4 %.
    (Re definition not given, so not checked.)

    Definitions
    -----------
    - G is the mass flux in ONE channel.
    - Bo = q / (G * h_lv) (Boiling number).
    - Re_ratio = Re_g / Re_l = (x * mu_l) / ((1 - x) * mu_v).
    - X_tt = ((1 - x) / x)^0.875 * (rho_v / rho_l)^0.5 * (mu_l / mu_v)^0.125
      (Lockhart-Martinelli parameter, turbulent-turbulent).
    - Saturated properties at p_sat: k_l, rho_l, rho_v, mu_l, mu_v, h_lv.
    - D_eq = geom['Dh'] = 2 * b (plate gap).
    - geom['beta'] is not used by this correlation.
    """
    # 1) Read inputs from cond and geom
    G = cond['G']           # [kg/m2/s] mass flux in one channel
    x = cond['x']           # [-] vapour quality
    q = cond['q']           # [W/m2] heat flux
    p_sat = cond['p']       # [Pa] saturation pressure
    D_eq = geom['Dh']       # [m] equivalent diameter (D_eq = 4a = 2b)

    # 2) Quality check and clamp (avoid division by zero)
    if not np.isfinite(x):
        raise ValueError("Lee-Kang-Kim boiling: cond['x'] must be finite.")
    x = min(max(x, 1e-4), 1.0 - 1e-4)

    # 3) Saturated properties at p_sat
    s = _sat_props(AS, p_sat)
    k_l = s["k_l"]          # [W/m/K] liquid thermal conductivity
    rho_l = s["rho_l"]      # [kg/m3] liquid density
    rho_v = s["rho_v"]      # [kg/m3] vapour density
    mu_l = s["mu_l"]        # [Pa s] liquid dynamic viscosity
    mu_v = s["mu_v"]        # [Pa s] vapour dynamic viscosity
    h_lv = s["h_lv"]        # [J/kg] latent heat

    # 4) Vapour viscosity must be available (finite)
    if not np.isfinite(mu_v):
        raise ValueError("Lee-Kang-Kim boiling: vapour viscosity not available "
                         "for this fluid at the given pressure.")

    # 5) Boiling number with guard against q = 0
    Bo = q / (G * h_lv)     # [-] Boiling number
    Bo = max(Bo, 1e-6)      # guard: q = 0 would give h = 0

    # 6) Ratio of vapour-alone to liquid-alone Reynolds numbers
    # Re_ratio = Re_g / Re_l = (x / (1 - x)) * (mu_l / mu_v)
    Re_ratio = (x / (1.0 - x)) * (mu_l / mu_v)

    # 7) Lockhart-Martinelli parameter (turbulent-turbulent, paper exponents)
    # X_tt = ((1 - x) / x)^0.875 * (rho_v / rho_l)^0.5 * (mu_l / mu_v)^0.125
    X_tt = ((1.0 - x) / x) ** 0.875 * (rho_v / rho_l) ** 0.5 * (mu_l / mu_v) ** 0.125

    # 8) Two branches depending on Re_ratio
    if Re_ratio < 9.0:
        # Equation for Re_ratio < 9
        h = (98.7 * (k_l / D_eq) * Re_ratio ** -0.0848 *
             Bo ** -0.0597 * X_tt ** 0.0973)
    else:
        # Equation for Re_ratio >= 9
        h = (234.9 * (k_l / D_eq) * Re_ratio ** -0.576 *
             Bo ** -0.275 * X_tt ** 0.66)

    # 9) Validity range warnings (list all parameters that are outside)
    warnings_list = []
    if not (14.5 <= G <= 33.6):
        warnings_list.append(f"G = {G:.1f} kg/m2/s (range 14.5-33.6)")
    if not (15000 <= q <= 30000):
        warnings_list.append(f"q = {q:.0f} W/m2 (range 15000-30000)")
    if not (0.09 <= x <= 0.6):
        warnings_list.append(f"x = {x:.3f} (range 0.09-0.6)")
    if not (2.3 <= Re_ratio <= 32.1):
        warnings_list.append(f"Re_ratio = {Re_ratio:.2f} (range 2.3-32.1)")
    if not (0.00019 <= Bo <= 0.001):
        warnings_list.append(f"Bo = {Bo:.6f} (range 0.00019-0.001)")
    if not (0.028 <= X_tt <= 0.3):
        warnings_list.append(f"X_tt = {X_tt:.4f} (range 0.028-0.3)")

    if warnings_list:
        msg = ("Lee-Kang-Kim boiling: outside validated range: " +
               "; ".join(warnings_list) + ". Fluid: water.")
        warnings.warn(msg, ValidityRangeWarning)

    return h

def BPHE_huang_sheer_boiling_HTC(AS, geom, cond):
    """
    Huang, Sheer & Bailey-McEwan (2012) flow boiling heat transfer coefficient
    for plate heat exchangers. International Journal of Refrigeration 35(2),
    325-335. DOI: 10.1016/j.ijrefrig.2011.11.002.
    Also: Huang, J. (2011) Performance analysis of plate heat exchangers used
    as refrigerant evaporators, PhD thesis, University of the Witwatersrand.

    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.boiling_plate.h_boiling_Huang_Sheer, which cites the
    references above. As noted by ht, the printed equation of the 2012 paper
    misses a square on the bubble diameter in the group
    (h_lv d_o^2 / alpha_l^2); the corrected, dimensionless form is used here.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        CoolProp state object set to the working fluid. Must be already
        initialised with the fluid name and backend.
    geom : dict
        Geometry dictionary (not used by this correlation; present for API
        compatibility). Keys read: none.
    cond : dict
        Operating conditions dictionary with keys:
        - 'q' : local heat flux [W/m2]
        - 'p' : saturation pressure [Pa]

    Returns
    -------
    h : float
        Flow boiling heat transfer coefficient [W/m2/K].

    Validity (data ranges from the paper)
    -------------------------------------
    - Heat flux q: 1850 - 10750 W/m2
    - Saturation temperature T_sat: 1.9 - 13.04 degC (275.05 - 286.19 K)
    - Fluids tested: mainly R134a and R507A (222 points), two points for
      ammonia and R12
    - Chevron angles: 28 - 60 deg
    - Mass flux G: 5.6 - 52.25 kg/m2/s
    - Quality x: 0.21 - 0.95

    The correlation depends only on fluid saturated properties, heat flux q,
    and saturation temperature T_sat. It does NOT use geometry, mass flux G,
    or quality x. Because T_sat appears explicitly, the correlation should be
    used with caution for fluids or temperatures far from the experimental data.

    Definitions
    -----------
    - d_o : bubble departure diameter [m]
    - alpha_l : thermal diffusivity of saturated liquid [m2/s]
    - Properties at saturation: k_l, rho_l, rho_v, cp_l, Pr_l, h_lv, sigma
    """
    # 1) Inputs from cond dictionary
    q = cond['q']                 # local heat flux [W/m2]
    p_sat = cond['p']             # saturation pressure [Pa]

    # Fixed parameters used in the correlation development
    contact_angle = 35.0          # contact angle [deg] (used as a number in degrees)
    g = 9.80665                   # gravitational acceleration [m/s2]

    # 2) Guard on heat flux (same style as Cooper)
    if not (np.isfinite(q) and q > 0.0):
        raise ValueError(f"Huang-Sheer boiling: heat flux q = {q} must be finite and > 0.")

    # 3) Saturated properties at p_sat
    s = _sat_props(AS, p_sat)
    k_l = s["k_l"]                # liquid thermal conductivity [W/m/K]
    rho_l = s["rho_l"]            # liquid density [kg/m3]
    rho_v = s["rho_v"]            # vapour density [kg/m3]
    cp_l = s["cp_l"]              # liquid specific heat [J/kg/K]
    Pr_l = s["Pr_l"]              # liquid Prandtl number [-]
    h_lv = s["h_lv"]              # latent heat [J/kg]
    sigma = s["sigma"]            # surface tension [N/m]

    # 4) Saturation temperature [K]
    # Do this AFTER _sat_props because _sat_props may return cached values
    # without touching AS. We need AS in a known state for T().
    AS.update(CP.PQ_INPUTS, p_sat, 0)
    T_sat = AS.T()                # saturation temperature [K]

    # 5) Bubble departure diameter d_o [m]
    # Equation from the paper (corrected form as used in ht library)
    d_o = 0.0146 * contact_angle * (2.0 * sigma / (g * (rho_l - rho_v))) ** 0.5

    # 6) Thermal diffusivity of the liquid alpha_l [m2/s]
    alpha_l = k_l / (rho_l * cp_l)

    # 7) Heat transfer coefficient h [W/m2/K]
    # Write the correlation in several steps with intermediate variables
    # for readability (one operation per line).
    group1 = k_l / d_o
    group2 = (q * d_o) / (k_l * T_sat)
    group3 = (h_lv * d_o ** 2) / (alpha_l ** 2)

    h = (1.87e-3
         * group1
         * (group2 ** 0.56)
         * (group3 ** 0.31)
         * (Pr_l ** 0.33))

    # 8) Validity warnings (do not stop, just warn)
    if not (1850.0 <= q <= 10750.0):
        warnings.warn(
            f"Huang-Sheer boiling: heat flux q = {q:.1f} W/m2 outside the "
            f"validated range [1850, 10750] W/m2: extrapolated.", ValidityRangeWarning
        )
    if not (275.05 <= T_sat <= 286.19):
        warnings.warn(
            f"Huang-Sheer boiling: saturation temperature T_sat = {T_sat:.2f} K "
            f"({T_sat - 273.15:.2f} degC) outside the validated range "
            f"[275.05, 286.19] K: extrapolated. "
            f"Correlation should be used with caution for fluids or temperatures "
            f"far from the experimental data (mainly R134a, R507A).", ValidityRangeWarning
        )

    return h

#%% TWO-PHASE : CONDENSING

def BPHE_shah_condensing_HTC(AS, geom, cond):
    """
    Shah (2021) condensation in plate heat exchangers, forced-convection branch from
    Longo et al. (2015). Checked against the paper 2026-11-21.

    Reference: Shah, M.M. (2021) Heat transfer during condensation in corrugated plate heat
    exchangers, Int. J. Refrigeration 127, 180-193, Eqs. (18), (19), (21), (36), section 3:
    h = larger of h_grav (36) and h_fc (21) if Re_eq < 1600, h = h_fc if Re_eq >= 1600.

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality of the cell) [-],
           p (saturation pressure) [Pa]
    geom : Dh (hydraulic diameter D_hyd = 2b, as given by the BPHE geometry) [m],
           phi (surface enlargement factor, NOT the chevron angle) [-], default 1.17
    Returns h [W/m2/K] on the ACTUAL (enlarged) heat transfer area.

    Definitions: Re_eq = G_eq D_hyd / mu_f, G_eq = G [(1 - x) + x (rho_f/rho_g)^0.5],
    Re_LS = G (1 - x) D_hyd / mu_f, D_hyd = 2b (paper: "the diameter used is D_hyd"). The paper's
    h is based on the PROJECTED area (factor phi in Eqs. 21, 36); the labothappy BPHE area is the
    actual area (phi = A_h / A_projected), so h_paper / phi is returned (= Eqs. 21, 36 without phi).
    Validity: chevron 30-75 deg (MAD 19.5 % for 40-65 deg), 18 fluids, p_r 0.0083-0.8,
    G 2.3-165 kg/m2/s, q 2.5-93.5 kW/m2, x 0-1.
    (2026-11-21: Dh / phi was used as diameter and the phi factor was kept with the actual area:
    h_fc +28 %, h_grav +23 % at phi = 1.17. Corrected.)
    """
    G = cond['G']
    x = cond['x']
    p_sat = cond['p']

    # D_hyd = 2b (Shah 2021: all correlations use D_hyd)
    Dh = geom['Dh']
    
    if not 2.3 <= G <= 165:
        warnings.warn("Shah condensation: mass flux outside the validated range "
                      "(2.3-165 kg/m2/s): extrapolated.", ValidityRangeWarning)

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

    # Eq. (21) divided by phi: actual-area basis
    h_fc = 1.875 * (k_l / Dh) * Re_eq**0.445 * Pr_l**(1/3)
    if Re_eq >= 1600:
        return h_fc
    # Eq. (36) divided by phi: actual-area basis
    h_grav = (1.32 * Re_l**(-1/3)
              * (rho_l * (rho_l - rho_v) * g * k_l**3 / mu_l**2) ** (1/3))
    return max(h_grav, h_fc)

def BPHE_han_condensing_HTC(AS, geom, cond):
    """
    Han, Lee & Kim (2003), J. Korean Phys. Soc. 43(1), 66-73 (to be checked).
    Condensation in brazed plate HX, Nu = Ge1 Re_eq^Ge2 Pr^(1/3).

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality) [-], p [Pa]
    geom : Dh (hydraulic diameter) [m], Lambda (corrugation pitch) [m],
           beta (chevron angle from the main flow direction, rad or deg;
           Han's own angle is measured from the horizontal axis: pi/2 - beta)
    Returns h [W/m2/K].
    """
    G = cond['G']    
    x = cond['x']
    p_sat = cond['p']
    
    Dh = geom['Dh']
    
    theta = np.radians(chevron_deg(geom["beta"]))   # from the flow direction = Han's pi/2 - beta
    p_co = geom["Lambda"]

    if not 0.0 < theta <= np.pi / 2:
        raise ValueError(f"Han condensation: angle {theta} rad from the flow direction "
                         "outside ]0, pi/2].")
    if not 13 <= G <= 34 or not 45 <= np.degrees(theta) <= 70:
        warnings.warn("Han condensation: outside the range of the data (G 13-34 kg/m2/s, "
                      "plates of 45/55/70 deg from the flow direction): extrapolated.",
                      ValidityRangeWarning)
    
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
            warnings.warn("Longo condensation: Re_eq > 3000, linear law extrapolated.",
                          ValidityRangeWarning)
        Nu = ((75 - 60) / (3000 - 1750) * (Re_eq - 1750) + 60) * Pr_l ** (1 / 3)

    return Nu * k_l / Dh