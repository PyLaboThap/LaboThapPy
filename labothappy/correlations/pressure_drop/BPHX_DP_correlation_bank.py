# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:51:41 2026

@author: basil
"""

#%% IMPORTS
import numpy as np
import warnings
from labothappy.correlations.convection.correlation_warnings import ValidityRangeWarning
from labothappy.correlations.convection.correlation_warnings import FluidWarning

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

def _fanning_tube(Re):
    """
    Smooth circular tube Fanning friction factor as used in ACHP's
    Lockhart-Martinelli implementation (laminar 16/Re, turbulent Blasius-type
    0.046 Re^-0.2, linear blend for 1000 < Re < 2000).
    """
    if Re < 1000:
        return 16 / Re
    elif Re > 2000:
        return 0.046 / Re ** 0.2
    else:
        w = (Re - 1000) / (2000 - 1000)
        return (1 - w) * 16 / Re + w * 0.046 / Re ** 0.2

#%% SINGLE PHASE : Friction DP

def BPHX_sinnott_SP_DP(AS, geom, cond):
    """
    Coulson & Richardson's Chemical Engineering Design Vol. 6
    
    R. K. Sinnott (1983)
    
    Eq. 12.78 + Eq. 12.79

    Friction term only (port loss: see plate_port_DP).
    geom['L'] must be the total path length (plate length x number of passes).
    """
    # 1) Condition Extraction
    G = cond['G']
    rho = cond['D']
    mu = cond['mu']

    # 2) Geometry Parameters
    n_passes = geom.get('n_passes', 1)
    D_port = geom.get('D_port', 0)
    Dh = geom['Dh']
    L = (geom['L']-2*D_port)*n_passes
        
    # 3) Reynolds number + Speed
    Re = G*Dh/mu
    u = G/rho
    
    # 4) Friction factor
    f = 0.6*Re**(-0.3)
    
    # 5) DP friction
    DP_f = 8 * f * (L/Dh) * rho * u**2 / 2

    return DP_f

def BPHX_kumar_SP_DP(AS, geom, cond):
    """
    Kumar, H., The Plate Heat Exchanger: Construction and Design,
    IChemE Symp. Ser. 86, 1275-1288, 1984.

        dP_f = 4 f (L/Dh) G^2/(2 rho)
        f    = C2 / Re^p   (Fanning, C2 and p depend on angle and Re)

    Friction term only (port loss: see plate_port_DP).
    Angle: same convention as kumar_plate_HTC (geom['beta']); values between
    tabulated angles are interpolated linearly, outside they are clamped.
    geom['L'] must be the total path length (plate length x number of passes).
    """
    def _kumar_f_row(rows, Re):
        for Re_max, C2, p in rows:
            if Re <= Re_max:
                return C2 / Re**p
    
    # Kumar (1984) friction factor f = C2 / Re^p (Fanning)
    # angle [deg] -> [(Re_max, C2, p), ...]; the last row has Re_max = inf
    _KUMAR_F = [
        (30, [(10,  50.0,  1.0),   (100, 19.4,  0.589), (np.inf, 2.99,  0.183)]),
        (45, [(15,  47.0,  1.0),   (300, 18.29, 0.652), (np.inf, 1.441, 0.206)]),
        (50, [(20,  34.0,  1.0),   (300, 11.25, 0.631), (np.inf, 0.772, 0.161)]),
        (60, [(40,  24.0,  1.0),   (400, 3.24,  0.457), (np.inf, 0.76,  0.215)]),
        (65, [(50,  24.0,  1.0),   (500, 2.8,   0.451), (np.inf, 0.639, 0.213)]),
    ]
    
    G, rho, mu = cond['G'], cond['D'], cond['mu']
    
    n_passes = geom.get('n_passes', 1)
    D_port = geom.get('D_port', 0)
    Dh = geom['Dh']
    L = (geom['L']-2*D_port)*n_passes
        
    beta = chevron_deg(geom['beta'])

    Re = G * Dh / mu

    # Validity: data on graph from Re = 0.1 to 10000, chevron angles 30 to 65 deg (fluids.friction.
    # friction_plate_Kumar notes); beta is already in Kumar's convention here (converted by BPHX_DP)
    if not 0.1 <= Re <= 10000 or not 30 <= beta <= 65:
        warnings.warn(f"Kumar DP: outside the range of the data (0.1 <= Re <= 10000, 30-65 deg): "
                      f"Re = {Re:.3g}, angle = {beta:.1f} deg, extrapolated or clamped.",
                      ValidityRangeWarning)

    # Friction factor: evaluated at each tabulated angle, then interpolated
    angles = [a for a, _ in _KUMAR_F]
    f_tab = [_kumar_f_row(rows, Re) for _, rows in _KUMAR_F]
    f = np.interp(beta, angles, f_tab)

    # Channel friction
    DP_f = 4 * f * (L / Dh) * G**2 / (2 * rho)

    return DP_f

def BPHX_thonon_SP_DP(AS, geom, cond):
    """
    Thonon, B., Design Method for Plate Evaporators and Condensers,
    BHR Group Conf. Ser. 18, 37-47, 1995.

        dP_f = 4 f (L/Dh) G^2/(2 rho)
        f    = C2 / Re^p   (Fanning-type, C2 and p depend on angle and Re)

    Friction term only (port loss: see plate_port_DP).
    Valid for 50 <= Re <= 15000. Angle convention as in thonon_plate_HTC;
    values between tabulated angles are interpolated linearly, outside clamped.
    geom['L'] must be the total path length (plate length x number of passes).
    """
    # Thonon (1995) friction factor f = C2 / Re^p
    # angle [deg] -> (Re_switch, (C2, p) for Re <= Re_switch, (C2, p) for Re > Re_switch)
    _THONON_F = [
        (30, 160,  (45.57, 0.670), (0.370,  0.172)),   # see note on 0.370
        (45, 200,  (18.19, 0.682), (0.6857, 0.172)),
        (60, 550,  (26.34, 0.830), (0.572,  0.217)),
        (75, 1000, (28.21, 0.900), (0.872,  0.392)),
    ]
    
    G, rho, mu = cond['G'], cond['D'], cond['mu']
    
    n_passes = geom.get('n_passes', 1)
    D_port = geom.get('D_port', 0)
    Dh = geom['Dh']
    L = (geom['L']-2*D_port)*n_passes
        
    beta = chevron_deg(geom['beta'])

    Re = G * Dh / mu
    if not 50 <= Re <= 15000:
        warnings.warn(f"Re = {Re:.1f} outside the Thonon range [50, 15000]: extrapolated.",
                      ValidityRangeWarning)

    # Friction factor at each tabulated angle, then interpolation in angle
    angles = [a for a, *_ in _THONON_F]
    f_tab = []
    for _, Re_sw, low, high in _THONON_F:
        C2, p = low if Re <= Re_sw else high
        f_tab.append(C2 / Re**p)
    f = np.interp(beta, angles, f_tab)

    # Channel friction
    DP_f = 4 * f * (L / Dh) * G**2 / (2 * rho)

    return DP_f

def BPHX_muley_manglik_SP_DP(AS, geom, cond):
    """
    Muley, A. & Manglik, R.M. (1999). J. Heat Transfer 121(1):110-117.

        f    = C3(beta) D4(phi) Re^-q(beta)      (Fanning)
        dP_f = 4 f (L/Dh) G^2/(2 rho)

    Friction term only (port loss: see plate_port_DP).
    Validité : Re >= 1000, 30 <= beta <= 60 deg, 1 <= phi <= 1.5.
    beta mesuré depuis la direction de l'écoulement (comme pour le Nu).
    geom['L'] = longueur totale du chemin (longueur de plaque x passes).
    geom['phi'] optionnel (défaut 1.17).
    """
    G, rho, mu = cond['G'], cond['D'], cond['mu']
    
    n_passes = geom.get('n_passes', 1)
    D_port = geom.get('D_port', 0)
    Dh = geom['Dh']
    L = (geom['L']-2*D_port)*n_passes
 
    beta = chevron_deg(geom['beta'])
    phi = geom.get('phi', 1.17)

    Re = G * Dh / mu
    if Re < 1000 or not 30 <= beta <= 60:
        warnings.warn(f"Muley-Manglik DP: Re = {Re:.0f}, beta = {beta:.0f} deg "
                      "hors du domaine (Re >= 1000, 30-60 deg): extrapolé.", ValidityRangeWarning)

    C3 = 2.917 - 0.1277*beta + 2.016e-3*beta**2
    D4 = 5.474 - 19.02*phi + 18.93*phi**2 - 5.341*phi**3
    q = 0.2 + 0.0577*np.sin(2*np.pi*beta/90 + 2.1)
    f = C3 * D4 * Re**(-q)

    DP_f = 4 * f * (L / Dh) * G**2 / (2 * rho)

    return DP_f

def BPHX_martin_SP_DP(AS, geom, cond):
    """
    Martin, H. (1996), Chem. Eng. Process. 35, 301-310; VDI Heat Atlas N6.
    Frictional pressure drop over the full flow length. Port losses are NOT
    included (see plate_port_DP).
    geom['beta'] is read through chevron_deg and must follow Martin's
    convention (from the main flow direction).
    """
    G, rho, mu = cond['G'], cond['D'], cond['mu']
   
    n_passes = geom.get('n_passes', 1)
    D_port = geom.get('D_port', 0)
    Dh = geom['Dh']/geom['phi']
    L = (geom['L']-2*D_port)*n_passes
    
    beta = np.radians(chevron_deg(geom['beta']))

    Re = G * Dh / mu
    # Validity: experimental data Re 200-10000, chevron angles 0-80 deg (fluids.friction.
    # friction_plate_Martin_1999 notes)
    if not 200 <= Re <= 10000 or not 0 <= np.degrees(beta) <= 80:
        warnings.warn(f"Martin DP: outside the range of the data (200 <= Re <= 10000, 0-80 deg): "
                      f"Re = {Re:.3g}, angle = {np.degrees(beta):.1f} deg, extrapolated.",
                      ValidityRangeWarning)
    zeta = _zeta_martin(Re, beta)
    dp = zeta * (L / Dh) * G**2 / (2 * rho)

    return dp

#%% TWO-PHASE : BOILING + CONDENSING

def BPHX_lockhart_martinelli_DP(AS, geom, cond):
    """
    Two-phase FRICTIONAL pressure drop with the Lockhart-Martinelli separated
    flow model, in the analytical form of Chisholm (1967).

    Original validity: adiabatic two-phase flow in horizontal, smooth circular
    tubes (Lockhart & Martinelli, 1949). Applying it to a plate heat exchanger
    is an extrapolation.

    cond : G (total mass flux in ONE channel) [kg/m2/s], x (mean quality of the cell) [-],
           p (saturation pressure) [Pa]
    geom : Dh [m] (2*b/enlargement_factor for a PHE),
           L (TOTAL flow path length: plate length x number of passes) [m],
           beta (chevron angle from the main flow direction, rad or deg),
           optional: lm_C (Chisholm constant; None -> chosen from the laminar/turbulent
                     regime of each phase: 5, 10, 12 or 20),
                     lm_p_f (default True: Martin's PHE friction factor for each phase,
                     C not refitted for plates; False: smooth tube friction factors, original method)

    Returns dp [Pa]: frictional two-phase pressure drop over the FULL flow length L
    (acceleration, gravity and port losses NOT included, see plate_port_DP).
    """
    
    G = cond['G']
    if not np.isfinite(cond['x']):
        raise ValueError("Lockhart-Martinelli DP: cond['x'] must be finite (two-phase cell expected).")
    x = min(max(float(cond['x']), 0.0), 1.0)
    p_sat = cond['p']

    Dh, L = geom['Dh'], geom['L']
    C = geom.get('lm_C')
    alpha = np.radians(chevron_deg(geom['beta'])) if geom.get('lm_p_f', True) else None

    # Saturated properties
    s = _sat_props(AS, p_sat)
    rho_l, rho_v, mu_l, mu_v = s["rho_l"], s["rho_v"], s["mu_l"], s["mu_v"]
    if not np.isfinite(mu_v):
        raise ValueError("Lockhart-Martinelli DP needs the vapour viscosity, not available for this fluid.")

    # 1. Reynolds number of each phase flowing ALONE at its actual flow rate
    Re_l = G * (1 - x) * Dh / mu_l
    Re_v = G * x * Dh / mu_v

    # 2. Single-phase frictional pressure gradients [Pa/m]
    def dpdz_single(Re, G_phase, rho):
        if G_phase <= 0:
            return 0.0
        if alpha is None:
            f_F = _fanning_tube(Re)                      # Fanning
            return 2 * f_F * G_phase**2 / (rho * Dh)
        zeta = _zeta_martin(Re, alpha)                   # Darcy
        return zeta * G_phase**2 / (2 * rho * Dh)

    dpdz_l = dpdz_single(Re_l, G * (1 - x), rho_l)
    dpdz_v = dpdz_single(Re_v, G * x, rho_v)

    # 3. Chisholm constant from the flow regime of each phase
    #    (1500 used as transition Re to ensure continuity, as in ACHP)
    if C is None:
        if Re_l > 1500 and Re_v > 1500:
            C = 20.0    # turbulent liquid - turbulent vapour
        elif Re_l < 1500 and Re_v > 1500:
            C = 12.0    # laminar liquid - turbulent vapour
        elif Re_l > 1500 and Re_v < 1500:
            C = 10.0    # turbulent liquid - laminar vapour
        else:
            C = 5.0     # laminar liquid - laminar vapour

    # 4. Two-phase gradient (Chisholm form, identical to phi_l^2*dpdz_l = phi_v^2*dpdz_v)
    dpdz = dpdz_l + C * np.sqrt(dpdz_l * dpdz_v) + dpdz_v

    return dpdz * L

#%% TWO-PHASE : BOILING 

from labothappy.correlations.convection.BPHE_HTC_correlation_bank import _sat_props

def BPHX_amalfi_boiling_DP(AS, geom, cond):
    """
    Amalfi, Vakili-Farahani & Thome (2016), international journal of refrigeration 61 (2016) 185 – 203
    Flow boiling and frictional pressure gradients in plate heat exchangers. Part 2: Comparison of literature methods to database and new prediction methods

    f_tp is a FANNING friction factor: dP = 4 f_tp (L/Dh) G^2 / (2 rho_m) = 2 f_tp (L/Dh) G^2 / rho_m
    (rho_m: homogeneous density). Returns the frictional dP [Pa] over the length L.
    """
    beta_max_deg=70.0
    g = 9.81
    
    p_sat = cond['p']
    x = cond['x']
    G = cond['G']

    if not np.isfinite(x):
        raise ValueError("Amalfi boiling: cond['x'] must be finite (two-phase cell expected).")
    x = min(max(x, 1e-4), 1.0 - 1e-4)
    
    beta = chevron_deg(geom['beta'])
    Dh, L = geom['Dh'], geom['L']
    
    # Saturated properties
    s = _sat_props(AS, p_sat)
    rho_l = s["rho_l"]
    rho_v = s["rho_v"]
    sigma = s["sigma"]

    beta_star = beta / beta_max_deg
    rho_star = rho_l / rho_v
    rho_m = 1.0 / (x / rho_v + (1.0 - x) / rho_l)       # homogeneous density
    Bd = (rho_l - rho_v) * g * Dh**2 / sigma

    We = G**2 * Dh / (rho_m * sigma)

    # Validity: ranges of the pressure drop database, Amalfi (2016) EPFL thesis 6856, Table 7.1:
    # Bd 2.40-49.1, rho* 19.1-1350, We_m 0.0267-150, beta* 0.429-0.929 (warning outside; beta* is
    # compared at the 3 digits of the table: 0.429 = 30/70 deg)
    if not 2.40 <= Bd <= 49.1:
        warnings.warn(f"Amalfi DP: Bd = {Bd:.3g} outside 2.40-49.1 (database): extrapolated.",
                      ValidityRangeWarning)
    if not 19.1 <= rho_star <= 1350:
        warnings.warn(f"Amalfi DP: rho_l/rho_v = {rho_star:.3g} outside 19.1-1350 (database): "
                      "extrapolated.", ValidityRangeWarning)
    if not 0.0267 <= We <= 150:
        warnings.warn(f"Amalfi DP: We_m = {We:.3g} outside 0.0267-150 (database): extrapolated.",
                      ValidityRangeWarning)
    if not 0.429 <= round(beta_star, 3) <= 0.929:
        warnings.warn(f"Amalfi DP: beta* = {beta_star:.3g} outside 0.429-0.929 (database): "
                      "extrapolated.", ValidityRangeWarning)
    
    C = 2.125*beta_star**9.993 + 0.955
    f_tp = C * 15.698 * We**-0.475 * Bd**0.255 * rho_star**-0.571 # fanning friction factor

    # Fanning notation: dP = 4 f (L/Dh) G^2 / (2 rho_m)
    return 4.0 * f_tp * (L / Dh) * G**2 / (2.0 * rho_m)


#%% TWO-PHASE : CONDENSING 

def BPHX_han_condensing_DP(AS, geom, cond):
    """
    Han, Lee & Kim (2003), J. Korean Phys. Soc. 43(1), 66-73 (to be checked).
    Frictional pressure drop of a condensing flow in a brazed plate heat exchanger.
    Port losses are NOT included (see plate_port_DP).

        f    = Ge3 Re_eq^Ge4
        Ge3  = 3521.1  (p_co/Dh)^4.17   theta^-7.75
        Ge4  = -1.024  (p_co/Dh)^0.0925 theta^-1.3
        dP_f = f (L/Dh) G^2 / rho_m

    cond : G (mass flux in ONE channel) [kg/m2/s], x (mean quality of the cell) [-], p [Pa]
    geom : Dh [m], Lambda (corrugation pitch) [m],
           beta (chevron angle from the main flow direction, rad or deg),
           L (TOTAL flow path length: plate length x number of passes) [m]
    Returns dP [Pa] over the FULL flow length.

    Differences with the previous version (to be checked against the paper):
      - Re_eq = G_eq Dh / mu_l with G_eq = G [(1-x) + x sqrt(rho_l/rho_v)] (equivalent Reynolds number,
        as in the Han condensation HTC); before: G Dh / mu.
      - rho_m = homogeneous density 1 / (x/rho_v + (1-x)/rho_l), and G^2 / rho_m; before: G^2 (rho_v + rho_l),
        which is not a pressure (kg^3 m^-7 s^-2).
      - L is the total path length, as in the other plate DP functions; before: L_v N_cp, where N_cp (number
        of channels per pass) multiplied the pressure drop by the number of channels.
      - theta = angle from the flow direction (Han's pi/2 - beta), as in the Han condensation HTC.
    """
    G = cond['G']
    x = min(max(cond['x'], 0.0), 1.0)
    p_sat = cond['p']

    Dh, L, p_co = geom['Dh'], geom['L'], geom['Lambda']
    theta = np.radians(chevron_deg(geom['beta']))            # from the flow direction = Han's pi/2 - beta
    if not 0.0 < theta <= np.pi / 2:
        raise ValueError(f"Han condensation DP: angle {theta} rad from the flow direction "
                         "outside ]0, pi/2].")
    if not 13 <= G <= 34 or not 45 <= np.degrees(theta) <= 70:
        warnings.warn("Han condensation DP: outside the range of the Han data (G 13-34 kg/m2/s, "
                      "45-70 deg from the flow direction): extrapolated.", ValidityRangeWarning)

    # Saturated properties
    s = _sat_props(AS, p_sat)
    rho_l, rho_v, mu_l = s["rho_l"], s["rho_v"], s["mu_l"]

    # Equivalent Reynolds number
    G_eq = G * ((1 - x) + x * np.sqrt(rho_l / rho_v))
    Re_eq = G_eq * Dh / mu_l

    # Friction factor
    Ge3 = 3521.1 * (p_co / Dh) ** 4.17 * theta ** (-7.75)
    Ge4 = -1.024 * (p_co / Dh) ** 0.0925 * theta ** (-1.3)
    f = Ge3 * Re_eq ** Ge4

    # Frictional pressure drop over the full flow length
    rho_m = 1.0 / (x / rho_v + (1.0 - x) / rho_l)
    return f * (L / Dh) * G**2 / rho_m



#%% TWO-PHASE : PLATE CORRELATIONS OF THE NCTU GROUP (Kuo 2005, Yan-Lio-Lin 1999, Hsieh-Lin 2002)
# Added 2026-11-26 (drafts written by gemini-flash-lite-latest, reviewed; equations checked
# against the papers, see LOGS reports/plate_papers_check_2026-11-26.md)


def BPHX_kuo_condensing_DP(AS, geom, cond):
    """
    Frictional pressure drop of a condensing flow in a vertical plate heat exchanger.
    
    Parameters:
    - AS: CoolProp abstract state object [-]
    - geom: dictionary containing geometry parameters
        - 'Dh': hydraulic diameter [m]
        - 'L': total flow path length [m]
    - cond: dictionary containing operating conditions
        - 'G': mass flux in ONE channel [kg/(m^2 s)]
        - 'x': mean quality of the cell [-]
        - 'q': heat flux [W/m^2]
        - 'p': saturation pressure [Pa]

    Returns:
    - dP: frictional pressure drop over the full length L [Pa]

    Reference:
    Kuo, W.S., Lie, Y.M., Hsieh, Y.Y., Lin, T.F. (2005) Condensation heat transfer 
    and pressure drop of refrigerant R-410A flow in a vertical plate heat exchanger. 
    Int. J. Heat Mass Transfer 48, 5205-5220. doi:10.1016/j.ijheatmasstransfer.2005.07.023.
    Equations (12), (15), (17), (21), (23), and (24).

    Validity range:
    - Mass flux G: 50 - 150 kg/(m^2 s)
    - Heat flux q: 5000 - 20000 W/m^2
    - Quality x: 0.10 - 0.80
    - Fluid: R-410A (FluidWarning otherwise)
    - State only (not checked): p_sat 1.44 - 1.95 MPa, chevron angle 60 deg
    - Mean deviation: 23 % (f_tp)

    Definitions:
    - Fanning friction factor (f_tp): defined via Eq. (15),
      dP_f = 2 * f_tp * G^2 * v_m * L / Dh.
    - Reynolds number (Re_eq): based on equivalent mass flux G_eq and hydraulic diameter Dh,
      using liquid dynamic viscosity mu_l.
    - Mass flux G: defined for one channel.
    - Quality x: cell mean quality.
    - Heat flux q: absolute heat flux.
    - Properties: saturated properties evaluated at saturation pressure p_sat.

    Equations checked against the paper (LaboThapPy, 2026-11-26).
    """
    # 1) Inputs extraction and quality clamp (Common rules)
    G = cond['G']
    x_raw = cond['x']
    if not np.isfinite(x_raw):
        raise ValueError("Kuo condensing DP: cond['x'] must be finite.")
    x = min(max(x_raw, 1e-4), 1.0 - 1e-4)
    
    q = abs(cond['q'])
    p_sat = cond['p']
    
    Dh = geom['Dh']
    L = geom['L']

    # 2) Saturated properties
    s = _sat_props(AS, p_sat)
    rho_l = s["rho_l"]
    rho_v = s["rho_v"]
    mu_l = s["mu_l"]
    h_lv = s["h_lv"]

    # 3) Equivalent mass flux and Reynolds number (Eq. 24 and Eq. 23)
    term_sqrt = np.sqrt(rho_l / rho_v)
    G_eq = G * ((1.0 - x) + x * term_sqrt)
    Re_eq = G_eq * Dh / mu_l

    # 4) Boiling number (Eq. 21) with guard
    Bo = q / (G * h_lv)
    Bo = max(Bo, 1e-6)

    # 5) Fanning friction factor (Eq. 17)
    f_tp = 21500.0 * (Re_eq ** (-1.14)) * (Bo ** (-0.085))

    # 6) Homogeneous specific volume (Eq. 12)
    v_m = (x / rho_v) + ((1.0 - x) / rho_l)

    # 7) Frictional pressure drop over the full length L (Eq. 15)
    dP = 2.0 * f_tp * (G ** 2) * v_m * (L / Dh)

    # Validity range checks and warnings
    out = []
    if G < 50.0 or G > 150.0:
        out.append(f"G = {G:.4g} kg/m2/s (50-150)")
    if q < 5000.0 or q > 20000.0:
        out.append(f"q = {q:.4g} W/m2 (5000-20000)")
    if x_raw < 0.10 or x_raw > 0.80:
        out.append(f"x = {x_raw:.4g} (0.10-0.80)")

    if out:
        warnings.warn(
            "BPHX_kuo_condensing_DP: outside the data of the paper: " + ", ".join(out) + ".",
            ValidityRangeWarning
        )

    fluid_name = AS.fluid_names()[0].upper().replace('-', '')
    if fluid_name != 'R410A':
        warnings.warn(
            f"BPHX_kuo_condensing_DP: fluid {AS.fluid_names()[0]}, "
            "correlation established for R410A.",
            FluidWarning
        )

    return dP


def BPHX_yan_lio_lin_condensing_DP(AS, geom, cond):
    """
    Yan, Y.-Y., Lio, H.-C., Lin, T.-F. (1999) Condensation heat transfer and pressure drop of
    refrigerant R-134a in a plate heat exchanger. Int. J. Heat Mass Transfer 42, 993-1006.
    doi:10.1016/S0017-9310(98)00217-8

    Frictional pressure drop of a condensing flow in a plate heat exchanger.
    Port losses, gravity, and acceleration terms are NOT included.

    Parameters
    ----------
    AS : CoolProp abstract state or equivalent
         CoolProp fluid state object.
    geom : dict
         Geometry dictionary containing:
         - 'Dh' : Hydraulic diameter [m] (Dh = 2 * b, where b is the channel gap).
         - 'L' : Total flow path length [m].
    cond : dict
         Condition dictionary containing:
         - 'G' : Mass flux in ONE channel [kg/m2/s].
         - 'x' : Mean quality of the cell [-].
         - 'q' : Heat flux [W/m2] (absolute value).
         - 'p' : Saturation pressure [Pa].

    Returns
    -------
    dP : float
         Frictional pressure drop over the full flow length L [Pa].

    Validity range
    --------------
    - Mass flux G : 60 to 120 kg/m2/s
    - Heat flux q : 10000 to 18000 W/m2
    - Fluid : R134a (FluidWarning otherwise)
    - State only (not checked): saturation pressure 0.7 to 0.9 MPa, Re_eq about 2500 to
      15000 (read on Fig. 11), chevron angle 60 deg

    Definitions
    -----------
    - Fanning friction factor f_tp = dP_f * Dh / (2 * G^2 * v_m * L)
    - Reynolds numbers based on hydraulic diameter Dh and liquid dynamic viscosity mu_l
    - G is the mass flux in one channel
    - Quality is the cell mean quality
    - Properties are saturated at the saturation pressure p
    - p_crit is read from the AbstractState independently of saturation property calls

    Equations checked against the paper (LaboThapPy, 2026-11-26).
    """
    # 1) Inputs extraction and quality check
    G = cond['G']
    x_raw = cond['x']
    q = abs(cond['q'])
    p_sat = cond['p']
    Dh = geom['Dh']
    L = geom['L']

    # Quality check: verify it is finite, then clamp
    if not np.isfinite(x_raw):
        raise ValueError("Yan Lio Lin condensing DP: cond['x'] must be finite.")
    x = min(max(x_raw, 1e-4), 1.0 - 1e-4)

    # 2) Critical pressure (independent of the state of AS)
    p_crit = AS.p_critical()

    # Saturated properties
    s = _sat_props(AS, p_sat)
    rho_l = s["rho_l"]
    rho_v = s["rho_v"]
    mu_l = s["mu_l"]
    h_lv = s["h_lv"]

    # 3) Equivalent mass flux and Reynolds number
    term_sqrt = np.sqrt(rho_l / rho_v)
    G_eq = G * ((1.0 - x) + (x * term_sqrt))
    Re_eq = G_eq * Dh / mu_l

    # 4) All-liquid Reynolds number
    Re_l = G * Dh / mu_l

    # 5) Boiling number with guard
    Bo = q / (G * h_lv)
    Bo = max(Bo, 1e-6)

    # 6) Fanning friction factor (paper Fig. 11)
    ratio_p = p_sat / p_crit
    f_tp = 94.75 * (Re_eq ** (-0.0467)) * (Re_l ** (-0.4)) * (Bo ** 0.5) * (ratio_p ** 0.8)

    # 7) Homogeneous specific volume
    v_m = (x / rho_v) + ((1.0 - x) / rho_l)

    # 8) Frictional pressure drop over the full length L [Pa]
    dP = 2.0 * f_tp * (G ** 2) * v_m * L / Dh

    # Validity range checks and warnings
    out = []
    if G < 60.0 or G > 120.0:
        out.append(f"G = {G:.4g} kg/m2/s (60-120)")
    if q < 10000.0 or q > 18000.0:
        out.append(f"q = {q:.4g} W/m2 (10000-18000)")

    if out:
        warnings.warn(
            "BPHX_yan_lio_lin_condensing_DP: outside the data of the paper: "
            + ", ".join(out) + ".",
            ValidityRangeWarning
        )

    # Fluid check
    fluid_name = AS.fluid_names()[0].upper().replace('-', '')
    if fluid_name != 'R134A':
        warnings.warn(
            f"BPHX_yan_lio_lin_condensing_DP: fluid {AS.fluid_names()[0]}, "
            "correlation established for R134a.",
            FluidWarning
        )

    return dP


def BPHX_hsieh_lin_boiling_DP(AS, geom, cond):
    """
    Frictional pressure drop of a saturated boiling flow in a plate heat exchanger.
    
    Parameters
    ----------
    AS : CoolProp abstract state or fluid object
        Working fluid state.
    geom : dict
        Plate geometry dictionary with keys:
        - 'Dh' : hydraulic diameter [m] (Dh = 2 * b, where b is channel gap)
        - 'L'  : total flow path length (plate length x number of passes) [m]
    cond : dict
        Operating conditions dictionary with keys:
        - 'G'   : mass flux in ONE channel [kg/m^2/s]
        - 'x'   : mean quality of the cell [-]
        - 'q'   : heat flux [W/m^2] (the sign depends on the side)
        - 'p'   : saturation pressure [Pa]

    Returns
    -------
    float
        Frictional pressure drop over the full length L [Pa].

    Reference
    ---------
    Hsieh, Y.Y., Lin, T.F. (2002) Saturated flow boiling heat transfer and pressure drop
    of refrigerant R-410A in a vertical plate heat exchanger.
    Int. J. Heat Mass Transfer 45, 1033-1044.
    doi:10.1016/S0017-9310(01)00219-8
    Equations (22), (26), (27), and (28).

    Validity
    --------
    - Mass flux G : 50 to 125 kg/m^2/s
    - Heat flux q : 5000 to 35000 W/m^2
    - Fluid : R-410A (FluidWarning otherwise)
    - State only (not checked): saturation pressure 1.08 to 1.44 MPa, chevron angle 60 deg
    - Accuracy : 85 % of the data within 20 %.

    Definitions
    -----------
    - Fanning friction factor f_tp defined by Eq. (22): dP = 2 * f_tp * G**2 * v_m * L / Dh
    - Reynolds number Re_eq based on equivalent mass flux G_eq and liquid viscosity mu_l
    - G is the mass flux in one channel [kg/m^2/s]
    - Quality x is the cell mean quality [-]
    - Heat flux q is the absolute heat flux [W/m^2]
    - Properties are saturated at pressure p

    Equations checked against the paper (LaboThapPy, 2026-11-26).
    """
    # 1) Inputs extraction and quality check / clamp
    G = cond['G']
    x_input = cond['x']
    if not np.isfinite(x_input):
        raise ValueError("Hsieh-Lin boiling DP: cond['x'] must be finite "
                         "(two-phase cell expected).")
    
    x = min(max(x_input, 1e-4), 1.0 - 1e-4)
    q = abs(cond['q'])
    p_sat = cond['p']

    Dh = geom['Dh']
    L = geom['L']

    # 2) Saturated properties at p_sat
    s = _sat_props(AS, p_sat)
    rho_l = s["rho_l"]
    rho_v = s["rho_v"]
    mu_l = s["mu_l"]

    # 3) Equivalent mass flux Eq. (28) and Reynolds number Eq. (27)
    ratio_density_sqrt = np.sqrt(rho_l / rho_v)
    term_quality_effect = (1.0 - x) + (x * ratio_density_sqrt)
    G_eq = G * term_quality_effect
    Re_eq = G_eq * Dh / mu_l

    # 4) Fanning friction factor Eq. (26)
    f_tp = 61000.0 * (Re_eq ** (-1.25))

    # 5) Homogeneous specific volume (definition of v_m used in Eq. 22)
    v_m = (x / rho_v) + ((1.0 - x) / rho_l)

    # 6) Frictional pressure drop Eq. (22) over the full length L
    dP = 2.0 * f_tp * (G ** 2) * v_m * (L / Dh)

    # Validity range checks and warnings
    out = []
    if G < 50.0 or G > 125.0:
        out.append(f"G = {G:.4g} kg/m2/s (50-125)")
    if q < 5000.0 or q > 35000.0:
        out.append(f"q = {q:.4g} W/m2 (5000-35000)")

    if out:
        warnings.warn(
            "BPHX_hsieh_lin_boiling_DP: outside the data of the paper: " + ", ".join(out) + ".",
            ValidityRangeWarning
        )

    # Fluid check
    fluid_name = AS.fluid_names()[0].upper().replace('-', '')
    if fluid_name != 'R410A':
        warnings.warn(
            f"BPHX_hsieh_lin_boiling_DP: fluid {AS.fluid_names()[0]}, "
            "correlation established for R410A.",
            FluidWarning
        )

    return dP
