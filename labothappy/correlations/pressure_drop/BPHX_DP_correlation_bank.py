# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:51:41 2026

@author: basil
"""

#%% IMPORTS
import numpy as np
import warnings

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
        warnings.warn(f"Re = {Re:.1f} outside the Thonon range [50, 15000]: extrapolated.")

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
                      "hors du domaine (Re >= 1000, 30-60 deg): extrapolé.")

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
    Dh = geom['Dh']*geom['phi']
    L = (geom['L']-2*D_port)*n_passes
    
    beta = np.radians(chevron_deg(geom['beta']))

    Re = G * Dh / mu
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
                     lm_plate_friction (default True: Martin's PHE friction factor for each phase,
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
    alpha = np.radians(chevron_deg(geom['beta'])) if geom.get('lm_plate_friction', True) else None

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

#%% TWO-PHASE : CONDENSING 

from labothappy.correlations.convection.BPHE_HTC_correlation_bank import _sat_props

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
    geom : Dh [m], p_co (corrugation pitch) [m],
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

    Dh, L, p_co = geom['Dh'], geom['L'], geom['p_co']
    theta = np.radians(chevron_deg(geom['beta']))            # from the flow direction = Han's pi/2 - beta
    if not 0.0 < theta <= np.pi / 2:
        raise ValueError(f"Han condensation DP: angle {theta} rad from the flow direction "
                         "outside ]0, pi/2].")
    if not 13 <= G <= 34 or not 45 <= np.degrees(theta) <= 70:
        warnings.warn("Han condensation DP: outside the range of the Han data (G 13-34 kg/m2/s, "
                      "45-70 deg from the flow direction): extrapolated.")

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

