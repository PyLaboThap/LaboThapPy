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

#%% SINGLE PHASE : Friction DP

def sinnott_plate_DP(AS, geom, cond):
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

def kumar_plate_DP(AS, geom, cond):
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

def thonon_plate_DP(AS, geom, cond):
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

def muley_manglik_plate_DP(AS, geom, cond):
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

def martin_plate_DP(AS, geom, cond):
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
    Dh = geom['Dh']
    L = (geom['L']-2*D_port)*n_passes
    
    beta = np.radians(chevron_deg(geom['beta']))

    Re = G * Dh / mu
    zeta = _zeta_martin(Re, beta)
    dp = zeta * (L / Dh) * G**2 / (2 * rho)

    return dp

#%% TWO-PHASE : GENERIC


#%% TWO-PHASE : BOILING 


#%% TWO-PHASE : CONDENSING 
