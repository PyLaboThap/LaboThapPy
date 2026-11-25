# -*- coding: utf-8 -*-
"""
Frictional pressure drop in curved (helically / spirally coiled) channels, (AS, geom, cond) formalism.
Prefix HCHE_..._DP. Every function returns the pressure drop over the FULL flow length geom['L'] [Pa]
(friction + gravity); no entrance / exit / bend losses, no acceleration.

geom : Dh (or D, D_t_i)  hydraulic diameter [m];  R_curv (or D_coil, P_coil) radius of curvature [m];
       L developed length [m];  K roughness [m] (default 0);  gamma [deg] (inclination, default 0)
cond : G [kg/m2/s] mass flux of ONE channel, mu [Pa s], rho or D [kg/m3]

Curved-channel friction (Darcy), curv = Dh / (2 R), Re = G Dh / mu, De = Re sqrt(curv):
    Re <  Re_cr = 20000 curv^0.32 (Ito)   f = f_s,lam * max(1, ratio_lam)
        ratio_lam = 0.1033 sqrt(De) / [sqrt(1 + 1.729/De) - sqrt(1.729/De)]^3   for De >= 11.6, else 1
        (Ito laminar, the same factor as in the HTC bank; f_s,lam = 64/Re, or fRe of the duct)
    Re >= Re_cr                           f = f_s,Churchill * max(1, [Re curv^2]^0.05)
        (Ito turbulent, established for 0.034 < Re curv^2 < 300)
The two branches are continuous at Re_cr (checked for curv ~ 0.02-0.05). These relations are written from the
literature, not from your previous spiral_ducts_DP.py: compare them with it before relying on the numbers.
"""
import warnings
import numpy as np

from labothappy.correlations.pressure_drop.pipe_DP_correlation_bank import friction_factor_churchill
from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import _D, _R_curv, _L, _need_cond

G_GRAVITY = 9.81


def curved_friction_factor(Re, curv, K, D, f_lam_Re=64.0):
    """Darcy friction factor of a curved channel (see module docstring). f_lam_Re = f Re of the straight duct
    in laminar flow (64 circular, 56.9 square, ...)."""
    De = Re * np.sqrt(curv)
    Re_cr = 20000.0 * curv ** 0.32
    if Re < Re_cr:
        ratio = 1.0
        if De >= 11.6:
            x = 1.729 / De
            ratio = max(1.0, 0.1033 * np.sqrt(De) / (np.sqrt(1.0 + x) - np.sqrt(x)) ** 3)
        return f_lam_Re / Re * ratio
    return friction_factor_churchill(K, D, Re) * max(1.0, (Re * curv ** 2) ** 0.05)


def _curved_sp_dp(geom, cond, f_lam_Re=64.0):
    _need_cond(cond, 'G', 'mu')
    D, L, R = _D(geom), _L(geom), _R_curv(geom)
    rho = cond['rho'] if 'rho' in cond else cond['D']
    G, Re = cond['G'], cond['G'] * D / cond['mu']
    f = curved_friction_factor(Re, D / (2.0 * R), geom.get('K', 0.0), D, f_lam_Re)
    return f * (L / D) * G ** 2 / (2.0 * rho) + G_GRAVITY * rho * L * np.sin(np.radians(geom.get('gamma', 0.0)))


def curved_pipe_ito_SP_DP(AS, geom, cond):
    """Single phase in a curved channel (Ito laminar / turbulent, see module docstring)."""
    return _curved_sp_dp(geom, cond)


#%% ---------------------------------------------------------------------------------------------------------
# Basile's banks (spiral_ducts_DP.py), (AS, geom, cond) version.
# ---------------------------------------------------------------------------------------------------------
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props, _clamp_x


def _dp_fanning(f, G, rho, L, D):
    """DP = 4 f (L/D) rho v^2 / 2   (f Fanning)."""
    return 4.0 * f * (L / D) * G ** 2 / (2.0 * rho)


def curved_pipe_sigalotti_SP_DP(AS, geom, cond):
    """Sigalotti, Alvarado-Rodriguez & Rendon (2023), Fluids 8(12) 308 (review), after White / Ito.
    FANNING f; gamma = curv = D/(2R); De = Re sqrt(curv).
    Re < Re_cr = 12730 curv^0.2 (Eq. 12, 0.001165 < curv < 0.0667): laminar, White (1929) Eq. (16)
        (11.6 < De < 2000, curv <= 0.066), 16/Re for De <= 11.6 (straight-pipe limit);
    Re >= Re_cr: turbulent, White (1932) Eq. (17) f = 0.08 Re^-1/4 + 0.012 sqrt(curv), valid
        1.5e4 < Re < 1e5 (the paper prints "FRe =": typo, it is F; F Re would give f ~ 1e-6).
    Warnings outside these ranges.
    (2026-11-21: Re_cr <= Re < 1.5e5 used Adler's LAMINAR law 16 * 0.1064 sqrt(De) / Re (f 20-70 % low,
    x3 jump at 1.5e5) and White Eq. 17 only above 1.5e5: corrected, Eq. 17 used above Re_cr.)"""
    _need_cond(cond, 'G', 'mu')
    D, L, R = _D(geom), _L(geom), _R_curv(geom)
    rho = cond['rho'] if 'rho' in cond else cond['D']
    G = cond['G']
    Re = G * D / cond['mu']
    curv = D / (2.0 * R)
    De = Re * np.sqrt(curv)
    if not 0.001165 < curv < 0.0667:
        warnings.warn(f"Sigalotti DP: curv = {curv:.4g} outside 0.001165-0.0667 (Re_cr, Eq. 12).")
    if Re < 12730.0 * curv ** 0.2:
        # Laminar: White (1929), Sigalotti Eq. (16)
        if curv > 0.066 or De > 2000:
            warnings.warn(f"Sigalotti DP laminar: De = {De:.0f}, curv = {curv:.4g} outside "
                          "De < 2000, curv <= 0.066 (White, Eq. 16).")
        if curv <= 0.066 and De > 11.6:
            f = 16.0 / (1.0 - (1.0 - (11.6 / De) ** 0.45) ** (1.0 / 0.45)) / Re
        else:
            f = 16.0 / Re
    else:
        # Turbulent: White (1932), Sigalotti Eq. (17), FANNING, valid 1.5e4 < Re < 1e5
        if Re < 1.5e4 or Re > 1.0e5:
            warnings.warn(f"Sigalotti DP turbulent: Re = {Re:.3g} outside 1.5e4-1e5 (White, Eq. 17).")
        f = 0.08 * Re ** -0.25 + 0.012 * np.sqrt(curv)
    return _dp_fanning(f, G, rho, L, D) + G_GRAVITY * rho * L * np.sin(np.radians(geom.get('gamma', 0.0)))


def _ito_turb_lo(G, D, mu, rho, L, curv):
    """Liquid-alone DP, Ito turbulent (Darcy): f = 0.304 Re^-0.25 + 0.029 sqrt(curv)."""
    Re = G * D / mu
    return (0.304 * Re ** -0.25 + 0.029 * np.sqrt(curv)) * (L / D) * G ** 2 / (2.0 * rho)


def curved_pipe_cioncolini_boiling_DP(AS, geom, cond, C=None):
    """Cioncolini, Santini & Ricotti (2008), Exp. Therm. Fluid Sci. 32, 1301-1312, Eq. (20) (also in
    Fsadni & Whitty 2016): Lockhart-Martinelli phi_l^2 = (1 + C/Xtt + 1/Xtt^2) with the Tarasova heating
    correction (1 + 0.0044 (q/G)^0.7) (dimensional: q in W/m2, G in kg/m2/s), liquid alone with Ito.
    This is an existing straight-pipe correlation that the authors found to fit their coil data (not
    their own fit); Tarasova was developed for 4.9-19.6 MPa, G 501-2577, q 101-1671 kW/m2.
    C = 20 if Re_lo > 2000 and Re_go > 2000, C = 10 if Re_lo > 2000 and Re_go < 2000 (paper);
    Re_lo = G d / mu_l, Re_go = G d / mu_v (total flow as liquid / as gas). C may be forced by argument.
    cond: G, p, x, q [W/m2] (heating; default 0 = no heating term; q < 0 -> 0).
    Paper range: water, d 4.03 / 4.98 mm, D/d 26.1-93.3, p 120-660 kPa, G 290-690, q 50-440 kW/m2.
    (2026-11-21: C = 10 case added.)"""
    _need_cond(cond, 'G', 'p', 'x')
    D, L, R = _D(geom), _L(geom), _R_curv(geom)
    G, x = cond['G'], _clamp_x(cond['x'], 1e-6)
    s = _sat_props(AS, cond['p'], cond)
    if C is None:
        # Chisholm constant of the paper (flow regime of liquid-only and gas-only flows)
        Re_lo = G * D / s['mu_l']
        Re_go = G * D / s['mu_v']
        if Re_lo <= 2000:
            warnings.warn(f"Cioncolini 2008: Re_lo = {Re_lo:.0f} <= 2000, C not given by the paper "
                          "(C = 20 used).")
        if Re_go < 2000:
            C = 10.0
        else:
            C = 20.0
    inv_Xtt = (x / (1 - x)) ** 0.9 * (s['rho_l'] / s['rho_v']) ** 0.5 * (s['mu_v'] / s['mu_l']) ** 0.1
    q = max(cond.get('q', 0.0), 0.0)
    phi_l2 = (1 + C * inv_Xtt + inv_Xtt ** 2) * (1 + 0.0044 * (q / G) ** 0.7)
    return phi_l2 * _ito_turb_lo((1 - x) * G, D, s['mu_l'], s['rho_l'], L, D / (2 * R))


def curved_pipe_ruffell_boiling_DP(AS, geom, cond):
    """Ruffell (1974), after Fsadni & Whitty (2016): phi_lo^2 = (v_m / v_l)(1 + F), liquid only with Ito.
    Steam-water, d = 10.7-18.6 mm, 6-18 MPa, G = 300-1800: far outside for most refrigerant cases."""
    _need_cond(cond, 'G', 'p', 'x')
    D, L, R = _D(geom), _L(geom), _R_curv(geom)
    G, x = cond['G'], min(max(cond['x'], 0.0), 1.0)
    s = _sat_props(AS, cond['p'], cond)
    v_m = x / s['rho_v'] + (1 - x) / s['rho_l']
    y, g = 2.0 * R / (100.0 * D), G / 1000.0
    F = (np.sin(1.16 * g)
         * (0.875 - 0.314 * y - 0.74 * g * (0.152 - 0.07 * y) - x * (0.155 * g + 0.7 - 0.19 * y))
         * (1 - 12 * (x - 0.3) * (x - 0.4) * (x - 0.5) * (x - 0.6)))
    phi_lo2 = v_m * s['rho_l'] * (1 + F)
    return phi_lo2 * _ito_turb_lo(G, D, s['mu_l'], s['rho_l'], L, D / (2 * R))


def curved_pipe_fsadni_condensing_DP(AS, geom, cond):
    """Your curved_pipe_condensing_DP: Cioncolini if curv <= 0.04 and x < 0.9, else Ruffell.
    The hard switch at x = 0.9 is replaced by a linear blend over 0.85 < x < 0.95: the jump (~25 %) made the
    pressure-drop fixed-point loop of the model cycle (not converged after 30 iterations)."""
    if D_over_2R(geom) > 0.04:
        return curved_pipe_ruffell_boiling_DP(AS, geom, cond)
    w = min(max((cond['x'] - 0.85) / 0.10, 0.0), 1.0)
    if w >= 1.0:
        return curved_pipe_ruffell_boiling_DP(AS, geom, cond)
    dp = curved_pipe_cioncolini_boiling_DP(AS, geom, cond)
    return dp if w <= 0.0 else (1.0 - w) * dp + w * curved_pipe_ruffell_boiling_DP(AS, geom, cond)


def D_over_2R(geom):
    return _D(geom) / (2.0 * _R_curv(geom))

