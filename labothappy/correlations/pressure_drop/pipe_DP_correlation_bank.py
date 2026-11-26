# -*- coding: utf-8 -*-
"""
In-pipe (tube-side) frictional pressure drop, in the (AS, geom, cond) formalism.
Correlations and friction factors from the labothappy pipe pressure drop library (E. Neven), written here
without the AbstractState-based interface: the fluid state comes in `cond`.

Each function returns the pressure drop over the FULL flow length geom['L'] [Pa]: friction + gravity
(geom['gamma'], default 0, positive = upward flow). No acceleration (it depends on the quality change of the
cell and cannot be scaled with the length fraction of the cell), no return bends or head losses.

cond : G (mass flux in ONE tube) [kg/m2/s] (or m_dot + geom n_t, n_passes), mu [Pa s], rho or D [kg/m3]
       (single phase); p [Pa], x [-] (two-phase; optional cond['mu_l'], cond['mu_v'], cond['sigma'] overrides
       for fluids without a CoolProp transport / surface tension model)
geom : D_t_i (or D_t_o and t_t), L (total tube flow path: straight length x passes) [m],
       K (absolute roughness) [m], default 0; gamma [deg] (inclination), default 0
"""

import math
import warnings
import numpy as np

from labothappy.correlations.convection.pipe_HTC_correlation_bank import (
    _tube_ID, _tube_G, _sat_props, _need_sigma)
from labothappy.correlations.convection.correlation_warnings import ValidityRangeWarning

EPS = 1e-12
G_GRAVITY = 9.81  # m/s^2

#%% FRICTION FACTORS (Darcy)

def friction_factor_churchill(K, d_hyd, Re):
    """Churchill (1977): one formula for laminar, transition and turbulent flow, with roughness."""
    Re = max(EPS, Re)
    laminar_term = (8.0 / Re)**12
    k_rel = K / max(EPS, d_hyd)
    ln_arg = 1.0 / ((7.0 / Re)**0.9 + 0.27 * k_rel)
    A = (2.457 * math.log(max(EPS, ln_arg)))**16
    B = (37530.0 / Re)**16
    return 8.0 * (laminar_term + 1.0 / (A + B)**1.5)**(1.0 / 12.0)


def friction_factor_swamee_jain(K, d_hyd, Re):
    """Swamee & Jain (1976): explicit Colebrook-White, turbulent flow (Re >= 5000).
    Validity: 5e3 <= Re <= 1e8, 1e-6 <= K/d <= 5e-2 (fluids.friction.Swamee_Jain_1976 notes);
    warning outside (the smooth limit K/d < 1e-6 is not flagged: the formula tends to the smooth
    pipe value)."""
    if Re < 5e3 or Re > 1e8:
        warnings.warn(f"Swamee-Jain friction factor: Re = {Re:.3g} outside 5e3-1e8: extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    if K / max(EPS, d_hyd) > 5e-2:
        warnings.warn(f"Swamee-Jain friction factor: K/d = {K / d_hyd:.3g} above 5e-2: "
                      "extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    Re = max(EPS, Re)
    log_arg = K / max(EPS, d_hyd) / 3.7 + 5.74 / Re**0.9
    return 0.25 / (math.log10(max(EPS, log_arg))**2)


def friction_factor_haaland(K, d_hyd, Re):
    """Haaland (1983): explicit Colebrook-White, 4e3 <= Re <= 1e8, 1e-6 <= K/d <= 0.05
    (also fluids.friction.Haaland notes); warning outside (the smooth limit K/d < 1e-6 is not
    flagged: the formula tends to the smooth pipe value)."""
    if Re < 4e3 or Re > 1e8:
        warnings.warn(f"Haaland friction factor: Re = {Re:.3g} outside 4e3-1e8: extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    if K / d_hyd > 0.05:
        warnings.warn(f"Haaland friction factor: K/d = {K / d_hyd:.3g} above 0.05: extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    return (1.0 / (-1.8 * np.log10((K / (3.75 * d_hyd))**1.11 + 6.9 / Re)))**2


def friction_factor_konakov(Re):
    """Konakov: f = (1.8 log10 Re - 1.5)^-2, smooth pipes, 4e3 <= Re <= 3e6 (warning outside)."""
    if Re < 4e3 or Re > 3e6:
        warnings.warn(f"Konakov friction factor: Re = {Re:.3g} outside 4e3-3e6: extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    return (1.8 * math.log10(Re) - 1.5)**(-2)


def friction_factor_petukhov(Re):
    """Petukhov (1970): f = (0.79 ln Re - 1.64)^-2, smooth pipes, 3e3 <= Re <= 5e6 (warning
    outside)."""
    if Re < 3e3 or Re > 5e6:
        warnings.warn(f"Petukhov friction factor: Re = {Re:.3g} outside 3e3-5e6: extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    return (0.79 * np.log(Re) - 1.64)**(-2)


def friction_factor_cheng_CO2(G, d_hyd, P, h, mu):
    """Cheng et al. (2024): supercritical CO2 in a horizontal tube (liquid-like / vapor-like / two-phase-like
    by the pseudo-vapor quality). D = 10 mm, G 496.7-1346.2 kg/m2/s, P 7.53-23.51 MPa, CO2 only.
    Needs labothappy.correlations.properties.supercritical (pseudo_boiling_bounds, pseudo_critical_temperature_CO2).
    Warning outside G 496.7-1346.2 kg/m2/s, P 7.53-23.51 MPa (ranges above), and for a diameter
    outside 8-12 mm (single tested diameter 10 mm; 8-12 mm = window of the labothappy automatic
    choice, select_pipe_DP)."""
    if G < 496.7 or G > 1346.2:
        warnings.warn(f"Cheng CO2 friction factor: G = {G:.4g} kg/m2/s outside 496.7-1346.2: "
                      "extrapolated.",
                      ValidityRangeWarning, stacklevel=2)
    if P < 7.53e6 or P > 23.51e6:
        warnings.warn(f"Cheng CO2 friction factor: P = {P / 1e6:.3g} MPa outside 7.53-23.51 MPa: "
                      "extrapolated.", ValidityRangeWarning, stacklevel=2)
    if d_hyd < 8e-3 or d_hyd > 12e-3:
        warnings.warn(f"Cheng CO2 friction factor: D = {d_hyd * 1e3:.3g} mm, tested diameter "
                      "10 mm: "
                      "extrapolated.", ValidityRangeWarning, stacklevel=2)
    import CoolProp.CoolProp as CP
    from labothappy.correlations.properties.supercritical import (
        pseudo_boiling_bounds, pseudo_critical_temperature_CO2)
    AS = CP.AbstractState('HEOS', 'CO2')
    AS.update(CP.HmassP_INPUTS, h, P)
    rho = AS.rhomass()
    T_pc, T_minus, T_plus = pseudo_boiling_bounds(P, 'CO2', T_pc_func=pseudo_critical_temperature_CO2)
    AS.update(CP.PT_INPUTS, P, T_minus)
    h_LL = AS.hmass()
    AS.update(CP.PT_INPUTS, P, T_plus)
    h_LV = AS.hmass()
    x_pb = (h - h_LL) / (h_LV - h_LL)
    Re = G * d_hyd / mu
    if x_pb < 0 or x_pb > 1:                          # vapor-like / liquid-like
        return 1.5 * (1.82 * np.log10(Re) - 1.64)**(-2)
    Fr_VL = G**2 * x_pb / (rho**2 * G_GRAVITY * d_hyd)  # two-phase-like
    return 4.379 * Re**(-0.388) * Fr_VL**(-0.0167)

#%% HELPERS

def _gravity(geom, rho):
    return G_GRAVITY * rho * geom['L'] * math.sin(math.radians(geom.get('gamma', 0.0)))


def _friction_dp(f, G, rho, D, L):
    return f * (L / D) * G**2 / (2.0 * rho)


def _sp_dp(AS, geom, cond, friction):
    """friction(K, D, Re, G, cond) -> Darcy f."""
    D, L, K = _tube_ID(geom), geom['L'], geom.get('K', 0.0)
    G = _tube_G(geom, cond)
    rho = cond['rho'] if 'rho' in cond else cond['D']
    Re = G * D / cond['mu']
    f = friction(K, D, Re, G, cond)
    return _friction_dp(f, G, rho, D, L) + _gravity(geom, rho)

#%% SINGLE PHASE

def pipe_churchill_SP_DP(AS, geom, cond):
    """Darcy-Weisbach, Churchill friction factor (all Re, roughness geom['K'])."""
    return _sp_dp(AS, geom, cond, lambda K, D, Re, G, c: friction_factor_churchill(K, D, Re))


def pipe_swamee_jain_SP_DP(AS, geom, cond):
    """Darcy-Weisbach, Swamee-Jain friction factor (turbulent, roughness)."""
    return _sp_dp(AS, geom, cond, lambda K, D, Re, G, c: friction_factor_swamee_jain(K, D, Re))


def pipe_haaland_SP_DP(AS, geom, cond):
    """Darcy-Weisbach, Haaland friction factor (turbulent, roughness)."""
    return _sp_dp(AS, geom, cond, lambda K, D, Re, G, c: friction_factor_haaland(K, D, Re))


def _smooth_only(name, geom):
    if geom.get('K', 0.0) > 0:
        warnings.warn(f"K={geom['K']} was given but the {name!r} correlation is smooth-pipe-only: roughness ignored.",
                      ValidityRangeWarning,
                      stacklevel=3)


def pipe_konakov_SP_DP(AS, geom, cond):
    """Darcy-Weisbach, Konakov friction factor (turbulent, smooth pipes)."""
    _smooth_only('Konakov', geom)
    return _sp_dp(AS, geom, cond, lambda K, D, Re, G, c: friction_factor_konakov(Re))


def pipe_petukhov_SP_DP(AS, geom, cond):
    """Darcy-Weisbach, Petukhov friction factor (turbulent, smooth pipes)."""
    _smooth_only('Petukhov', geom)
    return _sp_dp(AS, geom, cond, lambda K, D, Re, G, c: friction_factor_petukhov(Re))


def pipe_cheng_CO2_SP_DP(AS, geom, cond):
    """Darcy-Weisbach, Cheng et al. friction factor for supercritical CO2. cond : G, p, h, mu, rho."""
    if AS.fluid_names()[0] not in ('CarbonDioxide', 'CO2'):
        raise ValueError(f"Cheng CO2 friction factor is specific to CO2, got {AS.fluid_names()[0]!r}.")
    return _sp_dp(AS, geom, cond,
                  lambda K, D, Re, G, c: friction_factor_cheng_CO2(G, D, c['p'], c['h'], c['mu']))

#%% TWO-PHASE (EVAPORATION / CONDENSATION)

def _two_phase_dp(AS, geom, cond, model):
    """model(G, x, s, D, L, K, cond) -> frictional dP; gravity added with the homogeneous density."""
    D, L, K = _tube_ID(geom), geom['L'], geom.get('K', 0.0)
    G, x = _tube_G(geom, cond), min(max(cond['x'], 1e-6), 1 - 1e-6)
    s = _sat_props(AS, cond['p'], cond)
    dP = model(G, x, s, D, L, K, cond)
    rho_h = 1.0 / (x / s['rho_v'] + (1 - x) / s['rho_l'])
    return dP + _gravity(geom, rho_h)


def _friedel(G, x, s, D, L, K, cond):
    rho_l, rho_v, mu_l, mu_v, sigma = s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v'], s['sigma']
    rho_h = 1.0 / (x / rho_v + (1 - x) / rho_l)
    f_l = friction_factor_churchill(K, D, G * D / mu_l)
    f_v = friction_factor_churchill(K, D, G * D / mu_v)
    dP_lo = _friction_dp(f_l, G, rho_l, D, L)
    Fr = G**2 / (G_GRAVITY * D * rho_h**2)
    We = G**2 * D / (sigma * rho_h)
    F = x**0.78 * (1 - x)**0.224
    mu_ratio = mu_v / mu_l
    H = (rho_l / rho_v)**0.91 * mu_ratio**0.19 * (1 - mu_ratio)**0.7
    E = (1 - x)**2 + x**2 * (rho_l * f_v) / (rho_v * f_l)
    Phi_l2 = E + (3.24 * F * H) / max(EPS, Fr**0.0454 * We**0.035)
    return dP_lo * Phi_l2


def _pipe_friedel_dp(AS, geom, cond):
    """Friedel (1979), separated flow, recommended for mu_l / mu_v < 1000: dP = Phi_lo^2 dP_lo.
    Churchill friction factors, homogeneous density in Fr and We. Needs the surface tension
    (cond['sigma'] if CoolProp has none).
    Validity (fluids.two_phase.Friedel notes): works poorly when mu_l / mu_v > 1000; tested on
    diameters down to 4 mm. Warning outside."""
    def model(G, x, s, D, L, K, c):
        _need_sigma(AS, s)
        if s['mu_l'] / s['mu_v'] > 1000:
            warnings.warn(f"Friedel two-phase DP: mu_l / mu_v = {s['mu_l'] / s['mu_v']:.0f} > 1000 "
                          "(correlation known to work poorly there).", ValidityRangeWarning,
                          stacklevel=4)
        if D < 4e-3:
            warnings.warn(f"Friedel two-phase DP: D = {D * 1e3:.2f} mm below 4 mm (smallest tested "
                          "diameter): extrapolated.", ValidityRangeWarning, stacklevel=4)
        return _friedel(G, x, s, D, L, K, c)
    return _two_phase_dp(AS, geom, cond, model)


def _msh(G, x, s, D, L, K, cond):
    rho_l, rho_v, mu_l, mu_v = s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v']
    dP_lo = _friction_dp(friction_factor_churchill(K, D, G * D / mu_l), G, rho_l, D, L)
    dP_vo = _friction_dp(friction_factor_churchill(K, D, G * D / mu_v), G, rho_v, D, L)
    G_MSH = dP_lo + 2 * (dP_vo - dP_lo) * x
    return G_MSH * (1 - x)**(1 / 3) + dP_vo * x**3


def _pipe_msh_dp(AS, geom, cond):
    """Muller-Steinhagen & Heck (1986): empirical interpolation between liquid-only and vapor-only pressure drops
    (9300 points, D 4-392 mm). Churchill friction factors are used here instead of Swamee-Jain, which is
    turbulent only (the liquid-only flow of a tube can be laminar). No surface tension needed.
    Validity: tube diameters of the data bank 4-392 mm (Muller-Steinhagen & Heck 1986, Chem. Eng.
    Process. 20, 297-308, Appendix A, Tables A-1 and A-2); warning outside."""
    # Diameter range of the data bank (Appendix A of the paper)
    D = _tube_ID(geom)
    if D < 4e-3 or D > 392e-3:
        warnings.warn(f"MSH two-phase DP: D = {D * 1e3:.3g} mm outside 4-392 mm (data bank): "
                      "extrapolated.", ValidityRangeWarning, stacklevel=3)
    return _two_phase_dp(AS, geom, cond, _msh)


def _pipe_choi_dp(AS, geom, cond):
    """Choi, Kedzierski & Domanski (2001), evaporation / condensation in smooth and micro-fin tubes (modified
    Pierre): needs the quality change of the pipe, cond['x_out'] or cond['h_in'] and cond['h_out'], with cond['x']
    the inlet quality (cond['x_in'] if present). Written for a whole pipe: the quality change is that of the
    cell while L is the full flow length, so use it with a global evaluation (Correlation_Global), not cell by cell.
    Raises ValueError if the quality does not change.
    Validity: NIST data of Choi et al. (2001) and Choi et al. (1999, NISTIR 6333, Tables 2 and 4):
    D_h 4.3-11.1 mm, G 88-620 kg/m2/s, p 0.36-2.5 MPa, refrigerants R12, R22, R32, R125, R134a,
    R407C, R410A, R32/R134a; warning outside the D_h, G and p ranges."""
    def model(G, x, s, D, L, K, c):
        # Ranges of the NIST data (Choi 2001 p. 3-4; NISTIR 6333 Tables 2 and 4)
        out = []
        if D < 4.3e-3 or D > 11.1e-3:
            out.append(f"D_h = {D * 1e3:.3g} mm (4.3-11.1)")
        if G < 88.0 or G > 620.0:
            out.append(f"G = {G:.4g} kg/m2/s (88-620)")
        if c['p'] < 0.36e6 or c['p'] > 2.5e6:
            out.append(f"p = {c['p'] / 1e6:.3g} MPa (0.36-2.5)")
        if out:
            warnings.warn("Choi two-phase DP: outside the NIST data: " + ", ".join(out) + ".",
                          ValidityRangeWarning, stacklevel=4)
        x_su = c.get('x_in', c['x'])
        if 'x_out' in c:
            x_ex = c['x_out']
        elif 'h_in' in c and 'h_out' in c:
            x_ex = x_su + (c['h_out'] - c['h_in']) / s['h_lv']
        else:
            raise KeyError("Choi pressure drop: cond needs 'x_out' (or 'h_in' and 'h_out').")
        x_ex = min(max(x_ex, 0.0), 1.0)
        
        # if abs(x_ex - x_su) < 1e-9:
        #     raise ValueError("Choi pressure drop: the quality does not change over the pipe (K_f = 0).")
            
        rho_su = 1.0 / (x_su / s['rho_v'] + (1 - x_su) / s['rho_l'])
        rho_ex = 1.0 / (x_ex / s['rho_v'] + (1 - x_ex) / s['rho_l'])
        Re_fo = G * D / s['mu_l']
        K_f = abs(x_ex - x_su) * s['h_lv'] / (G_GRAVITY * L)
        f_N = 0.00506 * Re_fo**(-0.0951) * K_f**0.1554
        return f_N * L * (1.0 / rho_ex + 1.0 / rho_su) / D * G**2
    return _two_phase_dp(AS, geom, cond, model)

def _pipe_homogeneous_dp(AS, geom, cond):
    """Homogeneous model: McAdams mixture viscosity, homogeneous density, Churchill friction factor."""
    def model(G, x, s, D, L, K, c):
        rho_h = 1.0 / (x / s['rho_v'] + (1 - x) / s['rho_l'])
        mu_h = 1.0 / (x / s['mu_v'] + (1 - x) / s['mu_l'])
        return _friction_dp(friction_factor_churchill(K, D, G * D / mu_h), G, rho_h, D, L)
    return _two_phase_dp(AS, geom, cond, model)


def _pipe_lockhart_martinelli_dp(AS, geom, cond):
    """Lockhart-Martinelli, Chisholm form: dP = Phi_l^2 dP_l, Phi_l^2 = 1 + C / X + 1 / X^2,
    X^2 = dP_l / dP_v (each phase alone). C = 20 (tt), 12 (vt), 10 (tv), 5 (vv); Re limit 2000.
    cond['lm_C'] imposes C. Reference: Chisholm (1967); Lockhart & Martinelli (1949).
    Validity: developed for horizontal flow (fluids.two_phase.Lockhart_Martinelli notes): warning
    if the tube is inclined (geom['gamma'] not 0).
    Original data: air with benzene, kerosene, water and oils, pipe diameters 1.5-25.8 mm
    (Xu et al. 2012, Nucl. Eng. Des., Table 1): warning outside this diameter range."""
    if geom.get('gamma', 0.0) != 0.0:
        warnings.warn(f"Lockhart-Martinelli two-phase DP: inclined tube "
                      f"(gamma = {geom['gamma']} deg), "
                      "the correlation was developed for horizontal flow.", ValidityRangeWarning,
                      stacklevel=3)
    # Diameter range of the original data (Xu et al. 2012, Table 1)
    D_check = _tube_ID(geom)
    if D_check < 1.5e-3 or D_check > 25.8e-3:
        warnings.warn(f"Lockhart-Martinelli two-phase DP: D = {D_check * 1e3:.3g} mm outside "
                      "1.5-25.8 mm (original data): extrapolated.", ValidityRangeWarning,
                      stacklevel=3)

    def model(G, x, s, D, L, K, c):
        G_l, G_v = G * (1 - x), G * x
        dP_l = _friction_dp(friction_factor_churchill(K, D, G_l * D / s['mu_l']), G_l, s['rho_l'], D, L)
        dP_v = _friction_dp(friction_factor_churchill(K, D, G_v * D / s['mu_v']), G_v, s['rho_v'], D, L)
        X = np.sqrt(dP_l / dP_v)
        if 'lm_C' in c:
            C = c['lm_C']
        else:
            tl, tv = G_l * D / s['mu_l'] >= 2000, G_v * D / s['mu_v'] >= 2000
            C = {(True, True): 20.0, (False, True): 12.0, (True, False): 10.0, (False, False): 5.0}[(tl, tv)]
        return (1.0 + C / X + 1.0 / X**2) * dP_l
    return _two_phase_dp(AS, geom, cond, model)


def pipe_xu_fang_condensing_DP(AS, geom, cond):
    """
    Xu & Fang (2013) two-phase frictional pressure drop for CONDENSING flow in tubes.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state (used for saturated properties at cond['p']).
    geom : dict
        Geometry dictionary with keys:
        - 'D_t_i' or ('D_t_o', 't_t') : tube inner diameter [m]
        - 'L' : total flow length (straight length x passes) [m]
        - 'K' : absolute roughness [m], default 0.0
        - 'gamma' : inclination angle [deg], default 0.0 (positive = upward flow)
    cond : dict
        Operating conditions with keys:
        - 'G' : mass flux in ONE tube [kg/m2/s] (or 'm_dot' with geom 'n_t', 'n_passes')
        - 'p' : saturation pressure [Pa]
        - 'x' : vapour quality [-] (clipped to [1e-6, 1-1e-6])
        - optional 'sigma' : surface tension [N/m] if CoolProp has no model

    Returns
    -------
    dp : float
        Pressure drop [Pa] over geom['L'] = frictional + gravity (geom['gamma']),
        no acceleration term (consistent with other pipe DP functions in this file).

    Reference
    ---------
    Xu, Y., Fang, X. (2013). A new correlation of two-phase frictional pressure drop
    for condensing flow in pipes. Nuclear Engineering and Design, 263, 87-96.
    doi:10.1016/j.nucengdes.2013.04.017

    Implementation follows the independent implementation in the fluids library
    (C. Bell), fluids.two_phase.Xu_Fang, which cites the reference above.

    Validity
    --------
    Condensing flow, database of the paper (Xu & Fang 2013, Sec. 2, Table 1, Table 2): 525 points,
    R134a, NH3, R410A, R22, R236ea, R125, R32, R600a, R290; D_h 0.1-10.07 mm, G 20-800 kg/m2/s,
    q 2-55.3 kW/m2, x 0-1; fitted parameters Y = sqrt(dp_vo / dp_lo) 2.05-7.29,
    Fr_tp 0.28-95354, We_tp 5.31-12270 (homogeneous density). Warning outside the D, G, Y, Fr_tp
    and We_tp ranges (q is not an input of this function).

    Notes
    -----
    - Diameter in Reynolds numbers is the tube inner diameter.
    - Saturated liquid/vapour properties at cond['p'] are used.
    - Quality is the local quality cond['x'].
    - Churchill friction factor (Darcy) is used for single-phase friction factors,
      whereas the reference library uses its own Colebrook solution.
    """
    # 1. Read geometry and flow conditions
    D = _tube_ID(geom)                    # tube inner diameter [m]
    L = geom['L']                         # total flow length [m]
    K = geom.get('K', 0.0)                # absolute roughness [m]
    G = _tube_G(geom, cond)               # mass flux in one tube [kg/m2/s]
    x = min(max(cond['x'], 1e-6), 1 - 1e-6)  # quality, clipped
    p = cond['p']                         # saturation pressure [Pa]

    # 2. Saturated properties at pressure p
    s = _sat_props(AS, p, cond)           # rho_l, rho_v, mu_l, mu_v, (sigma)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    mu_v = s['mu_v']

    # Surface tension is needed for Weber number
    _need_sigma(AS, s)
    sigma = s['sigma']

    # 3. Liquid-only flow
    Re_lo = G * D / mu_l
    f_lo = friction_factor_churchill(K, D, Re_lo)
    dp_lo = _friction_dp(f_lo, G, rho_l, D, L)

    # 4. Vapour-only flow
    Re_vo = G * D / mu_v
    f_vo = friction_factor_churchill(K, D, Re_vo)
    dp_vo = _friction_dp(f_vo, G, rho_v, D, L)

    # 5. Ratio of vapour-only to liquid-only pressure drop
    Y2 = dp_vo / dp_lo

    # 6. Homogeneous density and velocity
    rho_h = 1.0 / (x / rho_v + (1.0 - x) / rho_l)
    v_h = G / rho_h

    # 7. Froude number (squared form, as in fluids library) and Weber number
    Fr = v_h**2 / (G_GRAVITY * D)
    We = v_h**2 * D * rho_h / sigma

    # Range of the database (Xu & Fang 2013, Table 1 and text after Eq. 3, Table 2)
    Y = Y2**0.5
    out = []
    if D < 0.1e-3 or D > 10.07e-3:
        out.append(f"D = {D * 1e3:.3g} mm (0.1-10.07)")
    if G < 20.0 or G > 800.0:
        out.append(f"G = {G:.4g} kg/m2/s (20-800)")
    if Y < 2.05 or Y > 7.29:
        out.append(f"Y = {Y:.3g} (2.05-7.29)")
    if Fr < 0.28 or Fr > 95354.0:
        out.append(f"Fr_tp = {Fr:.3g} (0.28-95354)")
    if We < 5.31 or We > 12270.0:
        out.append(f"We_tp = {We:.3g} (5.31-12270)")
    if out:
        warnings.warn("Xu-Fang condensing DP: outside the database: " + ", ".join(out) + ".",
                      ValidityRangeWarning, stacklevel=2)

    # 8. Two-phase multiplier phi_lo^2 (Xu & Fang 2013, as implemented in fluids)
    term_vapour = Y2 * x**3
    term_liquid = (1.0 - x**2.59)**0.632
    bracket = (1.0
               + 2.0 * x**1.17 * (Y2 - 1.0)
               + 0.00775 * x**-0.475 * Fr**0.535 * We**0.188)
    phi_lo2 = term_vapour + term_liquid * bracket

    # 9. Frictional pressure drop + gravity term with homogeneous density
    dp_friction = phi_lo2 * dp_lo
    return dp_friction + _gravity(geom, rho_h)


def pipe_chisholm_DP(AS, geom, cond):
    """Two-phase frictional pressure drop in tubes, Chisholm (1973) B-coefficient method
    (evaporating flow; used for evaporation and condensation).

    Implementation follows the independent implementation in the fluids library (C. Bell),
    fluids.two_phase.Chisholm (without the optional 1978 roughness correction), which cites
    the reference below (also in Thome, Wolverine Engineering Data Book III, 2004).

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state object.
    geom : dict
        Geometry dictionary with keys:
        - 'D_t_i' [m] (or 'D_t_o' [m] and 't_t' [m])
        - 'L' [m]: full flow length
        - 'K' [m]: absolute roughness (default 0)
        - 'gamma' [deg]: inclination angle (default 0, positive = upward)
    cond : dict
        Condition dictionary with keys:
        - 'p' [Pa]: pressure
        - 'x' [-]: vapor quality
        - 'G' [kg/m2/s]: mass flux in ONE tube (or 'm_dot' [kg/s] + geom n_t, n_passes)

    Returns
    -------
    float
        Total pressure drop over geom['L'] [Pa] = friction + gravity (geom['gamma']),
        no acceleration.

    Reference
    ---------
    Chisholm, D. (1973) Pressure gradients due to friction during the flow of evaporating
    two-phase mixtures in smooth tubes and channels. International Journal of Heat and Mass
    Transfer 16(2), 347-358. doi:10.1016/0017-9310(73)90063-X

    Notes: Darcy friction factors from friction_factor_churchill (labothappy choice, as Friedel
    and MSH of this file); the fluids library uses its own Colebrook solution. Saturated
    properties at cond['p']; liquid-only / vapour-only flows with the total mass flux.

    Validity: derived for TURBULENT flow in SMOOTH tubes (Blasius exponent n, Chisholm 1973,
    abstract and Sec. 7); no numerical data range in the paper. Warning when the liquid-only flow
    is laminar (Re_lo < 2300).
    """
    # 1. Read tube inner diameter [m], full length [m], and absolute roughness [m]
    D = _tube_ID(geom)
    L = geom['L']
    K = geom.get('K', 0.0)

    # Calculate mass flux in one tube [kg/m2/s]
    G = _tube_G(geom, cond)

    # Read and bound vapor quality [-] between 1e-6 and 1 - 1e-6
    x = cond['x']
    x = min(max(x, 1e-6), 1.0 - 1e-6)

    # Read pressure [Pa]
    p = cond['p']

    # 2. Get saturated fluid properties: densities [kg/m3] and dynamic viscosities [Pa s]
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    mu_v = s['mu_v']

    # 3. Liquid-only flow pressure drop calculation
    Re_lo = G * D / mu_l
    f_lo = friction_factor_churchill(K, D, Re_lo)
    dp_lo = _friction_dp(f_lo, G, rho_l, D, L)

    # Chisholm (1973) is derived for turbulent flow: warning if the liquid-only flow is laminar
    if Re_lo < 2300.0:
        warnings.warn(f"Chisholm two-phase DP: Re_lo = {Re_lo:.0f} < 2300 (laminar liquid-only "
                      "flow, outside the turbulent derivation).", ValidityRangeWarning,
                      stacklevel=2)

    # 4. Vapour-only flow pressure drop calculation
    Re_vo = G * D / mu_v
    f_vo = friction_factor_churchill(K, D, Re_vo)
    dp_vo = _friction_dp(f_vo, G, rho_v, D, L)

    # 5. Calculate Gamma parameter and Blasius exponent
    Gamma = (dp_vo / dp_lo) ** 0.5
    n = 0.25

    # 6. Calculate coefficient B based on Gamma and mass flux G (Chisholm 1973)
    if Gamma <= 9.5:
        if G <= 500.0:
            B = 4.8
        elif 500.0 < G < 1900.0:
            B = 2400.0 / G
        else:
            B = 55.0 / (G ** 0.5)
    elif Gamma <= 28.0:
        if G <= 600.0:
            B = 520.0 / (Gamma * (G ** 0.5))
        else:
            B = 21.0 / Gamma
    else:
        B = 15000.0 / ((Gamma ** 2) * (G ** 0.5))

    # 7. Calculate two-phase multiplier phi2
    exponent = (2.0 - n) / 2.0
    bracket = B * (x ** exponent) * ((1.0 - x) ** exponent) + (x ** (2.0 - n))
    phi2 = 1.0 + ((Gamma ** 2) - 1.0) * bracket

    # 8. Calculate frictional pressure drop
    dp_friction = phi2 * dp_lo

    # Calculate homogeneous density [kg/m3] for the gravity term
    rho_h = 1.0 / (x / rho_v + (1.0 - x) / rho_l)

    # Return total pressure drop: frictional + gravity
    return dp_friction + _gravity(geom, rho_h)


def pipe_gronnerud_boiling_DP(AS, geom, cond):
    """
    Gronnerud (1972) two-phase frictional pressure drop for evaporating refrigerants in tubes.

    Implementation follows the independent implementation in the fluids library
    (C. Bell), fluids.two_phase.Gronnerud, which cites the reference below as presented
    in the ASHRAE Handbook Fundamentals (2013) and Thome, Wolverine Engineering Data Book III
    (2004).

    Parameters
    ----------
    AS : CoolProp.AbstractState
        AbstractState for the fluid (used only for saturation properties via _sat_props).
    geom : dict
        Geometry dictionary with keys:
        - 'D_t_i' or ('D_t_o', 't_t') : tube inner diameter [m]
        - 'L' : total flow length (straight length x passes) [m]
        - 'K' : absolute roughness [m], default 0.0
        - 'gamma' : inclination angle [deg], default 0.0 (positive = upward flow)
    cond : dict
        Operating conditions dictionary with keys:
        - 'G' : mass flux in ONE tube [kg/m2/s] (or 'm_dot' with geom['n_t'], geom['n_passes'])
        - 'x' : vapor quality [-] (clipped to [1e-6, 1-1e-6])
        - 'p' : saturation pressure [Pa]

    Returns
    -------
    dp : float
        Total pressure drop [Pa] over geom['L'] = frictional (Gronnerud) + gravity (geom['gamma']),
        no acceleration term (consistent with the other pipe DP functions in this file).

    Reference
    ---------
    Gronnerud, R. (1972) Investigation of liquid hold-up, flow resistance and heat transfer in
    circulation type evaporators. 4. Two-phase flow resistance in boiling refrigerants. Proc.
    Freudenstadt Meeting, IIR, 1972-1, 127-138.

    Validity
    --------
    Developed for evaporators; 0 < x < 1 (fluids library notes). The original range (IIR 1972) is
    not available. Tested by Ould Didi, Kattan & Thome (2002, Int. J. Refrig. 25, 935-947):
    D 10.92 and 12.0 mm, G 100-500 kg/m2/s, x 0.04-1, R134a, R123, R402A, R404A, R502.
    Warning outside D 10-18 mm and G 100-500 kg/m2/s. The 10-18 mm diameter window is a
    LaboThapPy choice (Basile, 2026-11-26), wider than the two tested diameters: no source.

    Notes
    -----
    - Single-phase friction factor uses Churchill (1977) via friction_factor_churchill (Darcy),
      whereas the reference fluids library uses its own Colebrook solution.
    - Surface tension is not required by this correlation.
    - Homogeneous density is used for the gravity term: rho_h = 1 / (x/rho_v + (1-x)/rho_l).
    """
    # --- 1. Read geometry and conditions (common to all pipe DP functions) ---
    D = _tube_ID(geom)                     # tube inner diameter [m]
    L = geom['L']                          # total flow length [m]
    K = geom.get('K', 0.0)                 # absolute roughness [m]
    G = _tube_G(geom, cond)                # mass flux in one tube [kg/m2/s]
    x = min(max(cond['x'], 1e-6), 1 - 1e-6)  # vapor quality [-], clipped
    p = cond['p']                          # saturation pressure [Pa]

    # Validity: tested G range (Ould Didi 2002), D window 10-18 mm (LaboThapPy choice)
    out = []
    if D < 10e-3 or D > 18e-3:
        out.append(f"D = {D * 1e3:.3g} mm (10-18)")
    if G < 100.0 or G > 500.0:
        out.append(f"G = {G:.4g} kg/m2/s (100-500)")
    if out:
        warnings.warn("Gronnerud two-phase DP: outside the tested range: " + ", ".join(out) + ".",
                      ValidityRangeWarning, stacklevel=2)

    # --- 2. Saturated properties at pressure p ---
    s = _sat_props(AS, p, cond)            # dict with rho_l, rho_v, mu_l, mu_v (sigma not needed)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    mu_v = s['mu_v']

    # --- 3. Liquid-only Froude number (Gronnerud 1972, as implemented in fluids) ---
    v_lo = G / rho_l                       # liquid-only velocity [m/s]
    Fr_lo = v_lo**2 / (G_GRAVITY * D)      # liquid-only Froude number [-]

    # --- 4. Froude function f_Fr (Gronnerud 1972, as implemented in fluids) ---
    if Fr_lo >= 1.0:
        f_Fr = 1.0
    else:
        log_term = np.log(1.0 / Fr_lo)
        f_Fr = Fr_lo**0.3 + 0.0055 * log_term**2

    # --- 5. Friction gradient factor dpdz_Fr (Gronnerud 1972, as implemented in fluids) ---
    dpdz_Fr = f_Fr * (x + 4.0 * (x**1.8 - x**10 * f_Fr**0.5))

    # --- 6. Two-phase multiplier phi_gd (Gronnerud 1972, as implemented in fluids) ---
    property_ratio = (rho_l / rho_v) / (mu_l / mu_v)**0.25
    phi_gd = 1.0 + dpdz_Fr * (property_ratio - 1.0)

    # --- 7. Liquid-only pressure drop dp_lo (Gronnerud 1972, as implemented in fluids) ---
    Re_lo = G * D / mu_l                   # liquid-only Reynolds number [-]
    f_lo = friction_factor_churchill(K, D, Re_lo)  # Darcy friction factor [-]
    dp_lo = _friction_dp(f_lo, G, rho_l, D, L)     # liquid-only frictional DP [Pa]

    # --- 8. Two-phase frictional DP + gravity (Gronnerud 1972, as implemented in fluids) ---
    dp_friction = phi_gd * dp_lo           # two-phase frictional DP [Pa]

    # Homogeneous density for gravity term
    rho_h = 1.0 / (x / rho_v + (1.0 - x) / rho_l)
    dp_gravity = _gravity(geom, rho_h)     # gravity DP [Pa]

    return dp_friction + dp_gravity


def pipe_kim_mudawar_condensing_DP(AS, geom, cond):
    """
    Kim & Mudawar (2012) universal approach for two-phase frictional pressure drop in
    adiabatic and condensing mini/micro-channel flows.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state (used for critical pressure and saturation properties).
    geom : dict
        Geometry dictionary with keys:
        - 'D_t_i' or ('D_t_o', 't_t') : tube inner diameter [m]
        - 'L' : total flow length [m]
        - 'gamma' : inclination angle [deg], default 0 (positive = upward)
    cond : dict
        Flow conditions with keys:
        - 'G' : mass flux in one tube [kg/m2/s] (or 'm_dot' with 'n_t', 'n_passes')
        - 'x' : vapour quality [-]
        - 'p' : pressure [Pa]

    Returns
    -------
    dP : float
        Total pressure drop over geom['L'] [Pa] = frictional + gravity (geom['gamma']),
        no acceleration term (consistent with other pipe DP functions in this file).

    Reference
    ---------
    Kim, S.-M., Mudawar, I. (2012) Universal approach to predicting two-phase frictional
    pressure drop for adiabatic and condensing mini/micro-channel flows.
    International Journal of Heat and Mass Transfer 55(11-12), 3246-3261.
    doi:10.1016/j.ijheatmasstransfer.2012.02.047

    Implementation follows the independent implementation in the fluids library
    (C. Bell), fluids.two_phase.Kim_Mudawar, which cites the reference above
    (also Kim & Mudawar 2014, IJHMT 77, 74-97, doi:10.1016/j.ijheatmasstransfer.2014.04.035).

    Validity
    --------
    Hydraulic diameter 0.0695-6.22 mm; mass flux 4-8528 kg/m2/s; reduced pressure
    0.0052-0.91 (p / p_critical). 7115 data points from 36 sources; fluids: air, CO2,
    N2, water, ethanol, R12, R22, R134a, R236ea, R245fa, R404A, R407C, propane,
    methane, ammonia.

    Notes
    -----
    This model uses its OWN Darcy friction factor (from the paper), roughness is NOT used:
    Re < 2000: f = 64 / Re; 2000 <= Re < 20000: f = 0.316 * Re**-0.25; Re >= 20000:
    f = 0.184 * Re**-0.2. The branch logic is written inline for liquid and vapour phases.
    """
    # 1. Read geometry and flow conditions
    D = _tube_ID(geom)                    # hydraulic diameter [m]
    L = geom['L']                         # total flow length [m]
    G = _tube_G(geom, cond)               # mass flux in one tube [kg/m2/s]
    x = min(max(cond['x'], 1e-6), 1 - 1e-6)  # quality, clipped to avoid 0/1
    p = cond['p']                         # pressure [Pa]

    # 2. Saturated properties at pressure p
    s = _sat_props(AS, p, cond)           # rho_l, rho_v, mu_l, mu_v
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    mu_v = s['mu_v']

    # Surface tension (needed for Suratman number)
    _need_sigma(AS, s)
    sigma = s['sigma']

    # 3. Validity checks (warn if outside ranges)
    # Hydraulic diameter: 0.0695-6.22 mm
    if D < 0.0695e-3 or D > 6.22e-3:
        warnings.warn(
            f"Kim-Mudawar condensing: D_h = {D*1000:.2f} mm outside validity range "
            f"[0.0695, 6.22] mm", ValidityRangeWarning,
            stacklevel=2
        )
    # Mass flux: 4-8528 kg/m2/s
    if G < 4 or G > 8528:
        warnings.warn(
            f"Kim-Mudawar condensing: G = {G:.1f} kg/m2/s outside validity range "
            f"[4, 8528] kg/m2/s", ValidityRangeWarning,
            stacklevel=2
        )
    # Reduced pressure: 0.0052-0.91
    p_crit = AS.p_critical()              # critical pressure [Pa]
    p_red = p / p_crit
    if p_red < 0.0052 or p_red > 0.91:
        warnings.warn(
            f"Kim-Mudawar condensing: reduced pressure p/p_crit = {p_red:.4f} "
            f"outside validity range [0.0052, 0.91]", ValidityRangeWarning,
            stacklevel=2
        )

    # 4. Liquid flowing alone (Kim & Mudawar 2012)
    G_l = G * (1 - x)                     # liquid mass flux [kg/m2/s]
    Re_l = G_l * D / mu_l                 # liquid Reynolds number [-]

    # Liquid friction factor (paper's own correlation, no roughness)
    if Re_l < 2000:
        f_l = 64.0 / Re_l
    elif Re_l < 20000:
        f_l = 0.316 * Re_l ** -0.25
    else:
        f_l = 0.184 * Re_l ** -0.2

    dp_l = _friction_dp(f_l, G_l, rho_l, D, L)  # liquid-only pressure drop [Pa]

    # 5. Vapour flowing alone (Kim & Mudawar 2012)
    G_v = G * x                           # vapour mass flux [kg/m2/s]
    Re_v = G_v * D / mu_v                 # vapour Reynolds number [-]

    # Vapour friction factor (paper's own correlation, no roughness)
    if Re_v < 2000:
        f_v = 64.0 / Re_v
    elif Re_v < 20000:
        f_v = 0.316 * Re_v ** -0.25
    else:
        f_v = 0.184 * Re_v ** -0.2

    dp_v = _friction_dp(f_v, G_v, rho_v, D, L)  # vapour-only pressure drop [Pa]

    # 6. Liquid-only Reynolds number (Kim & Mudawar 2012)
    Re_lo = G * D / mu_l                  # based on total mass flux

    # 7. Suratman number of the vapour (Kim & Mudawar 2012)
    Su_vo = rho_v * sigma * D / (mu_v ** 2)

    # 8. Martinelli parameter X and transition Reynolds number (Kim & Mudawar 2012)
    X = (dp_l / dp_v) ** 0.5
    Re_c = 2000.0                         # transition Reynolds number

    # 9. Chisholm parameter C (Kim & Mudawar 2012) - four regimes
    rho_ratio = rho_l / rho_v

    if Re_l < Re_c and Re_v < Re_c:
        # laminar liquid, laminar vapour
        C = 3.5e-5 * Re_lo ** 0.44 * Su_vo ** 0.50 * rho_ratio ** 0.48
    elif Re_l < Re_c and Re_v >= Re_c:
        # laminar liquid, turbulent vapour
        C = 0.0015 * Re_lo ** 0.59 * Su_vo ** 0.19 * rho_ratio ** 0.36
    elif Re_l >= Re_c and Re_v < Re_c:
        # turbulent liquid, laminar vapour
        C = 8.7e-4 * Re_lo ** 0.17 * Su_vo ** 0.50 * rho_ratio ** 0.14
    else:
        # turbulent liquid, turbulent vapour
        C = 0.39 * Re_lo ** 0.03 * Su_vo ** 0.10 * rho_ratio ** 0.35

    # 10. Two-phase multiplier and frictional pressure drop (Kim & Mudawar 2012)
    phi_l2 = 1.0 + C / X + 1.0 / (X ** 2)
    dp_friction = phi_l2 * dp_l

    # 11. Homogeneous density for gravity term
    rho_h = 1.0 / (x / rho_v + (1 - x) / rho_l)

    # 12. Total pressure drop = friction + gravity (no acceleration)
    return dp_friction + _gravity(geom, rho_h)


#%% PUBLIC NAMES: <phase> TAG, as the pipe HTC (same implementation for evaporation and condensation)

def pipe_friedel_boiling_DP(AS, geom, cond):
    """Friedel two-phase frictional pressure drop, evaporation (see _pipe_friedel_dp)."""
    return _pipe_friedel_dp(AS, geom, cond)


def pipe_friedel_condensing_DP(AS, geom, cond):
    """Friedel two-phase frictional pressure drop, condensation (see _pipe_friedel_dp)."""
    return _pipe_friedel_dp(AS, geom, cond)


def pipe_msh_boiling_DP(AS, geom, cond):
    """Msh two-phase frictional pressure drop, evaporation (see _pipe_msh_dp)."""
    return _pipe_msh_dp(AS, geom, cond)


def pipe_msh_condensing_DP(AS, geom, cond):
    """Msh two-phase frictional pressure drop, condensation (see _pipe_msh_dp)."""
    return _pipe_msh_dp(AS, geom, cond)


def pipe_choi_boiling_DP(AS, geom, cond):
    """Choi two-phase frictional pressure drop, evaporation (see _pipe_choi_dp)."""
    return _pipe_choi_dp(AS, geom, cond)


def pipe_choi_condensing_DP(AS, geom, cond):
    """Choi two-phase frictional pressure drop, condensation (see _pipe_choi_dp)."""
    return _pipe_choi_dp(AS, geom, cond)


def pipe_homogeneous_boiling_DP(AS, geom, cond):
    """Homogeneous two-phase frictional pressure drop, evaporation (see _pipe_homogeneous_dp)."""
    return _pipe_homogeneous_dp(AS, geom, cond)


def pipe_homogeneous_condensing_DP(AS, geom, cond):
    """Homogeneous two-phase frictional pressure drop, condensation (see _pipe_homogeneous_dp)."""
    return _pipe_homogeneous_dp(AS, geom, cond)


def pipe_lockhart_martinelli_boiling_DP(AS, geom, cond):
    """Lockhart Martinelli two-phase frictional pressure drop, evaporation (see _pipe_lockhart_martinelli_dp)."""
    return _pipe_lockhart_martinelli_dp(AS, geom, cond)


def pipe_lockhart_martinelli_condensing_DP(AS, geom, cond):
    """Lockhart Martinelli two-phase frictional pressure drop, condensation (see _pipe_lockhart_martinelli_dp)."""
    return _pipe_lockhart_martinelli_dp(AS, geom, cond)