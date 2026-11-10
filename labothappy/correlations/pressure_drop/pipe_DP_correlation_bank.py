# -*- coding: utf-8 -*-
"""
In-pipe (tube-side) frictional pressure drop, in the (AS, geom, cond) formalism.
Correlations and friction factors from the labothappy pipe pressure drop library (E. Neven), written here
without the AbstractState-based interface: the fluid state comes in `cond`.

Each function returns the pressure drop over the FULL flow length geom['L'] [Pa]: friction + gravity
(geom['theta'], default 0, positive = upward flow). No acceleration (it depends on the quality change of the
cell and cannot be scaled with the length fraction of the cell), no return bends or head losses.

cond : G (mass flux in ONE tube) [kg/m2/s] (or m_dot + geom n_tubes, n_passes), mu [Pa s], rho or D [kg/m3]
       (single phase); p [Pa], x [-] (two-phase; optional cond['mu_l'], cond['mu_v'], cond['sigma'] overrides
       for fluids without a CoolProp transport / surface tension model)
geom : Tube_ID (or Tube_OD and Tube_t), L (total tube flow path: straight length x passes) [m],
       K (absolute roughness) [m], default 0; theta [deg], default 0
"""

import math
import warnings
import numpy as np

from labothappy.correlations.convection.pipe_HTC_correlation_bank import (
    _tube_ID, _tube_G, _sat_props, _need_sigma)

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
    """Swamee & Jain (1976): explicit Colebrook-White, turbulent flow (Re >= 5000)."""
    Re = max(EPS, Re)
    log_arg = K / max(EPS, d_hyd) / 3.7 + 5.74 / Re**0.9
    return 0.25 / (math.log10(max(EPS, log_arg))**2)


def friction_factor_haaland(K, d_hyd, Re):
    """Haaland (1983): explicit Colebrook-White, 4e3 <= Re <= 1e8, 1e-6 <= K/d <= 0.05."""
    return (1.0 / (-1.8 * np.log10((K / (3.75 * d_hyd))**1.11 + 6.9 / Re)))**2


def friction_factor_konakov(Re):
    """Konakov: f = (1.8 log10 Re - 1.5)^-2, smooth pipes, 4e3 <= Re <= 3e6."""
    return (1.8 * math.log10(Re) - 1.5)**(-2)


def friction_factor_petukhov(Re):
    """Petukhov (1970): f = (0.79 ln Re - 1.64)^-2, smooth pipes, 3e3 <= Re <= 5e6."""
    return (0.79 * np.log(Re) - 1.64)**(-2)


def friction_factor_cheng_CO2(G, d_hyd, P, h, mu):
    """Cheng et al. (2024): supercritical CO2 in a horizontal tube (liquid-like / vapor-like / two-phase-like
    by the pseudo-vapor quality). D = 10 mm, G 496.7-1346.2 kg/m2/s, P 7.53-23.51 MPa, CO2 only.
    Needs labothappy.correlations.properties.supercritical (pseudo_boiling_bounds, pseudo_critical_temperature_CO2)."""
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
    return G_GRAVITY * rho * geom['L'] * math.sin(math.radians(geom.get('theta', 0.0)))


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
    (cond['sigma'] if CoolProp has none)."""
    def model(G, x, s, D, L, K, c):
        _need_sigma(AS, s)
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
    turbulent only (the liquid-only flow of a tube can be laminar). No surface tension needed."""
    return _two_phase_dp(AS, geom, cond, _msh)


def _pipe_choi_dp(AS, geom, cond):
    """Choi, Kedzierski & Domanski (2001), evaporation / condensation in smooth and micro-fin tubes (modified
    Pierre): needs the quality change of the pipe, cond['x_out'] or cond['h_in'] and cond['h_out'], with cond['x']
    the inlet quality (cond['x_in'] if present). Written for a whole pipe: the quality change is that of the
    cell while L is the full flow length, so use it with a global evaluation (Correlation_Global), not cell by cell.
    Raises ValueError if the quality does not change."""
    def model(G, x, s, D, L, K, c):
        x_su = c.get('x_in', c['x'])
        if 'x_out' in c:
            x_ex = c['x_out']
        elif 'h_in' in c and 'h_out' in c:
            x_ex = x_su + (c['h_out'] - c['h_in']) / s['h_lv']
        else:
            raise KeyError("Choi pressure drop: cond needs 'x_out' (or 'h_in' and 'h_out').")
        x_ex = min(max(x_ex, 0.0), 1.0)
        if abs(x_ex - x_su) < 1e-9:
            raise ValueError("Choi pressure drop: the quality does not change over the pipe (K_f = 0).")
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
    cond['lm_C'] imposes C. Reference: Chisholm (1967); Lockhart & Martinelli (1949)."""
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