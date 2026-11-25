# -*- coding: utf-8 -*-
"""
Shell-side pressure drop of shell & tube heat exchangers, in the (AS, geom, cond) formalism
of the STHE shell HTC correlations (same geometry helpers: coherent HTC and DP of a side).

Each function returns the pressure drop of the whole shell side [Pa].

cond : m_dot (flow of the whole side) [kg/s], mu [Pa s], rho or D (density) [kg/m3],
       mu_w (optional, mu_w = mu if absent)
geom : see each function; n_B (number of baffles; default round(L_t / B) - 1);
       n_parallel (optional) divides the flow rate.
"""

import warnings

import numpy as np

from labothappy.correlations.convection.STHE_HTC_correlation_bank import (
    s_max_kern, s_max, s_L, d_h, bell_delaware_coefs)
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props

#%% HELPERS

def _rho(cond):
    return cond['rho'] if 'rho' in cond else cond['D']

def _n_baffles(geom):
    if 'n_B' in geom:
        return geom['n_B']
    return np.round(geom['L_t'] / geom['B']) - 1

#%% KERN

def STHE_kern_shell_DP(AS, geom, cond):
    """Kern shell-side pressure drop [Pa].
        DP = f G^2 Ds (N_b + 1) / (2 rho De phi_s),  f = exp(0.576 - 0.19 ln Re),  phi_s = (mu / mu_w)^0.14

    cond : m_dot, mu, rho (or D), mu_w (optional)
    geom : D_t_o, PR, D_s, B, theta, n_B (or L_t)
    Reference: Process Heat Transfer - D. Q. Kern"""
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)
    rho, mu = _rho(cond), cond['mu']
    mu_w = cond.get('mu_w', mu)

    S_T = s_max_kern(geom['D_t_o'], geom['PR'], geom['D_s'],
                     geom['B'], geom['theta'])          # m^2
    D_hydro = d_h(geom['D_t_o'], geom['PR'], geom['theta'])

    G = m_dot / S_T                                                         # kg/(m^2 s)
    Re = G * D_hydro / mu
    f = np.exp(0.576 - 0.19 * np.log(Re))                                   # [-]
    phi_s = (mu / mu_w)**0.14

    return (f * G**2 * geom['D_s'] * (_n_baffles(geom) + 1)) / (2 * rho * D_hydro * phi_s)

#%% DONOHUE

def STHE_donohue_shell_DP(AS, geom, cond):
    """Donohue-type shell-side pressure drop [Pa]: Kern's form with the Reynolds number of the
    Donohue HTC (geometric mean of the cross-flow and window velocities, on the tube diameter), no
    viscosity correction.

    cond : m_dot, mu, rho (or D)
    geom : D_t_o, PR, D_s, B, theta, BC [%], n_t,
           n_B (or L_t)
    Reference: Process Heat Transfer - D. Q. Kern"""
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)
    rho, mu = _rho(cond), cond['mu']
    d_o = geom['D_t_o']

    S_T = s_max_kern(d_o, geom['PR'], geom['D_s'],
                     geom['B'], geom['theta'])          # m^2 (as in the Donohue HTC)
    S_L = s_L(geom['BC'], geom['D_s'], geom['n_t'], d_o)   # m^2
    D_hydro = d_h(d_o, geom['PR'], geom['theta'])

    # Re = rho V_R d_o / mu, V_R = sqrt(V_L V_T), V = m_dot / (rho S): the density cancels
    Re = m_dot * d_o / (mu * np.sqrt(S_T * S_L))
    f = np.exp(0.576 - 0.19 * np.log(Re))
    G = m_dot / S_T

    return (f * G**2 * geom['D_s'] * (_n_baffles(geom) + 1)) / (2 * rho * D_hydro)

#%% BELL-DELAWARE

def STHE_bell_delaware_shell_DP(AS, geom, cond):
    """Bell-Delaware shell-side pressure drop [Pa] (excluding the nozzles): ideal tube bank crossflow,
    baffle window, and the leakage (R_l), bundle bypass (R_b) and unequal baffle spacing (R_s) corrections.

    cond : m_dot, mu, rho (or D), mu_w (optional), cp not needed
    geom : D_t_o, PR, D_s, D_OTL, B, theta, BC [%], delta_TB,
           delta_BS, n_t, N_strips, B_i, B_o, n_B (or L_t)

    The leakage clearances are computed as in STHE_bell_delaware_shell_HTC (same bundle model).

    TO CHECK against the source (kept from the previous version, not verified):
      - the end zones are counted once (DP_PE), not twice
      - the window wetted perimeter uses theta_DS * D_s (not theta_DS * D_s / 2)

    References: https://cheguide.com/heat_exchanger_rating.html"""
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)
    rho, mu = _rho(cond), cond['mu']
    mu_w = cond.get('mu_w', mu)
    p = geom
    N_b = _n_baffles(p)

    "1) Ideal tube bank"

    S_m = s_max(p['D_t_o'], p['PR'], p['D_s'], p['D_OTL'],
                p['B'], p['theta'])                      # m^2
    G_shell = m_dot / S_m
    Re = G_shell * p['D_t_o'] / mu                                         # on the tube diameter, as the HTC

    a_vec, b_vec = bell_delaware_coefs(p['theta'], Re)
    b = b_vec[2] / (1 + 0.14 * Re**b_vec[3])
    f_i = b_vec[0] * (1.33 / p['PR'])**b * Re**b_vec[1]

    P_T = p['D_t_o'] * p['PR']
    if p['theta'] == 0 or p['theta'] == 90:  # square inline
        P_p = P_T
    elif p['theta'] == 60:                 # triangular staggered
        P_p = P_T * (3**0.5 / 2)
    elif p['theta'] == 45:                 # square staggered
        P_p = P_T / (2**0.5)
    else:
        raise ValueError(f"Unknown theta {p['theta']!r}: 0 or 90, 45 or 60.")

    # Effective tube rows crossed between two baffles
    N_tcc = (p['D_s'] / P_p) * (1 - 2 * p['BC'] / 100)

    # Friction viscosity correction (mu_w / mu)^0.14 (the inverse of the heat transfer one)
    DP_ideal = 2 * f_i * (G_shell**2 / rho) * (mu_w / mu)**0.14 * N_tcc

    "2) Baffle leakage : R_l"

    D_CTL = p['D_OTL'] - p['D_t_o']
    theta_CTL = 2 * np.arccos((p['D_s'] * (1 - 2 * p['BC'] / 100)) / D_CTL)
    F_w = (theta_CTL - np.sin(theta_CTL)) / (2 * np.pi)

    theta_DS = 2 * np.arccos(1 - 2 * p['BC'] / 100)
    # Clearance areas (tube-to-baffle and shell-to-baffle), as in the Bell-Delaware HTC
    S_tb = (np.pi / 4) * ((p['D_t_o'] + p['delta_TB'])**2 - p['D_t_o']**2) * p['n_t'] * (1 - F_w)
    S_sb = 0.00436 * p['D_s'] * p['delta_BS'] * (180 / np.pi) * (2 * np.pi - theta_DS)

    r_L = (S_sb + S_tb) / S_m
    r_S = S_sb / (S_sb + S_tb)
    p_exp = 0.8 - 0.15 * (1 + r_S)
    R_l = np.exp(-1.33 * (1 + r_S) * r_L**p_exp)

    "3) Window section"

    SWG = (p['D_s']**2 / 8) * (theta_DS - np.sin(theta_DS))
    SWT = p['n_t'] * F_w * (np.pi * p['D_t_o']**2 / 4)
    SW = SWG - SWT
    GW = m_dot / (S_m * SW)**0.5
    DW = 4 * SW / (np.pi * p['D_t_o'] * p['n_t'] * F_w + theta_DS * p['D_s'])

    N_tcw = (0.8 / P_p) * (p['D_s'] * (p['BC'] / 100) - ((p['D_s'] - (p['D_OTL'] - p['D_t_o'])) / 2))

    if Re >= 100:       # turbulent
        DP_w = N_b * R_l * (2 + 0.6 * N_tcw) * GW**2 / (2 * rho)
    else:               # laminar (viscosity mu, not cp)
        a_coef = (26 * GW * mu) / rho
        b_coef = N_tcw / (P_T - p['D_t_o'])
        c_coef = p['B'] / DW**2
        d_coef = GW**2 / rho
        DP_w = N_b * R_l * (a_coef * (b_coef + c_coef) + d_coef)

    "4) Bundle bypass : R_b"

    r_ss = p['N_strips'] / N_tcc
    if r_ss < 0.5:
        C_r = 3.7 if Re >= 100 else 4.5
        S_b = p['B'] * (p['D_s'] - p['D_OTL'] - p['D_t_o'] / 2)
        R_b = np.exp(-C_r * (S_b / S_m) * (1 - (2 * r_ss)**(1 / 3)))
    else:
        R_b = 1.0

    "5) Unequal baffle spacing : R_s"

    n = 0.2 if Re >= 100 else 1.0
    R_s = 0.5 * ((p['B'] / p['B_i'])**(2 - n)
                 + (p['B'] / p['B_o'])**(2 - n))

    "6) Total"

    DP_PC = (N_b - 1) * DP_ideal * R_l * R_b                 # central baffle spaces
    DP_PE = DP_ideal * (1 + N_tcw / N_tcc) * R_b * R_s       # entrance & exit baffle spaces
    return DP_PE + DP_PC + DP_w

#%% TWO-PHASE : ISHIHARA

def STHE_ishihara_two_phase_DP(AS, geom, cond):
    """
    Frictional two-phase pressure drop [Pa] of the shell side (whole shell) for
    two-phase cross flow over a horizontal tube bundle using the Ishihara
    correlation.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Fluid property provider (moved to saturation by _sat_props).
    geom : dict
        Geometry dictionary containing:
        - 'D_t_o' : Outer diameter of the tubes [m]
        - 'PR' : Ratio of tube pitch to tube OD [-]
        - 'D_s' : Inner diameter of the shell [m]
        - 'B' : Baffle spacing [m]
        - 'theta' : Tube layout angle (0, 45 or 90 = square, others = triangular) [deg]
        - 'n_parallel' : Optional, number of parallel streams (default 1) [-]
        - Additional keys required by STHE_kern_shell_DP (e.g., 'n_B' or 'L_t').
    cond : dict
        Condition dictionary containing:
        - 'm_dot' : Total mass flow rate of the whole side [kg/s]
        - 'p' : Pressure [Pa]
        - 'x' : Vapour quality [-]

    Returns
    -------
    float
        Frictional two-phase pressure drop [Pa] (positive, no momentum or gravity terms).

    References
    ----------
    Ishihara, K., Palen, J.W., Taborek, J. (1980), Critical review of correlations
    for predicting two-phase flow pressure drop across tube banks, Heat Transfer
    Engineering 1(3), 23-33; as given by Consolini, Robinson, Thome (2006),
    Heat Transfer Engineering 27(3), doi 10.1080/01457630500453832,
    Eqs. 14a, 14b, 16a, 16b (C = 8).

    Validity Range
    --------------
    - Mass flux G >= 42 kg/(m^2 s) (emits a warning if lower, as the method
      underpredicts at low mass flux).
    - Vapour quality x between 0 and 1 (clipped with a warning if out of bounds).
    - TODO: Validity beyond the specified mass flux threshold is unknown.
    """
    # Step 1: Quality check and clipping
    x = cond['x']
    if x < 0.0 or x > 1.0:
        warnings.warn(
            f"Vapour quality x = {x} is out of the [0, 1] range. "
            "It has been clipped."
        )
        x = np.clip(x, 0.0, 1.0)

    # Shell-side mass flux calculation
    n_parallel = geom.get('n_parallel', 1)
    m_dot_shell = cond['m_dot'] / n_parallel
    S = s_max_kern(
        geom['D_t_o'],
        geom['PR'],
        geom['D_s'],
        geom['B'],
        geom['theta']
    )
    G = m_dot_shell / S

    # Warning if mass flux is below the recommended threshold
    if G < 42.0:
        warnings.warn(
            f"Mass flux G = {G:.2f} kg/(m^2 s) is below 42 kg/(m^2 s). "
            "The Ishihara method may underpredict pressure drop at low mass flux."
        )

    # Step 2: Saturated properties at pressure p
    sat_props = _sat_props(AS, cond['p'], cond)
    rho_l = sat_props['rho_l']
    rho_v = sat_props['rho_v']
    mu_l = sat_props['mu_l']
    mu_v = sat_props['mu_v']

    # Step 4: Special cases for single-phase extremes
    if x <= 0.0:
        cond_l = dict(cond)
        cond_l['m_dot'] = cond['m_dot']
        cond_l['mu'] = mu_l
        cond_l['rho'] = rho_l
        cond_l['mu_w'] = mu_l
        return STHE_kern_shell_DP(AS, geom, cond_l)

    if x >= 1.0:
        cond_g = dict(cond)
        cond_g['m_dot'] = cond['m_dot']
        cond_g['mu'] = mu_v
        cond_g['rho'] = rho_v
        cond_g['mu_w'] = mu_v
        return STHE_kern_shell_DP(AS, geom, cond_g)

    # Step 3: Single-phase pressure drops using the Kern function
    # Liquid-only single-phase pressure drop (Consolini Eqs. 15a/15b)
    cond_l = dict(cond)
    cond_l['m_dot'] = cond['m_dot'] * (1.0 - x)
    cond_l['mu'] = mu_l
    cond_l['rho'] = rho_l
    cond_l['mu_w'] = mu_l
    DP_l = STHE_kern_shell_DP(AS, geom, cond_l)

    # Vapour-only single-phase pressure drop
    cond_g = dict(cond)
    cond_g['m_dot'] = cond['m_dot'] * x
    cond_g['mu'] = mu_v
    cond_g['rho'] = rho_v
    cond_g['mu_w'] = mu_v
    DP_g = STHE_kern_shell_DP(AS, geom, cond_g)

    # Step 5: Martinelli parameter for turbulent-turbulent flow (Consolini et al.)
    term_x = (1.0 - x) / x
    term_rho = rho_v / rho_l
    term_mu = mu_l / mu_v
    X_tt = (term_x**0.9) * (term_rho**0.5) * (term_mu**0.1)

    # Step 6: Liquid Reynolds number on tube outer diameter
    Re_l = G * geom['D_t_o'] * (1.0 - x) / mu_l

    # Step 7: Two-phase multiplier and total frictional pressure drop
    if Re_l > 2000.0:
        # Re_l > 2000: liquid multiplier (Consolini Eqs. 16a and 14a)
        Phi_l2 = 1.0 + (8.0 / X_tt) + (1.0 / (X_tt**2))
        DP = Phi_l2 * DP_l
    else:
        # Re_l <= 2000: vapour multiplier (Consolini Eqs. 16b and 14b)
        Phi_g2 = 1.0 + (8.0 * X_tt) + (X_tt**2)
        DP = Phi_g2 * DP_g

    return float(DP)
