# -*- coding: utf-8 -*-
"""
Shell-side pressure drop of shell & tube heat exchangers, in the (AS, geom, cond) formalism
of the STHE shell HTC correlations (same geometry helpers: coherent HTC and DP of a side).

Each function returns the pressure drop of the whole shell side [Pa].

cond : m_dot (flow of the whole side) [kg/s], mu [Pa s], rho or D (density) [kg/m3],
       mu_w (optional, mu_w = mu if absent)
geom : see each function; cross_passes (number of baffles; default round(Tube_L / central_spacing) - 1);
       n_parallel (optional) divides the flow rate.
"""

import numpy as np

from labothappy.correlations.convection.STHE_HTC_correlation_bank import (
    s_max_kern, s_max, s_L, d_h, bell_delaware_coefs)

#%% HELPERS

def _rho(cond):
    return cond['rho'] if 'rho' in cond else cond['D']

def _n_baffles(geom):
    if 'cross_passes' in geom:
        return geom['cross_passes']
    return np.round(geom['Tube_L'] / geom['central_spacing']) - 1

#%% KERN

def STHE_kern_shell_DP(AS, geom, cond):
    """Kern shell-side pressure drop [Pa].
        DP = f G^2 Ds (N_b + 1) / (2 rho De phi_s),  f = exp(0.576 - 0.19 ln Re),  phi_s = (mu / mu_w)^0.14

    cond : m_dot, mu, rho (or D), mu_w (optional)
    geom : Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout, cross_passes (or Tube_L)
    Reference: Process Heat Transfer - D. Q. Kern"""
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)
    rho, mu = _rho(cond), cond['mu']
    mu_w = cond.get('mu_w', mu)

    S_T = s_max_kern(geom['Tube_OD'], geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])          # m^2
    D_hydro = d_h(geom['Tube_OD'], geom['pitch_ratio'], geom['tube_layout'])

    G = m_dot / S_T                                                         # kg/(m^2 s)
    Re = G * D_hydro / mu
    f = np.exp(0.576 - 0.19 * np.log(Re))                                   # [-]
    phi_s = (mu / mu_w)**0.14

    return (f * G**2 * geom['Shell_ID'] * (_n_baffles(geom) + 1)) / (2 * rho * D_hydro * phi_s)

#%% DONOHUE

def STHE_donohue_shell_DP(AS, geom, cond):
    """Donohue-type shell-side pressure drop [Pa]: Kern's form with the Reynolds number of the
    Donohue HTC (geometric mean of the cross-flow and window velocities, on the tube diameter), no
    viscosity correction.

    cond : m_dot, mu, rho (or D)
    geom : Tube_OD, pitch_ratio, Shell_ID, central_spacing, tube_layout, Baffle_cut [%], n_tubes,
           cross_passes (or Tube_L)
    Reference: Process Heat Transfer - D. Q. Kern"""
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)
    rho, mu = _rho(cond), cond['mu']
    d_o = geom['Tube_OD']

    S_T = s_max_kern(d_o, geom['pitch_ratio'], geom['Shell_ID'],
                     geom['central_spacing'], geom['tube_layout'])          # m^2 (as in the Donohue HTC)
    S_L = s_L(geom['Baffle_cut'], geom['Shell_ID'], geom['n_tubes'], d_o)   # m^2
    D_hydro = d_h(d_o, geom['pitch_ratio'], geom['tube_layout'])

    # Re = rho V_R d_o / mu, V_R = sqrt(V_L V_T), V = m_dot / (rho S): the density cancels
    Re = m_dot * d_o / (mu * np.sqrt(S_T * S_L))
    f = np.exp(0.576 - 0.19 * np.log(Re))
    G = m_dot / S_T

    return (f * G**2 * geom['Shell_ID'] * (_n_baffles(geom) + 1)) / (2 * rho * D_hydro)

#%% BELL-DELAWARE

def STHE_bell_delaware_shell_DP(AS, geom, cond):
    """Bell-Delaware shell-side pressure drop [Pa] (excluding the nozzles): ideal tube bank crossflow,
    baffle window, and the leakage (R_l), bundle bypass (R_b) and unequal baffle spacing (R_s) corrections.

    cond : m_dot, mu, rho (or D), mu_w (optional), cp not needed
    geom : Tube_OD, pitch_ratio, Shell_ID, D_OTL, central_spacing, tube_layout, Baffle_cut [%], clear_TB,
           clear_BS, n_tubes, N_strips, inlet_spacing, outlet_spacing, cross_passes (or Tube_L)

    The leakage clearances are computed as in STHE_bell_delaware_shell_HTC (same bundle model).

    TO CHECK against the source (kept from the previous version, not verified):
      - the end zones are counted once (DP_PE), not twice
      - the window wetted perimeter uses theta_DS * Shell_ID (not theta_DS * Shell_ID / 2)

    References: https://cheguide.com/heat_exchanger_rating.html"""
    m_dot = cond['m_dot'] / geom.get('n_parallel', 1)
    rho, mu = _rho(cond), cond['mu']
    mu_w = cond.get('mu_w', mu)
    p = geom
    N_b = _n_baffles(p)

    "1) Ideal tube bank"

    S_m = s_max(p['Tube_OD'], p['pitch_ratio'], p['Shell_ID'], p['D_OTL'],
                p['central_spacing'], p['tube_layout'])                      # m^2
    G_shell = m_dot / S_m
    Re = G_shell * p['Tube_OD'] / mu                                         # on the tube diameter, as the HTC

    a_vec, b_vec = bell_delaware_coefs(p['tube_layout'], Re)
    b = b_vec[2] / (1 + 0.14 * Re**b_vec[3])
    f_i = b_vec[0] * (1.33 / p['pitch_ratio'])**b * Re**b_vec[1]

    P_T = p['Tube_OD'] * p['pitch_ratio']
    if p['tube_layout'] == 0:                    # square inline
        P_p = P_T
    elif p['tube_layout'] == 60:                 # triangular staggered
        P_p = P_T * (3**0.5 / 2)
    elif p['tube_layout'] == 45:                 # square staggered
        P_p = P_T / (2**0.5)
    else:
        raise ValueError(f"Unknown tube_layout {p['tube_layout']!r}: 0, 45 or 60.")

    # Effective tube rows crossed between two baffles
    N_tcc = (p['Shell_ID'] / P_p) * (1 - 2 * p['Baffle_cut'] / 100)

    # Friction viscosity correction (mu_w / mu)^0.14 (the inverse of the heat transfer one)
    DP_ideal = 2 * f_i * (G_shell**2 / rho) * (mu_w / mu)**0.14 * N_tcc

    "2) Baffle leakage : R_l"

    D_CTL = p['D_OTL'] - p['Tube_OD']
    theta_CTL = 2 * np.arccos((p['Shell_ID'] * (1 - 2 * p['Baffle_cut'] / 100)) / D_CTL)
    F_w = (theta_CTL - np.sin(theta_CTL)) / (2 * np.pi)

    theta_DS = 2 * np.arccos(1 - 2 * p['Baffle_cut'] / 100)
    # Clearance areas (tube-to-baffle and shell-to-baffle), as in the Bell-Delaware HTC
    S_tb = (np.pi / 4) * ((p['Tube_OD'] + p['clear_TB'])**2 - p['Tube_OD']**2) * p['n_tubes'] * (1 - F_w)
    S_sb = 0.00436 * p['Shell_ID'] * p['clear_BS'] * (180 / np.pi) * (2 * np.pi - theta_DS)

    r_L = (S_sb + S_tb) / S_m
    r_S = S_sb / (S_sb + S_tb)
    p_exp = 0.8 - 0.15 * (1 + r_S)
    R_l = np.exp(-1.33 * (1 + r_S) * r_L**p_exp)

    "3) Window section"

    SWG = (p['Shell_ID']**2 / 8) * (theta_DS - np.sin(theta_DS))
    SWT = p['n_tubes'] * F_w * (np.pi * p['Tube_OD']**2 / 4)
    SW = SWG - SWT
    GW = m_dot / (S_m * SW)**0.5
    DW = 4 * SW / (np.pi * p['Tube_OD'] * p['n_tubes'] * F_w + theta_DS * p['Shell_ID'])

    N_tcw = (0.8 / P_p) * (p['Shell_ID'] * (p['Baffle_cut'] / 100) - ((p['Shell_ID'] - (p['D_OTL'] - p['Tube_OD'])) / 2))

    if Re >= 100:       # turbulent
        DP_w = N_b * R_l * (2 + 0.6 * N_tcw) * GW**2 / (2 * rho)
    else:               # laminar (viscosity mu, not cp)
        a_coef = (26 * GW * mu) / rho
        b_coef = N_tcw / (P_T - p['Tube_OD'])
        c_coef = p['central_spacing'] / DW**2
        d_coef = GW**2 / rho
        DP_w = N_b * R_l * (a_coef * (b_coef + c_coef) + d_coef)

    "4) Bundle bypass : R_b"

    r_ss = p['N_strips'] / N_tcc
    if r_ss < 0.5:
        C_r = 3.7 if Re >= 100 else 4.5
        S_b = p['central_spacing'] * (p['Shell_ID'] - p['D_OTL'] - p['Tube_OD'] / 2)
        R_b = np.exp(-C_r * (S_b / S_m) * (1 - (2 * r_ss)**(1 / 3)))
    else:
        R_b = 1.0

    "5) Unequal baffle spacing : R_s"

    n = 0.2 if Re >= 100 else 1.0
    R_s = 0.5 * ((p['central_spacing'] / p['inlet_spacing'])**(2 - n)
                 + (p['central_spacing'] / p['outlet_spacing'])**(2 - n))

    "6) Total"

    DP_PC = (N_b - 1) * DP_ideal * R_l * R_b                 # central baffle spaces
    DP_PE = DP_ideal * (1 + N_tcw / N_tcc) * R_b * R_s       # entrance & exit baffle spaces
    return DP_PE + DP_PC + DP_w