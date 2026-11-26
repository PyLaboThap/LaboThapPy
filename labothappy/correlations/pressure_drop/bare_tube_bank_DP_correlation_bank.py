# -*- coding: utf-8 -*-
"""
Pressure drop on the OUTER side of banks of bare (plain) tubes in crossflow, (AS, geom, cond)
formalism. Replaces the pressure drop functions of the old convection/tube_bank_htc.py
(euler_coeff, euler_coeff_stag, DP_external_flow_*_bank, tube_bank_DP), 2026-11-26.

Returns the pressure drop [Pa] of the WHOLE bank (all rows crossed by the flow), friction only.

geom : D_t_o, L_t, n_t, n_rows, n_passes, P_V (S_T), P_H (S_L), theta, n_parallel (optional)
       (same names as the TFHE functions)
cond : m_dot (WHOLE outer side) [kg/s], rho [kg/m3], mu [Pa s] (bulk properties of the cell)

Function (draft written by gemini-flash-lite-latest, reviewed and corrected by Claude):
  tube_bank_zukauskas_SP_DP    single phase, Zukauskas Euler-number polynomials   'Zukauskas'
  tube_bank_ishihara_two_phase_DP  two-phase, Ishihara (1980) on the Zukauskas DP  'Ishihara'
"""
import warnings
import numpy as np
from labothappy.correlations.convection.correlation_warnings import ValidityRangeWarning
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props


def tube_bank_zukauskas_SP_DP(AS, geom, cond):
    """
    Single-phase pressure drop of a bank of bare tubes in crossflow (outer side), Zukauskas
    Euler-number polynomials, inline and staggered layouts.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        - 'D_t_o' : tube outer diameter [m]
        - 'L_t' : tube length [m]
        - 'n_t' : total number of tubes [-]
        - 'n_rows' : rows per pass [-]
        - 'n_passes' : number of passes [-]
        - 'P_V' : transverse pitch S_T (normal to the flow) [m]
        - 'P_H' : longitudinal pitch S_L (parallel to the flow) [m]
        - 'theta' : tube layout angle [deg] (0 or 90 = inline; 30, 45 or 60 = staggered)
        - 'n_parallel' : optional, number of identical banks in parallel [-] (default 1)
    cond : dict
        - 'm_dot' : flow rate of the WHOLE outer side [kg/s]
        - 'rho' : bulk density of the cell [kg/m3]
        - 'mu' : bulk dynamic viscosity of the cell [Pa s]

    Returns
    -------
    dP : float
        Pressure drop of the whole bank [Pa] = N_L Eu rho V_max^2 / 2 (Eu per row). No entrance /
        exit losses, no gravity, no acceleration.

    Reference
    ---------
    Beale, S. "Tube banks, crossflow over", Thermopedia,
    doi:10.1615/AtoZ.t.tube_banks_crossflow_over (https://www.thermopedia.com/content/1211/),
    based on Zukauskas, A. (1972) Heat transfer from tubes in crossflow, Advances in Heat
    Transfer 8, 93-160, and Zukauskas & Ulinskas; also Incropera et al., Fundamentals of Heat
    and Mass Transfer (tube banks).
    Coefficients of Thermopedia Tables 1 and 2 (from the Heat Exchanger Design Handbook, 1983),
    Eq. (7) Eu = sum c_i / Re^i; checked against the saved Thermopedia page (2026-11-26): equal
    to the previous labothappy tube_bank_htc.py except c3 of in-line a = 1.25, Re 3-2e3
    (-0.286e3, the previous code had -0.286e4).

    Validity
    --------
    Inline square banks a = S_T/D = 1.25, 1.5, 2.0; staggered banks b = S_L/D = 1.25, 1.5, 2.0,
    2.5 (Table 2 is headed 'b' and indexed by b here, as in the previous code; NOTE: Figure 3 of
    the same page labels the same curves with a = S_T/D: to be clarified); Re_max ranges of each
    table (3 or 7 or 100
    up to 1e6 or 2e6). Outside: nearest table / nearest Re range and ONE ValidityRangeWarning
    (the previous code returned 0 silently); non-square inline banks: correction factor not
    applied (warning).

    Definitions
    -----------
    Re_max = rho V_max D / mu, V_max = velocity in the minimum flow section (Thermopedia /
    Incropera); Eu = c0 + c1/Re + c2/Re^2 + c3/Re^3 + c4/Re^4 per row; bulk properties.
    """
    out = []

    # 1. Read geometry and conditions
    D = geom['D_t_o']
    L_t = geom['L_t']
    n_t = geom['n_t']
    n_rows = geom['n_rows']
    n_passes = geom['n_passes']
    S_T = geom['P_V']
    S_L = geom['P_H']
    theta = geom['theta']
    n_parallel = geom.get('n_parallel', 1)

    m_dot_total = cond['m_dot']
    rho = cond['rho']
    mu = cond['mu']

    # Rows crossed by the flow, tubes per row, flow rate of one bank
    N_L = n_rows * n_passes
    tubes_per_row = n_t / N_L
    m_dot = m_dot_total / n_parallel

    # 2. Layout from theta and pitch ratios
    if theta == 0 or theta == 90:
        layout = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        layout = 'staggered'
    else:
        raise ValueError(f"tube_bank_zukauskas_SP_DP: theta = {theta!r} deg, expected 0 or 90 "
                         "(inline), 30, 45 or 60 (staggered).")
    a = S_T / D
    b = S_L / D

    # 3. Frontal velocity, maximum velocity V_max and Reynolds number
    u_front = m_dot / (rho * tubes_per_row * S_T * L_t)
    if layout == 'inline':
        V_max = S_T / (S_T - D) * u_front
    else:
        S_D = np.sqrt(S_L**2 + (S_T / 2)**2)
        if S_D < (S_T + D) / 2:
            V_max = S_T / (2 * (S_D - D)) * u_front
        else:
            V_max = S_T / (S_T - D) * u_front
    Re_max = rho * V_max * D / mu

    # 4. Tabulated pitch ratio nearest to the actual one
    if layout == 'inline':
        actual_ratio = a
        tabulated = [1.25, 1.5, 2.0]
        ratio_name = "a = S_T/D"
    else:
        actual_ratio = b
        tabulated = [1.25, 1.5, 2.0, 2.5]
        ratio_name = "b = S_L/D"
    pitch_val = tabulated[0]
    for value in tabulated:
        if abs(value - actual_ratio) < abs(pitch_val - actual_ratio):
            pitch_val = value
    if abs(actual_ratio - pitch_val) > 0.01:
        out.append(f"{ratio_name} = {actual_ratio:.3g} not tabulated, table of {pitch_val} used")
    if layout == 'inline' and abs(a - b) > 0.01:
        out.append("non-square inline bank: correction factor not applied")

    # 5. Euler number coefficients [c0, c1, c2, c3, c4] (nearest Re range if outside)
    if layout == 'inline':
        if pitch_val == 1.25:
            if Re_max <= 2e3:
                # c3 = -0.286e3 (Thermopedia Table 1; the old tube_bank_htc.py had -0.286e4,
                # which gives Eu < 0 at Re = 3)
                coeffs = [0.272, 0.207e3, 0.102e3, -0.286e3, 0]
                Re_min_table, Re_max_table = 3, 2e6
            else:
                coeffs = [0.267, 0.249e4, -0.927e7, 0.1e11, 0]
                Re_min_table, Re_max_table = 3, 2e6
        elif pitch_val == 1.5:
            if Re_max <= 2e3:
                coeffs = [0.263, 0.867e2, -0.202, 0, 0]
                Re_min_table, Re_max_table = 3, 2e6
            else:
                coeffs = [0.235, 0.197e4, -0.124e8, 0.312e11, -0.274e14]
                Re_min_table, Re_max_table = 3, 2e6
        else:  # 2.0
            if Re_max <= 800:
                coeffs = [0.188, 0.566e2, -0.646e3, 0.601e4, -0.183e5]
                Re_min_table, Re_max_table = 7, 2e6
            else:
                coeffs = [0.247, -0.595, 0.15, -0.137, 0.396]
                Re_min_table, Re_max_table = 7, 2e6
    else:  # staggered
        if pitch_val == 1.25:
            if Re_max <= 1e3:
                coeffs = [0.795, 0.247e3, 0.335e3, -0.155e4, 0.241e4]
            else:
                coeffs = [0.245, 0.339e4, -0.984e7, 0.132e11, -0.599e13]
            Re_min_table, Re_max_table = 3, 2e6
        elif pitch_val == 1.5:
            if Re_max <= 1e3:
                coeffs = [0.683, 0.111e3, -0.973e2, -0.426e3, 0.574e3]
            else:
                coeffs = [0.203, 0.248e4, -0.758e7, 0.104e11, -0.482e13]
            Re_min_table, Re_max_table = 3, 2e6
        elif pitch_val == 2.0:
            if Re_max <= 1e2:
                coeffs = [0.713, 0.448e2, -0.126e3, -0.582e3, 0]
            elif Re_max <= 1e4:
                coeffs = [0.343, 0.303e3, -0.717e5, 0.88e7, -0.38e9]
            else:
                coeffs = [0.162, 0.181e4, 0.792e8, -0.165e13, 0.872e16]
            Re_min_table, Re_max_table = 7, 1e6
        else:  # 2.5
            if Re_max <= 5e3:
                coeffs = [0.33, 0.989e2, -0.148e5, 0.192e7, -0.862e8]
            else:
                coeffs = [0.119, 0.498e4, -0.507e8, 0.251e12, -0.463e15]
            Re_min_table, Re_max_table = 1e2, 2e6
    if Re_max < Re_min_table or Re_max > Re_max_table:
        out.append(f"Re_max = {Re_max:.3g} outside {Re_min_table:g}-{Re_max_table:g} of the table")

    # 6. One warning for all the out-of-range items
    if out:
        warnings.warn(
            "tube_bank_zukauskas_SP_DP: outside the range of the correlation: "
            + ", ".join(out) + ".",
            ValidityRangeWarning
        )

    # 7. Euler number per row: Eu = c0 + c1/Re + c2/Re^2 + c3/Re^3 + c4/Re^4
    c0 = coeffs[0]
    c1 = coeffs[1]
    c2 = coeffs[2]
    c3 = coeffs[3]
    c4 = coeffs[4]
    Eu = c0 + c1 / Re_max + c2 / Re_max**2 + c3 / Re_max**3 + c4 / Re_max**4

    # 8. Pressure drop over the whole bank [Pa]
    dP = N_L * Eu * rho * V_max**2 / 2

    return dP


def tube_bank_ishihara_two_phase_DP(AS, geom, cond):
    """
    Tube-bank version of STHE_ishihara_two_phase_DP: frictional two-phase pressure drop of the
    WHOLE bank for two-phase crossflow over a bank of horizontal tubes, Ishihara, Palen & Taborek
    (1980) as given by Consolini, Robinson & Thome (2006) Eqs. 14a, 14b, 16a, 16b (C = 8).

    Parameters
    ----------
    AS : CoolProp AbstractState
        Fluid property provider (moved to saturation by _sat_props).
    geom : dict
        Geometry dictionary containing:
        - 'D_t_o' : Outer diameter of the tubes [m]
        - 'L_t' : Tube length [m]
        - 'n_t' : Total number of tubes [-]
        - 'n_rows' : Rows per pass [-]
        - 'n_passes' : Number of passes [-]
        - 'P_V' : Transverse pitch S_T (normal to the flow) [m]
        - 'P_H' : Longitudinal pitch S_L (parallel to the flow) [m]
        - 'theta' : Tube layout angle [deg] (0 or 90 = inline; 30, 45 or 60 = staggered)
        - 'n_parallel' : Optional, number of parallel identical banks [-] (default 1)
    cond : dict
        Condition dictionary containing:
        - 'm_dot' : Total mass flow rate of the whole bank [kg/s]
        - 'p' : Pressure [Pa]
        - 'x' : Vapour quality [-]

    Returns
    -------
    float
        Frictional two-phase pressure drop of the whole bank [Pa] (positive, no momentum or
        gravity terms).

    References
    ----------
    Ishihara, K., Palen, J.W., Taborek, J. (1980), Critical review of correlations
    for predicting two-phase flow pressure drop across tube banks, Heat Transfer
    Engineering 1(3), 23-33; as given by Consolini, Robinson, Thome (2006),
    Heat Transfer Engineering 27(3), doi 10.1080/01457630500453832,
    Eqs. 14a, 14b, 16a, 16b (C = 8).

    Validity Range
    --------------
    - Mass flux G >= 42 kg/(m^2 s) in the minimum flow section (emits a warning if lower).
    - Vapour quality x between 0 and 1 (clipped with a warning if out of bounds).

    Definitions
    -----------
    G = mass flux in the minimum flow section of the bank. Single-phase pressure drops are
    computed using tube_bank_zukauskas_SP_DP. Bulk properties from saturation.
    """
    # Step 1: Quality check and clipping
    x = cond['x']
    if x < 0.0 or x > 1.0:
        warnings.warn(
            f"tube_bank_ishihara_two_phase_DP: vapour quality x = {x} is out of the [0, 1] range. "
            "It has been clipped.", ValidityRangeWarning
        )
        x = np.clip(x, 0.0, 1.0)

    # Step 2: Read geometry and compute layout angle
    D = geom['D_t_o']
    L_t = geom['L_t']
    n_t = geom['n_t']
    n_rows = geom['n_rows']
    n_passes = geom['n_passes']
    S_T = geom['P_V']
    S_L = geom['P_H']
    theta = geom['theta']
    n_parallel = geom.get('n_parallel', 1)

    if theta == 0 or theta == 90:
        layout = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        layout = 'staggered'
    else:
        raise ValueError(
            f"tube_bank_ishihara_two_phase_DP: theta = {theta!r} deg, expected 0 or 90 "
            "(inline), 30, 45 or 60 (staggered)."
        )

    # Step 3: Mass flux in the minimum flow section of the bank
    N_L = n_rows * n_passes
    tubes_per_row = n_t / N_L
    m_dot_bank = cond['m_dot'] / n_parallel
    G_front = m_dot_bank / (tubes_per_row * S_T * L_t)

    if layout == 'inline':
        G = G_front * S_T / (S_T - D)
    else:
        S_D = np.sqrt(S_L**2 + (S_T / 2)**2)
        if S_D < (S_T + D) / 2:
            G = G_front * S_T / (2 * (S_D - D))
        else:
            G = G_front * S_T / (S_T - D)

    # Warning if mass flux is below the recommended threshold
    if G < 42.0:
        warnings.warn(
            f"tube_bank_ishihara_two_phase_DP: mass flux G = {G:.2f} kg/(m^2 s) is below 42. "
            "The Ishihara method may underpredict pressure drop at low mass flux.",
            ValidityRangeWarning
        )

    # Step 4: Saturated properties at pressure p
    sat_props = _sat_props(AS, cond['p'], cond)
    rho_l = sat_props['rho_l']
    rho_v = sat_props['rho_v']
    mu_l = sat_props['mu_l']
    mu_v = sat_props['mu_v']

    # Step 5: Special cases for single-phase extremes
    if x <= 0.0:
        cond_l = dict(cond)
        cond_l['m_dot'] = cond['m_dot']
        cond_l['mu'] = mu_l
        cond_l['rho'] = rho_l
        return tube_bank_zukauskas_SP_DP(AS, geom, cond_l)

    if x >= 1.0:
        cond_g = dict(cond)
        cond_g['m_dot'] = cond['m_dot']
        cond_g['mu'] = mu_v
        cond_g['rho'] = rho_v
        return tube_bank_zukauskas_SP_DP(AS, geom, cond_g)

    # Step 6: Single-phase pressure drops using tube_bank_zukauskas_SP_DP
    # Liquid-only single-phase pressure drop (Consolini Eqs. 15a/15b)
    cond_l = dict(cond)
    cond_l['m_dot'] = cond['m_dot'] * (1.0 - x)
    cond_l['mu'] = mu_l
    cond_l['rho'] = rho_l
    DP_l = tube_bank_zukauskas_SP_DP(AS, geom, cond_l)

    # Vapour-only single-phase pressure drop
    cond_g = dict(cond)
    cond_g['m_dot'] = cond['m_dot'] * x
    cond_g['mu'] = mu_v
    cond_g['rho'] = rho_v
    DP_g = tube_bank_zukauskas_SP_DP(AS, geom, cond_g)

    # Step 7: Martinelli parameter for turbulent-turbulent flow (Consolini et al.)
    term_x = (1.0 - x) / x
    term_rho = rho_v / rho_l
    term_mu = mu_l / mu_v
    X_tt = (term_x**0.9) * (term_rho**0.5) * (term_mu**0.1)

    # Step 8: Liquid Reynolds number on tube outer diameter
    Re_l = G * D * (1.0 - x) / mu_l

    # Step 9: Two-phase multiplier and total frictional pressure drop
    if Re_l > 2000.0:
        # Re_l > 2000: liquid multiplier (Consolini Eqs. 16a and 14a)
        Phi_l2 = 1.0 + (8.0 / X_tt) + (1.0 / (X_tt**2))
        DP = Phi_l2 * DP_l
    else:
        # Re_l <= 2000: vapour multiplier (Consolini Eqs. 16b and 14b)
        Phi_g2 = 1.0 + (8.0 * X_tt) + (X_tt**2)
        DP = Phi_g2 * DP_g

    return float(DP)
