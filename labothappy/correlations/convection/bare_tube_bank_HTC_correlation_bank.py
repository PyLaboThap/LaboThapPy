# -*- coding: utf-8 -*-
"""
Heat transfer on the OUTER side of banks of bare (plain) tubes in crossflow, (AS, geom, cond)
formalism. Replaces the old convection/tube_bank_htc.py (2026-11-26).

Same geometry names as the TFHE functions:
geom : D_t_o (tube outer diameter) [m], L_t (tube length) [m], n_t (total number of tubes) [-],
       n_rows (rows per pass) and n_passes [-] (rows crossed by the flow N_L = n_rows * n_passes),
       P_V (transverse pitch S_T, normal to the flow) [m], P_H (longitudinal pitch S_L) [m],
       theta (layout angle [deg]: 90 or 0 = inline, 30 / 45 / 60 = staggered),
       n_parallel (optional, default 1)
cond : m_dot (flow rate of the WHOLE outer side) [kg/s], bulk properties rho, mu, k, Pr;
       Pr_w (or T_w and p) for the wall Prandtl number

Functions (drafts written by gemini-flash-lite-latest, reviewed and corrected by Claude):
  tube_bank_zukauskas_SP_HTC        single phase, Zukauskas (Thermopedia / Incropera),
                                    registered as 'Zukauskas'
  tube_bank_esdu_73031_SP_HTC       single phase, ESDU 73031 (as in ht / Hewitt 1994),
                                    registered as 'ESDU_73031'
  tube_nusselt_film_condensing_HTC  film condensation on ONE horizontal tube (Nusselt),
                                    registered as 'Nusselt_single_tube'
  tube_bank_shah_boiling_HTC        boiling on the bank, Shah (2017), registered as 'Shah'
  tube_bank_rose_kern_condensing_HTC condensation on the bank, Rose / Nusselt x Kern inundation,
                                    registered as 'Rose_Kern'
  tube_lienhard_eichhorn_CHF        CRITICAL heat flux in crossflow boiling
                                    (not an HTC, not registered)
"""
import warnings
import numpy as np
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props, _need_sigma
from labothappy.correlations.convection.BPHE_HTC_correlation_bank import BPHE_cooper_boiling_HTC
from labothappy.correlations.convection.correlation_warnings import (
    AssumedParameterWarning, ValidityRangeWarning)


def tube_bank_zukauskas_SP_HTC(AS, geom, cond):
    """
    Single-phase heat transfer coefficient of a bank of bare tubes in crossflow (outer side),
    Zukauskas method, inline and staggered layouts.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Used only to compute the wall Prandtl number from cond['T_w'] and cond['p'] when
        cond['Pr_w'] is not given (AS is moved to (p, T_w)). May be None if Pr_w is given.
    geom : dict
        - 'D_t_o' : tube outer diameter [m]
        - 'L_t' : tube length [m]
        - 'n_t' : total number of tubes [-]
        - 'n_rows' : rows per pass [-]
        - 'n_passes' : number of passes [-]
        - 'P_V' : transverse pitch S_T (normal to the flow) [m]
        - 'P_H' : longitudinal pitch S_L (parallel to the flow) [m]
        - 'theta' : tube layout angle [deg] (90 or 0 = inline; 30, 45 or 60 = staggered)
        - 'n_parallel' : optional, number of identical banks in parallel [-] (default 1)
    cond : dict
        - 'm_dot' : flow rate of the whole outer side [kg/s]
        - 'rho' : bulk density [kg/m3]
        - 'mu' : bulk dynamic viscosity [Pa s]
        - 'k' : bulk thermal conductivity [W/m/K]
        - 'Pr' : bulk Prandtl number [-]
        - 'Pr_w' : optional, Prandtl number at the wall temperature [-]
        - 'T_w', 'p' : optional, wall temperature [K] and pressure [Pa] (used when Pr_w is missing)

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K], referred to the outer area of the bare tubes.

    Reference
    ---------
    Beale, S. "Tube banks, crossflow over", Thermopedia,
    doi:10.1615/AtoZ.t.tube_banks_crossflow_over (https://www.thermopedia.com/content/1211/),
    based on Zukauskas, A. (1972) Heat transfer from tubes in crossflow, Advances in Heat
    Transfer 8, 93-160, and Zukauskas & Ulinskas; also Incropera et al., Fundamentals of Heat
    and Mass Transfer (tube banks). 100 <= Re_max < 1000: single isolated cylinder,
    Churchill & Bernstein (1977).
    Coefficients as in the previous labothappy tube_bank_htc.py, checked against
    Thermopedia by B. Chaudoir (2026-11-26).

    Validity
    --------
    10 <= Re_max <= 2e6, 0.5 <= Pr <= 500, 1.25 <= S_T/D <= 2.5 (as in the previous code);
    inline banks with S_T/S_L <= 0.7: correlation not established. One ValidityRangeWarning
    outside, the value is computed anyway.

    Definitions
    -----------
    Re_max = rho V_max D / mu, V_max = velocity in the minimum flow section (Thermopedia /
    Incropera definitions, the same V_max for the regime choice and for Nu); bulk properties of
    the cell (mean temperature), Pr_w at the wall temperature; row correction C_2 for N_L < 20.
    Corrections vs the previous code: one V_max for the staggered regime choice and Nu; no crash
    for inline S_T/S_L <= 0.7 (warning instead); C_2 defined for 16 < N_L < 20 (0.99).
    """
    # 1. Read geometry and operating conditions
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
    k = cond['k']
    Pr = cond['Pr']

    # 2. Layout from theta
    if theta == 0 or theta == 90:
        layout = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        layout = 'staggered'
    else:
        raise ValueError(f"tube_bank_zukauskas_SP_HTC: theta = {theta} deg, expected 0 or 90 "
                         "(inline), 30, 45 or 60 (staggered).")

    # 3. Number of rows crossed by the flow, tubes per row, pitch ratio
    N_L = n_rows * n_passes
    tubes_per_row = n_t / N_L
    ratio_T_L = S_T / S_L

    # 4. Frontal (approach) velocity
    m_dot_path = m_dot_total / n_parallel
    u_front = m_dot_path / (rho * tubes_per_row * S_T * L_t)

    # 5. Maximum velocity V_max (Thermopedia / Incropera)
    if layout == 'inline':
        V_max = S_T / (S_T - D) * u_front
    else:
        S_D = np.sqrt(S_L**2 + (S_T / 2.0)**2)
        if S_D < (S_T + D) / 2.0:
            V_max = S_T / (2.0 * (S_D - D)) * u_front
        else:
            V_max = S_T / (S_T - D) * u_front

    # 6. Reynolds number based on V_max and the tube outer diameter
    Re_max = rho * V_max * D / mu

    # 7. Validity range: one warning
    out = []
    if not (10 <= Re_max <= 2e6):
        out.append(f"Re_max = {Re_max:.1f} (10-2e6)")
    if not (0.5 <= Pr <= 500):
        out.append(f"Pr = {Pr:.2f} (0.5-500)")
    if not (1.25 <= S_T / D <= 2.5):
        out.append(f"S_T/D = {S_T / D:.2f} (1.25-2.5)")
    if layout == 'inline' and ratio_T_L <= 0.7:
        out.append(f"S_T/S_L = {ratio_T_L:.2f} <= 0.7 (inline banks with S_T/S_L < 0.7 are "
                   "inefficient, correlation not established)")
    if out:
        warnings.warn(
            "tube_bank_zukauskas_SP_HTC: outside the range of the correlation: "
            + ", ".join(out) + ".",
            ValidityRangeWarning
        )

    # 8. Wall Prandtl number
    Pr_w = None
    if 'Pr_w' in cond:
        Pr_w = cond['Pr_w']
    elif 'T_w' in cond and 'p' in cond:
        try:
            import CoolProp.CoolProp as CP
            AS.update(CP.PT_INPUTS, cond['p'], cond['T_w'])
            Pr_w = AS.Prandtl()
        except Exception:
            Pr_w = None
    if Pr_w is None:
        Pr_w = Pr
        warnings.warn(
            "tube_bank_zukauskas_SP_HTC: wall Prandtl number unknown, (Pr/Pr_w)^0.25 = 1.",
            AssumedParameterWarning
        )

    # 9. Nusselt number
    if 100 <= Re_max < 1000:
        # Single isolated cylinder, Churchill & Bernstein (1977)
        term_re_cb = (Re_max / 282000.0)**(5.0 / 8.0)
        term_bracket_cb = 1.0 + term_re_cb
        term_power_cb = term_bracket_cb**(4.0 / 5.0)
        term_pr_cb = (0.4 / Pr)**(2.0 / 3.0)
        term_denom_cb = (1.0 + term_pr_cb)**(1.0 / 4.0)
        term_num_cb = 0.62 * np.sqrt(Re_max) * (Pr**(1.0 / 3.0))
        Nu = 0.3 + (term_num_cb / term_denom_cb) * term_power_cb
    else:
        # Zukauskas constants C and m (Thermopedia / Incropera tables)
        if layout == 'inline':
            if Re_max < 100:
                C = 0.80
                m = 0.40
            elif Re_max < 2e5:
                C = 0.27
                m = 0.63
            else:
                C = 0.021
                m = 0.84
        else:
            if Re_max < 100:
                C = 0.90
                m = 0.40
            elif Re_max < 2e5:
                if ratio_T_L < 2.0:
                    C = 0.35 * (ratio_T_L**(1.0 / 5.0))
                    m = 0.60
                else:
                    C = 0.40
                    m = 0.60
            else:
                C = 0.021
                m = 0.84

        # Row correction factor C_2 (1 for 20 rows or more)
        N_L_values = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 13.0, 16.0]
        if layout == 'inline':
            C_2_values = [0.70, 0.80, 0.86, 0.90, 0.92, 0.95, 0.97, 0.98, 0.99]
        else:
            C_2_values = [0.64, 0.76, 0.84, 0.89, 0.92, 0.95, 0.97, 0.98, 0.99]
        if N_L >= 20.0:
            C_2 = 1.0
        else:
            C_2 = np.interp(N_L, N_L_values, C_2_values)

        # Zukauskas equation: Nu = C_2 C Re_max^m Pr^0.36 (Pr/Pr_w)^0.25
        term_re_zuk = Re_max**m
        term_pr_zuk = Pr**0.36
        term_prop_zuk = (Pr / Pr_w)**0.25
        Nu = C_2 * C * term_re_zuk * term_pr_zuk * term_prop_zuk

    # 10. Heat transfer coefficient referred to the outer area of the bare tubes
    h = Nu * k / D

    return h


def tube_nusselt_film_condensing_HTC(AS, geom, cond):
    """
    Laminar film condensation on the outside of ONE horizontal tube (Nusselt). No row
    (inundation) correction: for tube bundles use the STHE correlation 'Rose_Kern'.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Fluid state (moved by the saturated property calls).
    geom : dict
        - 'D_t_o' : tube outer diameter [m]
    cond : dict
        - 'p' : saturation pressure [Pa]
        - 'T_w' : tube outer wall temperature [K]
        - optional 'mu_l', 'k_l' overrides of the liquid transport properties (see _sat_props)

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K], referred to the outer area of the tube.

    Reference
    ---------
    Nusselt, W. (1916) Die Oberflaechenkondensation des Wasserdampfes, Z. VDI 60, 541-546 and
    569-575 (pages not checked); form for a horizontal tube with the modified latent heat as in
    Incropera et al., Fundamentals of Heat and Mass Transfer (film condensation on a horizontal
    tube, C = 0.729; h'_fg = h_fg (1 + 0.68 Ja)).
    Replaces ext_tube_film_condens of the previous tube_bank_htc.py, which returned the Nusselt
    number instead of h (k_l / D missing) and used the uncorrected latent heat.

    Validity and definitions
    ------------------------
    Laminar film, quiescent (or slow) vapour, single tube. No numerical range check.
    Properties: saturated liquid and vapour at p (the textbook evaluates the liquid properties at
    the film temperature (T_sat + T_w) / 2). T_sat = saturation temperature at p.
    """
    # 1. Inputs: D [m], p_sat [Pa], T_w [K]
    D = geom['D_t_o']
    p_sat = cond['p']
    T_w = cond['T_w']

    # 2. Saturated properties at p_sat
    s = _sat_props(AS, p_sat, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    k_l = s['k_l']
    cp_l = s['cp_l']
    h_lv = s['h_lv']
    T_sat = s['T_sat']

    # 3. Wall subcooling: the wall must be below the saturation temperature
    dT = T_sat - T_w
    if dT <= 0:
        raise ValueError(
            f"tube_nusselt_film_condensing_HTC: no condensation, the wall (T_w = {T_w:.2f} K) "
            f"is not below the saturation temperature (T_sat = {T_sat:.2f} K)."
        )

    # 4. Jakob number and modified latent heat h'_fg = h_fg (1 + 0.68 Ja)
    Ja = cp_l * dT / h_lv
    h_lv_mod = h_lv * (1.0 + 0.68 * Ja)

    # 5. Nusselt correlation for a horizontal tube, C = 0.729
    g = 9.81
    term = g * rho_l * (rho_l - rho_v) * (k_l ** 3) * h_lv_mod / (mu_l * dT * D)
    h = 0.729 * (term ** 0.25)

    # 6. Heat transfer coefficient [W/m2/K] (outer area of the tube)
    return h


def tube_lienhard_eichhorn_CHF(AS, geom, cond):
    """
    CRITICAL heat flux (maximum heat flux) of forced-convection boiling of a liquid in
    crossflow over ONE cylinder, Lienhard & Eichhorn (1976).
    This is NOT a heat transfer coefficient and it is not registered in the dispatcher: it
    gives the heat flux that must not be exceeded (burnout).

    Parameters
    ----------
    AS : CoolProp AbstractState
        Fluid state (moved by the saturated property calls).
    geom : dict
        - 'D_t_o' : tube outer diameter [m]
    cond : dict
        - 'p' : saturation pressure [Pa]
        - 'V' : approach velocity of the liquid crossflow [m/s]
        - optional 'sigma' : surface tension [N/m] if CoolProp has none

    Returns
    -------
    q_max : float
        Critical heat flux [W/m2].

    Reference
    ---------
    Lienhard, J.H., Eichhorn, R. (1976) Peak boiling heat flux on cylinders in a cross flow,
    Int. J. Heat Mass Transfer 19, 1135-1142 (pages not checked); as given in Incropera et al.,
    Fundamentals of Heat and Mass Transfer (forced convection boiling, external flow).
    Equations as transcribed in the previous labothappy tube_bank_htc.py (ext_tube_conv_boil),
    which divided q_max by (T_tube - T_sat) and returned it as a heat transfer coefficient.
    TODO: check the low / high velocity criterion against the textbook.

    Definitions
    -----------
    We_D = rho_v V^2 D / sigma; saturated properties at p.
    """
    # 1. Inputs
    D = geom['D_t_o']
    p_sat = cond['p']
    V = cond['V']

    # 2. Saturated properties and surface tension
    s = _sat_props(AS, p_sat, cond)
    _need_sigma(AS, s)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    h_lv = s['h_lv']
    sigma = s['sigma']

    # 3. Weber number based on the vapour density and the approach velocity
    We_D = rho_v * (V ** 2) * D / sigma

    # 4. Low-velocity coefficient
    coef_low = (1.0 / np.pi) * (1.0 + (4.0 / We_D) ** (1.0 / 3.0))

    # 5. Regime and critical heat flux
    limit_value = (0.275 / np.pi) * (rho_l / rho_v) ** 0.5 + 1.0
    if coef_low >= limit_value:
        # Low velocity region
        q_max = rho_v * V * h_lv * coef_low
    else:
        # High velocity region
        term_1 = (rho_l / rho_v) ** (3.0 / 4.0) / (169.0 * np.pi)
        term_2 = (rho_l / rho_v) ** 0.5 / (19.2 * np.pi * We_D ** (1.0 / 3.0))
        q_max = rho_v * V * h_lv * (term_1 + term_2)

    # 6. Critical heat flux [W/m2]
    return q_max


def tube_bank_esdu_73031_SP_HTC(AS, geom, cond):
    """
    Single-phase heat transfer coefficient of a bank of bare tubes in crossflow (outer side),
    ESDU 73031 method, inline and staggered layouts.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Used only to compute the wall Prandtl number from cond['T_w'] and cond['p'] when
        cond['Pr_w'] is not given (AS is moved to (p, T_w)). May be None if Pr_w is given.
    geom : dict
        - 'D_t_o' : tube outer diameter [m]
        - 'L_t' : tube length [m]
        - 'n_t' : total number of tubes [-]
        - 'n_rows' : rows per pass [-]
        - 'n_passes' : number of passes [-]
        - 'P_V' : transverse pitch S_T (normal to the flow) [m]
        - 'P_H' : longitudinal pitch S_L (parallel to the flow) [m]
        - 'theta' : tube layout angle [deg] (90 or 0 = inline; 30, 45 or 60 = staggered)
        - 'n_parallel' : optional, number of identical banks in parallel [-] (default 1)
    cond : dict
        - 'm_dot' : flow rate of the whole outer side [kg/s]
        - 'rho' : bulk density [kg/m3]
        - 'mu' : bulk dynamic viscosity [Pa s]
        - 'k' : bulk thermal conductivity [W/m/K]
        - 'Pr' : bulk Prandtl number [-]
        - 'Pr_w' : optional, Prandtl number at the wall temperature [-]
        - 'T_w', 'p' : optional, wall temperature [K] and pressure [Pa] (used when Pr_w is missing)

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K], referred to the outer area of the bare tubes.

    Reference
    ---------
    "Convective heat transfer during crossflow of fluids over plain tube banks", ESDU 73031 (1973),
    as shown in Hewitt, G.F., Shires, G.L., Bott, T.R. (1994) Process Heat Transfer, CRC Press.
    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.conv_tube_bank.Nu_ESDU_73031, which cites the references above (layout taken from
    geom['theta'] here, not guessed from the pitch ratio as in ht).

    Validity
    --------
    10 <= Re_max <= 2e6; S_T/D 1.2-4 (inline) or 1-4 (staggered); N_L >= 3 (for N_L < 3 the row
    factor of 3 rows is used).
    Accuracy about 15 %; F1 exponent 0.26 for heating and cooling, liquids and gases (as ht).
    One ValidityRangeWarning outside, the value is computed anyway.

    Definitions
    -----------
    Re_max = rho V_max D / mu, V_max = velocity in the minimum flow section (Thermopedia /
    Incropera definitions, the same V_max for the regime choice and for Nu); bulk properties of
    the cell (mean temperature), Pr_w at the wall temperature; straight banks only (no inclination
    factor).
    """
    # 1. Read geometry and operating conditions
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
    k = cond['k']
    Pr = cond['Pr']

    # 2. Layout from theta
    if theta == 0 or theta == 90:
        layout = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        layout = 'staggered'
    else:
        raise ValueError(f"tube_bank_esdu_73031_SP_HTC: theta = {theta} deg, expected 0 or 90 "
                         "(inline), 30, 45 or 60 (staggered).")

    # 3. Number of rows crossed by the flow, tubes per row
    N_L = n_rows * n_passes
    tubes_per_row = n_t / N_L

    # 4. Frontal (approach) velocity
    m_dot_path = m_dot_total / n_parallel
    u_front = m_dot_path / (rho * tubes_per_row * S_T * L_t)

    # 5. Maximum velocity V_max (Thermopedia / Incropera definitions)
    if layout == 'inline':
        V_max = S_T / (S_T - D) * u_front
    else:
        S_D = np.sqrt(S_L**2 + (S_T / 2.0)**2)
        if S_D < (S_T + D) / 2.0:
            V_max = S_T / (2.0 * (S_D - D)) * u_front
        else:
            V_max = S_T / (S_T - D) * u_front

    # 6. Reynolds number based on V_max and the tube outer diameter
    Re_max = rho * V_max * D / mu

    # 7. Wall Prandtl number
    Pr_w = None
    if 'Pr_w' in cond:
        Pr_w = cond['Pr_w']
    elif 'T_w' in cond and 'p' in cond:
        try:
            import CoolProp.CoolProp as CP
            AS.update(CP.PT_INPUTS, cond['p'], cond['T_w'])
            Pr_w = AS.Prandtl()
        except Exception:
            Pr_w = None
    if Pr_w is None:
        Pr_w = Pr
        warnings.warn(
            "tube_bank_esdu_73031_SP_HTC: wall Prandtl number unknown, Pr_w = Pr assumed.",
            AssumedParameterWarning
        )

    # 8. Validity range: one warning list
    out = []
    if not (10 <= Re_max <= 2e6):
        out.append(f"Re_max = {Re_max:.1f} (10-2e6)")
    if layout == 'inline':
        if not (1.2 <= S_T / D <= 4.0):
            out.append(f"S_T/D = {S_T / D:.2f} (1.2-4 inline)")
    else:
        if not (1.0 <= S_T / D <= 4.0):
            out.append(f"S_T/D = {S_T / D:.2f} (1-4 staggered)")
    if N_L < 3:
        out.append(f"N_L = {N_L} < 3: row factor of 3 rows used")
    if out:
        warnings.warn(
            "tube_bank_esdu_73031_SP_HTC: outside the range of the correlation: "
            + ", ".join(out) + ".",
            ValidityRangeWarning
        )

    # 9. Constants a and m depending on layout and Re_max
    if layout == 'inline':
        if Re_max <= 300.0:
            a = 0.742
            m = 0.431
        elif Re_max <= 2e5:
            a = 0.211
            m = 0.651
        else:
            a = 0.116
            m = 0.700
    else:
        if Re_max <= 300.0:
            a = 1.309
            m = 0.360
        elif Re_max <= 2e5:
            a = 0.273
            m = 0.635
        else:
            a = 0.124
            m = 0.700

    # 10. Property correction factor F1
    F1 = (Pr / Pr_w)**0.26

    # 11. Tube row correction factor F2
    if N_L >= 10.0:
        F2 = 1.0
    else:
        if layout == 'inline':
            F2_values = [0.8479, 0.8957, 0.9306, 0.9551, 0.9724, 0.9839, 0.9902]
        else:
            F2_values = [0.8593, 0.8984, 0.9268, 0.9482, 0.9650, 0.9777, 0.9868]

        if N_L <= 2.0:
            F2 = F2_values[0]
        else:
            index = int(N_L) - 3
            F2 = F2_values[index]

    # 12. Nusselt number and heat transfer coefficient
    term_re = Re_max**m
    term_pr = Pr**0.34
    Nu = a * term_re * term_pr * F1 * F2
    h = Nu * k / D

    return h


def tube_bank_shah_boiling_HTC(AS, geom, cond):
    """
    Local boiling heat transfer coefficient on the outer surface of tubes in a horizontal
    plain-tube bank (crossflow boiling) [W/m2/K], tube-bank version of
    STHE_shah_bundle_boiling_HTC (same equations; mass flux in the minimum flow section of the
    bank instead of Kern's shell cross-flow area).

    Parameters
    ----------
    AS : CoolProp abstract state object
        Fluid property provider.
    geom : dict
        - 'D_t_o' : tube outer diameter [m]
        - 'L_t' : tube length [m]
        - 'n_t' : total number of tubes [-]
        - 'n_rows' : rows per pass [-]
        - 'n_passes' : number of passes [-]
        - 'P_V' : transverse pitch S_T (normal to the flow) [m]
        - 'P_H' : longitudinal pitch S_L (parallel to the flow) [m]
        - 'theta' : tube layout angle [deg] (90 or 0 = inline; 30, 45 or 60 = staggered)
        - 'n_parallel' : optional, number of parallel banks [-] (default 1)
    cond : dict
        - 'm_dot' : total mass flow rate on the outer side [kg/s]
        - 'p' : saturation pressure [Pa]
        - 'x' : local (cell mean) vapour quality [-]
        - 'q' : local heat flux on the outer tube surface [W/m^2]

    Returns
    -------
    h : float
        Local boiling heat transfer coefficient on the nominal outer area [W/(m^2 K)].

    Reference
    ---------
    Shah, M.M. (2017), A correlation for heat transfer during boiling on bundles of horizontal
    plain and enhanced tubes, International Journal of Refrigeration 78, doi 10.1016/j.ijrefrig.
    2017.03.010, Eqs. 13-20 (modified correlation, section 3.2).
    Cooper, M.G. (1984), Advances in Heat Transfer 16, 157-239 (h_cooper).

    Validity Range (Shah 2017 Table 1)
    ---------------------------------
    p_red : 0.005 - 0.2866
    G : 0.15 - 1391 kg/(m^2 s)
    q : 1e3 - 1e6 W/m^2
    x : 0 - 0.98
    D : 3e-3 - 25.4e-3 m
    S_T / D : 1.17 - 2.0
    Re_L : 7 - 52117
    Database: 2173 data points, 12 fluids (water, ammonia, R-11, R-12, R-22, R-113,
    R-134a, R-236fa, R-410A, R-507A, isobutane, n-pentane; Shah 2017 Table 1).
    Mean absolute deviation 15.2 %.

    Definitions
    -----------
    Re_L = G * D / mu_l (all mass as liquid, tube OD).
    Bo and Fr are defined with G (mass flux in the minimum flow section of the tube bank).
    All fluid properties are evaluated at saturation temperature corresponding to pressure p.
    q is the heat flux on the outer tube surface.
    x is the local (cell mean) vapour quality.
    h_cooper is the Cooper (1984) original pool boiling heat transfer coefficient using M^-0.5
    and roughness term through the existing library function BPHE_cooper_boiling_HTC.

    """
    # 1. Read geometry parameters
    D = geom['D_t_o']
    L_t = geom['L_t']
    n_t = geom['n_t']
    n_rows = geom['n_rows']
    n_passes = geom['n_passes']
    S_T = geom['P_V']
    S_L = geom['P_H']
    theta = geom['theta']
    n_parallel = geom.get('n_parallel', 1)

    # 2. Layout verification from theta
    if theta == 0 or theta == 90:
        layout = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        layout = 'staggered'
    else:
        raise ValueError(
            f"tube_bank_shah_boiling_HTC: theta = {theta} deg, expected 0 or 90 "
            "(inline), 30, 45 or 60 (staggered)."
        )

    # 3. Number of rows crossed and tubes per row
    N_L = n_rows * n_passes
    tubes_per_row = n_t / N_L

    # 4. Frontal mass flux
    m_dot_bank = cond['m_dot'] / n_parallel
    G_front = m_dot_bank / (tubes_per_row * S_T * L_t)

    # 5. Mass flux in the minimum flow section of the bank
    if layout == 'inline':
        G = G_front * S_T / (S_T - D)
    else:
        S_D = np.sqrt(S_L**2 + (S_T / 2.0)**2)
        if S_D < (S_T + D) / 2.0:
            G = G_front * S_T / (2.0 * (S_D - D))
        else:
            G = G_front * S_T / (S_T - D)

    # 6. Read operating conditions
    x = cond['x']
    g = 9.81

    # Check quality bounds
    if x < 0.0 or x > 1.0:
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: quality x = {x} outside [0, 1]. Clipped.",
            ValidityRangeWarning
        )
        x = np.clip(x, 0.0, 1.0)

    # 7. Validate heat flux
    q = abs(cond['q'])
    if not np.isfinite(q) or q <= 0.0:
        raise ValueError(
            f"tube_bank_shah_boiling_HTC: heat flux q = {q} must be finite and > 0."
        )

    # 8. Compute Cooper pool boiling heat transfer coefficient
    cond_cooper = dict(cond)
    cond_cooper['q'] = q
    h_cooper = BPHE_cooper_boiling_HTC(AS, geom, cond_cooper)

    # 9. Obtain saturation properties and critical pressure
    p_critical = AS.p_critical()
    p_red = cond['p'] / p_critical
    sat_props = _sat_props(AS, cond['p'], cond)
    rho_l = sat_props['rho_l']
    mu_l = sat_props['mu_l']
    k_l = sat_props['k_l']
    Pr_l = sat_props['Pr_l']
    h_lv = sat_props['h_lv']

    # 10. Dimensionless numbers: Boiling number, Froude number, and modified Y parameter (Eq. 13)
    Bo = q / (G * h_lv)
    Fr = (G**2) / ((rho_l**2) * g * D)
    Y = Bo * (Fr**0.3)

    # 11. Plain tube factor
    F_pb = 1.0

    # Liquid-only Reynolds number (Eq. 20), also used in validity check
    Re_L = G * D / mu_l

    # 12. High-Y branch (Eq. 14)
    if Y > 0.008:
        h = F_pb * h_cooper

    # 13. Low-Y branch (Y <= 0.008): maximum of Eq. 14, Eq. 20, Eq. 15, and Eq. 16
    else:
        # Eq. 14: h_14
        h_14 = F_pb * h_cooper

        # Liquid-only coefficient (Eq. 20)
        h_LT = 0.21 * (Re_L**0.62) * (Pr_l**0.4) * k_l / D

        # Eq. 15: phi0 components (Eqs. 17, 18, 19)
        phi_17 = 443.0 * (Bo**0.65) * F_pb
        phi_18 = 31.0 * (Bo**0.33) * F_pb
        phi_19 = 1.0
        phi0 = max(phi_17, phi_18, phi_19)
        h_15 = phi0 * h_LT

        # Eq. 16 (only if x > 0)
        if x > 0.0:
            Z = (((1.0 - x) / x)**0.8) * (p_red**0.4)
            phi = 2.3 / ((Z**0.08) * (Fr**0.22))
            h_16 = phi * h_LT
            h = max(h_14, h_15, h_16)
        else:
            h = max(h_14, h_15)

    # 14. Validity checks (Shah 2017 Table 1)
    pitch_ratio = S_T / D
    if not (0.005 <= p_red <= 0.2866):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: reduced pressure p_red = {p_red:.4f} "
            "is outside validity range [0.005, 0.2866].",
            ValidityRangeWarning
        )
    if not (0.15 <= G <= 1391.0):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: mass flux G = {G:.2f} kg/(m^2 s) "
            "is outside validity range [0.15, 1391].",
            ValidityRangeWarning
        )
    if not (1e3 <= q <= 1e6):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: heat flux q = {q:.1f} W/m^2 "
            "is outside validity range [1e3, 1e6].",
            ValidityRangeWarning
        )
    if not (0.0 <= x <= 0.98):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: quality x = {x:.3f} "
            "is outside validity range [0, 0.98].",
            ValidityRangeWarning
        )
    if not (3e-3 <= D <= 25.4e-3):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: tube OD D = {D:.5f} m "
            "is outside validity range [3e-3, 25.4e-3].",
            ValidityRangeWarning
        )
    if not (1.17 <= pitch_ratio <= 2.0):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: pitch ratio S_T/D = {pitch_ratio:.2f} "
            "is outside validity range [1.17, 2.0].",
            ValidityRangeWarning
        )
    if not (7.0 <= Re_L <= 52117.0):
        warnings.warn(
            f"tube_bank_shah_boiling_HTC: liquid Reynolds number Re_L = {Re_L:.1f} "
            "is outside validity range [7, 52117].",
            ValidityRangeWarning
        )

    # Return final heat transfer coefficient
    return h


def tube_bank_rose_kern_condensing_HTC(AS, geom, cond):
    """
    Mean condensation heat transfer coefficient on the outer surface of a bank of horizontal tubes
    (tube-bank version of STHE_rose_kern_condensing_HTC, same equations; mass flux in the minimum
    flow section of the bank instead of Kern's shell cross-flow area).

    Parameters
    ----------
    AS : CoolProp AbstractState
        CoolProp abstract state object for property evaluation.
    geom : dict
        Geometry dictionary containing:
        - D_t_o : Tube outside diameter [m]
        - L_t : Tube length [m]
        - n_t : Total number of tubes [-]
        - n_rows : Rows per pass [-]
        - n_passes : Number of passes [-]
        - P_V : Transverse pitch S_T (normal to the flow) [m]
        - P_H : Longitudinal pitch S_L (parallel to the flow) [m]
        - theta : Tube layout angle [deg] (90 or 0 = inline; 30, 45 or 60 = staggered)
        - n_parallel : Optional, number of parallel paths [-] (default 1)
    cond : dict
        Condition dictionary containing:
        - m_dot : Flow rate of the whole outer side [kg/s]
        - p : Pressure [Pa]
        - x : Vapour quality [-]
        - T_w : Wall temperature [K]

    Returns
    -------
    h : float
        Mean condensation heat transfer coefficient [W/(m^2 K)] based on the
        outer tube surface area (pi * D * L), bundle-averaged (inundation included).

    Reference
    ---------
    - Bonneau, Josset, Melot, Auvity (2019), Synthesis of heat transfer correlations
      for condensation outside horizontal tube bundles, Nuclear Engineering and Design 349,
      doi 10.1016/j.nucengdes.2019.04.005, section 4 (recommendations); equation numbers
      below are those of Bonneau 2019.
    - Rose, J.W. (1984), Int. J. Heat Mass Transfer 27(1), 39-47.
    - Fujii et al. (1979), as given by Bonneau 2019 Eq. 70.
    - Nusselt, W. (1916), VDI-Z. 60, 541-546, 569-575.
    - Kern, D.Q. (1958), inundation factor N^(-1/6), as given by Bonneau 2019 Eq. 110.

    Validity
    --------
    - Fluids: steam and refrigerants (Bonneau 2019 database).
    - Rose single-tube correlations: F < 1e3 (gravity-controlled domain for F > 1e3).

    Definitions
    -----------
    - Nu and Re evaluated on the tube outer diameter (D_t_o = D).
    - Liquid properties evaluated at saturation temperature via _sat_props.
    - dT = T_sat - T_w.
    - U = G * x / rho_v, where G is the mass flux in the minimum flow section of the bank
      (labothappy assumption).
    - Inundation (labothappy assumption): vapour flows downward through the bank, number of tubes
      in a vertical column N = N_L = n_rows * n_passes, Kern factor F_inund = N**(-1/6).
    """
    # 1) Retrieve geometry and flow conditions
    m_dot_total = cond['m_dot']
    n_parallel = geom.get('n_parallel', 1)
    m_dot_bank = m_dot_total / n_parallel

    D_t_o = geom['D_t_o']
    L_t = geom['L_t']
    n_t = geom['n_t']
    n_rows = geom['n_rows']
    n_passes = geom['n_passes']
    S_T = geom['P_V']
    S_L = geom['P_H']
    theta = geom['theta']

    # 2) Layout validation from theta
    if theta == 0 or theta == 90:
        layout = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        layout = 'staggered'
    else:
        raise ValueError(
            f"tube_bank_rose_kern_condensing_HTC: theta = {theta} deg, expected 0 or 90 "
            "(inline), 30, 45 or 60 (staggered)."
        )

    # 3) Number of rows crossed by the flow, tubes per row
    N_L = n_rows * n_passes
    tubes_per_row = n_t / N_L

    # 4) Frontal (approach) mass flux
    G_front = m_dot_bank / (tubes_per_row * S_T * L_t)

    # 5) Maximum mass flux G in the minimum flow section of the bank
    if layout == 'inline':
        G = G_front * S_T / (S_T - D_t_o)
    else:
        S_D = np.sqrt(S_L**2 + (S_T / 2.0)**2)
        if S_D < (S_T + D_t_o) / 2.0:
            G = G_front * S_T / (2.0 * (S_D - D_t_o))
        else:
            G = G_front * S_T / (S_T - D_t_o)

    x = cond['x']
    if x < 0.0 or x > 1.0:
        warnings.warn(
            f"tube_bank_rose_kern_condensing_HTC: vapour quality x = {x} outside [0, 1], "
            "clipping to bounds.",
            ValidityRangeWarning
        )
        x = np.clip(x, 0.0, 1.0)

    g = 9.81

    # 6) Obtain saturation properties at pressure p
    s = _sat_props(AS, cond['p'], cond)
    rho_l = s['rho_l']
    mu_l = s['mu_l']
    k_l = s['k_l']
    cp_l = s['cp_l']
    Pr_l = s['Pr_l']
    rho_v = s['rho_v']
    mu_v = s['mu_v']
    h_lv = s['h_lv']
    T_sat = s['T_sat']

    # 7) Temperature difference and check
    dT = T_sat - cond['T_w']
    if dT <= 0.0:
        raise ValueError(
            f"tube_bank_rose_kern_condensing_HTC: no condensation possible, dT = {dT} K <= 0."
        )

    # 8) Vapour velocity approaching the tubes
    U = G * x / rho_v

    # 9) Dimensionless groups (Bonneau nomenclature p. 2, liquid properties)
    Ja = cp_l * dT / h_lv
    H = Ja / Pr_l
    R = np.sqrt(rho_l * mu_l / (rho_v * mu_v))
    P = (Pr_l / Ja) * (rho_v / rho_l)
    Re_t = rho_l * U * D_t_o / mu_l
    Fr = U**2 / (g * D_t_o)

    if U <= 0.0:
        F = np.inf
    else:
        F = Pr_l / (Fr * Ja)

    # 10) Single-tube Nusselt number Nu_1 = h_1 * D / k_l
    D = D_t_o
    if U <= 0.0 or F > 1000.0:
        # a) Nusselt Eq. 12 (gravity-controlled)
        Nu_1 = 0.728 * (rho_l**2 * g * D**3 * h_lv / (mu_l * k_l * dT))**0.25
    elif P > 10.0:
        # b) Fujii Eq. 70
        chi = 0.90 * (1.0 + 1.0 / (R * H))**(1.0 / 3.0)
        Nu_1 = (chi**4.0 + 0.276 * F)**0.25 * Re_t**0.5
    elif P <= F / 8.0:
        # c) Rose Eq. 90
        numerator = 0.9 * (1.0 + 1.0 / (R * H))**(1.0 / 3.0) + 0.728 * F**0.5
        denominator = (1.0 + 3.44 * F**0.5 + F)**0.25
        Nu_1 = (numerator / denominator) * Re_t**0.5
    else:
        # d) Rose Eq. 94
        pressure_term = (1.0 + 1.81 * P)**0.209
        numerator = (0.636 * pressure_term * (1.0 + 1.0 / (R * H))**(1.0 / 3.0)
                     + 0.728 * F**0.5)
        denominator = (1.0 + 3.51 * F**0.53 + F)**0.25
        Nu_1 = (numerator / denominator) * Re_t**0.5

    h_1 = Nu_1 * k_l / D

    # 11) Inundation factor (Kern Eq. 110 with labothappy column count N = N_L)
    F_inund = N_L**(-1.0 / 6.0)

    # 12) Mean heat transfer coefficient over the bundle
    h = h_1 * F_inund

    return h
