# -*- coding: utf-8 -*-
"""
Single-phase heat transfer correlations of Printed Circuit Heat Exchangers (semi-circular channels),
(AS, geom, cond) formalism, prefix PCHE_..._HTC.

geom : D_c (channel diameter), L_c (developed channel length), alpha (bending angle of a zigzag
       channel [deg]); Dh (optional, default pi*D_c/(2+pi))
cond : G (mass flux in ONE channel), mu, k, Pr, rho, T (bulk), m_dot is not needed.
"""
import warnings
import numpy as np


def _dh(geom):
    return geom.get('Dh', np.pi * geom['D_c'] / (2 + np.pi))


def _need(cond, *keys):
    miss = [k for k in keys if k not in cond]
    if miss:
        raise KeyError(f"PCHE HTC: cond is missing {miss} (single-phase bulk properties are required).")


def PCHE_straight_HTC(AS, geom, cond):
    """Straight channel: Stephan-Preusser (Re < 2300) / Gnielinski (Re >= 2300).
    Ref: S.W. Lee et al. (2022), thermal-hydraulic performance of PCHE for SFR recuperators."""
    _need(cond, 'G', 'mu', 'k', 'Pr')
    Dh, L = _dh(geom), geom['L_c']
    Re, Pr = cond['G'] * Dh / cond['mu'], cond['Pr']
    if Re < 2300:
        x = Re * Pr * Dh / L
        Nu = 3.657 + 0.0677 * x ** 1.33 / (1 + 0.1 * Pr * (Re * Dh / L) ** 0.3)
    else:
        f = (1.82 * np.log10(Re) - 1.64) ** -2
        Nu = (f / 8 * (Re - 1000) * Pr / (1 + 12.7 * np.sqrt(f / 8) * (Pr ** (2 / 3) - 1))) \
            * (1 + (Dh / L) ** (2 / 3))
    return Nu * cond['k'] / Dh


def _bennett_chen_inputs(geom, cond):
    """Common inputs of the Bennett & Chen correlations.
    m_dot : mass flow of ONE channel [kg/s]; Re : Reynolds number [local value of the cell, the paper uses
    the inlet value]; theta : bend angle of the paper [rad], the INTERIOR angle (pi = straight channel),
    so theta = pi - alpha with geom['alpha'] the deviation from the flow axis [deg] (Lee convention);
    Lf/Dh : developed channel length over hydraulic diameter."""
    _need(cond, 'G', 'mu')
    Dh, D_c = _dh(geom), geom['D_c']
    m_dot = cond['G'] * np.pi * D_c ** 2 / 8
    Re = cond['G'] * Dh / cond['mu']
    theta = np.pi - np.radians(geom['alpha'])
    return m_dot, Re, theta, geom['L_c'] / Dh


def PCHE_bennett_chen_hot_HTC(AS, geom, cond):
    """Zigzag channel, HOT-side Nusselt number of Bennett & Chen (2020), Thermal Science and Engineering
    Progress 17, 100383 (doi 10.1016/j.tsep.2019.100383), Eq. (22), Table 5 (also Table 3 of Lee et al. 2022,
    where the sign of a1 and the '+' between the terms are lost):
        Nu = a1 Re^a2 m^a3 theta^a4 (Lf/Dh)^a5 + a6 m^a7 T_in^a8
    Supercritical CO2, CFD fit (r2 = 0.64): +-30 % at best, mean error 25 % for Re >= 1200, only
    'a rough estimate'. Validity (Table 3): m 5e-5 .. 5e-4 kg/s per channel, bend angle 84-170 deg,
    T_h,in 500-900 K, p_h,in 2-10 MPa, Dh ~ 0.6-1.2 mm.
    WARNING: with the published coefficients the two terms are almost equal (~1e4) and the ratio
    term2/term1 is < 1 for every Re >= ~10 in the validity range, so Nu comes out NEGATIVE for all
    admissible inputs (checked on m = 5e-5 .. 5e-4 kg/s, T = 500-900 K, Re 1500-16000): the printed
    Table 5 very probably contains a typo. This function raises ValueError in that case; do not use it
    until the coefficients are checked. The cold-side formula behaves normally (Nu ~ 15-60).
    cond : G, mu, k, T_in (inlet temperature of the side [K], default cond['T'])."""
    _need(cond, 'k')
    m, Re, th, LD = _bennett_chen_inputs(geom, cond)
    T_in = cond.get('T_in', cond['T'])
    Nu = (-2.2631e6 * Re ** 1.0515e-3 * m ** 6.5405e-1 * th ** -2.0665e-4 * LD ** 8.0885e-5
          + 2.2967e6 * m ** 6.5429e-1 * T_in ** -1.3319e-3)
    if Nu <= 0:
        raise ValueError(f"PCHE_bennett_chen_hot_HTC: Nu = {Nu:.3g} <= 0, outside the range of the fit "
                         "(Re >= 1200 and m in 5e-5 .. 5e-4 kg/s per channel).")
    return Nu * cond['k'] / _dh(geom)


def PCHE_bennett_chen_cold_HTC(AS, geom, cond):
    """COLD-side counterpart, Eq. (26), Table 6 (p in MPa, assumed from Table 3; T in K):
        Nu = a1 Re^a2 m^a3 theta^a4 (Lf/Dh)^a5 + a6 m^a7 p_c,in^a8 + a9 T_c,in^a10 p_h,in^a11
    Validity (Table 3): T_c,in 305-450 K, p_c,in 7.5-20 MPa, p_h,in 2-10 MPa, m 5e-5 .. 5e-4 kg/s.
    cond : G, mu, k, T_in, p_in [Pa] of the cold side and p_other_in [Pa] (inlet pressure of the hot side)."""
    _need(cond, 'k', 'p_other_in')
    m, Re, th, LD = _bennett_chen_inputs(geom, cond)
    T_in, p_c, p_h = cond.get('T_in', cond['T']), cond.get('p_in', cond['p']) / 1e6, cond['p_other_in'] / 1e6
    Nu = (3.7047e-12 * Re ** 3.0848 * m ** -1.7767 * th ** -2.6196e-1 * LD ** -2.9133
          + 3.0520e3 * m ** 5.6786e-1 * p_c ** -8.7924e-2
          + 1.5024e1 * T_in ** -7.1049 * p_h ** -1.0171e1)
    return Nu * cond['k'] / _dh(geom)


def PCHE_bennett_chen_HTC(AS, geom, cond):
    """Hot formula for the hot side, cold formula for the cold side (cond['heating'] = cold side)."""
    return (PCHE_bennett_chen_cold_HTC if cond.get('heating') else PCHE_bennett_chen_hot_HTC)(AS, geom, cond)


def PCHE_lee_zigzag_HTC(AS, geom, cond):
    """Zigzag semicircular channel of a PCHE: heat transfer coefficient h [W/m2/K].

    Re <= 2000 (laminar / transition), Lee et al. (2024), Eq. (3):
        Nu = 0.0479 Re^0.7089 (h/p)^0.30 Pr^0.33
        h/p = zigzag amplitude / pitch. Table I of the paper (alpha 10/20/30/40 deg) gives
        h/p = tan(alpha) / 2 to 4 digits, so h/p is computed from the bending angle.
        Validity: 300 <= Re <= 2000, 10 <= alpha <= 40 deg, water (CFD, Transition SST),
        D_c = 1.5 mm; +-15 % vs CFD, +-20 % vs Kim & No (2011) and Chen et al. (2016) data.
        Ref: Y. Lee, S. Lee, H.B. Park, K.W. Seo, D. Lee, Thermal-hydraulic performance analysis
        of zigzag channel PCHE according to bending angle, Trans. Korean Nuclear Society Spring
        Meeting, Jeju, May 9-10, 2024.
    Re > 2000 (turbulent), Jiang et al. (2018), Eq. (23a), "high-angle" (40 deg) channel:
        Nu = 0.0845 Re^0.721 Pr^(1/3)
        Validity: 3000 <= Re <= 20600 (fit of Heatric data), sCO2, D_c = 2 mm, alpha = 40 deg.
        Ref: Y. Jiang, E. Liese, S.E. Zitney, D. Bhattacharyya, Design and dynamic modeling of
        printed circuit heat exchangers for supercritical carbon dioxide Brayton power cycles
        (2018).

    Definitions: Re = G Dh / mu with Dh = pi D_c / (2 + pi) (Jiang Eq. 8; Lee does not state it),
    bulk properties of the cell; Nu on Dh.
    geom : alpha [deg] (bending angle from the main flow direction), D_c [m] (or Dh [m])
    cond : G [kg/m2/s] (one channel), mu [Pa s], k [W/m/K], Pr [-]
    Returns h [W/m2/K]. Outside the validity ranges: warning, value extrapolated."""
    _need(cond, 'G', 'mu', 'k', 'Pr')
    Dh = _dh(geom)
    alpha = geom['alpha']
    Pr = cond['Pr']

    # Reynolds number on the hydraulic diameter of the semicircular channel
    Re = cond['G'] * Dh / cond['mu']

    if Re <= 2000:
        # A zigzag needs 0 < alpha < 90 deg (tan(alpha) > 0)
        if not 0 < alpha < 90:
            raise ValueError(f"PCHE_lee_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")
        # Lee et al. (2024): range of the CFD data
        if not 10 <= alpha <= 40:
            warnings.warn(f"PCHE_lee_zigzag_HTC: bending angle {alpha} deg outside the range "
                          "10-40 deg of Lee et al. (2024): extrapolated.")
        if Re < 300:
            warnings.warn(f"PCHE_lee_zigzag_HTC: Re = {Re:.0f} below 300, lower limit of "
                          "Lee et al. (2024): extrapolated.")
        # Amplitude over pitch of the zigzag, from the bending angle (Table I: h/p = tan(alpha)/2)
        h_over_p = np.tan(np.radians(alpha)) / 2
        # Lee et al. (2024), Eq. (3)
        Nu = 0.0479 * Re ** 0.7089 * h_over_p ** 0.30 * Pr ** 0.33
    else:
        # Jiang et al. (2018): fit of Heatric data, high-angle (40 deg) channel only
        if Re < 3000 or Re > 20600:
            warnings.warn(f"PCHE_lee_zigzag_HTC: Re = {Re:.0f} outside 3000-20600, range of "
                          "Jiang et al. (2018) Eq. (23a): extrapolated.")
        if alpha != 40:
            warnings.warn(f"PCHE_lee_zigzag_HTC: Jiang et al. (2018) Eq. (23a) is for a 40 deg "
                          f"channel, used here with alpha = {alpha} deg.")
        # Jiang et al. (2018), Eq. (23a)
        Nu = 0.0845 * Re ** 0.721 * Pr ** (1 / 3)
    return Nu * cond['k'] / Dh


def PCHE_jiang_zigzag_HTC(AS, geom, cond):
    """Zigzag semicircular channel of a PCHE: single-phase TURBULENT heat transfer coefficient
    h [W/m2/K] from Jiang et al. (2018). Two fits are provided: a "low-angle" (30 deg) and a
    "high-angle" (40 deg) channel. The 30 deg fit is used for alpha <= 35 deg, the 40 deg fit above
    (same rule as the pressure-drop function PCHE_lee_jiang_zigzag_SP_DP).

    Reference:
        Y. Jiang, E. Liese, S.E. Zitney, D. Bhattacharyya, Design and dynamic modeling of printed
        circuit heat exchangers for supercritical carbon dioxide Brayton power cycles (2018),
        Eqs. (22a), (23a), (8), Fig. 3. Preprint, journal not stated. The fits are regressions of
        Heatric data (Le Pierres et al., 2011).

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : bending angle of the zigzag from the main flow direction (0 < alpha < 90).
            'D_c' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi) (Jiang Eq. 8).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        3000 <= Re <= 20600 (range of the Heatric data used for the fit).
        alpha = 30 deg (Eq. 22a) or 40 deg (Eq. 23a) only; other angles are extrapolated.
        Fluid: supercritical CO2. Channel diameter D_c = 2 mm (Heatric data).
        Fully developed turbulent flow assumed. No laminar branch: for Re < 3000 use 'Lee'.

    Definitions:
        Re = G * Dh / mu  with Dh = pi * D_c / (2 + pi) (Jiang Eq. 8).
        Nu = h * Dh / k   on the same Dh.
        Bulk properties of the cell; G is the mass flux of one channel.
        AS is not used (common interface).
    """
    # Required bulk properties
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Hydraulic diameter (Jiang Eq. 8) and bending angle
    Dh = _dh(geom)
    alpha = geom['alpha']

    # A zigzag channel requires 0 < alpha < 90 deg
    if not 0 < alpha < 90:
        raise ValueError(f"PCHE_jiang_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")

    # Reynolds number on the hydraulic diameter
    Re = cond['G'] * Dh / cond['mu']
    Pr = cond['Pr']

    # Validity warning for Reynolds number (Heatric data range)
    if Re < 3000 or Re > 20600:
        warnings.warn(f"PCHE_jiang_zigzag_HTC: Re = {Re:.0f} outside 3000-20600, "
                      "range of Jiang et al. (2018) fit: extrapolated.")

    # Validity warning for bending angle (fits exist for 30 and 40 deg only)
    if alpha != 30 and alpha != 40:
        warnings.warn("PCHE_jiang_zigzag_HTC: Jiang et al. (2018) gives fits for 30 and 40 deg "
                      "only: the 30 deg fit is used for alpha <= 35 deg, the 40 deg fit above.")

    # Choice of correlation: low-angle (30 deg) for alpha <= 35, high-angle (40 deg) above
    if alpha <= 35:
        # Jiang et al. (2018), Eq. (22a): low-angle (30 deg) channel
        Nu = 0.0176 * Re ** 0.809 * Pr ** (1 / 3)
    else:
        # Jiang et al. (2018), Eq. (23a): high-angle (40 deg) channel
        Nu = 0.0845 * Re ** 0.721 * Pr ** (1 / 3)

    # Heat transfer coefficient h = Nu * k / Dh
    h = Nu * cond['k'] / Dh
    return h