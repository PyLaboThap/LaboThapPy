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
# Saturated properties for the two-phase PCHE correlations (Kwon et al. 2018)
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props


def _dh(geom):
    # Dh given (any channel section): use it; D_c is only read when Dh is missing (semicircle)
    if geom.get('Dh') is not None:
        return geom['Dh']
    return np.pi * geom['D_ch'] / (2 + np.pi)


def _need(cond, *keys):
    miss = [k for k in keys if k not in cond]
    if miss:
        raise KeyError(f"PCHE HTC: cond is missing {miss} (single-phase bulk properties are required).")


def PCHE_straight_HTC(AS, geom, cond):
    """Straight channel: Stephan-Preusser (Re < 2300) / Gnielinski (Re >= 2300).
    Ref: S.W. Lee et al. (2022), thermal-hydraulic performance of PCHE for SFR recuperators."""
    _need(cond, 'G', 'mu', 'k', 'Pr')
    Dh, L = _dh(geom), geom['L_ch']
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
    the inlet value]; theta : bend angle of the paper [rad], the INTERIOR angle at the apex of the zigzag
    (Fig. 1 of the paper; pi = straight channel). Each segment deviates by alpha from the main flow
    direction (geom['alpha'] [deg], Lee convention), so theta = pi - 2 alpha
    (Table 3: theta 84-170 deg <-> alpha 5-48 deg). Corrected 2026-11-23 (A2-bis, was pi - alpha);
    Lf/Dh : developed channel length over hydraulic diameter."""
    _need(cond, 'G', 'mu')
    Dh, D_c = _dh(geom), geom['D_ch']
    m_dot = cond['G'] * np.pi * D_c ** 2 / 8
    Re = cond['G'] * Dh / cond['mu']
    theta = np.pi - 2 * np.radians(geom['alpha'])
    return m_dot, Re, theta, geom['L_ch'] / Dh


def PCHE_bennett_chen_hot_HTC(AS, geom, cond):
    """Zigzag channel, HOT-side Nusselt number of Bennett & Chen (2020), Thermal Science and Engineering
    Progress 17, 100383 (doi 10.1016/j.tsep.2019.100383), Eq. (22), Table 5 (also Table 3 of Lee et al. 2022,
    where the sign of a1 and the '+' between the terms are lost):
        Nu = a1 Re^a2 m^a3 theta^a4 (Lf/Dh)^a5 + a6 m^a7 T_r^a8
    T_r = T_h,in / T_critical is the REDUCED inlet temperature: stated only in K. Bennett's PhD
    dissertation (UNLV 2019, doi 10.34917/18608592, p. 187, Eq. 80, Table 30 = Table 5 of the paper).
    With T in K (as the paper seems to imply) Nu is negative everywhere. Checked 2026-11-23 (A2-ter):
    with T_r the printed coefficients give 17.65 for the thesis Ishizuka case (thesis: 17.63) and
    reproduce the thesis values of the 50 CFD cases (median 2 %).
    m = mass flow of ONE hot channel [kg/s] (the thesis CFD cases have two hot channels per cold one).
    Supercritical CO2, CFD fit (r2 = 0.64): +-30 % at best, mean error 25 % for Re >= 1200, only
    'a rough estimate'. Validity (Table 3): Re >= 1200, bend angle 84-170 deg, T_h,in 500-900 K,
    p_h,in 2-10 MPa, Dh ~ 0.6-1.2 mm, m about 2.5e-5 .. 2.5e-4 kg/s per hot channel.
    The two terms nearly cancel (~1e4 each): a ValueError is kept as a guard if Nu <= 0.
    cond : G, mu, k, T_in (inlet temperature of the side [K], default cond['T']).
    AS : CoolProp AbstractState of the fluid (T_critical)."""
    _need(cond, 'k')
    m, Re, th, LD = _bennett_chen_inputs(geom, cond)
    # Inlet temperature: cond['T'] is read only when 'T_in' is missing
    if 'T_in' in cond:
        T_in = cond['T_in']
    else:
        T_in = cond['T']
    # Reduced inlet temperature (thesis p. 187)
    T_r = T_in / AS.T_critical()
    Nu = (-2.2631e6 * Re ** 1.0515e-3 * m ** 6.5405e-1 * th ** -2.0665e-4 * LD ** 8.0885e-5
          + 2.2967e6 * m ** 6.5429e-1 * T_r ** -1.3319e-3)
    if Nu <= 0:
        raise ValueError(f"PCHE_bennett_chen_hot_HTC: Nu = {Nu:.3g} <= 0, outside the range of the fit "
                         "(Re >= 1200 and m in 5e-5 .. 5e-4 kg/s per channel).")
    return Nu * cond['k'] / _dh(geom)


def PCHE_bennett_chen_cold_HTC(AS, geom, cond):
    """COLD-side counterpart, Eq. (26), Table 6 (= thesis Eq. 84, Table 31):
        Nu = a1 Re^a2 m^a3 theta^a4 (Lf/Dh)^a5 + a6 m^a7 p_r,c^a8 + a9 T_r,c^a10 p_h,in^a11
    Reduced temperature and pressure (K. Bennett PhD dissertation, UNLV 2019, p. 195): T_r,c = T_c,in /
    T_critical, p_r,c = p_c,in / p_critical. p_h,in (hot inlet pressure) kept in MPa: labothappy reading,
    reducing it too breaks the agreement with the thesis values; the third term is then negligible.
    Checked 2026-11-23 (A2-ter): thesis values of the 50 CFD cases, median ratio 0.996.
    m = mass flow of ONE cold channel [kg/s].
    Validity (Table 3): T_c,in 305-450 K, p_c,in 7.5-20 MPa, p_h,in 2-10 MPa, m 5e-5 .. 5e-4 kg/s.
    cond : G, mu, k, T_in, p_in [Pa] of the cold side and p_other_in [Pa] (inlet pressure of the hot side).
    AS : CoolProp AbstractState of the fluid (T_critical, p_critical)."""
    _need(cond, 'k', 'p_other_in')
    m, Re, th, LD = _bennett_chen_inputs(geom, cond)
    # Inlet temperature and pressure: cond['T'] / cond['p'] are read only when T_in / p_in are missing
    if 'T_in' in cond:
        T_in = cond['T_in']
    else:
        T_in = cond['T']
    if 'p_in' in cond:
        p_in = cond['p_in']
    else:
        p_in = cond['p']
    # Reduced cold inlet temperature and pressure (thesis p. 195); hot inlet pressure in MPa
    T_r = T_in / AS.T_critical()
    p_r = p_in / AS.p_critical()
    p_h = cond['p_other_in'] / 1e6
    Nu = (3.7047e-12 * Re ** 3.0848 * m ** -1.7767 * th ** -2.6196e-1 * LD ** -2.9133
          + 3.0520e3 * m ** 5.6786e-1 * p_r ** -8.7924e-2
          + 1.5024e1 * T_r ** -7.1049 * p_h ** -1.0171e1)
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
    2000 < Re < 3000: linear blend in Re between Lee Eq. (3) at Re = 2000 and Jiang Eq. (23a) at
        Re = 3000 (labothappy choice, as the 'Lee_Jiang' pressure drop; no warning).
    Re >= 3000 (turbulent), Jiang et al. (2018), Eq. (23a), "high-angle" (40 deg) channel:
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
    elif Re < 3000:
        # 2000 < Re < 3000: between the two ranges (Lee up to 2000, Jiang from 3000).
        # Linear blend in Re between Lee Eq. (3) at Re = 2000 and Jiang Eq. (23a) at Re = 3000,
        # same rule as the PCHE pressure drop 'Lee_Jiang' (2026-11-23, proposed_changes E19:
        # without it Nu jumped x2.9 at Re = 2000). No range warning inside the blend.
        if not 0 < alpha < 90:
            raise ValueError(f"PCHE_lee_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")
        h_over_p = np.tan(np.radians(alpha)) / 2
        # End point of the laminar branch: Lee et al. (2024), Eq. (3) at Re = 2000
        Nu_lee_2000 = 0.0479 * 2000.0 ** 0.7089 * h_over_p ** 0.30 * Pr ** 0.33
        # Start point of the turbulent branch: Jiang et al. (2018), Eq. (23a) at Re = 3000
        Nu_jiang_3000 = 0.0845 * 3000.0 ** 0.721 * Pr ** (1 / 3)
        # Weight of the turbulent end: 0 at Re = 2000, 1 at Re = 3000
        weight = (Re - 2000.0) / 1000.0
        Nu = Nu_lee_2000 + weight * (Nu_jiang_3000 - Nu_lee_2000)
    else:
        # Jiang et al. (2018): fit of Heatric data, high-angle (40 deg) channel only
        if Re > 20600:
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
    h [W/m2/K] from Jiang et al. (2018). The paper gives a "low-angle" (30 deg) and a
    "high-angle" (40 deg) fit. Since 2026-11-23 (E21) the 40 deg fit is used for EVERY angle: the
    30 deg fit is ~2 x below the independent sources at 32.5-37 deg (LOGS report
    e16_bennett_check_2026-11-23.md). NOTE: the pressure-drop function PCHE_lee_jiang_zigzag_SP_DP
    still switches to the 30 deg friction fit for alpha <= 35 deg.

    Reference:
        Y. Jiang, E. Liese, S.E. Zitney, D. Bhattacharyya, Design and dynamic modeling of printed
        circuit heat exchangers for supercritical carbon dioxide Brayton power cycles (2018),
        Eqs. (22a), (23a), (8), Fig. 3. Preprint, journal not stated. The fits are regressions of
        Heatric data (Le Pierres et al., 2011).

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : bending angle of the zigzag from the main flow direction (0 < alpha < 90).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
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
        Eq. (23a) fitted at alpha = 40 deg; other angles are extrapolated (warning).
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

    # Validity warning for bending angle (Eq. 23a is the 40 deg fit)
    if alpha != 40:
        warnings.warn("PCHE_jiang_zigzag_HTC: Jiang et al. (2018) gives fits for 30 and 40 deg "
                      f"only; the 40 deg fit (Eq. 23a) is used for every angle (alpha = {alpha} deg).")

    # Jiang et al. (2018), Eq. (23a): high-angle (40 deg) channel, used for EVERY angle.
    # 2026-11-23 (E21, Basile's choice): the low-angle fit Eq. (22a), Nu = 0.0176 Re^0.809 Pr^(1/3),
    # is no longer used: at 32.5-37 deg it gives ~0.5 x every independent source (Kim SG 2016,
    # Meshram 2016, Katz 2021) and ~0.9 x a straight channel; Eq. (23a) is within 0.7-1.4 x of them.
    Nu = 0.0845 * Re ** 0.721 * Pr ** (1 / 3)

    # Heat transfer coefficient h = Nu * k / Dh
    h = Nu * cond['k'] / Dh
    return h


def PCHE_ngo_zigzag_HTC(AS, geom, cond):
    """Single-phase turbulent heat transfer coefficient in a zigzag channel of a PCHE
    from Ngo et al. (2007), zigzag-fin correlation (Eq. 29).

    Reference:
        T.L. Ngo, Y. Kato, K. Nikitin, T. Ishizuka (2007), Heat transfer and pressure drop
        correlations of microchannel heat exchangers with S-shaped and zigzag fins for carbon
        dioxide cycles, Experimental Thermal and Fluid Science 32 (2007) 560-570,
        doi:10.1016/j.expthermflusci.2007.06.006, Eq. (29) (zigzag fins) and Eqs. (12), (17).

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : bending angle of the zigzag from the main flow direction
                            (same as Ngo's fin angle theta, Fig. 3b).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        3500 < Re < 22000, 0.75 < Pr < 2.2, supercritical CO2 (experiments, 6-12 MPa),
        single geometry: fin angle 52 deg, Dh = 1.09 mm, fin gap 1.31 mm, fin depth 0.94 mm,
        fin width 0.8 mm, fin pitch 7.565 mm (x) and 3.426 mm (y), 316L stainless steel; U
        reproduced with a standard deviation of +-3.0 %. No laminar branch: for Re < 3500 use
        'Lee'.
        Note: Ngo's channels are etched "zigzag fins" (rectangular-like cross-section
        1.31 x 0.94 mm), not semicircular channels.

    Definitions:
        Re = G * Dh / mu, Nu = h * Dh / k on the hydraulic diameter of the channel; G of one
        channel (mass flow over the free flow area); bulk properties of the cell (the paper used
        properties averaged along the exchanger, Eq. 17); returns h [W/m2/K] on the channel
        wall. AS is not used.
    """
    # Check that required bulk properties are present in the condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Get the hydraulic diameter and bending angle from the geometry dictionary
    Dh = _dh(geom)
    alpha = geom['alpha']

    # A zigzag needs 0 < alpha < 90 deg: otherwise raise ValueError
    if not 0 < alpha < 90:
        raise ValueError(f"PCHE_ngo_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")

    # Extract properties for the calculations
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Reynolds number (Ngo Eq. 17: Re = W Dh / (Ac mu), W/Ac = G of one channel)
    Re = G * Dh / mu

    # Validity warning for Reynolds number (Eq. 29 range: 3500 to 22000)
    if Re < 3500 or Re > 22000:
        warnings.warn(f"PCHE_ngo_zigzag_HTC: Re = {Re:.0f} outside 3500-22000, "
                      "range of Ngo et al. (2007) Eq. 29: extrapolated.")

    # Validity warning for Prandtl number (Eq. 29 range: 0.75 to 2.2)
    if Pr < 0.75 or Pr > 2.2:
        warnings.warn(f"PCHE_ngo_zigzag_HTC: Pr = {Pr:.2f} outside 0.75-2.2, "
                      "range of Ngo et al. (2007) Eq. 29: extrapolated.")

    # Validity warning for fin angle (Ngo et al. tested a single geometry with 52 deg)
    if alpha != 52:
        warnings.warn(f"PCHE_ngo_zigzag_HTC: Ngo et al. (2007) tested a single geometry, "
                      f"fin angle 52 deg (used here with alpha = {alpha} deg).")

    # Nusselt number correlation for zigzag fins (Ngo et al. 2007, Eq. 29)
    term_re = Re ** 0.629
    term_pr = Pr ** 0.317
    Nu = 0.1696 * term_re * term_pr

    # Heat transfer coefficient (Ngo Eq. 12)
    h = Nu * k / Dh

    return h


def PCHE_ngo_S_fin_HTC(AS, geom, cond):
    """Single-phase TURBULENT heat transfer coefficient in an S-shaped fin channel
    of a PCHE from Ngo et al. (2007), S-shaped fin correlation (Eq. 28).

    Reference:
        T.L. Ngo, Y. Kato, K. Nikitin, T. Ishizuka (2007), Heat transfer and pressure drop
        correlations of microchannel heat exchangers with S-shaped and zigzag fins for carbon
        dioxide cycles, Experimental Thermal and Fluid Science 32 (2007) 560-570,
        doi:10.1016/j.expthermflusci.2007.06.006, Eq. (28) (S-shaped fins) and Eqs. (12), (17).

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : fin angle from the main flow direction (0 < alpha < 90)
                            (Ngo's theta, Fig. 3a: straight middle part of the S fin).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        3500 < Re < 23000, 0.75 < Pr < 2.2, supercritical CO2 (experiments, 6-12 MPa),
        single geometry: fin angle 52 deg, Dh = 1.09 mm, fin gap 1.31 mm, fin depth 0.94 mm,
        fin width 0.8 mm, fin pitch 7.565 mm (x) and 3.426 mm (y), 316L stainless steel; U
        reproduced with a standard deviation of +-2.3 %. No laminar branch: no S-fin correlation
        exists for Re < 3500 in labothappy.
        Note: discontinuous S-shaped fins (rectangular-like passage 1.31 x 0.94 mm between fins),
        not semicircular channels.

    Definitions:
        Re = G * Dh / mu, Nu = h * Dh / k on the hydraulic diameter of the channel; G of one
        channel (mass flow over the free flow area); bulk properties of the cell (the paper used
        properties averaged along the exchanger, Eq. 17); returns h [W/m2/K] on the channel
        wall. AS is not used.
    """
    # Check that required bulk properties are present in the condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Get the hydraulic diameter and fin angle from the geometry dictionary
    Dh = _dh(geom)
    alpha = geom['alpha']

    # The fin angle must be 0 < alpha < 90 deg: otherwise raise ValueError
    if not 0 < alpha < 90:
        raise ValueError(f"PCHE_ngo_S_fin_HTC: fin angle {alpha} deg outside ]0, 90[.")

    # Extract properties for the calculations
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Reynolds number (Ngo Eq. 17: Re = W Dh / (Ac mu), W/Ac = G of one channel)
    Re = G * Dh / mu

    # Validity warning for Reynolds number (Eq. 28 range: 3500 to 23000)
    if Re < 3500 or Re > 23000:
        warnings.warn(f"PCHE_ngo_S_fin_HTC: Re = {Re:.0f} outside 3500-23000, "
                      "range of Ngo et al. (2007) Eq. 28: extrapolated.")

    # Validity warning for Prandtl number (Eq. 28 range: 0.75 to 2.2)
    if Pr < 0.75 or Pr > 2.2:
        warnings.warn(f"PCHE_ngo_S_fin_HTC: Pr = {Pr:.2f} outside 0.75-2.2, "
                      "range of Ngo et al. (2007) Eq. 28: extrapolated.")

    # Validity warning for fin angle (Ngo et al. tested a single geometry with 52 deg)
    if alpha != 52:
        warnings.warn(f"PCHE_ngo_S_fin_HTC: Ngo et al. (2007) tested a single geometry, "
                      f"fin angle 52 deg (used here with alpha = {alpha} deg).")

    # Nusselt number correlation for S-shaped fins (Ngo et al. 2007, Eq. 28)
    term_re = Re ** 0.593
    term_pr = Pr ** 0.430
    Nu = 0.1740 * term_re * term_pr

    # Heat transfer coefficient (Ngo Eq. 12)
    h = Nu * k / Dh

    return h


def PCHE_kim_no_2012_zigzag_HTC(AS, geom, cond):
    """Single-phase LAMINAR heat transfer coefficient in a zigzag semicircular channel of a PCHE
    from Kim & No (2012), Nusselt correlation Eq. (4) with the constants of Table 5.

    Reference:
        I.H. Kim, H.C. No (2012), Physical model development and optimal design of PCHE for
        intermediate heat exchangers in HTGRs, Nuclear Engineering and Design 243, 243-250,
        doi:10.1016/j.nucengdes.2011.11.020, Eqs. (2), (4), Tables 3 and 5.

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : bending angle of the zigzag from the main flow direction
                            (0 < alpha < 90, paper Eq. 9).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
            'P_zz' [m] : optional pitch of the zigzag (Kim & No Fig. 2).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        Laminar flow, Re up to about 2500 (range of the CFD data, Figs. 4-5; no explicit Re
        range in the paper). HELIUM only (FLUENT CFD; "the Prandtl number was nearly constant
        in our test ranges", so Nu has no Pr term: do not use for liquids). Angle 5-45 deg.
        Dh 0.922 or 1.222 mm, pitch 12.3 or 24.6 mm (three CFD geometries of Table 5).
        RMS error of each fit 0.4-3 % (Table 5). H1 boundary condition (constant axial heat
        flux, peripherally uniform wall temperature).

    Definitions:
        Nu = h * Dh / k on the hydraulic diameter; Re = G * Dh / mu; G of one channel;
        bulk properties of the cell (the paper: local pitch-averaged properties).
        alpha = angle of the zigzag from the main flow direction; P_zz as in Kim & No
        Fig. 2.
        Interpolation in angle and nearest-geometry choice: added in labothappy (not in the paper).
        Returns h [W/m2/K] on the channel wall. AS is not used.
    """
    # Required bulk properties
    _need(cond, 'G', 'mu', 'k')

    # Hydraulic diameter and bending angle
    Dh = _dh(geom)
    alpha = geom['alpha']

    # A zigzag needs 0 < alpha < 90 deg
    if not 0 < alpha < 90:
        raise ValueError(f"PCHE_kim_no_2012_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")

    # Zigzag pitch (optional, default from paper geometry (a))
    pitch = geom.get('P_zz')
    if pitch is None:
        pitch = 24.6e-3
        warnings.warn("PCHE_kim_no_2012_zigzag_HTC: geom['P_zz'] not given: "
                      "24.6 mm assumed")

    # Extract properties
    G = cond['G']
    mu = cond['mu']
    k = cond['k']

    # Reynolds number (Kim & No Eq. 2)
    Re = G * Dh / mu

    # Validity warning for Reynolds number (CFD range up to 2500)
    if Re > 2500:
        warnings.warn("PCHE_kim_no_2012_zigzag_HTC: Re = {:.0f} above the laminar CFD range of "
                      "Kim & No (2012) (Re up to 2500, Figs. 4-5)".format(Re))

    # Table 5 data: three geometries (a), (b), (c)
    # Angles [deg] (same for all three tables)
    angles = [5, 10, 15, 20, 25, 30, 35, 40, 45]

    # Geometry (a): Dh = 0.922 mm, pitch = 24.6 mm
    Dh_a = 0.922e-3
    pitch_a = 24.6e-3
    c_a = [0.00036, 0.0022, 0.00544, 0.00894, 0.01321, 0.01877, 0.02578, 0.03489, 0.04885]
    d_a = [1.18213, 0.99841, 0.91361, 0.86708, 0.82738, 0.78945, 0.75403, 0.71951, 0.67896]

    # Geometry (b): Dh = 0.922 mm, pitch = 12.3 mm
    Dh_b = 0.922e-3
    pitch_b = 12.3e-3
    c_b = [0.00061, 0.00251, 0.00607, 0.0118, 0.02165, 0.03006, 0.03845, 0.04829, 0.06368]
    d_b = [1.13917, 1.01844, 0.93106, 0.85964, 0.79051, 0.76007, 0.73793, 0.7159, 0.68516]

    # Geometry (c): Dh = 1.222 mm, pitch = 24.6 mm
    Dh_c = 1.222e-3
    pitch_c = 24.6e-3
    c_c = [0.00071, 0.00314, 0.0083, 0.01703, 0.02182, 0.0299, 0.04247, 0.05186, 0.05988]
    d_c = [1.10341, 0.96567, 0.86054, 0.79007, 0.77285, 0.7403, 0.70055, 0.68039, 0.66801]

    # Collect tables for distance calculation
    tables = [
        (Dh_a, pitch_a, c_a, d_a),
        (Dh_b, pitch_b, c_b, d_b),
        (Dh_c, pitch_c, c_c, d_c)
    ]

    # Choose the table with the minimum relative distance (rule of labothappy, not of the paper)
    min_dist = float('inf')
    chosen_table = None
    for Dh_t, pitch_t, c_t, d_t in tables:
        dist = ((Dh - Dh_t) / Dh_t) ** 2 + ((pitch - pitch_t) / pitch_t) ** 2
        if dist < min_dist:
            min_dist = dist
            chosen_table = (Dh_t, pitch_t, c_t, d_t)

    Dh_table, pitch_table, c_table, d_table = chosen_table

    # Warn if geometry differs by more than 10% from chosen table
    if abs(Dh - Dh_table) / Dh_table > 0.1 or abs(pitch - pitch_table) / pitch_table > 0.1:
        warnings.warn("PCHE_kim_no_2012_zigzag_HTC: geometry outside the three CFD geometries of "
                      "Kim & No (2012): nearest table used")

    # Angle handling: clamp to 5-45 deg range
    alpha_used = alpha
    if alpha < 5:
        warnings.warn("PCHE_kim_no_2012_zigzag_HTC: bending angle {:.1f} deg outside 5-45 deg, "
                      "value at the nearest table angle".format(alpha))
        alpha_used = 5
    elif alpha > 45:
        warnings.warn("PCHE_kim_no_2012_zigzag_HTC: bending angle {:.1f} deg outside 5-45 deg, "
                      "value at the nearest table angle".format(alpha))
        alpha_used = 45

    # Find neighbouring table angles for linear interpolation
    # angles list is sorted: [5, 10, 15, 20, 25, 30, 35, 40, 45]
    alpha_1 = None
    alpha_2 = None
    idx_1 = None
    idx_2 = None
    for i, ang in enumerate(angles):
        if ang <= alpha_used:
            alpha_1 = ang
            idx_1 = i
        if ang >= alpha_used and alpha_2 is None:
            alpha_2 = ang
            idx_2 = i

    # If alpha_used exactly matches a table angle, use that directly
    if alpha_1 == alpha_2:
        c_1 = c_table[idx_1]
        d_1 = d_table[idx_1]
        Nu = 4.089 + c_1 * Re ** d_1
    else:
        # Linear interpolation in angle between the two neighbouring table angles
        c_1 = c_table[idx_1]
        d_1 = d_table[idx_1]
        c_2 = c_table[idx_2]
        d_2 = d_table[idx_2]

        # Nusselt at lower angle (Kim & No Eq. 4: Nu = 4.089 + c * Re^d)
        Nu_1 = 4.089 + c_1 * Re ** d_1
        # Nusselt at upper angle
        Nu_2 = 4.089 + c_2 * Re ** d_2

        # Linear interpolation in angle
        Nu = Nu_1 + (Nu_2 - Nu_1) * (alpha_used - alpha_1) / (alpha_2 - alpha_1)

    # Heat transfer coefficient h = Nu * k / Dh
    h = Nu * k / Dh

    return h


def PCHE_kim_no_2011_zigzag_HTC(AS, geom, cond):
    """Single-phase laminar heat transfer coefficient in a zigzag semicircular channel
    of a printed circuit heat exchanger (PCHE) from Kim & No (2011), Eq. (22).

    Reference:
        I.H. Kim, H.C. No (2011), Thermal hydraulic performance analysis of a printed circuit heat
        exchanger using a helium-water test loop and numerical simulations, Applied Thermal
        Engineering 31, 4064-4073, doi:10.1016/j.applthermaleng.2011.08.012, Eq. (22).

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : bending angle of the zigzag from the main flow
                            direction (0 < alpha < 90).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        Laminar, Re up to 2500, helium and water (Pr range not stated in the paper),
        single geometry: angle 15 deg, channel diameter 1.51 mm, Dh 0.922 mm, pitch 24.6 mm;
        average error 3.6 %, maximum 26 %.

    Definitions:
        Re = G * Dh / mu and Nu = h * Dh / k on the hydraulic diameter; G of one channel;
        bulk properties of the cell (paper: local pitch-averaged properties).
        Returns h [W/m2/K] on the channel wall. AS is not used.
    """
    # Check that required bulk properties are present in the condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Get the hydraulic diameter and bending angle from the geometry dictionary
    Dh = _dh(geom)
    alpha = geom['alpha']

    # A zigzag needs 0 < alpha < 90 deg: otherwise raise ValueError
    if not 0 < alpha < 90:
        raise ValueError(f"PCHE_kim_no_2011_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")

    # Extract properties for the calculations
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Reynolds number on the hydraulic diameter
    Re = G * Dh / mu

    # Validity warning for Reynolds number (CFD range up to Re 2500)
    if Re > 2500:
        warnings.warn(f"PCHE_kim_no_2011_zigzag_HTC: Re = {Re:.0f} above the range of "
                      "Kim & No (2011) (CFD up to Re 2500): extrapolated.")

    # Validity warning for bending angle (tested a single geometry with 15 deg)
    if alpha != 15:
        warnings.warn(f"PCHE_kim_no_2011_zigzag_HTC: Kim & No (2011) studied a single geometry, "
                      f"angle 15 deg (Dh 0.922 mm, pitch 24.6 mm), used here with "
                      f"alpha = {alpha} deg.")

    # Nusselt number correlation (Kim & No 2011, Eq. 22)
    # 4.089 is the Nu of a straight semicircular duct, H1 condition
    term_re_pr = Re * (Pr ** 0.58)
    Nu = 4.089 + 0.00365 * term_re_pr

    # Heat transfer coefficient
    h = Nu * k / Dh

    return h


def PCHE_liu_straight_HTC(AS, geom, cond):
    """Single-phase TURBULENT heat transfer coefficient [W/m2/K] of supercritical CO2 COOLED in a
    STRAIGHT semicircular channel of a printed circuit heat exchanger (PCHE) from Liu et al. (2020),
    Nusselt correlation Eq. (19) with a wall / bulk specific heat ratio.

    Reference:
        S.-h. Liu, Y.-p. Huang, J.-f. Wang, et al. (2020), Experimental study of thermal-hydraulic
        performance of a printed circuit heat exchanger with straight channels, International
        Journal of Heat and Mass Transfer 160, 120109, doi:10.1016/j.ijheatmasstransfer.2020.120109,
        Eqs. (8), (19).

    Arguments:
        AS : CoolProp AbstractState of the fluid (moved to the wall state during evaluation).
        geom : dict with keys
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell and operating conditions
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity (bulk).
            'k'   [W/m/K]   : thermal conductivity (bulk).
            'Pr'  [-]       : Prandtl number (bulk).
            'cp'  [J/kg/K]  : specific heat capacity (bulk).
            'p'   [Pa]      : pressure.
            'T'   [K]       : bulk temperature.
            'T_w' [K]       : wall temperature.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        3600 < Re_b < 36500, supercritical CO2 cooled by water, pressure 7.53-11.97 MPa,
        bulk temperature 26.4-121.3 C, mass flux 50-780 kg/m2/s, Pr_b 0.9-10.6 (Table 3);
        straight semicircular channels, channel diameter 1.87 mm, D_eq 1.14 mm, length 500 mm;
        80.1 % of the data within +-20 %, 92.4 % within +-30 %.

    Definitions:
        Nu = h D_eq / k, Re = rho_b u_b D_eq / mu_b = G Dh / mu (bulk properties, subscript b);
        cp_w = specific heat at the wall temperature and the cell pressure, cp = bulk specific heat.
        G is the mass flux of one channel. Returns h [W/m2/K] on the channel wall.
    """
    import CoolProp.CoolProp as CP

    # Check that required properties are present in the condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr', 'cp', 'p', 'T', 'T_w')

    # Hydraulic diameter and cell properties (bulk), pressure and wall temperature
    Dh = _dh(geom)
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']
    cp = cond['cp']
    p = cond['p']
    T = cond['T']
    T_w = cond['T_w']

    # 1. Reynolds number (Liu Eq. 8, bulk properties, equivalent diameter)
    Re = G * Dh / mu

    # 2. Wall specific heat: update AS to wall state (pressure and wall temperature)
    AS.update(CP.PT_INPUTS, p, T_w)
    cp_w = AS.cpmass()

    # Validity warnings for Reynolds number (Eq. 19 range: 3600 to 36500)
    if Re < 3600 or Re > 36500:
        warnings.warn(f"PCHE_liu_straight_HTC: Re = {Re:.0f} outside 3600-36500, "
                      "range of Liu et al. (2020) Eq. 19: extrapolated.")

    # Validity warnings for Prandtl number (Table 3 range: 0.9 to 10.6)
    if Pr < 0.9 or Pr > 10.6:
        warnings.warn(f"PCHE_liu_straight_HTC: Pr = {Pr:.2f} outside 0.9-10.6, "
                      "range of Liu et al. (2020) Table 3: extrapolated.")

    # Validity warning for cooling mode (wall temperature should be lower than bulk temperature)
    if T_w > T:
        warnings.warn("PCHE_liu_straight_HTC: Liu et al. (2020) correlation is for sCO2 COOLING "
                      "(wall colder than the bulk).")

    # 3. Liu et al. (2020) Eq. (19) Nusselt number calculation (one factor per line)
    term_re = 0.1229 * Re ** 0.6021
    term_pr = Pr ** 0.3
    term_cp = (cp_w / cp) ** 0.1310
    Nu = term_re * term_pr * term_cp

    # 4. Heat transfer coefficient calculation
    h = Nu * k / Dh

    return h


def PCHE_kwon_straight_HTC(AS, geom, cond):
    """Single-phase TURBULENT heat transfer coefficient [W/m2/K] in a STRAIGHT semicircular PCHE
    channel, modified Dittus-Boelter correlation of Kwon et al. (2018), Eq. (5), fitted on
    nitrogen GAS at cryogenic / ambient temperature.

    Reference:
        D. Kwon, L. Jin, W. Jung, S. Jeong (2018), Experimental investigation of heat transfer
        coefficient of mini-channel PCHE (printed circuit heat exchanger), Cryogenics 92, 41-49,
        doi:10.1016/j.cryogenics.2018.03.011, Eqs. (1), (5).

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        8500 < Re < 17000, nitrogen gas on both sides (hot side 239-296 K, cold side 114-244 K,
        0.33-0.68 MPa, G 52-149 kg/m2/s, Tables 3 and 6); predicts the measured U within 5 %.
        Pr range not stated in the paper (nitrogen gas: Pr about 0.7-0.8); the Pr exponent 0.33
        was fixed by the authors, not fitted.
        ONE geometry only: STRAIGHT semicircular channels (short inlet / outlet bends), channel
        diameter 3 mm, hydraulic diameter 1.83 mm, pitch 4 mm, length 117 mm, 10 channels per
        layer, 1 hot + 1 cold layer, stainless steel, counter-flow. The h values were extracted
        from the overall U with 1/U = 1/h_h + 1/h_c (wall neglected).

    Definitions:
        Re = G * Dh / mu and Nu = h * Dh / k on the hydraulic diameter (1.83 mm in the paper);
        G of one channel; bulk properties of the cell (paper: REFPROP 9.1 properties at the mean
        temperature and mean pressure of the side). Returns h [W/m2/K] on the channel wall.
    """
    # Check that required bulk properties are present in the condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Get the hydraulic diameter and required properties
    Dh = _dh(geom)
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # 1. Reynolds number (Kwon Eq. 1)
    Re = G * Dh / mu

    # Validity warning for Reynolds number (range of Kwon et al. (2018) Eq. 5)
    if Re < 8500 or Re > 17000:
        warnings.warn(f"PCHE_kwon_straight_HTC: Re = {Re:.0f} outside 8500-17000, "
                      "range of Kwon et al. (2018) Eq. 5: extrapolated.")

    # 2. Kwon et al. (2018) Eq. (5) Nusselt number (one factor per line)
    term_re = 0.0274 * Re ** 0.78
    term_pr = Pr ** 0.33
    Nu = term_re * term_pr

    # 3. Heat transfer coefficient h = Nu * k / Dh
    h = Nu * k / Dh

    return h


def PCHE_kwon_condensing_HTC(AS, geom, cond):
    """CONDENSATION heat transfer coefficient [W/m2/K] in a STRAIGHT semicircular PCHE channel,
    Akers-type correlation of Kwon et al. (2018), Eqs. (13) and (15), fitted on CO2 condensing at
    about 0.7 MPa (T_sat about 224 K).

    Reference:
        D. Kwon, L. Jin, W. Jung, S. Jeong (2018), Experimental investigation of heat transfer
        coefficient of mini-channel PCHE (printed circuit heat exchanger), Cryogenics 92, 41-49,
        doi:10.1016/j.cryogenics.2018.03.011, Eqs. (13), (15).

    Arguments:
        AS : CoolProp AbstractState of the fluid (moved to saturation by _sat_props).
        geom : dict with keys
            'D_ch' [m] : channel diameter (used if 'Dh' not provided).
            'Dh'  [m] : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with operating conditions
            'G'   [kg/m2/s] : total mass flux of ONE channel (both phases).
            'p'   [Pa]      : pressure of the cell.
            'x'   [-]       : local (cell mean) vapour quality.

    Returns:
        h [W/m2/K] : condensation heat transfer coefficient on the channel wall.

    Validity:
        CO2 (hot side) condensed by cold helium gas, p 0.69-0.70 MPa (T_sat about 224 K),
        G 56-104 kg/m2/s, x 0.68-1 (Re_eq computed with the mean quality of the channel),
        3100 < Re_eq < 6200 (Table 10); predicts the 5 data points within 10 %, but the
        uncertainty of the measured condensation h is 24-29 % (Table 12). Only 5 data points,
        one fluid. The Pr exponent 0.33 was fixed by the authors, not fitted.
        ONE geometry only: STRAIGHT semicircular channels (short inlet / outlet bends), channel
        diameter 3 mm, hydraulic diameter 1.83 mm, pitch 4 mm, length 117 mm, 10 channels per
        layer, 1 hot + 1 cold layer, stainless steel, counter-flow. The h values were extracted
        from the overall U with 1/U = 1/h_h + 1/h_c (wall neglected).

    Definitions:
        Re_eq = G Dh / mu_l [(1 - x) + x (rho_l / rho_v)^0.5] (Akers equivalent all-liquid
        Reynolds number), Nu_l = h Dh / k_l; saturated liquid (l) and vapour (v) properties at the
        cell pressure; x = local (cell mean) vapour quality (the paper uses the mean quality of
        the whole channel); G of one channel (total mass flux, both phases).
        No heat flux or wall temperature needed. Returns h [W/m2/K] on the channel wall.
        (Independent check: with CoolProp, this Re_eq reproduces Table 10 within 1 %.)
    """
    _need(cond, 'G', 'p', 'x')

    # Critical pressure, read before _sat_props moves AS
    p_crit = AS.p_critical()

    Dh = _dh(geom)
    G = cond['G']
    p = cond['p']
    x = cond['x']

    # Keep the quality inside [0, 1]
    if x < 0.0:
        x = 0.0
    if x > 1.0:
        x = 1.0

    # Saturated properties at the cell pressure
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    k_l = s['k_l']
    Pr_l = s['Pr_l']

    # 1. Equivalent Reynolds number, Kwon Eq. (13) (Akers definition)
    density_ratio = rho_l / rho_v
    term_quality = (1.0 - x) + x * density_ratio ** 0.5
    Re_eq = G * Dh / mu_l * term_quality

    # Validity warnings (Table 10 range)
    if Re_eq < 3100 or Re_eq > 6200:
        warnings.warn(f"PCHE_kwon_condensing_HTC: Re_eq = {Re_eq:.0f} outside 3100-6200, "
                      "range of Kwon et al. (2018) Eq. (15): extrapolated.")
    if x < 0.68:
        warnings.warn(f"PCHE_kwon_condensing_HTC: x = {x:.2f} below 0.68, "
                      "the data cover only qualities 0.68-1.")
    p_red = p / p_crit
    if p_red < 0.08 or p_red > 0.11:
        warnings.warn(f"PCHE_kwon_condensing_HTC: reduced pressure p / p_crit = {p_red:.3f} "
                      "outside 0.08-0.11, range of Kwon et al. (2018) (CO2 at 0.7 MPa).")

    # 2. Kwon et al. (2018) Eq. (15) liquid Nusselt number (one factor per line)
    term_re = 1.8 * Re_eq ** 0.29
    term_pr = Pr_l ** 0.33
    Nu_l = term_re * term_pr

    # 3. Heat transfer coefficient, Nu_l = h Dh / k_l (Table 11)
    h = Nu_l * k_l / Dh

    return h


def PCHE_kwon_film_boiling_HTC(AS, geom, cond):
    """FILM FLOW BOILING heat transfer coefficient [W/m2/K] in a STRAIGHT semicircular PCHE
    channel, Sieder-Tate-type correlation of Kwon et al. (2018), Eqs. (7) and (11), fitted on
    liquid NITROGEN evaporated by nitrogen gas at about 290 K (very large wall superheat: film
    boiling, vapour layer on the wall). Explicit choice only: NOT a nucleate / convective
    boiling correlation.

    Reference:
        D. Kwon, L. Jin, W. Jung, S. Jeong (2018), Experimental investigation of heat transfer
        coefficient of mini-channel PCHE (printed circuit heat exchanger), Cryogenics 92, 41-49,
        doi:10.1016/j.cryogenics.2018.03.011, Eqs. (7), (11).

    Arguments:
        AS : CoolProp AbstractState of the fluid (moved to saturation, then to the wall state).
        geom : dict with keys
            'D_ch' [m] : channel diameter (used if 'Dh' not provided).
            'Dh'  [m] : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with operating conditions
            'G'   [kg/m2/s] : total mass flux of ONE channel.
            'p'   [Pa]      : pressure of the cell.
            'x'   [-]       : mean vapour quality of the cell.
            'T_w' [K]       : wall temperature estimate.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        Liquid nitrogen (cold side) at 0.49-0.71 MPa, G 27-60 kg/m2/s, inlet subcooled, outlet
        quality 0.25-0.61, heated by nitrogen gas at 285-296 K (film boiling regime assumed by
        the authors from the large temperature difference); 2100 < Re_tp < 2500 (Re_tp,m computed
        with HALF the outlet quality); predicts the 7 data points within 20 %; boiling h only
        81-122 W/m2/K. The wall temperature was not measured: the authors deduced it from
        q / A = h (T_w - T_sat) (Eq. 10). Exponents 0.33 (Pr) and 0.5 (viscosity ratio) fixed by
        the authors, 1.22e-5 and 1.89 fitted.
        ONE geometry only: STRAIGHT semicircular channels (short inlet / outlet bends), channel
        diameter 3 mm, hydraulic diameter 1.83 mm, pitch 4 mm, length 117 mm, 10 channels per
        layer, 1 hot + 1 cold layer, stainless steel, counter-flow.

    Definitions:
        Re_tp = G Dh / mu_v [x + (1 - x) rho_v / rho_l] (homogeneous two-phase Reynolds number on
        the vapour viscosity), Nu_v = h Dh / k_v; saturated vapour (v) and liquid (l) properties at
        the cell pressure; mu_w = vapour viscosity at the wall temperature T_w and the cell
        pressure; x = local (cell mean) quality; G of one channel (total mass flux).
        Returns h [W/m2/K] on the channel wall.
        (Independent check: with CoolProp, this Re_tp is 5-10 % above Table 8.)
    """
    import CoolProp.CoolProp as CP

    # Check that required properties are present in the condition dictionary
    _need(cond, 'G', 'p', 'x', 'T_w')

    # Hydraulic diameter, mass flux, pressure, quality and wall temperature
    Dh = _dh(geom)
    G = cond['G']
    p = cond['p']
    x = cond['x']
    T_w = cond['T_w']

    # Keep the quality inside [0, 1]
    if x < 0.0:
        x = 0.0
    if x > 1.0:
        x = 1.0

    # Saturated properties at the cell pressure
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_v = s['mu_v']
    k_v = s['k_v']
    Pr_v = s['Pr_v']
    T_sat = s['T_sat']

    # Vapour viscosity at the wall temperature (superheated vapour at the wall)
    if T_w > T_sat:
        AS.update(CP.PT_INPUTS, p, T_w)
        mu_w = AS.viscosity()
    else:
        warnings.warn("PCHE_kwon_film_boiling_HTC: no wall superheat (T_w <= T_sat): "
                      "film boiling correlation not applicable; mu_w = mu_v used.")
        mu_w = mu_v

    # 1. Homogeneous two-phase Reynolds number, Kwon Eq. (7)
    density_ratio = rho_v / rho_l
    Re_tp = G * Dh / mu_v * (x + (1 - x) * density_ratio)

    # Validity warning for Re_tp (narrow range of Eq. 11)
    if Re_tp < 2100 or Re_tp > 2500:
        warnings.warn(f"PCHE_kwon_film_boiling_HTC: Re_tp = {Re_tp:.0f} outside 2100-2500, "
                      "range of Kwon et al. (2018) Eq. 11: extrapolated.")

    # Validity warning for quality (data from 0 to 0.61)
    if x > 0.61:
        warnings.warn(f"PCHE_kwon_film_boiling_HTC: x = {x:.2f} above 0.61, "
                      "upper limit of Kwon et al. (2018) data: extrapolated.")

    # 2. Kwon et al. (2018) Eq. (11) vapour Nusselt number (one factor per line)
    term_re = 1.22e-5 * Re_tp ** 1.89
    term_pr = Pr_v ** 0.33
    term_mu = (mu_v / mu_w) ** 0.5
    Nu_v = term_re * term_pr * term_mu

    # 3. Heat transfer coefficient, Nu_v = h Dh / k_v (Table 9)
    h = Nu_v * k_v / Dh

    return h


def PCHE_katz_zigzag_HTC(AS, geom, cond):
    """Single-phase heat transfer coefficient [W/m2/K] in a ZIGZAG PCHE channel,
    Katz et al. (2021), Table 7 / Eq. (9): one correlation per fluid tested
    (sCO2 and helium fitted separately on the overall UA).

    Reference:
        A. Katz, S.R. Aakre, M.H. Anderson, et al. (2021), Experimental investigation of
        pressure drop and heat transfer in high temperature supercritical CO2 and helium in a
        printed-circuit heat exchanger, International Journal of Heat and Mass Transfer 171,
        121089, doi:10.1016/j.ijheatmasstransfer.2021.121089, Eqs. (9)-(12), Table 7.

    Arguments:
        AS : CoolProp AbstractState of the fluid (used for fluid name only).
        geom : dict with keys
            'alpha' [deg] : zigzag bending angle from the main flow direction
                            (labothappy convention; = fin angle / zigzag pitch angle).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        sCO2: Re 500-18000 (laminar to turbulent), Pr 0.66-1.25, up to 20 MPa, up to 550 C.
        Helium: Re 400-3200 (laminar / transition), Pr ~0.66 (absorbed in constant),
                up to 4 MPa, up to 550 C.
        102 heat transfer data points (Appendix A), all within +-15 %.
        The two Nu correlations were fitted jointly on the overall conductance with a
        continuity constraint at Re = 2300 (Eq. 12).
        ONE geometry only: 316L PCHE, zigzag pitch angle 37 deg (= labothappy alpha),
        wavelength 6.63 mm, corner radii 1.5 mm, semi-elliptical channel A_C 0.95 mm2,
        D_H = 4 A_C / P = 1 mm, 144 cold (sCO2) / 128 hot (helium) channels,
        developed length 98 mm (Table 3).

    Definitions:
        Nu = h * D_H / k, Re = G * D_H / mu, Pr = bulk Prandtl number.
        G is the mass flux of one channel; bulk properties of the cell.
        D_H is the hydraulic diameter of the channel (1 mm in the paper).
        AS is used only to identify the fluid (helium vs sCO2).
    """
    # Required bulk properties
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Hydraulic diameter and bending angle
    Dh = _dh(geom)
    alpha = geom['alpha']

    # Warn if geometry differs from the single tested geometry (alpha = 37 deg)
    if alpha != 37:
        warnings.warn(f"PCHE_katz_zigzag_HTC: geom['alpha'] = {alpha} deg, "
                      "correlation is for a single geometry with alpha = 37 deg.")

    # Extract properties
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Reynolds number on the hydraulic diameter
    Re = G * Dh / mu

    # Identify fluid (CoolProp name, lower case)
    fluid = AS.fluid_names()[0].lower()

    # Helium correlation: Nu = 0.8072 * Re^0.4359 (Pr absorbed in constant, Pr ~ 0.66)
    if fluid in ('helium',):
        # Validity warnings for helium
        if Re < 400 or Re > 3200:
            warnings.warn(f"PCHE_katz_zigzag_HTC: Re = {Re:.0f} outside 400-3200, "
                          "helium range of Katz et al. (2021): extrapolated.")
        if Pr < 0.6 or Pr > 0.7:
            warnings.warn(f"PCHE_katz_zigzag_HTC: Pr = {Pr:.2f} outside 0.6-0.7, "
                          "helium Prandtl term absorbed in the constant: extrapolated.")

        # Katz et al. (2021), Table 7, helium correlation (Eq. 9 for helium)
        Nu = 0.8072 * Re ** 0.4359

    # sCO2 correlation: Nu = 0.0261 * Re^0.8765 * Pr^(1/3)
    else:
        # Warn if fluid is not CO2 (the correlation was fitted on sCO2 data)
        if fluid not in ('carbondioxide', 'co2'):
            warnings.warn(f"PCHE_katz_zigzag_HTC: fluid = '{fluid}', "
                          "sCO2 correlation used (fitted on supercritical CO2 data).")

        # Validity warnings for sCO2
        if Re < 500 or Re > 18000:
            warnings.warn(f"PCHE_katz_zigzag_HTC: Re = {Re:.0f} outside 500-18000, "
                          "sCO2 range of Katz et al. (2021): extrapolated.")
        if Pr < 0.66 or Pr > 1.25:
            warnings.warn(f"PCHE_katz_zigzag_HTC: Pr = {Pr:.2f} outside 0.66-1.25, "
                          "sCO2 range of Katz et al. (2021): extrapolated.")

        # Katz et al. (2021), Table 7, sCO2 correlation (Eq. 9 for sCO2)
        Nu = 0.0261 * Re ** 0.8765 * Pr ** (1.0 / 3.0)

    # Heat transfer coefficient h = Nu * k / Dh
    h = Nu * k / Dh

    return h


def PCHE_kim_sg_zigzag_HTC(AS, geom, cond):
    """Single-phase sCO2 heat transfer coefficient in a zigzag PCHE channel,
    Kim S.G. et al. (2016), Table 12 (CFD-aided correlations for two fin angles).

    Reference:
        S.G. Kim, Y. Lee, Y. Ahn, J.I. Lee (2016), CFD aided approach to design
        printed circuit heat exchangers for supercritical CO2 Brayton cycle application,
        Annals of Nuclear Energy 92, 175-185, doi:10.1016/j.anucene.2016.01.019,
        Table 12.

    Arguments:
        AS : not used (common interface).
        geom : dict with keys
            'alpha' [deg] : zigzag bending angle from the main flow direction
                            (labothappy convention; = "fin angle" / "zigzag pitch angle"
                            of the papers).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        sCO2, 0.7 < Pr < 1.0.
        32.5 deg fit: 2000 < Re < 58000 (R^2 = 0.999).
        40 deg fit:   2000 < Re < 55000 (R^2 = 0.995).
        CFD (ANSYS CFX) validated on the Ishizuka et al. (2005) PCHE:
        semicircular channels width 1.9 / 1.8 mm, depth 0.9 mm, Dh 1.15 mm,
        fin angles 32.5 and 40 deg, sCO2 recuperator.
        Nu = h * Dh / k, channel-averaged CFD values, Re = G * Dh / mu.

    Definitions:
        Re = G * Dh / mu, Nu = h * Dh / k on the hydraulic diameter;
        G of one channel; bulk properties of the cell.
        The correlation has NO Prandtl number term (paper Table 12).
        Fit selection: alpha <= 36.25 deg (midpoint) -> 32.5 deg fit, else 40 deg fit.
        AS is not used (common interface).
    """
    # Required bulk properties
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Hydraulic diameter and bending angle
    Dh = _dh(geom)
    alpha = geom['alpha']

    # A zigzag needs 0 < alpha < 90 deg
    if not 0 < alpha < 90:
        raise ValueError(f"PCHE_kim_sg_zigzag_HTC: bending angle {alpha} deg outside ]0, 90[.")

    # Extract properties
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Reynolds number on the hydraulic diameter
    Re = G * Dh / mu

    # Midpoint between the two fitted angles (32.5 and 40 deg)
    midpoint = 36.25

    # Choose correlation based on bending angle
    if alpha <= midpoint:
        # 32.5 deg fit: Nu = 0.0292 * Re^0.8138
        fit_angle = 32.5
        coeff = 0.0292
        exp_re = 0.8138
        re_min = 2000
        re_max = 58000
    else:
        # 40 deg fit: Nu = 0.0188 * Re^0.8742
        fit_angle = 40.0
        coeff = 0.0188
        exp_re = 0.8742
        re_min = 2000
        re_max = 55000

    # Warn if alpha differs from the fit angle by more than 0.5 deg
    if abs(alpha - fit_angle) > 0.5:
        warnings.warn(f"PCHE_kim_sg_zigzag_HTC: fit for {fit_angle} deg used with "
                      f"alpha = {alpha} deg (difference > 0.5 deg).")

    # Warn if Reynolds number is outside the fit range
    if Re < re_min or Re > re_max:
        warnings.warn(f"PCHE_kim_sg_zigzag_HTC: Re = {Re:.0f} outside {re_min}-{re_max}, "
                      f"range of the {fit_angle} deg fit: extrapolated.")

    # Warn if Prandtl number is outside the validity range (0.7 < Pr < 1.0)
    if Pr < 0.7 or Pr > 1.0:
        warnings.warn(f"PCHE_kim_sg_zigzag_HTC: Pr = {Pr:.2f} outside 0.7-1.0, "
                      "validity range of Kim et al. (2016): extrapolated.")

    # Nusselt number from the chosen correlation (Table 12, no Pr term)
    Nu = coeff * Re ** exp_re

    # Heat transfer coefficient h = Nu * k / Dh
    h = Nu * k / Dh

    return h


def PCHE_baik_zigzag_HTC(AS, geom, cond):
    """Single-phase heat transfer coefficient [W/m2/K] in a zigzag printed circuit heat
    exchanger (PCHE) channel, Baik et al. (2017), Table 3: one correlation for water
    (laminar) and one for CO2 (turbulent, near the critical point).

    Reference:
        S. Baik, S.G. Kim, J. Lee, J.I. Lee (2017), Study on CO2-water printed circuit heat
        exchanger performance operating under various CO2 phases for S-CO2 power cycle
        application, Applied Thermal Engineering 113, 1536-1546,
        doi:10.1016/j.applthermaleng.2016.11.132, Table 3, Fig. 8.

    Arguments:
        AS : CoolProp AbstractState of the fluid (used to obtain the fluid name).
        geom : dict with keys
            'alpha' [deg] : bending angle of the zigzag from the main flow direction
                            (KAIST S-CO2PE precooler design: 32.5 deg).
            'D_ch' [m]   : channel diameter (optional, used to compute modified diameter).
            'Dh'  [m]   : optional hydraulic diameter (fallback if D_c not in geom).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        Water side: 50 < Re < 200 (laminar).
        CO2 side: 15000 < Re < 85000, Pr 2-33, 26-43 C, 7.3-8.6 MPa (liquid, supercritical
        and gas-like CO2 near the critical point).
        ONE geometry: KAIST S-CO2PE precooler (23 kW), zigzag 32.5 deg, wave length 9 mm,
        channel diameter 1.8 mm, corners rounded R 1 mm, SS316L; correlations from CFD
        (ANSYS CFX) validated on the tests, pitch-averaged (bend form losses included);
        headers excluded.

    Definitions:
        Baik uses a MODIFIED hydraulic diameter D_mod = sqrt(4 A / pi) (circle of the same
        flow area, 1.167 mm for their etched channel of 1.8 mm diameter, instead of
        4A/P = 0.973 mm) in BOTH Re = G * D_mod / mu and Nu = h * D_mod / k (and f).
        In this implementation: if 'D_ch' in geom, D_mod = geom['D_ch'] / sqrt(2.0)
        (ideal semicircle of diameter D_c: A = pi * D_c^2 / 8); otherwise geom['Dh'] is
        used with a warning that the modified diameter of Baik could not be computed.
        G is unchanged. Bulk properties of the cell; G of one channel.
    """
    # Read required bulk properties from condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Extract operating conditions
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Modified hydraulic diameter D_mod (Baik et al. 2017, Fig. 8: same flow area)
    if 'D_ch' in geom:
        D_mod = geom['D_ch'] / np.sqrt(2.0)
    else:
        D_mod = _dh(geom)
        warnings.warn(
            "PCHE_baik_zigzag_HTC: 'D_ch' not in geom, "
            "modified diameter of Baik could not be computed from D_c; using geom['Dh']."
        )

    # Reynolds number based on the modified hydraulic diameter
    Re = G * D_mod / mu

    # Warn if the bending angle differs from the tested geometry (32.5 deg)
    alpha = geom.get('alpha')
    if alpha != 32.5:
        warnings.warn(
            f"PCHE_baik_zigzag_HTC: bending angle alpha = {alpha} deg differs from "
            "the 32.5 deg geometry of Baik et al. (2017)."
        )

    # Fluid name from CoolProp AbstractState
    fluid = AS.fluid_names()[0].lower()

    # Select correlation based on the fluid name
    if fluid == 'water':
        # Water side correlation (Table 3, laminar water flow)
        Nu = 0.2829 * Re ** 0.6686

        # Validity check for water side Reynolds number
        if Re < 50 or Re > 200:
            warnings.warn(
                f"PCHE_baik_zigzag_HTC: Re = {Re:.1f} outside the water validity range "
                "50-200 of Baik et al. (2017): extrapolated."
            )

    else:
        # CO2 / other fluid branch (Table 3, turbulent near-critical CO2)
        Nu = 0.8405 * Re ** 0.5704 * Pr ** 1.08

        # Warn if fluid is not carbon dioxide
        if fluid not in ('carbondioxide',):
            warnings.warn(
                f"PCHE_baik_zigzag_HTC: fluid '{fluid}' is not 'carbondioxide' "
                "(correlations developed for CO2)."
            )

        # Validity checks for CO2 side
        if Re < 15000 or Re > 85000:
            warnings.warn(
                f"PCHE_baik_zigzag_HTC: Re = {Re:.0f} outside the CO2 validity range "
                "15000-85000 of Baik et al. (2017): extrapolated."
            )

        if Pr < 2 or Pr > 33:
            warnings.warn(
                f"PCHE_baik_zigzag_HTC: Pr = {Pr:.2f} outside the CO2 validity range "
                "2-33 of Baik et al. (2017): extrapolated."
            )

    # Heat transfer coefficient h = Nu * k / D_mod
    h = Nu * k / D_mod

    return h


def PCHE_saeed_zigzag_HTC(AS, geom, cond):
    """Single-phase heat transfer coefficient [W/m2/K] in a zigzag Printed Circuit Heat
    Exchanger (PCHE) channel, Saeed et al. (2020), Table 6.

    Reference:
        M. Saeed, A.S. Berrouk, M.S. Siddiqui, A.A. Awais (2020), Numerical investigation of
        thermal and hydraulic characteristics of sCO2-water printed circuit heat exchangers
        with zigzag channels, Energy Conversion and Management 224, 113375,
        doi:10.1016/j.enconman.2020.113375, Eqs. (7)-(10), Table 6.

    Arguments:
        AS : CoolProp AbstractState of the fluid (used to determine the fluid name).
        geom : dict with keys
            'alpha' [deg] : fin angle / zigzag bending angle from the main flow direction.
            'D_ch' [m]     : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]     : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number.

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        Table 6: sCO2 3000 <= Re <= 60000, 2 <= Pr <= 13 (the text of Section 3 says 5000-70000
        and 1-13: the table is used); water 50 <= Re <= 15000; about 94 % of the CFD points within
        +-10 %. CFD (72 periodic units, segment-averaged data), ONE geometry: fin angle 40 deg
        (= labothappy alpha), pitch along length 9 mm, transverse pitch 2.35 mm, Dh 1.106 mm,
        SS316L, sCO2 precooler cooled by water (Table 2).

    Definitions:
        Re = G * Dh / mu, Nu = h * Dh / k, segment-averaged. Bulk properties of the cell,
        G of one channel.
    """
    # Read required bulk properties and operating conditions
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Hydraulic diameter and bending angle
    Dh = _dh(geom)
    alpha = geom.get('alpha', 0.0)

    # Extract conditions
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # 1. Reynolds number based on hydraulic diameter
    Re = G * Dh / mu

    # Determine fluid name from AbstractState
    fluid = AS.fluid_names()[0].lower()

    # 2. Warn if geometry bending angle is not 40 deg
    if alpha != 40:
        warnings.warn(f"PCHE_saeed_zigzag_HTC: fin angle alpha = {alpha} deg (Saeed et al. 2020 "
                      "correlations were derived for 40 deg).")

    # 3. Fluid-dependent Nusselt number correlation
    if fluid == 'water':
        # Water correlation: Saeed et al. (2020), Table 6
        if Re < 50 or Re > 15000:
            warnings.warn(f"PCHE_saeed_zigzag_HTC: Re = {Re:.0f} outside the validity range "
                          "50-15000 for water (Saeed et al. 2020).")
        Nu = 0.08 * Re ** 0.77

    else:
        # sCO2 correlation (also used, with a warning, for any other fluid)
        if fluid not in ('carbondioxide', 'co2'):
            warnings.warn(f"PCHE_saeed_zigzag_HTC: fluid '{fluid}' not tested by Saeed et al. "
                          "(2020) (water or sCO2): sCO2 correlation used.")
        if Re < 3000 or Re > 60000:
            warnings.warn(f"PCHE_saeed_zigzag_HTC: Re = {Re:.0f} outside the validity range "
                          "3000-60000 for sCO2 (Saeed et al. 2020).")
        if Pr < 2 or Pr > 13:
            warnings.warn(f"PCHE_saeed_zigzag_HTC: Pr = {Pr:.2f} outside the validity range "
                          "2-13 for sCO2 (Saeed et al. 2020).")
        # Saeed et al. (2020) Table 6, sCO2
        Nu = 0.475 * Re ** 0.61 * Pr ** 0.17

    # 4. Heat transfer coefficient h = Nu * k / Dh
    h = Nu * k / Dh

    return h


def PCHE_chen_zigzag_HTC(AS, geom, cond):
    """Single-phase heat transfer coefficient [W/m2/K] of HELIUM in a ZIGZAG semicircular PCHE
    channel (laminar and transition), Chen et al. (2016), Eq. (12) (two branches).

    Reference:
        M. Chen, X. Sun, R.N. Christensen, I. Skavdahl, V. Utgikar, P. Sabharwall (2016),
        Pressure drop and heat transfer characteristics of a high-temperature printed circuit
        heat exchanger, Applied Thermal Engineering 108, 1409-1417,
        doi:10.1016/j.applthermaleng.2016.07.149, Section 4.2, Eq. (12).

    Arguments:
        AS : CoolProp AbstractState (used only for the fluid name).
        geom : dict with keys
            'alpha' [deg] : zigzag bending angle from the main flow direction
                            (labothappy convention).
            'D_ch' [m]   : channel diameter (used if 'Dh' not provided).
            'Dh'  [m]   : optional hydraulic diameter; default pi*D_c/(2+pi).
        cond : dict with bulk properties of the cell (single-phase)
            'G'   [kg/m2/s] : mass flux in ONE channel.
            'mu'  [Pa s]    : dynamic viscosity.
            'k'   [W/m/K]   : thermal conductivity.
            'Pr'  [-]       : Prandtl number (not used in the correlation; checked by _need).

    Returns:
        h [W/m2/K] : heat transfer coefficient on the channel wall.

    Validity:
        Helium on both sides, 1400 <= Re <= 3558, 82 data points (67 used for the fit);
        transition onset Re about 2200.
        ONE geometry only (Alloy 617, 8 hot + 8 cold plates, 11 channels each):
        semicircular channels nominal diameter about 2 mm, span-wise pitch 2.5 mm,
        zigzag angle 15 deg from the flow direction (= labothappy alpha),
        bend radius of curvature 4 mm; helium on both sides, up to 464 / 802 C and 2.7 MPa,
        flow up to 39 kg/h.

    Definitions:
        Re = G * Dh / mu, Nu = h * Dh / k on the hydraulic diameter; G of one channel;
        bulk properties of the cell. No Prandtl number dependence in the correlation.
        AS is used only to check that the fluid is helium.
    """
    # Check that required bulk properties are present in the condition dictionary
    _need(cond, 'G', 'mu', 'k', 'Pr')

    # Get the hydraulic diameter and bending angle from the geometry dictionary
    Dh = _dh(geom)
    alpha = geom['alpha']

    # Extract properties for the calculations
    G = cond['G']
    mu = cond['mu']
    k = cond['k']

    # 1. Reynolds number on the hydraulic diameter
    Re = G * Dh / mu

    # 2. Fluid name check (CoolProp name, lower case)
    fluid = AS.fluid_names()[0].lower()
    if fluid not in ('helium',):
        warnings.warn(f"PCHE_chen_zigzag_HTC: fluid is '{fluid}', correlation is for helium only.")

    # 3. Validity warning for Reynolds number (Chen et al. 2016 range: 1400 to 3558)
    if Re < 1400 or Re > 3558:
        warnings.warn(f"PCHE_chen_zigzag_HTC: Re = {Re:.0f} outside 1400-3558, "
                      "range of Chen et al. (2016) Eq. (12): extrapolated.")

    # 4. Validity warning for bending angle (tested geometry has alpha = 15 deg)
    if alpha != 15:
        warnings.warn(f"PCHE_chen_zigzag_HTC: Chen et al. (2016) tested a single geometry, "
                      f"zigzag angle 15 deg (used here with alpha = {alpha} deg).")

    # 5. Nusselt number from Chen et al. (2016), Eq. (12) - two branches
    if Re <= 2200:
        # Laminar branch: 1400 <= Re <= 2200
        Nu = 0.05516 * Re ** 0.69195
    else:
        # Transition branch: 2200 < Re <= 3558
        Nu = 0.09221 * Re ** 0.62507

    # 6. Heat transfer coefficient h = Nu * k / Dh
    h = Nu * k / Dh

    return h