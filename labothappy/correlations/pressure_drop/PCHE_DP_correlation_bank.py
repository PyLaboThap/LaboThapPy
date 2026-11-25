# -*- coding: utf-8 -*-
"""
PCHE (printed circuit heat exchanger) channel-side frictional pressure drop, in the
(AS, geom, cond) formalism of the other pressure drop banks.

Each function returns the FRICTIONAL pressure drop [Pa] over the full developed channel length
geom['L'] (no gravity, no acceleration, no ports / headers).

geom : Dh [m] (hydraulic diameter of the semicircular channel, pi D_c / (2 + pi)), L [m] (developed
       channel length), alpha [deg] (zigzag bending angle from the main flow direction)
cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho or D [kg/m3]
       (bulk properties of the cell)
"""

import warnings
import numpy as np
# Saturated properties for the two-phase PCHE pressure drop (Bae et al. 2019)
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props
from labothappy.correlations.convection.correlation_warnings import (
    AssumedParameterWarning, FluidWarning, ValidityRangeWarning)

#%% SINGLE PHASE

def PCHE_lee_jiang_zigzag_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop in a zigzag semicircular channel of a PCHE [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - Darcy f used in DP (Lee's Fanning f multiplied by 4).
    - Re = G * Dh / mu, with Dh = pi * D_c / (2 + pi) (Jiang Eq. 8).
    - Bulk properties of the cell.
    - G of one channel.
    - L = developed length (length along the zigzag, Jiang Eq. 1).
    - Returns frictional DP [Pa] over L.

    References:
    - Y. Lee, S. Lee, H.B. Park, K.W. Seo, D. Lee, Thermal-hydraulic performance analysis of zigzag
      channel PCHE according to bending angle, Trans. Korean Nuclear Society Spring Meeting, Jeju,
      May 9-10, 2024, Eq. (2). Friction factor is FANNING.
    - Y. Jiang, E. Liese, S.E. Zitney, D. Bhattacharyya, Design and dynamic modeling of printed
      circuit heat exchangers for supercritical carbon dioxide Brayton power cycles (2018),
      Eqs. (22b), (23b), (7). Friction factor is DARCY.

    Validity range and fluids tested:
    - Lee: 300 <= Re <= 2000, 10 <= alpha <= 40 deg, water (CFD), D_c = 1.5 mm, f +-20 % vs
      Kim & No (2011).
    - Jiang: 3000 <= Re <= 20600, alpha 30 and 40 deg, sCO2, D_c = 2 mm (fit of Heatric data).
    - Between Re 2000 and 3000: linear blend (no data).

    geom : Dh [m], L [m], alpha [deg] (bending angle of the zigzag from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # A zigzag needs 0 < alpha < 90 deg: otherwise raise ValueError
    if not (0 < alpha < 90):
        raise ValueError(f"PCHE_lee_jiang_zigzag_SP_DP: bending angle {alpha} deg outside ]0, 90[.")

    # Reynolds number on hydraulic diameter
    Re = G * Dh / mu

    # Amplitude over pitch of the zigzag, from bending angle (Table I of Lee)
    h_over_p = np.tan(np.radians(alpha)) / 2

    # Laminar friction, Lee et al. (2024), Eq. (2), used for Re <= 2000
    if Re <= 2000:
        # Warn if Re < 300 (lower limit of Lee)
        if Re < 300:
            warnings.warn(
                f"PCHE_lee_jiang_zigzag_SP_DP: Re = {Re:.0f} below 300, lower limit of "
                "Lee et al. (2024): extrapolated.", ValidityRangeWarning
            )

        # Warn if alpha outside 10-40 deg (range of Lee's CFD)
        if not (10 <= alpha <= 40):
            warnings.warn(
                f"PCHE_lee_jiang_zigzag_SP_DP: bending angle {alpha} deg outside the range "
                "10-40 deg of Lee et al. (2024): extrapolated.", ValidityRangeWarning
            )

        # Lee et al. (2024) Eq. (2) for Fanning friction factor
        f_fanning = 29.1955 * (Re ** -0.5223) * (h_over_p ** 0.864)

        # Convert Fanning to Darcy friction factor
        f_darcy = 4.0 * f_fanning

    # Turbulent friction, Jiang et al. (2018), used for Re >= 3000
    elif Re >= 3000:
        # Warn if Re > 20600 (upper limit of Jiang's fit)
        if Re > 20600:
            warnings.warn(
                f"PCHE_lee_jiang_zigzag_SP_DP: Re = {Re:.0f} above 20600, upper limit of "
                "Jiang et al. (2018): extrapolated.", ValidityRangeWarning
            )

        # Warn if alpha is not 30 or 40 deg
        if alpha != 30 and alpha != 40:
            warnings.warn(
                f"PCHE_lee_jiang_zigzag_SP_DP: Jiang et al. (2018) gives fits for 30 and 40 deg "
                f"channels only, used here with alpha = {alpha} deg.", ValidityRangeWarning
            )

        # Select low angle (<= 35 deg) or high angle (> 35 deg) correlation
        if alpha <= 35:
            # Jiang et al. (2018) Eq. (22b), low angle (30 deg) channel Darcy friction factor
            f_darcy = 0.3905 * (Re ** -0.0355)
        else:
            # Jiang et al. (2018) Eq. (23b), high angle (40 deg) channel Darcy friction factor
            f_darcy = 1.336 * (Re ** -0.1268)

    # Transition region 2000 < Re < 3000: linear interpolation in Re (no warnings)
    else:
        # Calculate Darcy friction factor at Re = 2000 using Lee et al. (2024)
        f_fanning_2000 = 29.1955 * (2000.0 ** -0.5223) * (h_over_p ** 0.864)
        f_low = 4.0 * f_fanning_2000

        # Calculate Darcy friction factor at Re = 3000 using Jiang et al. (2018)
        if alpha <= 35:
            f_high = 0.3905 * (3000.0 ** -0.0355)
        else:
            f_high = 1.336 * (3000.0 ** -0.1268)

        # Linear interpolation weight
        weight = (Re - 2000.0) / (3000.0 - 2000.0)

        # Interpolated Darcy friction factor
        f_darcy = f_low + weight * (f_high - f_low)

    # Frictional pressure drop over developed length L (Jiang Eq. 7, Darcy)
    term_friction = f_darcy * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_ngo_zigzag_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop in a zigzag channel of a PCHE [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - FANNING friction factor f of Eq. (31) (the "pressure-drop factor" of the paper is defined
      by Eq. (23), DP = 2 L W^2 / (Dh Ac^2 rho) f); Darcy = 4 x Fanning.
    - Re = G * Dh / mu, on the hydraulic diameter (Eq. 17).
    - Bulk properties of the cell.
    - Mass flux G of one channel.
    - L: "active flow length" in the paper; taken here as the developed channel length geom['L'].
      (L back-computed from the measured DP of the paper's Table 4 is 1.0-1.2 m, longer than the
      0.745 m block: L is not the axial block length.)
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no ports).

    Reference:
    - T.L. Ngo, Y. Kato, K. Nikitin, T. Ishizuka (2007), Heat transfer and
      pressure drop correlations of microchannel heat exchangers with S-shaped
      and zigzag fins for carbon dioxide cycles, Experimental Thermal and Fluid
      Science 32 (2007) 560-570, doi:10.1016/j.expthermflusci.2007.06.006,
      Eqs. (23), (31) (zigzag fins) and (17).

    Validity range and fluids tested:
    - 3500 < Re < 22000, supercritical CO2 (experiments). No Pr dependence of f found
      (Pr 0.75-2.2).
    - Single geometry tested: fin angle 52 deg, Dh = 1.09 mm, fin gap 1.31 mm,
      fin depth 0.94 mm, fin width 0.8 mm.
    - Measured DP reproduced with a standard deviation of +-13.5 %.
    - Inlet / outlet (distributor) losses neglected in the fit (paper, Eq. 23).
    - Note: Ngo's channels are etched "zigzag fins" (rectangular-like
      cross-section), not semicircular.

    geom : Dh [m], L [m], alpha [deg] (bending angle of the zigzag from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # A zigzag needs 0 < alpha < 90 deg: otherwise raise ValueError
    if not (0 < alpha < 90):
        raise ValueError(
            f"PCHE_ngo_zigzag_SP_DP: bending angle {alpha} deg outside ]0, 90[."
        )

    # Reynolds number on hydraulic diameter (Ngo Eq. 17)
    Re = G * Dh / mu

    # Check Reynolds number range (Ngo et al. Eq. 31 validity)
    if Re < 3500 or Re > 22000:
        warnings.warn(
            f"PCHE_ngo_zigzag_SP_DP: Re = {Re:.0f} outside validity range "
            "3500 - 22000 of Ngo et al. (2007): extrapolated.", ValidityRangeWarning
        )

    # Check fin angle against the single tested geometry
    if alpha != 52:
        warnings.warn(
            f"PCHE_ngo_zigzag_SP_DP: Ngo et al. (2007) tested a single geometry, "
            f"fin angle 52 deg, used here with alpha = {alpha} deg.", ValidityRangeWarning
        )

    # Fanning pressure-drop factor, Ngo et al. (2007) Eq. (31) for zigzag fins
    f_fanning = 0.1924 * (Re ** -0.091)

    # Convert Fanning friction factor to Darcy friction factor
    f_darcy = 4.0 * f_fanning

    # Frictional pressure drop, Ngo et al. (2007) Eq. (23)
    term_friction = f_darcy * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_ngo_S_fin_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop in an S-shaped fin channel of a PCHE [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - FANNING friction factor f of Eq. (30) (the "pressure-drop factor" of the paper is defined
      by Eq. (23), DP = 2 L W^2 / (Dh Ac^2 rho) f); Darcy = 4 x Fanning.
    - Re = G * Dh / mu, on the hydraulic diameter (Eq. 17).
    - Bulk properties of the cell.
    - Mass flux G of one channel.
    - L: "active flow length" in the paper; taken here as the flow-path length geom['L'].
      (L back-computed from the measured DP of the paper's Table 4 (S-fin exchanger) is
      1.0-1.2 m, and Table 1 area 0.5099 m2 / (96 channels x wetted perimeter 4.5 mm) gives
      1.18 m: both longer than the 0.745 m block, so L is not the axial block length.)
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no ports).

    Reference:
    - T.L. Ngo, Y. Kato, K. Nikitin, T. Ishizuka (2007), Heat transfer and
      pressure drop correlations of microchannel heat exchangers with S-shaped
      and zigzag fins for carbon dioxide cycles, Experimental Thermal and Fluid
      Science 32 (2007) 560-570, doi:10.1016/j.expthermflusci.2007.06.006,
      Eqs. (23), (30) (S-shaped fins) and (17).

    Validity range and fluids tested:
    - 3500 < Re < 23000, supercritical CO2 (experiments). No Pr dependence of f found
      (Pr 0.75-2.2).
    - Single geometry tested: fin angle 52 deg, Dh = 1.09 mm, fin gap 1.31 mm,
      fin depth 0.94 mm, fin width 0.8 mm.
    - Measured DP reproduced with a standard deviation of +-16.6 %.
    - Inlet / outlet (distributor) losses neglected in the fit (paper, Eq. 23).
    - Note: discontinuous S-shaped fins (rectangular-like passage between fins),
      not semicircular channels.

    geom : Dh [m], L [m], alpha [deg] (fin angle from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # The fin angle must be between 0 and 90 degrees: otherwise raise ValueError
    if not (0 < alpha < 90):
        raise ValueError(
            f"PCHE_ngo_S_fin_SP_DP: fin angle {alpha} deg outside ]0, 90[."
        )

    # Reynolds number on hydraulic diameter, Ngo et al. (2007) Eq. (17)
    Re = G * Dh / mu

    # Check Reynolds number range against Eq. (30) validity range
    if Re < 3500 or Re > 23000:
        warnings.warn(
            f"PCHE_ngo_S_fin_SP_DP: Re = {Re:.0f} outside validity range "
            "3500 - 23000 of Ngo et al. (2007): extrapolated.", ValidityRangeWarning
        )

    # Check fin angle against the single tested geometry
    if alpha != 52:
        warnings.warn(
            f"PCHE_ngo_S_fin_SP_DP: Ngo et al. (2007) tested a single geometry, "
            f"fin angle 52 deg, used here with alpha = {alpha} deg.", ValidityRangeWarning
        )

    # Fanning pressure-drop factor, Ngo et al. (2007) Eq. (30) for S-shaped fins
    f_fanning = 0.4545 * (Re ** -0.340)

    # Convert Fanning friction factor to Darcy friction factor
    f_darcy = 4.0 * f_fanning

    # Frictional pressure drop over flow-path length L, Ngo et al. (2007) Eq. (23)
    term_friction = f_darcy * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_kim_no_2012_zigzag_SP_DP(AS, geom, cond):
    """Single-phase LAMINAR frictional pressure drop in a zigzag semicircular channel of a PCHE
    [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - f is FANNING friction factor (paper Eq. 1: f = DP * Dh / (2 * L * rho * v^2), with L the
      channel length of a pitch); Darcy = 4 x Fanning.
    - Re = G * Dh / mu on the hydraulic diameter of the semicircular channel.
    - Bulk properties of the cell (the paper used local pitch-averaged properties).
    - Mass flux G of one channel.
    - alpha = angle of the zigzag from the main flow direction (paper Eq. 9: L = L_c * cos(alpha)).
    - P_zz: pitch of the zigzag as in Kim & No Fig. 2 (24.6 mm for the tested PCHE).
    - Interpolation in angle and nearest-geometry choice: added in labothappy (not in the paper).
    - Returns the FRICTIONAL pressure drop [Pa] over L (no gravity, no acceleration, no ports).

    Reference:
    - I.H. Kim, H.C. No (2012), Physical model development and optimal design of PCHE for
      intermediate heat exchangers in HTGRs, Nuclear Engineering and Design 243, 243-250,
      doi:10.1016/j.nucengdes.2011.11.020, Eqs. (1)-(3), Tables 3 and 4.

    Validity range and fluids tested:
    - Laminar flow, Re up to about 2500 (range of the CFD data, Figs. 4 and 5; the paper gives
      no explicit Re range).
    - Helium (FLUENT CFD, local pitch-averaged values; the model was validated against the helium
      tests of Kim et al. 2009, 15 deg, Dh 0.922 mm, pitch 24.6 mm).
    - Angle 5-45 deg, Dh 0.922 or 1.222 mm, pitch 12.3 or 24.6 mm.
    - RMS error of each fit 0.8-3 % (Table 4).

    geom : Dh [m], L [m] (developed channel length), alpha [deg] (bending angle of the zigzag
           from the main flow direction), P_zz [m] (optional)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Optional zigzag pitch, fallback to 24.6 mm with warning
    pitch = geom.get('P_zz')
    if pitch is None:
        pitch = 24.6e-3
        warnings.warn(
            "PCHE_kim_no_2012_zigzag_SP_DP: geom['P_zz'] not given: 24.6 mm assumed.",
            AssumedParameterWarning
        )

    # A zigzag needs 0 < alpha < 90 deg: otherwise raise ValueError
    if not (0 < alpha < 90):
        raise ValueError(
            f"PCHE_kim_no_2012_zigzag_SP_DP: bending angle {alpha} deg outside ]0, 90[."
        )

    # Reynolds number on hydraulic diameter, Kim & No Eq. (2)
    Re = G * Dh / mu

    # Warn if Re is above the laminar CFD range of Kim & No (2012)
    if Re > 2500:
        warnings.warn(
            f"PCHE_kim_no_2012_zigzag_SP_DP: Re = {Re:.0f} above the laminar CFD range of "
            "Kim & No (2012) (Re up to 2500, Figs. 4-5): extrapolated.", ValidityRangeWarning
        )

    # Table data from Kim & No (2012), Table 4
    angles = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]

    # Table (a): Dh = 0.922 mm, pitch = 24.6 mm
    Dh_a = 0.922e-3
    pitch_a = 24.6e-3
    a_a = [0.0019, 0.01775, 0.06455, 0.08918, 0.1442, 0.21115, 0.29342, 0.39158, 0.51539]
    b_a = [1.1187, 0.90795, 0.81021, 0.8136, 0.79529, 0.78615, 0.78118, 0.7799, 0.77944]

    # Table (b): Dh = 0.922 mm, pitch = 12.3 mm
    Dh_b = 0.922e-3
    pitch_b = 12.3e-3
    a_b = [0.00366, 0.02536, 0.0696, 0.12817, 0.19392, 0.27701, 0.37934, 0.49995, 0.65898]
    b_b = [1.07165, 0.91728, 0.85362, 0.83085, 0.82838, 0.82549, 0.82413, 0.82555, 0.82508]

    # Table (c): Dh = 1.222 mm, pitch = 24.6 mm
    Dh_c = 1.222e-3
    pitch_c = 24.6e-3
    a_c = [0.0034, 0.02342, 0.06677, 0.12748, 0.17458, 0.25418, 0.35159, 0.46869, 0.62339]
    b_c = [1.0502, 0.8863, 0.81258, 0.78473, 0.79345, 0.7867, 0.78015, 0.78205, 0.78214]

    # Relative distance to each table (rule of labothappy, not of the paper)
    dist_a = ((Dh - Dh_a) / Dh_a) ** 2 + ((pitch - pitch_a) / pitch_a) ** 2
    dist_b = ((Dh - Dh_b) / Dh_b) ** 2 + ((pitch - pitch_b) / pitch_b) ** 2
    dist_c = ((Dh - Dh_c) / Dh_c) ** 2 + ((pitch - pitch_c) / pitch_c) ** 2

    # Choose the table with the smallest distance
    if dist_a <= dist_b and dist_a <= dist_c:
        a_list = a_a
        b_list = b_a
        chosen_Dh = Dh_a
        chosen_pitch = pitch_a
    elif dist_b <= dist_a and dist_b <= dist_c:
        a_list = a_b
        b_list = b_b
        chosen_Dh = Dh_b
        chosen_pitch = pitch_b
    else:
        a_list = a_c
        b_list = b_c
        chosen_Dh = Dh_c
        chosen_pitch = pitch_c

    # Warn if the chosen table geometry deviates significantly from inputs
    if abs(Dh - chosen_Dh) / chosen_Dh > 0.1 or abs(pitch - chosen_pitch) / chosen_pitch > 0.1:
        warnings.warn(
            "PCHE_kim_no_2012_zigzag_SP_DP: geometry outside the three CFD geometries of "
            "Kim & No (2012): nearest table used.", ValidityRangeWarning
        )

    # Handle angle limits and warnings
    if alpha < 5.0:
        warnings.warn(
            f"PCHE_kim_no_2012_zigzag_SP_DP: angle {alpha} deg outside 5-45 deg, "
            "value at the nearest table angle (5 deg) used.", ValidityRangeWarning
        )
        alpha_used = 5.0
    elif alpha > 45.0:
        warnings.warn(
            f"PCHE_kim_no_2012_zigzag_SP_DP: angle {alpha} deg outside 5-45 deg, "
            "value at the nearest table angle (45 deg) used.", ValidityRangeWarning
        )
        alpha_used = 45.0
    else:
        alpha_used = alpha

    # Find the two neighbouring table angles for linear interpolation
    alpha_1 = angles[0]
    alpha_2 = angles[-1]
    index_1 = 0
    index_2 = len(angles) - 1

    for i in range(len(angles) - 1):
        if angles[i] <= alpha_used <= angles[i + 1]:
            alpha_1 = angles[i]
            alpha_2 = angles[i + 1]
            index_1 = i
            index_2 = i + 1
            break

    # Extract coefficients for the two neighbouring angles
    a_1 = a_list[index_1]
    b_1 = b_list[index_1]
    a_2 = a_list[index_2]
    b_2 = b_list[index_2]

    # Compute f*Re for the two angles using Kim & No Eq. (3)
    # Note: 15.78 is the f*Re of a straight semicircular duct, Table 3
    fRe_1 = 15.78 + a_1 * (Re ** b_1)
    fRe_2 = 15.78 + a_2 * (Re ** b_2)

    # Interpolate linearly in angle
    if alpha_1 == alpha_2:
        fRe = fRe_1
    else:
        fRe = fRe_1 + (fRe_2 - fRe_1) * (alpha_used - alpha_1) / (alpha_2 - alpha_1)

    # Fanning friction factor, f = (f*Re) / Re
    f = fRe / Re

    # Frictional pressure drop, Kim & No Eq. (1) / Eq. (16)
    term_friction = 4.0 * f * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_kim_no_2011_zigzag_SP_DP(AS, geom, cond):
    """Single-phase LAMINAR frictional pressure drop in a zigzag semicircular channel of a PCHE
    [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - FANNING friction factor f of Eq. (19); Darcy = 4 x Fanning.
    - Re = G * Dh / mu, on the hydraulic diameter.
    - Bulk properties of the cell (local pitch-averaged properties in the paper).
    - Mass flux G of one channel.
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no ports).
    - Note on viscosity correction: The viscosity correction of Eq. (21)
      (C.F = 10.939 - 11.014 nu_s / nu_b, fitted on the water side only) is NOT implemented
      here. According to the authors, Eq. (19) is conservative (it overestimates
      DP where the wall-bulk temperature difference is large).

    Reference:
    - I.H. Kim, H.C. No (2011), Thermal hydraulic performance analysis of a printed circuit heat
      exchanger using a helium-water test loop and numerical simulations, Applied Thermal
      Engineering 31, 4064-4073, doi:10.1016/j.applthermaleng.2011.08.012, Eq. (19).

    Validity range and fluids tested:
    - Laminar flow, Re up to 2500 (experiments helium 350 < Re < 1200 in Kim et al. 2009,
      CFD extended to 2500).
    - Helium and water.
    - Single geometry tested: angle 15 deg, channel diameter 1.51 mm, Dh 0.922 mm,
      pitch 24.6 mm, wavy channel length 765 mm.

    geom : Dh [m], L [m], alpha [deg] (bending angle of the zigzag from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # A zigzag needs 0 < alpha < 90 deg: otherwise raise ValueError
    if not (0 < alpha < 90):
        raise ValueError(
            f"PCHE_kim_no_2011_zigzag_SP_DP: bending angle {alpha} deg outside ]0, 90[."
        )

    # Reynolds number on hydraulic diameter
    Re = G * Dh / mu

    # Check Reynolds number range against Kim & No (2011) validity range
    if Re > 2500:
        warnings.warn(
            f"PCHE_kim_no_2011_zigzag_SP_DP: Re = {Re:.0f} above the range "
            "of Kim & No (2011) (CFD up to Re 2500).", ValidityRangeWarning
        )

    # Check bending angle against the single tested geometry
    if alpha != 15:
        warnings.warn(
            f"PCHE_kim_no_2011_zigzag_SP_DP: Kim & No (2011) studied a single geometry, "
            f"angle 15 deg (Dh 0.922 mm, pitch 24.6 mm), used here with alpha = {alpha} deg.",
            ValidityRangeWarning
        )

    # Fanning friction factor times Reynolds number, Kim & No (2011) Eq. (19)
    # (15.78 is the f*Re of a straight semicircular duct)
    fRe = 15.78 + 0.004868 * (Re ** 0.8416)

    # Fanning friction factor
    f_fanning = fRe / Re

    # Frictional pressure drop over developed length L
    term_friction = 4.0 * f_fanning * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_liu_straight_SP_DP(AS, geom, cond):
    """Single-phase TURBULENT frictional pressure drop in a STRAIGHT semicircular channel of a
    PCHE [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - f is FANNING friction factor (paper nomenclature); Darcy = 4 x Fanning.
    - Re on the equivalent (hydraulic) diameter, bulk properties; G of one channel.
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no ports).

    Reference:
    - S.-h. Liu, Y.-p. Huang, J.-f. Wang, et al. (2020), Experimental study of thermal-hydraulic
      performance of a printed circuit heat exchanger with straight channels, International
      Journal of Heat and Mass Transfer 160, 120109, doi:10.1016/j.ijheatmasstransfer.2020.120109,
      Eqs. (8), (13).

    Validity range and fluids tested:
    - 3600 < Re_b < 36500, supercritical CO2 cooled by water, p 7.53-11.97 MPa,
      bulk temperature 26.4-121.3 C, G 50-780 kg/m2/s, Pr_b 0.9-10.6 (Table 3);
      straight semicircular channels, channel diameter 1.87 mm, D_eq 1.14 mm, length 500 mm;
      100 % of the data within +-20 %.

    geom : Dh [m], L [m]
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Reynolds number on hydraulic diameter, Liu et al. (2020) Eq. (8)
    Re = G * Dh / mu

    # Check Reynolds number range against Eq. (13) validity range
    if Re < 3600 or Re > 36500:
        warnings.warn(
            f"PCHE_liu_straight_SP_DP: Re = {Re:.0f} outside validity range "
            "3600 - 36500 of Liu et al. (2020): extrapolated.", ValidityRangeWarning
        )

    # Fanning friction factor, Liu et al. (2020) Eq. (13)
    f_fanning = 0.05776 * (Re ** -0.2192)

    # Frictional pressure drop over length L
    term_friction = 4.0 * f_fanning * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_meshram_zigzag_SP_DP(AS, geom, cond):
    """Single-phase TURBULENT frictional pressure drop in a zigzag semicircular channel
    of a printed circuit heat exchanger (PCHE) [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - f is FANNING (Eq. 9); Darcy = 4 x Fanning.
    - Re = G * Dh / mu, on the hydraulic diameter, bulk properties.
    - G is the mass flux of ONE channel.
    - Returns the FRICTIONAL pressure drop [Pa] over L (no gravity, no acceleration, no ports).

    Reference:
    - A. Meshram, A.K. Jaiswal, S.D. Khivsara, J.D. Ortega, C. Ho, R. Bapat, P. Dutta (2016),
      Modeling and analysis of a printed circuit heat exchanger for supercritical CO2 power
      cycle applications, Applied Thermal Engineering 109, 861-870,
      doi:10.1016/j.applthermaleng.2016.05.033, Eqs. (9), (10), Table 5.

    Validity range and fluids tested:
    - 5000 < Re < 32000, supercritical CO2 (CFD, ANSYS Fluent with REFPROP).
    - Cold side of a recuperator, 400 K < T_b < 520 K.
    - Single geometry: channel diameter 2 mm, dh = 1.22 mm, zigzag bend angle 108 deg as
      defined in Meshram Fig. 2b (its relation to labothappy geom['alpha'] is not stated
      in the text: alpha is NOT used), linear pitch 12 mm.
    - Other rows of Table 5 (hot fluid, other temperature ranges, straight channels,
      pitch / angle fits) are not implemented.

    geom : Dh [m], L [m] (developed channel length)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3]),
           T [K] (optional, bulk temperature for validity check)
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Optional bulk temperature
    T_bulk = cond.get('T')

    # Reynolds number on hydraulic diameter (Meshram Eq. 10)
    Re = G * Dh / mu

    # Warn if Re is outside the CFD study validity range
    if Re < 5000 or Re > 32000:
        warnings.warn(
            f"PCHE_meshram_zigzag_SP_DP: Re = {Re:.0f} outside validity range "
            "5000 - 32000 of Meshram et al. (2016): extrapolated.", ValidityRangeWarning
        )

    # Warn if bulk temperature is outside the cold-fluid fit range (400 K to 520 K)
    if T_bulk is not None:
        if T_bulk < 400 or T_bulk > 520:
            warnings.warn(
                f"PCHE_meshram_zigzag_SP_DP: T = {T_bulk:.1f} K outside 400 - 520 K "
                "range of the Meshram et al. (2016) cold-fluid zigzag fit.", ValidityRangeWarning
            )

    # Fanning friction factor, Meshram et al. (2016) Table 5 (cold fluid) and Eq. (9)
    f_fanning = 0.869 * (Re ** -0.512) + 0.041

    # Convert Fanning friction factor to Darcy friction factor
    f_darcy = 4.0 * f_fanning

    # Frictional pressure drop over developed length L
    term_friction = f_darcy * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_katz_zigzag_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop in a zigzag channel of a PCHE [Pa].

    Definitions:
    - Fanning friction factor f is used in DP (Katz et al. Eq. 4:
      DP = 4 * f_fanning * (L / Dh) * G^2 / (2 * rho)).
    - Re = G * Dh / mu.
    - Bulk properties of the cell.
    - Mass flux G of one channel.
    - L = developed channel length (distance between pressure ports corrected for zigzag).
    - One fit per fluid: helium -> Eq. (14), any other fluid -> sCO2 Eq. (13) (warning if not CO2).
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no headers).

    Reference:
    - A. Katz, S.R. Aakre, M.H. Anderson, et al. (2021), Experimental investigation of pressure
      drop and heat transfer in high temperature supercritical CO2 and helium in a printed-circuit
      heat exchanger, International Journal of Heat and Mass Transfer 171, 121089,
      doi:10.1016/j.ijheatmasstransfer.2021.121089, Eqs. (4), (13), (14), Tables 5-6.

    Validity range and fluids tested:
    - Table 5: sCO2 20-200 C, 20 MPa, Re 500-16000 (13455 points); helium 20-550 C, 4 MPa,
      Re 300-3000 (8832 points); fits within +-10 %.
    - ONE geometry tested: 316L PCHE, zigzag pitch angle 37 deg, wavelength 6.63 mm, corner radii
      1.5 mm, semi-elliptical channel A_C 0.95 mm2, D_H = 1 mm, developed length 98 mm.

    AS : CoolProp abstract state (used to identify the fluid).
    geom : Dh [m], L [m], alpha [deg] (zigzag bending angle from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Read fluid name from AS
    fluid = AS.fluid_names()[0].lower()

    # Warn if alpha != 37 deg (geometry tested by Katz et al.)
    if alpha != 37:
        warnings.warn(
            f"PCHE_katz_zigzag_SP_DP: geometry angle {alpha} deg differs from the "
            "37 deg angle of Katz et al. (2021): extrapolated.", ValidityRangeWarning
        )

    # Reynolds number on the hydraulic diameter, Re = G Dh / mu
    Re = G * Dh / mu

    # Select correlation based on fluid name (Katz et al. Table 6 / Model 2)
    if fluid == 'helium':
        # Warn if Re is outside the helium validity range (300 to 3000)
        if Re < 300 or Re > 3000:
            warnings.warn(
                f"PCHE_katz_zigzag_SP_DP: Re = {Re:.0f} outside validity range "
                "300 - 3000 for helium of Katz et al. (2021): extrapolated.", ValidityRangeWarning
            )

        # Katz et al. (2021) Eq. (14), Model 2 for helium (term 1)
        term_one = 0.34 * (Re ** -0.229)

        # Katz et al. (2021) Eq. (14), Model 2 for helium (term 2)
        term_two = 5.532 / Re

        # Total Fanning friction factor for helium
        f_fanning = term_one + term_two

    else:
        # Warn if fluid is not carbondioxide (or recognized sCO2)
        if fluid != 'carbondioxide':
            warnings.warn(
                f"PCHE_katz_zigzag_SP_DP: fluid '{fluid}' not explicitly recognized "
                "as sCO2 or helium; sCO2 correlation applied.", FluidWarning
            )

        # Warn if Re is outside the sCO2 validity range (500 to 16000)
        if Re < 500 or Re > 16000:
            warnings.warn(
                f"PCHE_katz_zigzag_SP_DP: Re = {Re:.0f} outside validity range "
                "500 - 16000 for sCO2 of Katz et al. (2021): extrapolated.", ValidityRangeWarning
            )

        # Katz et al. (2021) Eq. (13), Model 2 for sCO2 (term 1)
        term_one = 0.09 * (Re ** -0.103)

        # Katz et al. (2021) Eq. (13), Model 2 for sCO2 (term 2)
        term_two = 28.593 / Re

        # Total Fanning friction factor for sCO2
        f_fanning = term_one + term_two

    # Frictional pressure drop over developed length L (Katz et al. Eq. 4, Fanning)
    term_friction = 4.0 * f_fanning * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_kim_sg_zigzag_SP_DP(AS, geom, cond):
    """Single-phase sCO2 frictional pressure drop in a zigzag PCHE channel [Pa].

    AS is not used: kept for the common (AS, geom, cond) interface.

    Definitions:
    - FANNING friction factor f used in Eq. (5):
      f = DP / (2 * rho * V^2) * (Dh / L).
    - Darcy friction factor is 4 * f_fanning, used in the pressure drop formula:
      DP = 4 * f_fanning * (L / Dh) * G^2 / (2 * rho).
    - Re = G * Dh / mu, on the hydraulic diameter.
    - Bulk properties of the cell.
    - Mass flux G of one channel.
    - L = developed channel length [m].
    - alpha = geom['alpha'] [deg] = zigzag bending angle from the main flow direction;
      alpha <= 36.25 deg (midpoint) -> 32.5 deg fit, else 40 deg fit.
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no headers).

    Reference:
    - S.G. Kim, Y. Lee, Y. Ahn, J.I. Lee (2016), CFD aided approach to design printed circuit heat
      exchangers for supercritical CO2 Brayton cycle application, Annals of Nuclear Energy 92,
      175-185, doi:10.1016/j.anucene.2016.01.019, Eq. (5), Table 12.

    Validity range and fluids tested:
    - sCO2, 0.7 < Pr < 1.0.
    - R^2 = 0.983 for 32.5 deg, R^2 = 0.908 for 40 deg.
    - CFD (ANSYS CFX) validated on the Ishizuka et al. (2005) PCHE: semicircular channels width
      1.9 / 1.8 mm, depth 0.9 mm, Dh 1.15 mm, fin angles 32.5 and 40 deg, sCO2 recuperator.
    - Fit validity ranges: 2000 < Re < 58000 for 32.5 deg; 2000 < Re < 55000 for 40 deg.

    geom : Dh [m], L [m], alpha [deg] (zigzag bending angle from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Reynolds number on hydraulic diameter
    Re = G * Dh / mu

    # Check which angle is closer (32.5 deg or 40 deg) and apply corresponding fit
    if alpha <= 36.25:
        # Warn if Re is outside the validity range for the 32.5 deg fit (2000 - 58000)
        if Re < 2000 or Re > 58000:
            warnings.warn(
                f"PCHE_kim_sg_zigzag_SP_DP: Re = {Re:.0f} outside validity range "
                "2000 - 58000 for 32.5 deg of Kim S.G. et al. (2016): extrapolated.",
                ValidityRangeWarning
            )

        # Warn if alpha is more than 0.5 deg away from 32.5 deg
        if abs(alpha - 32.5) > 0.5:
            warnings.warn(
                f"PCHE_kim_sg_zigzag_SP_DP: angle {alpha} deg is more than 0.5 deg "
                "away from the fit angle 32.5 deg of Kim S.G. et al. (2016).", ValidityRangeWarning
            )

        # Fanning friction factor for 32.5 deg, Kim S.G. et al. (2016) Table 12
        f_fanning = 0.2515 * (Re ** -0.2031)

    else:
        # Warn if Re is outside the validity range for the 40 deg fit (2000 - 55000)
        if Re < 2000 or Re > 55000:
            warnings.warn(
                f"PCHE_kim_sg_zigzag_SP_DP: Re = {Re:.0f} outside validity range "
                "2000 - 55000 for 40 deg of Kim S.G. et al. (2016): extrapolated.",
                ValidityRangeWarning
            )

        # Warn if alpha is more than 0.5 deg away from 40 deg
        if abs(alpha - 40.0) > 0.5:
            warnings.warn(
                f"PCHE_kim_sg_zigzag_SP_DP: angle {alpha} deg is more than 0.5 deg "
                "away from the fit angle 40.0 deg of Kim S.G. et al. (2016).", ValidityRangeWarning
            )

        # Fanning friction factor for 40 deg, Kim S.G. et al. (2016) Table 12
        f_fanning = 0.2881 * (Re ** -0.1322)

    # Convert Fanning friction factor to Darcy friction factor
    f_darcy = 4.0 * f_fanning

    # Frictional pressure drop over developed length L (Fanning definition of Eq. 5)
    term_friction = f_darcy * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_baik_zigzag_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop in a zigzag PCHE channel [Pa],
    Baik et al. (2017), Table 3: one fit for water (laminar), one for CO2 (turbulent).

    Definitions:
    - Fanning friction factor f is used in the pressure drop equation:
      DP = 4 * f_fanning * (L / D_mod) * G^2 / (2 * rho). (Bae et al. 2019 quote the Darcy
      form of the CO2 fit, 0.2992 Re^-0.19 = 4 x 0.0748.)
    - Modified hydraulic diameter D_mod = sqrt(4 * A / pi) (circle of the same
      flow area, Baik Fig. 8). If 'D_ch' in geom, D_mod = geom['D_ch'] / sqrt(2.0) (ideal
      semicircle); otherwise D_mod = geom['Dh'] with a warning.
    - Re = G * D_mod / mu, based on the modified diameter.
    - Bulk properties of the cell.
    - Mass flux G of one channel.
    - Developed channel length L [m].
    - Bending angle alpha [deg] (zigzag bending angle from the main flow direction).
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no headers).

    Reference:
    - S. Baik, S.G. Kim, J. Lee, J.I. Lee (2017), Study on CO2-water printed
      circuit heat exchanger performance operating under various CO2 phases for
      S-CO2 power cycle application, Applied Thermal Engineering 113, 1536-1546,
      doi:10.1016/j.applthermaleng.2016.11.132, Eqs. (4), (5), Table 3, Fig. 8.

    Validity range and fluids tested:
    - Water: Re between 50 and 200 (Baik Table 3).
    - Carbon dioxide (supercritical / CO2): Re between 15000 and 85000 (Baik Table 3).
    - Pitch-averaged f (zigzag bend losses included, headers excluded: Baik Table 4 gives
      header K 0.72 + 0.42 inlet, 1.00 + 0.05 outlet).
    - ONE geometry tested: KAIST S-CO2PE precooler (23 kW), zigzag 32.5 deg,
      wave length 9 mm, channel diameter 1.8 mm, corners rounded R 1 mm, SS316L.
      Correlations from CFD (ANSYS CFX) validated on tests.

    AS : CoolProp abstract state (used to identify fluid name via AS.fluid_names()[0].lower())
    geom : D_ch [m] (optional, semicircular channel diameter), Dh [m] (fallback hydraulic diameter),
           L [m] (developed channel length), alpha [deg] (zigzag bending angle)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read geometry inputs
    L = geom['L']
    alpha = geom['alpha']

    # Read condition inputs
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Determine modified hydraulic diameter D_mod (Baik Fig. 8)
    if 'D_ch' in geom:
        D_mod = geom['D_ch'] / np.sqrt(2.0)
    else:
        D_mod = geom['Dh']
        warnings.warn(
            "PCHE_baik_zigzag_SP_DP: 'D_ch' not found in geom; using geom['Dh'] "
            "as D_mod, but Baik (2017) uses a modified hydraulic diameter.", AssumedParameterWarning
        )

    # Check bending angle against the single geometry tested by Baik (32.5 deg)
    if alpha != 32.5:
        warnings.warn(
            f"PCHE_baik_zigzag_SP_DP: bending angle {alpha} deg differs from "
            "the 32.5 deg geometry of Baik et al. (2017).", ValidityRangeWarning
        )

    # Reynolds number based on the modified hydraulic diameter
    Re = G * D_mod / mu

    # Identify fluid name from AS
    fluid = AS.fluid_names()[0].lower()

    # Select correlation based on fluid name (Baik Table 3)
    if fluid == 'water':
        # Warn if Re is outside the water validity range (50 - 200)
        if Re < 50 or Re > 200:
            warnings.warn(
                f"PCHE_baik_zigzag_SP_DP: Re = {Re:.1f} outside the water validity range "
                "50 - 200 of Baik et al. (2017): extrapolated.", ValidityRangeWarning
            )

        # Fanning friction factor for water, Baik Table 3
        f_fanning = 6.9982 * (Re ** -0.766)

    else:
        # Check if fluid is carbondioxide, warn if not
        if fluid != 'carbondioxide':
            warnings.warn(
                f"PCHE_baik_zigzag_SP_DP: fluid '{fluid}' is not 'carbondioxide'; "
                "using the carbon dioxide correlation of Baik et al. (2017).", FluidWarning
            )

        # Warn if Re is outside the carbon dioxide validity range (15000 - 85000)
        if Re < 15000 or Re > 85000:
            warnings.warn(
                f"PCHE_baik_zigzag_SP_DP: Re = {Re:.1f} outside the carbon dioxide validity "
                "range 15000 - 85000 of Baik et al. (2017): extrapolated.", ValidityRangeWarning
            )

        # Fanning friction factor for carbon dioxide, Baik Table 3
        f_fanning = 0.0748 * (Re ** -0.19)

    # Frictional pressure drop over developed length L (Baik Eqs. 4-5, Fanning)
    term_friction = 4.0 * f_fanning * (L / D_mod)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_saeed_zigzag_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop in a zigzag channel of a PCHE [Pa],
    Saeed et al. (2020), Table 6: one fit for sCO2, one for water.

    Definitions:
    - f is a DARCY friction factor: Eq. (11) f = (dp/dz) 2 Dh / (rho V^2) (checked: their
      Fig. 16 plots Baik's Fanning fit multiplied by 4). DP = f_darcy (L / Dh) G^2 / (2 rho).
    - Re = G * Dh / mu, on the hydraulic diameter.
    - Bulk properties of the cell.
    - Mass flux G of one channel.
    - L = developed channel length [m].
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no headers).

    Reference:
    - M. Saeed, A.S. Berrouk, M.S. Siddiqui, A.A. Awais (2020), Numerical investigation of thermal
      and hydraulic characteristics of sCO2-water printed circuit heat exchangers with zigzag
      channels, Energy Conversion and Management 224, 113375, doi:10.1016/j.enconman.2020.113375,
      Eqs. (11)-(12), Table 6.

    Validity range and fluids tested:
    - sCO2: Re 3000 - 60000, Pr 2 - 13 (Table 6; the text of Section 3 says 5000-70000).
    - Water: Re 50 - 15000.
    - About 90% of the CFD points within +-10%.
    - CFD (72 periodic units, segment-averaged data), one geometry: fin angle 40 deg
      (= labothappy alpha), pitch along length 9 mm, transverse pitch 2.35 mm, Dh 1.106 mm,
      SS316L, sCO2 precooler cooled by water (Table 2).

    AS : CoolProp abstract state (used to determine fluid name).
    geom : Dh [m], L [m], alpha [deg] (zigzag bending angle from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read geometry parameters
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']

    # Read condition parameters
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Determine fluid name from CoolProp abstract state
    fluid = AS.fluid_names()[0].lower()

    # Warn if zigzag bending angle is not 40 deg (Saeed et al. geometry)
    if alpha != 40:
        warnings.warn(
            f"PCHE_saeed_zigzag_SP_DP: bending angle {alpha} deg differs from the "
            "40 deg geometry of Saeed et al. (2020).", ValidityRangeWarning
        )

    # Reynolds number on the hydraulic diameter, Re = G Dh / mu
    Re = G * Dh / mu

    # Select correlation based on fluid name
    if fluid == 'water':
        # Warn if Reynolds number is outside the water validity range (50 - 15000)
        if Re < 50 or Re > 15000:
            warnings.warn(
                f"PCHE_saeed_zigzag_SP_DP: Re = {Re:.0f} outside water validity range "
                "50 - 15000 of Saeed et al. (2020): extrapolated.", ValidityRangeWarning
            )

        # Darcy friction factor for water, Saeed et al. (2020) Table 6
        f_darcy = 3.04 * (Re ** -0.37)

    else:
        # Warn if fluid is not carbondioxide
        if fluid != 'carbondioxide':
            warnings.warn(
                f"PCHE_saeed_zigzag_SP_DP: fluid '{fluid}' not explicitly tested in "
                "Saeed et al. (2020) (tested: water, sCO2).", FluidWarning
            )

        # Warn if Reynolds number is outside the sCO2 validity range (3000 - 60000)
        if Re < 3000 or Re > 60000:
            warnings.warn(
                f"PCHE_saeed_zigzag_SP_DP: Re = {Re:.0f} outside sCO2 validity range "
                "3000 - 60000 of Saeed et al. (2020): extrapolated.", ValidityRangeWarning
            )

        # Darcy friction factor for sCO2, Saeed et al. (2020) Table 6
        f_darcy = 0.13 * (Re ** -0.044)

    # Frictional pressure drop over developed length L (Darcy, Saeed et al. 2020 Eq. 11)
    term_friction = f_darcy * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP


def PCHE_chen_zigzag_SP_DP(AS, geom, cond):
    """Single-phase frictional pressure drop of helium in a zigzag PCHE channel [Pa],
    Chen et al. (2016), Eq. (11).

    AS is used to check the fluid name ('helium').

    Definitions:
    - FANNING friction factor f used in DP (DP = 4 * f_fanning * (L / Dh) * G^2 / (2 * rho)).
    - Re = G * Dh / mu, with Dh the hydraulic diameter of the semicircular channel [m].
    - Bulk properties of the cell.
    - Mass flux G of ONE channel [kg/m2/s].
    - L = developed channel length [m].
    - geom['alpha'] [deg] = zigzag bending angle from the main flow direction.
    - Returns the frictional pressure drop [Pa] over L (no gravity, no acceleration, no headers).

    Reference:
    - M. Chen, X. Sun, R.N. Christensen, I. Skavdahl, V. Utgikar, P. Sabharwall (2016),
      Pressure drop and heat transfer characteristics of a high-temperature printed circuit
      heat exchanger, Applied Thermal Engineering 108, 1409-1417,
      doi:10.1016/j.applthermaleng.2016.07.149, Eqs. (8), (11).

    Validity range and fluids tested:
    - 164 data points (72 laminar, 92 transition), 98 % within +-10 %.
    - ONE geometry (Alloy 617, 8 hot + 8 cold plates, 11 channels each): semicircular channels
      nominal diameter about 2 mm, span-wise pitch 2.5 mm, zigzag angle 15 deg from the flow
      direction (alpha = 15 deg), bend radius of curvature 4 mm.
    - Helium on both sides, up to 464 / 802 C and 2.7 MPa, flow up to 39 kg/h.
    - Re validity range: 1400 <= Re <= 3558 (laminar for Re <= 2200,
      transition for 2200 < Re <= 3558).

    geom : Dh [m], L [m], alpha [deg] (zigzag bending angle from the main flow direction)
    cond : G [kg/m2/s] (mass flux in ONE channel), mu [Pa s], rho [kg/m3] (or D [kg/m3])
    """
    # Read inputs from geometry dictionary
    Dh = geom['Dh']
    L = geom['L']
    alpha = geom['alpha']

    # Read inputs from condition dictionary
    G = cond['G']
    mu = cond['mu']

    # Read density, fallback to 'D' if 'rho' is not in cond
    if 'rho' in cond:
        rho = cond['rho']
    else:
        rho = cond['D']

    # Get fluid name from AS and check if it is helium
    fluid = AS.fluid_names()[0].lower()
    if fluid != "helium":
        warnings.warn(
            f"PCHE_chen_zigzag_SP_DP: fluid is '{fluid}', but Chen et al. (2016) "
            "correlation was developed and tested for helium.", FluidWarning
        )

    # Check zigzag bending angle
    if alpha != 15:
        warnings.warn(
            f"PCHE_chen_zigzag_SP_DP: alpha = {alpha} deg, but Chen et al. (2016) "
            "tested strictly alpha = 15 deg.", ValidityRangeWarning
        )

    # Reynolds number on the hydraulic diameter, Re = G Dh / mu
    Re = G * Dh / mu

    # Check Reynolds number validity range (1400 <= Re <= 3558)
    if Re < 1400 or Re > 3558:
        warnings.warn(
            f"PCHE_chen_zigzag_SP_DP: Re = {Re:.1f} outside the validity range "
            "1400 - 3558 of Chen et al. (2016): extrapolated.", ValidityRangeWarning
        )

    # Fanning friction factor, Chen et al. (2016) Eq. (11)
    if Re <= 2200:
        # Laminar regime (1400 <= Re <= 2200)
        f_fanning = 17.639 * (Re ** -0.8861)
    else:
        # Transition regime (2200 < Re <= 3558): constant, mean of 92 points
        f_fanning = 0.019044

    # Frictional pressure drop over developed length L (Fanning definition, Eq. 8)
    term_friction = 4.0 * f_fanning * (L / Dh)
    term_dynamic_head = (G ** 2) / (2.0 * rho)
    DP = term_friction * term_dynamic_head

    return DP

#%% TWO PHASE

def PCHE_bae_two_phase_DP(AS, geom, cond):
    """Two-phase frictional pressure drop of CO2 near the critical point in a zigzag PCHE
    channel [Pa], Bae et al. (2019): Baik's gas-only Darcy friction factor times homogeneous
    (HEM) two-phase multiplier with the Cicchitti viscosity and form-loss correction factor.
    Used for condensation and evaporation cells.

    Definitions:
    - Darcy friction factor f used in DP (explicitly stated: f_go = 0.2992 Re^-0.19,
      which is 4 x Baik's Fanning 0.0748).
    - Modified hydraulic diameter D_mod = sqrt(4 * A / pi) (circle of the same flow area).
      If 'D_ch' in geom: D_mod = geom['D_ch'] / sqrt(2.0) (ideal semicircle of diameter D_c:
      A = pi * D_c^2 / 8); otherwise use geom['Dh'] and warn that the modified diameter
      could not be computed. G is unchanged (same area). Bae does not restate the diameter:
      Baik's is assumed (same KAIST PCHE).
    - Re_go = G * D_mod / mu_v, on the vapour viscosity and total mass flux G.
    - Saturated properties from _sat_props(AS, p, cond).
    - Homogeneous density rho_h = 1 / (x / rho_v + (1 - x) / rho_l).
    - Cicchitti viscosity mu_h = x * mu_v + (1 - x) * mu_l (Bae Eq. 9).
    - HEM multiplier phi2_go = (rho_v / rho_h) * (mu_h / mu_v) ** 0.19 (Bae Eq. 12).
    - Form-loss correction = max(1.0, -20.3 * p / p_crit + 19.9) (Bae Eq. 16).
    - Returns the FRICTIONAL pressure drop [Pa] over L (frictional + bend form losses; the HEM
      acceleration term, Bae Eq. 14, is NOT included).

    References:
    - S.J. Bae, J. Kwon, S.G. Kim, I. Son, J.I. Lee (2019), Condensation heat transfer and
      multi-phase pressure drop of CO2 near the critical point in a printed circuit heat
      exchanger, International Journal of Heat and Mass Transfer 129, 1206-1221,
      doi:10.1016/j.ijheatmasstransfer.2018.10.055, Table 4, Eqs. (9), (12), (16).
    - Friction factor of S. Baik, S.G. Kim, J. Lee, J.I. Lee (2017), Study on CO2-water
      printed circuit heat exchanger performance operating under various CO2 phases for S-CO2
      power cycle application, Applied Thermal Engineering 113, 1536-1546,
      doi:10.1016/j.applthermaleng.2016.11.132, Table 3.

    Validity range and fluids tested:
    - CO2 26.2-29.6 C, 6.62-6.98 MPa, x 0.77-1, about 1 kg/s in the Baik PCHE (ONE geometry
      (KAIST S-CO2PE precooler, 23 kW): zigzag 32.5 deg, wave length 9 mm, channel d 1.8 mm,
      corners rounded R 1 mm, SS316L; correlations from CFD (ANSYS CFX) validated on the
      tests, pitch-averaged, bend form losses included; headers excluded).
    - Within about 5 % of the measured core DP with the correction (Fig. 21), 10-35 % without.
      The correction is 1 above p / p_crit = 0.931.

    AS : CoolProp abstract state of the fluid (moved to saturation by _sat_props).
    geom : L [m] (developed channel length), D_ch [m] (optional, channel diameter),
           Dh [m] (fallback hydraulic diameter)
    cond : G [kg/m2/s] (mass flux in ONE channel), p [Pa] (pressure),
           x [-] (vapour quality)
    """
    # Check required condition keys
    required_cond_keys = ['G', 'p', 'x']
    for key in required_cond_keys:
        if key not in cond:
            raise KeyError(f"PCHE_bae_two_phase_DP: missing required condition key '{key}'.")

    # Read geometry and conditions
    L = geom['L']
    G = cond['G']
    p = cond['p']
    x = cond['x']

    # Read fluid name and warn if not carbondioxide
    fluid = AS.fluid_names()[0].lower()
    if fluid != 'carbondioxide':
        warnings.warn(
            f"PCHE_bae_two_phase_DP: fluid is '{fluid}', but correlation was developed "
            "for carbon dioxide.", FluidWarning
        )

    # Keep quality x in [0, 1]
    if x < 0.0:
        x = 0.0
    if x > 1.0:
        x = 1.0

    # Warn if quality is outside the validity range x 0.77 - 1
    if x < 0.77:
        warnings.warn(
            f"PCHE_bae_two_phase_DP: quality x = {x:.3f} outside validity range "
            "0.77 - 1.0 of Bae et al. (2019): extrapolated.", ValidityRangeWarning
        )

    # Read critical pressure first (before _sat_props moves AS)
    p_crit = AS.p_critical()

    # Warn if pressure ratio p / p_crit is outside 0.89 - 0.95 (Bae data 6.62 - 6.98 MPa)
    pressure_ratio = p / p_crit
    if pressure_ratio < 0.89 or pressure_ratio > 0.95:
        warnings.warn(
            f"PCHE_bae_two_phase_DP: pressure ratio p/p_crit = {pressure_ratio:.3f} "
            "outside validity range 0.89 - 0.95 of Bae et al. (2019): extrapolated.",
            ValidityRangeWarning
        )

    # Compute modified hydraulic diameter (Baik's modified diameter D = sqrt(4 * A / pi))
    if 'D_ch' in geom:
        D_c = geom['D_ch']
        D_mod = D_c / np.sqrt(2.0)
    else:
        D_mod = geom['Dh']
        warnings.warn(
            "PCHE_bae_two_phase_DP: 'D_ch' not found in geom. Used geom['Dh'] as "
            "D_mod, but Baik's modified diameter could not be precisely computed.",
            AssumedParameterWarning
        )

    # Read saturated properties using library helper
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    mu_v = s['mu_v']

    # 1. Gas-only Reynolds number
    Re_go = G * D_mod / mu_v

    # Warn if Re_go is outside the range 15000 - 85000 of Baik's CO2 fit
    if Re_go < 15000 or Re_go > 85000:
        warnings.warn(
            f"PCHE_bae_two_phase_DP: Re_go = {Re_go:.0f} outside validity range "
            "15000 - 85000 of Baik's CO2 fit: extrapolated.", ValidityRangeWarning
        )

    # 2. Baik gas-only Darcy friction factor as quoted by Bae (Table 4)
    f_go = 0.2992 * (Re_go ** -0.19)

    # 3. Gas-only frictional pressure drop
    DP_go = f_go * (L / D_mod) * (G ** 2) / (2.0 * rho_v)

    # 4. Homogeneous density and Cicchitti viscosity (Bae Eq. 9)
    rho_h = 1.0 / (x / rho_v + (1.0 - x) / rho_l)
    mu_h = x * mu_v + (1.0 - x) * mu_l

    # 5. HEM multiplier (Bae Eq. 12)
    phi2_go = (rho_v / rho_h) * ((mu_h / mu_v) ** 0.19)

    # 6. Form-loss correction factor (Bae Eq. 16)
    correction = -20.3 * (p / p_crit) + 19.9
    if correction < 1.0:
        correction = 1.0

    # 7. Total two-phase frictional pressure drop (excluding HEM acceleration term)
    DP = DP_go * phi2_go * correction

    return DP
