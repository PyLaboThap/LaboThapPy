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
                "Lee et al. (2024): extrapolated."
            )

        # Warn if alpha outside 10-40 deg (range of Lee's CFD)
        if not (10 <= alpha <= 40):
            warnings.warn(
                f"PCHE_lee_jiang_zigzag_SP_DP: bending angle {alpha} deg outside the range "
                "10-40 deg of Lee et al. (2024): extrapolated."
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
                "Jiang et al. (2018): extrapolated."
            )

        # Warn if alpha is not 30 or 40 deg
        if alpha != 30 and alpha != 40:
            warnings.warn(
                f"PCHE_lee_jiang_zigzag_SP_DP: Jiang et al. (2018) gives fits for 30 and 40 deg "
                f"channels only, used here with alpha = {alpha} deg."
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
