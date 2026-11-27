# -*- coding: utf-8 -*-
"""
Finned-side heat transfer correlations of tube & fins heat exchangers (TFHE), (AS, geom, cond) formalism.

Source : Handbook for transversely finned tube heat exchanger design,
         E. Pis'mennyi, G. Polupan, I. Carvajal-Mariscal, F. Sanchez-Silva, I. Pioro.
         Staggered bank of transversely finned tubes.

Both functions return the REDUCED heat transfer coefficient of the finned side
    h_rdc = (A_f/A * eta_f * psi_eta + A_t/A) * h_c       [W/m2/K]
i.e. referred to the TOTAL outer area geom['A'] (tube + fins): the fin efficiency is already inside,
so UA = 1 / (1/(h_rdc A) + ...) with A = total outer area.

geom : D_t_o, t_t, L_t, D_f_o, t_f, f_density, k_f, n_t, n_rows, n_passes, w (bundle width),
       PR (annular) or P_V, P_H (square), n_parallel (optional)
cond : m_dot (WHOLE finned side), rho, mu, k, Pr   (bulk properties of the gas)
"""
import warnings
import numpy as np
from scipy.special import i0, i1, k0, k1      # Bessel functions (Kern & Kraus fin efficiency)
from labothappy.correlations.convection.correlation_warnings import ValidityRangeWarning

_NEEDED = ('D_t_o', 't_t', 'L_t', 'D_f_o', 't_f', 'f_density', 'k_f', 'n_t', 'n_rows',
           'n_passes', 'w')


def _need(geom, cond, extra=()):
    miss = [k for k in _NEEDED + tuple(extra) if k not in geom]
    if miss:
        raise KeyError(f"TFHE HTC: geom is missing {miss}.")
    miss = [k for k in ('m_dot', 'rho', 'mu', 'k', 'Pr') if k not in cond]
    if miss:
        raise KeyError(f"TFHE HTC: cond is missing {miss} (single-phase bulk properties of the gas).")


def _fin_efficiency_straight_rect(t_fin, L_fin, h_c, k_f, D_fin, D_tube):
    """Efficiency eta_f of an annular fin treated as a straight rectangular fin of corrected length
    (Pis'mennyi), and the correction factor psi_eta. mL = sqrt(2 h/(k t)) L'."""
    if h_c < 0:
        raise ValueError(f"TFHE: negative heat transfer coefficient {h_c}.")
    m = np.sqrt(2 * h_c / (k_f * t_fin))
    L_p = L_fin * (1 + (0.191 + 0.054 * (D_fin / D_tube)) * np.log(D_fin / D_tube))
    mL = m * L_p
    eta_f = 1.0 if mL == 0 else np.tanh(mL) / mL
    psi_eta = 1 - 0.016 * (D_fin / D_tube - 1) * (1 + np.tanh(2 * mL - 1))
    return eta_f, psi_eta


def _bank_h(geom, cond, area_fractions, Psi_f, sigma_1, sigma_2, pitch_for_diag):
    """Common part: gas velocity in the narrowest section, h_c, fin efficiency, reduced coefficient."""
    D, t = geom['D_t_o'], geom['t_f']
    L, N_rows = geom['L_t'], geom['n_rows'] * geom['n_passes']
    n_f = L * geom['f_density'] - 1
    s = (L - n_f * t) / (n_f + 1)                                   # fin spacing
    Fin_L = (geom['D_f_o'] - D) / 2
    n_tpr = geom['n_t'] / N_rows                                      # tubes per row

    D_fin_c = D + 2 * Fin_L * t / s                                       # conventional fin diameter
    diag = np.sqrt(2) * pitch_for_diag                                    # staggered bank
    psi_c = (pitch_for_diag - D_fin_c) / (diag - D_fin_c)
    S_free = geom['L_t'] * geom['w'] - n_tpr * D_fin_c * L
    S_flow = S_free * (2 / psi_c) if psi_c > 2 else S_free
    u = cond['m_dot'] / geom.get('n_parallel', 1) / cond['rho'] / S_flow
    nu = cond['mu'] / cond['rho']                                         # KINEMATIC viscosity

    X = sigma_1 / sigma_2 - 1.26 / Psi_f - 2
    n = 0.7 + 0.08 * np.tanh(X) + 0.005 * Psi_f
    C_q = (1.36 - np.tanh(X)) * (1.1 / (Psi_f + 8) - 0.014)
    if N_rows > 8:
        C_z = 1.0
    elif sigma_1 / sigma_2 > 2:
        C_z = 3.5 * N_rows ** 0.03 - 2.72
    else:
        C_z = 3.15 * N_rows ** 0.05 - 2.5
    h_c = 1.13 * C_z * C_q * (cond['k'] / D) * (u * D / nu) ** n * cond['Pr'] ** 0.33

    eta_f, psi_eta = _fin_efficiency_straight_rect(t, Fin_L, h_c, geom['k_f'], geom['D_f_o'], D)
    Af_A, At_A = area_fractions
    return (Af_A * eta_f * psi_eta + At_A) * h_c


def TFHE_pismennyi_annular_HTC(AS, geom, cond):
    """Annular (round) fins, staggered bank. sigma_1 = sigma_2 = geom['PR'], pitch = PR * D_t_o."""
    _need(geom, cond, ('PR',))
    D, t, L = geom['D_t_o'], geom['t_f'], geom['L_t']
    n_f = L * geom['f_density'] - 1
    s = (L - n_f * t) / (n_f + 1)
    # areas of ONE tube (the fractions do not depend on the number of tubes)
    A_tip = 2 * np.pi * (geom['D_f_o'] / 2) * t * n_f
    A_faces = 2 * np.pi * ((geom['D_f_o'] / 2) ** 2 - (D / 2) ** 2) * n_f
    A_tube = 2 * np.pi * (D / 2) * (L - n_f * t)
    A_tot = A_tip + A_faces + A_tube
    Psi_f = (geom['D_f_o'] ** 2 - D ** 2 + 2 * geom['D_f_o'] * t) / (2 * D * s) + 1 - t / s
    sigma = geom['PR']
    return _bank_h(geom, cond, ((A_tip + A_faces) / A_tot, A_tube / A_tot), Psi_f, sigma, sigma,
                   geom.get('P', sigma * D))


def TFHE_pismennyi_square_HTC(AS, geom, cond):
    """Square (plate) fins, staggered bank. sigma_1 = P_V / D_t_o, sigma_2 = P_H / D_t_o."""
    _need(geom, cond, ('P_V', 'P_H'))
    D, t, L = geom['D_t_o'], geom['t_f'], geom['L_t']
    n_f = L * geom['f_density'] - 1
    s = (L - n_f * t) / (n_f + 1)
    A_r = 2 * (geom['D_f_o'] ** 2 - 0.785 * D ** 2 + 2 * geom['D_f_o'] * t) * (L / s)   # fins of ONE tube
    A_t = np.pi * D * (L - n_f * t)                                                    # bare tube of ONE tube
    A_tot = A_r + A_t
    Psi_f = 2 * (geom['D_f_o'] ** 2 - 0.785 * D ** 2 + 2 * geom['D_f_o'] * t) / (np.pi * D * s) + 1 - t / s
    return _bank_h(geom, cond, (A_r / A_tot, A_t / A_tot), Psi_f, geom['P_V'] / D, geom['P_H'] / D,
                   np.sqrt(geom['P_V'] * geom['P_H']))


def TFHE_pismennyi_HTC(AS, geom, cond):
    """Pis'mennyi correlation routed on geom['fin_type'] : 'Annular' (round fins) or 'Square' (plate fins),
    case insensitive."""
    ft = str(geom.get('fin_type', '')).strip().lower()
    if ft == 'annular':
        return TFHE_pismennyi_annular_HTC(AS, geom, cond)
    if ft == 'square':
        return TFHE_pismennyi_square_HTC(AS, geom, cond)
    raise ValueError(f"TFHE HTC: geom['fin_type'] = {geom.get('fin_type')!r}: expected 'Annular' or 'Square'.")

#%% CIRCULAR (ANNULAR) FINS: Briggs-Young, ESDU high fin, Ganguli-VDI (Kern & Kraus fin efficiency)

def TFHE_briggs_young_HTC(AS, geom, cond):
    """
    Gas-side heat transfer coefficient of a staggered (triangular pitch) bank of tubes with
    circular (annular) fins, Briggs & Young (1963). Air coolers.

    Definitions: Re = rho V_max D_t_o / mu, V_max = gas velocity in the minimum flow area A_min
    (min of the transverse and diagonal gaps, fins included, as ht / fluids); bulk gas
    properties; Nu on D_t_o. Fin efficiency of Kern & Kraus (1972) (exact annular fin, Bessel
    functions). Returns the REDUCED h referred to the TOTAL outer area (tube + fins), fin
    efficiency included (labothappy convention, as the Pis'mennyi functions)
    = ht h_bare_tube_basis / A_increase (ht returns h on the bare-tube area).

    Implementation follows ht.air_cooler.h_Briggs_Young (C. Bell), which cites the reference
    below; areas and minimum flow area as fluids.geometry.AirCooledExchanger; fin efficiency
    of Kern & Kraus (1972) as ht.core.fin_efficiency_Kern_Kraus.

    Validity (ht notes; warning outside): 1000 < Re < 8000, 11.13 mm < D_t_o < 40.89 mm,
    1.42 mm < fin height < 16.57 mm, 0.33 mm < fin thickness < 2.02 mm,
    1.30 mm < fin pitch < 4.06 mm, 24.49 mm < transverse pitch < 111 mm; triangular
    (staggered) banks; air.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        Dictionary of geometry parameters with SI units:
        - 'fin_type' : str, fin type ('Annular')
        - 'D_t_o' : float, tube outer (root) diameter [m]
        - 'D_f_o' : float, fin outer diameter [m]
        - 't_f' : float, fin thickness [m]
        - 'f_density' : float, number of fins per metre of tube [1/m]
        - 'L_t' : float, tube length [m]
        - 'k_f' : float, fin thermal conductivity [W/m/K]
        - 'n_rows' : float or int, number of tube rows [-]
        - 'n_passes' : float or int, number of tube passes [-]
        - 'n_t' : float or int, total number of tubes [-]
        - 'P_V' : float, transverse pitch (normal to gas flow) [m]
        - 'P_H' : float, longitudinal pitch (parallel to gas flow) [m]
        - 'n_parallel' : float or int, optional, number of parallel paths [-]
        - 'theta' : float, optional, layout angle [deg]: 0/90 inline, 30/45/60 staggered (default 45)
    cond : dict
        Dictionary of operating conditions with SI units:
        - 'm_dot' : float, total gas mass flow rate through the finned side [kg/s]
        - 'rho' : float, bulk gas density [kg/m3]
        - 'mu' : float, bulk gas dynamic viscosity [Pa s]
        - 'k' : float, bulk gas thermal conductivity [W/m/K]
        - 'Pr' : float, bulk gas Prandtl number [-]

    Returns
    -------
    h_reduced : float
        Reduced heat transfer coefficient referred to the total outer area [W/m2/K].

    Reference
    ---------
    Briggs, D.E., Young, E.H. (1963) Convection heat transfer and pressure drop of air flowing
    across triangular banks of finned tubes. Chemical Engineering Progress Symposium Series
    59(41) (Heat Transfer - Houston).
    """
    # Step 1. Check the fin shape
    fin_type = str(geom.get('fin_type', 'Annular')).strip().lower()
    if fin_type != 'annular':
        raise ValueError(
            f"TFHE Briggs-Young: correlation is for circular (annular) fins only. "
            f"Received geom['fin_type'] = {geom.get('fin_type')!r}."
        )

    # Tube arrangement check
    # Tube arrangement from the layout angle theta [deg] (default 45, i.e. staggered):
    # 0 or 90 = inline, 30, 45 or 60 = staggered
    theta = geom.get('theta', 45)
    if theta == 0 or theta == 90:
        tube_arrangement = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        tube_arrangement = 'staggered'
    else:
        raise ValueError(f"TFHE: theta = {theta!r} deg: expected 0 or 90 (inline), "
                         "30, 45 or 60 (staggered).")
    if tube_arrangement == 'inline':
        warnings.warn(
            "Briggs-Young: correlation for staggered (triangular pitch) banks, used here for "
            "an inline bank.",
            ValidityRangeWarning
        )

    # Step 2. Inputs (SI units)
    D = geom['D_t_o']
    D_fin = geom['D_f_o']
    t_fin = geom['t_f']
    fin_density = geom['f_density']
    L_tube = geom['L_t']
    k_f = geom['k_f']
    tube_rows = geom['n_rows'] * geom['n_passes']
    tubes_per_row = geom['n_t'] / tube_rows
    pitch_normal = geom['P_V']
    pitch_parallel = geom['P_H']
    n_parallel = geom.get('n_parallel', 1)
    m_dot = cond['m_dot'] / n_parallel
    rho = cond['rho']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Step 3. Fin dimensions and areas of ONE tube
    fin_height = (D_fin - D) / 2
    fin_interval = 1.0 / fin_density
    bare_length = fin_interval - t_fin
    A_tube_showing = np.pi * D * L_tube * (1.0 - t_fin / fin_interval)
    A_per_fin = 2.0 * np.pi / 4.0 * (D_fin**2 - D**2) + np.pi * D_fin * t_fin
    A_fin = fin_density * L_tube * A_per_fin
    A_total = A_tube_showing + A_fin

    # Step 4. Minimum flow area and maximum velocity
    pitch_diagonal = np.sqrt(pitch_parallel**2 + (pitch_normal / 2.0)**2)
    blockage = D + 2.0 * fin_density * fin_height * t_fin
    A_diagonal = 2.0 * tubes_per_row * L_tube * (pitch_diagonal - blockage)
    A_normal = tubes_per_row * L_tube * (pitch_normal - blockage)
    A_min = min(A_diagonal, A_normal)
    V_max = m_dot / (A_min * rho)
    Re = rho * V_max * D / mu

    # Validity checks (ranges of the ht notes): warning, then computed anyway
    if not (1000.0 < Re < 8000.0):
        warnings.warn(f"Briggs-Young: Re = {Re:.0f} outside the validity range (1000, 8000).",
                      ValidityRangeWarning)
    if not (0.01113 < D < 0.04089):
        warnings.warn(f"Briggs-Young: D_t_o = {D * 1e3:.2f} mm outside the validity range "
                      "(11.13 mm, 40.89 mm).", ValidityRangeWarning)
    if not (0.00142 < fin_height < 0.01657):
        warnings.warn(f"Briggs-Young: fin height = {fin_height * 1e3:.2f} mm outside the "
                      "validity range (1.42 mm, 16.57 mm).", ValidityRangeWarning)
    if not (0.00033 < t_fin < 0.00202):
        warnings.warn(f"Briggs-Young: fin thickness = {t_fin * 1e3:.3f} mm outside the "
                      "validity range (0.33 mm, 2.02 mm).", ValidityRangeWarning)
    if not (0.00130 < fin_interval < 0.00406):
        warnings.warn(f"Briggs-Young: fin pitch = {fin_interval * 1e3:.2f} mm outside the "
                      "validity range (1.30 mm, 4.06 mm).", ValidityRangeWarning)
    if not (0.02449 < pitch_normal < 0.111):
        warnings.warn(f"Briggs-Young: transverse pitch = {pitch_normal * 1e3:.2f} mm outside "
                      "the validity range (24.49 mm, 111 mm).", ValidityRangeWarning)

    # Nusselt number, Briggs & Young (1963) (as ht h_Briggs_Young; Pr exponent 1/3)
    term_re = Re**0.681
    term_pr = Pr**(1.0 / 3.0)
    term_geom1 = (bare_length / fin_height)**0.2
    term_geom2 = (bare_length / t_fin)**0.1134
    Nu = 0.134 * term_re * term_pr * term_geom1 * term_geom2

    # Heat transfer coefficient of the finned surface, before the fin efficiency
    h = Nu * k / D

    # Step 5. Fin efficiency (Kern & Kraus, 1972) and reduced HTC
    m_fin = np.sqrt(2.0 * h / (k_f * t_fin))
    r_e = D_fin / 2.0
    r_o = D / 2.0
    numerator = i1(m_fin * r_e) * k1(m_fin * r_o) - k1(m_fin * r_e) * i1(m_fin * r_o)
    denominator = i0(m_fin * r_o) * k1(m_fin * r_e) + i1(m_fin * r_e) * k0(m_fin * r_o)
    eta_fin = 2.0 * r_o / (m_fin * (r_e**2 - r_o**2)) * numerator / denominator
    h_reduced = (eta_fin * A_fin + A_tube_showing) / A_total * h

    return h_reduced


def TFHE_esdu_high_fin_HTC(AS, geom, cond):
    """
    Gas-side heat transfer coefficient of a staggered bank of high-fin tubes
    (circular fins) using the ESDU method.

    Nu = 0.242 Re^0.658 (bare length / fin height)^0.297 (P_normal / P_parallel)^-0.091
         Pr^(1/3) F1 F2
    F2 = row correction: 0.76, 0.84, 0.92 for 1, 2, 3 rows, 1 for 4 rows or more.
    F1 = (Pr / Pr_wall)^0.26 (wall property correction) is NOT applied (F1 = 1): the wall
    Prandtl number is not available in cond.

    Definitions: Re = rho V_max D_t_o / mu, V_max = gas velocity in the minimum flow area A_min
    (min of the transverse and diagonal gaps, fins included, as ht / fluids); bulk gas
    properties; Nu on D_t_o. Fin efficiency of Kern & Kraus (1972) (exact annular fin, Bessel
    functions). Returns the REDUCED h referred to the TOTAL outer area (tube + fins), fin
    efficiency included (labothappy convention, as the Pis'mennyi functions).

    Implementation follows ht.air_cooler.h_ESDU_high_fin (C. Bell), which cites the
    reference below; areas and minimum flow area as fluids.geometry.AirCooledExchanger;
    fin efficiency of Kern & Kraus (1972) as ht.core.fin_efficiency_Kern_Kraus.

    Validity (data of ESDU 86022 as listed in ht; warning outside): 4 to 11 fins per inch,
    D_t_o 3/8 to 2 inch, fin height 1/3 to 5/8 inch, D_f_o / D_t_o 1.2 to 2.4,
    5000 <= Re <= 50000; staggered banks; gases.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        Geometry dictionary containing:
        - 'fin_type' : str, fin type (must be 'annular')
        - 'D_t_o' : float, tube outer (root) diameter [m]
        - 'D_f_o' : float, fin outer diameter [m]
        - 't_f' : float, fin thickness [m]
        - 'f_density' : float, fin density [1/m]
        - 'L_t' : float, tube length [m]
        - 'k_f' : float, fin thermal conductivity [W/m/K]
        - 'n_rows' : float, number of tube rows [-]
        - 'n_passes' : float, number of tube passes [-]
        - 'n_t' : float, total number of tubes [-]
        - 'P_V' : float, transverse pitch (90 deg to gas flow) [m]
        - 'P_H' : float, longitudinal pitch (along gas flow) [m]
        - 'n_parallel' : float (optional), number of parallel banks [-]
        - 'theta' : float, optional, layout angle [deg]: 0/90 inline, 30/45/60 staggered (default 45)
    cond : dict
        Conditions dictionary containing:
        - 'm_dot' : float, total gas mass flow rate through the finned side [kg/s]
        - 'rho' : float, bulk gas density [kg/m3]
        - 'mu' : float, bulk gas dynamic viscosity [Pa s]
        - 'k' : float, bulk gas thermal conductivity [W/m/K]
        - 'Pr' : float, bulk gas Prandtl number [-]

    Returns
    -------
    h_reduced : float
        Reduced heat transfer coefficient referred to the total outer area
        (tube + fins), fin efficiency included [W/m2/K].
        (= ht h_bare_tube_basis / A_increase, where ht returns h on the bare-tube area).

    Reference
    ---------
    ESDU 86022 (1986) High-fin staggered tube banks: heat transfer and pressure drop
    for turbulent single phase gas flow; as given in Hewitt, G.F., Shires, G.L.,
    Bott, T.R. (1994) Process Heat Transfer, CRC Press.
    """
    # Step 1. Check the fin shape
    fin_type = str(geom.get('fin_type', 'Annular')).strip().lower()
    if fin_type != 'annular':
        raise ValueError(
            "TFHE: correlation is for circular (annular) fins only; "
            f"received geom['fin_type'] = {geom.get('fin_type')!r}."
        )

    # Step 2. Inputs (SI units)
    D = geom['D_t_o']
    D_fin = geom['D_f_o']
    t_fin = geom['t_f']
    fin_density = geom['f_density']
    L_tube = geom['L_t']
    k_f = geom['k_f']

    tube_rows = geom['n_rows'] * geom['n_passes']
    tubes_per_row = geom['n_t'] / tube_rows
    pitch_normal = geom['P_V']
    pitch_parallel = geom['P_H']
    n_parallel = geom.get('n_parallel', 1)

    m_dot = cond['m_dot'] / n_parallel
    rho = cond['rho']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Validity checks (warnings.warn for each limit exceeded, compute anyway)
    # Fin density range: 4 to 11 fins per inch -> 157.48 to 433.07 [1/m]
    if not (157.48 <= fin_density <= 433.07):
        warnings.warn(
            f"Fin density {fin_density} [1/m] is outside the validity range "
            "[157.48, 433.07] [1/m] (ESDU 86022).", ValidityRangeWarning
        )

    # Tube outer diameter range: 3/8 to 2 inch -> 9.525 mm to 50.8 mm
    if not (0.009525 <= D <= 0.0508):
        warnings.warn(
            f"Tube outer diameter {D} [m] is outside the validity range "
            "[0.009525, 0.0508] [m] (ESDU 86022).", ValidityRangeWarning
        )

    fin_height = (D_fin - D) / 2
    # Fin height range: 1/3 to 5/8 inch -> 8.467 mm to 15.875 mm
    if not (0.008467 <= fin_height <= 0.015875):
        warnings.warn(
            f"Fin height {fin_height} [m] is outside the validity range "
            "[0.008467, 0.015875] [m] (ESDU 86022).", ValidityRangeWarning
        )

    fin_ratio = D_fin / D
    # Fin tip to root diameter ratio range: 1.2 to 2.4
    if not (1.2 <= fin_ratio <= 2.4):
        warnings.warn(
            f"Fin tip to root diameter ratio {fin_ratio} is outside the "
            "validity range [1.2, 2.4] (ESDU 86022).", ValidityRangeWarning
        )

    # Tube arrangement check: staggered banks only
    # Tube arrangement from the layout angle theta [deg] (default 45, i.e. staggered):
    # 0 or 90 = inline, 30, 45 or 60 = staggered
    theta = geom.get('theta', 45)
    if theta == 0 or theta == 90:
        tube_arrang = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        tube_arrang = 'staggered'
    else:
        raise ValueError(f"TFHE: theta = {theta!r} deg: expected 0 or 90 (inline), "
                         "30, 45 or 60 (staggered).")
    if tube_arrang == 'inline':
        warnings.warn(
            "Tube arrangement is 'inline', but ESDU 86022 is valid for "
            "staggered banks only.", ValidityRangeWarning
        )

    # Step 3. Fin dimensions and areas of ONE tube (fin_height computed above)
    fin_interval = 1 / fin_density
    bare_length = fin_interval - t_fin
    A_tube_showing = np.pi * D * L_tube * (1 - t_fin / fin_interval)
    A_per_fin = 2 * np.pi / 4 * (D_fin**2 - D**2) + np.pi * D_fin * t_fin
    A_fin = fin_density * L_tube * A_per_fin
    A_total = A_tube_showing + A_fin

    # Step 4. Minimum flow area, maximum velocity, and Reynolds number
    pitch_diagonal = np.sqrt(pitch_parallel**2 + (pitch_normal / 2)**2)
    blockage = D + 2 * fin_density * fin_height * t_fin
    A_diagonal = 2 * tubes_per_row * L_tube * (pitch_diagonal - blockage)
    A_normal = tubes_per_row * L_tube * (pitch_normal - blockage)
    A_min = min(A_diagonal, A_normal)
    V_max = m_dot / (A_min * rho)
    Re = rho * V_max * D / mu

    # Reynolds number validity range: 5000 <= Re <= 50000
    if not (5000 <= Re <= 50000):
        warnings.warn(
            f"Reynolds number Re = {Re} is outside the validity range "
            "[5000, 50000] (ESDU 86022).", ValidityRangeWarning
        )

    # ESDU 86022 (Hewitt et al. 1994)
    # 1. Row correction factor F2
    if tube_rows < 2:
        F2 = 0.76
    elif tube_rows < 3:
        F2 = 0.84
    elif tube_rows < 4:
        F2 = 0.92
    else:
        F2 = 1.0

    # 2. Nusselt number
    term_re = Re**0.658
    term_geom = (bare_length / fin_height)**0.297
    term_pitch = (pitch_normal / pitch_parallel)**-0.091
    term_pr = Pr**(1 / 3)
    Nu = 0.242 * term_re * term_geom * term_pitch * term_pr * F2

    # 3. The wall-property correction F1 = (Pr / Pr_wall)^0.26 is NOT applied (F1 = 1:
    # the wall Prandtl number is not available in cond).

    # 4. Heat transfer coefficient of the finned surface, before the fin efficiency
    h = Nu * k / D

    # Step 5. Fin efficiency, Kern & Kraus (1972), and reduced h
    m_fin = np.sqrt(2 * h / (k_f * t_fin))
    r_e = D_fin / 2
    r_o = D / 2
    numerator = i1(m_fin * r_e) * k1(m_fin * r_o) - k1(m_fin * r_e) * i1(m_fin * r_o)
    denominator = i0(m_fin * r_o) * k1(m_fin * r_e) + i1(m_fin * r_e) * k0(m_fin * r_o)
    eta_fin = 2 * r_o / (m_fin * (r_e**2 - r_o**2)) * numerator / denominator
    h_reduced = (eta_fin * A_fin + A_tube_showing) / A_total * h

    return h_reduced


def TFHE_ganguli_vdi_HTC(AS, geom, cond):
    """Gas-side heat transfer coefficient of an in-line or staggered bank of tubes with
    circular (annular) fins, using the Ganguli et al. (1985) correlation with the VDI
    Heat Atlas (2010) modifications.

    Nu = C Re^0.6 Pr^(1/3) A_increase^-0.15, A_increase = total outer area / bare tube area.
    C (VDI): inline 0.22 (0.2 below 4 rows); staggered 0.38 (0.2, 0.33, 0.36 for 1, 2, 3 rows).

    Definitions: Re = rho V_max D_t_o / mu, V_max = gas velocity in the minimum flow area A_min
    (min of the transverse and diagonal gaps, fins included, as ht / fluids); bulk gas
    properties; Nu on D_t_o. Fin efficiency of Kern & Kraus (1972) (exact annular fin, Bessel
    functions). Returns the REDUCED h referred to the TOTAL outer area (tube + fins), fin
    efficiency included (labothappy convention, as the Pis'mennyi functions).

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        Geometry dictionary containing:
        - 'fin_type' : str, must be 'annular' (case insensitive)
        - 'D_t_o' : float, tube outer (root) diameter [m]
        - 'D_f_o' : float, fin outer diameter [m]
        - 't_f' : float, fin thickness [m]
        - 'f_density' : float, number of fins per metre of tube [1/m]
        - 'L_t' : float, tube length [m]
        - 'k_f' : float, fin thermal conductivity [W/m/K]
        - 'n_rows' : float, number of tube rows [ - ]
        - 'n_passes' : float, number of tube passes [ - ]
        - 'n_t' : float, total number of tubes [ - ]
        - 'P_V' : float, transverse pitch (90 deg to gas flow) [m]
        - 'P_H' : float, longitudinal pitch (along gas flow) [m]
        - 'theta' : float, layout angle [deg]: 0/90 inline, 30/45/60 staggered (note: ht guesses
          arrangement from pitch ratios close to 1, while labothappy reads geom['theta'])
        - 'n_parallel' : float, optional, number of parallel passes [ - ]
    cond : dict
        Condition dictionary containing:
        - 'm_dot' : float, whole finned side gas mass flow rate [kg/s]
        - 'rho' : float, bulk gas density [kg/m3]
        - 'mu' : float, bulk gas dynamic viscosity [Pa s]
        - 'k' : float, bulk gas thermal conductivity [W/m/K]
        - 'Pr' : float, bulk gas Prandtl number [ - ]

    Returns
    -------
    h_reduced : float
        Reduced heat transfer coefficient referred to the total outer area
        (tube + fins), fin efficiency included [W/m2/K].
        Calculated as: h_reduced = ht h_bare_tube_basis / A_increase,
        where ht returns h on the bare-tube area.

    Reference
    ---------
    Ganguli, A., Tung, S.S., Taborek, J. (1985) Parametric study of air-cooled heat
    exchanger finned tube geometry. AIChE Symposium Series 81, 122-128; with the
    modifications of the VDI Heat Atlas, 2nd edition (2010).

    Implementation
    --------------
    Implementation follows ht.air_cooler.h_Ganguli_VDI (C. Bell), which cites the
    reference above; areas and minimum flow area as fluids.geometry.AirCooledExchanger;
    fin efficiency of Kern & Kraus (1972) as ht.core.fin_efficiency_Kern_Kraus.

    Validity
    --------
    Validity range: VDI Heat Atlas (2010) M1 (Eqs. 15-18): 1e3 < Re_d < 1e5,
    5 <= A/A_t0 <= 30 (accuracy +-10 % to +-25 %), industrial data; warning outside. Air.
    """
    # Step 1. Check the fin shape
    fin_type = str(geom.get('fin_type', 'Annular')).strip().lower()
    if fin_type != 'annular':
        raise ValueError(
            f"TFHE Ganguli VDI HTC: correlation is for circular (annular) fins only. "
            f"Received geom['fin_type'] = {geom.get('fin_type')!r}."
        )

    # Step 2. Inputs (SI units)
    D = geom['D_t_o']
    D_fin = geom['D_f_o']
    t_fin = geom['t_f']
    fin_density = geom['f_density']
    L_tube = geom['L_t']
    k_f = geom['k_f']
    tube_rows = geom['n_rows'] * geom['n_passes']
    tubes_per_row = geom['n_t'] / tube_rows
    pitch_normal = geom['P_V']
    pitch_parallel = geom['P_H']
    n_parallel = geom.get('n_parallel', 1)
    m_dot = cond['m_dot'] / n_parallel
    rho = cond['rho']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Read arrangement from geometry
    # Tube arrangement from the layout angle theta [deg] (default 45, i.e. staggered):
    # 0 or 90 = inline, 30, 45 or 60 = staggered
    theta = geom.get('theta', 45)
    if theta == 0 or theta == 90:
        arrangement = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        arrangement = 'staggered'
    else:
        raise ValueError(f"TFHE: theta = {theta!r} deg: expected 0 or 90 (inline), "
                         "30, 45 or 60 (staggered).")
    if arrangement not in ('inline', 'staggered'):
        raise ValueError(
            f"TFHE Ganguli VDI HTC: inline or staggered bank expected, "
            f"got theta = {theta!r} deg."
        )

    # Step 3. Fin dimensions and areas of ONE tube
    fin_height = (D_fin - D) / 2
    fin_interval = 1 / fin_density
    A_bare_tube = np.pi * D * L_tube
    A_tube_showing = np.pi * D * L_tube * (1 - t_fin / fin_interval)
    A_per_fin = 2 * np.pi / 4 * (D_fin**2 - D**2) + np.pi * D_fin * t_fin
    A_fin = fin_density * L_tube * A_per_fin
    A_total = A_tube_showing + A_fin
    A_increase = A_total / A_bare_tube

    # Step 4. Minimum flow area and maximum velocity
    pitch_diagonal = np.sqrt(pitch_parallel**2 + (pitch_normal / 2)**2)
    blockage = D + 2 * fin_density * fin_height * t_fin
    A_diagonal = 2 * tubes_per_row * L_tube * (pitch_diagonal - blockage)
    A_normal = tubes_per_row * L_tube * (pitch_normal - blockage)
    A_min = min(A_diagonal, A_normal)
    V_max = m_dot / (A_min * rho)
    Re = rho * V_max * D / mu

    # Step 5. VDI coefficient based on arrangement and number of rows
    if arrangement == 'inline':
        if tube_rows < 4:
            coefficient = 0.2
        else:
            coefficient = 0.22
    else:
        if tube_rows < 1.5:
            coefficient = 0.2
        elif tube_rows < 2.5:
            coefficient = 0.33
        elif tube_rows < 3.5:
            coefficient = 0.36
        else:
            coefficient = 0.38

    # Nusselt number correlation (Ganguli et al. / VDI)
    Nu = coefficient * Re**0.6 * Pr**(1 / 3) * A_increase**-0.15

    # Heat transfer coefficient of the finned surface, before the fin efficiency
    h = Nu * k / D

    # Step 6. Fin efficiency, Kern & Kraus (1972), and reduced h
    m_fin = np.sqrt(2 * h / (k_f * t_fin))
    r_e = D_fin / 2
    r_o = D / 2

    numerator = i1(m_fin * r_e) * k1(m_fin * r_o) - k1(m_fin * r_e) * i1(m_fin * r_o)
    denominator = i0(m_fin * r_o) * k1(m_fin * r_e) + i1(m_fin * r_e) * k0(m_fin * r_o)
    eta_fin = 2 * r_o / (m_fin * (r_e**2 - r_o**2)) * numerator / denominator

    h_reduced = (eta_fin * A_fin + A_tube_showing) / A_total * h

    # Range of validity: VDI Heat Atlas (2010) M1, Eqs. (15)-(16) "fitted within +-10 % to +-25 %
    # for 1e3 < Re_d < 1e5 and 5 <= A/A_t0 <= 30" (Re_d on d_0 and the velocity in the smallest
    # cross section, as here). Range check added 2026-11-27 (warning only).
    out = []
    if Re < 1e3 or Re > 1e5:
        out.append(f"Re_d = {Re:.4g} outside 1e3-1e5")
    if A_increase < 5 or A_increase > 30:
        out.append(f"A/A_t0 = {A_increase:.3g} outside 5-30")
    if out:
        warnings.warn("Ganguli VDI (VDI Heat Atlas M1): outside the range of validity: "
                      + ", ".join(out) + ".", ValidityRangeWarning)

    return h_reduced


#%% CONTINUOUS PLAIN PLATE FINS (fin_type 'Square'): Wang, Chi, Chang (2000), Schmidt fin efficiency

def TFHE_wang_plain_fin_HTC(AS, geom, cond):
    """
    Gas-side heat transfer coefficient of a staggered bank of round tubes with
    continuous plain plate fins (Wang, Chi, Chang 2000). The correlation includes
    the fin-collar contact conductance. Returns the REDUCED heat transfer coefficient
    referred to the TOTAL outer area (fins + exposed tube), surface efficiency included
    (labothappy convention, as the Pis'mennyi functions): h_reduced = eta_o * h.

    Definitions: Re_Dc = rho * u_max * D_c / mu, u_max = gas velocity in the minimum
    free flow area A_c, D_c = D_t_o + 2 * t_f (fin collar diameter). Bulk gas
    properties. j = Colburn factor = h * Pr^(2/3) / (rho * u_max * cp).
    Fin efficiency: Schmidt (1949) equivalent circular fin, staggered layout (ACHP then
    replaces it by a Perrotin-Clodic (2006) form: NOT used here).
    Independent implementation: tested against ACHP PlainFins (h before efficiency)
    and Schmidt by hand.

    Equations as implemented in ACHP (Purdue), ACHP/FinCorrelations.py, function PlainFins.

    Validity: unknown (TODO, paper not available); air-side coils with small tubes
    (D_c about 6-13 mm).

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        Dictionary of geometry parameters with SI units:
        - 'fin_type' : str, fin type (must be 'Square' for continuous plate fins)
        - 'D_t_o' : float, tube outer diameter [m]
        - 't_f' : float, fin thickness [m]
        - 'f_density' : float, number of fins per metre of tube [1/m]
        - 'L_t' : float, tube length [m]
        - 'k_f' : float, fin thermal conductivity [W/m/K]
        - 'n_rows' : float or int, number of tube rows [-]
        - 'n_passes' : float or int, number of tube passes [-]
        - 'n_t' : float or int, total number of tubes [-]
        - 'P_V' : float, transverse pitch (normal to gas flow) [m]
        - 'P_H' : float, longitudinal pitch (along gas flow) [m]
        - 'n_parallel' : float or int, optional, number of parallel paths [-]
        - 'theta' : float, optional, layout angle [deg]: 0/90 inline, 30/45/60 staggered (default 45)
    cond : dict
        Dictionary of operating conditions with SI units (bulk gas properties):
        - 'm_dot' : float, total gas mass flow rate through the finned side [kg/s]
        - 'rho' : float, bulk gas density [kg/m3]
        - 'mu' : float, bulk gas dynamic viscosity [Pa s]
        - 'k' : float, bulk gas thermal conductivity [W/m/K]
        - 'Pr' : float, bulk gas Prandtl number [-]

    Returns
    -------
    h_reduced : float
        Reduced heat transfer coefficient referred to the total outer area
        (fins + exposed tube), surface efficiency included [W/m2/K].

    Reference
    ---------
    Wang, C.-C., Chi, K.-Y., Chang, C.-J. (2000) Heat transfer and friction
    characteristics of plain fin-and-tube heat exchangers, part II: Correlation.
    International Journal of Heat and Mass Transfer (2000). Volume / pages / DOI:
    TODO: check (paper not available).

    Fin efficiency: Schmidt, T.E. (1949) Heat transfer calculations for extended
    surfaces. Refrigerating Engineering 57, 351-357 (as cited by ACHP; not checked).
    """
    # Step 1. Check fin type (must be continuous plate fins, 'Square' in labothappy)
    fin_type = str(geom.get('fin_type', 'Square')).strip().lower()
    if fin_type != 'square':
        warnings.warn(
            "Wang plain fin: correlation for continuous plate fins (fin_type 'Square'). "
            f"Received geom['fin_type'] = {geom.get('fin_type')!r}.",
            ValidityRangeWarning
        )

    # Step 2. Check tube arrangement (staggered banks only)
    # Tube arrangement from the layout angle theta [deg] (default 45, i.e. staggered):
    # 0 or 90 = inline, 30, 45 or 60 = staggered
    theta = geom.get('theta', 45)
    if theta == 0 or theta == 90:
        tube_arrangement = 'inline'
    elif theta == 30 or theta == 45 or theta == 60:
        tube_arrangement = 'staggered'
    else:
        raise ValueError(f"TFHE: theta = {theta!r} deg: expected 0 or 90 (inline), "
                         "30, 45 or 60 (staggered).")
    if tube_arrangement == 'inline':
        warnings.warn(
            "Wang plain fin: staggered banks only. "
            f"Received geom['theta'] = {theta!r} deg (inline).",
            ValidityRangeWarning
        )

    # Step 3. Read geometry parameters (SI units)
    D_o = geom['D_t_o']              # tube outer diameter [m]
    t_fin = geom['t_f']              # fin thickness [m]
    f_density = geom['f_density']      # fins per metre [1/m]
    L_tube = geom['L_t']            # tube length [m]
    k_f = geom['k_f']              # fin thermal conductivity [W/m/K]
    tube_rows = geom['n_rows'] * geom['n_passes']  # number of tube rows crossed by gas [-]
    tubes_per_row = geom['n_t'] / tube_rows     # tubes per row [-]
    P_t = geom['P_V']              # transverse pitch [m]
    P_l = geom['P_H']              # longitudinal pitch [m]
    n_parallel = geom.get('n_parallel', 1)          # number of parallel paths [-]

    # Step 4. Read operating conditions (SI units)
    m_dot_total = cond['m_dot']        # total mass flow rate [kg/s]
    rho = cond['rho']                  # bulk gas density [kg/m3]
    mu = cond['mu']                    # bulk gas dynamic viscosity [Pa s]
    k = cond['k']                      # bulk gas thermal conductivity [W/m/K]
    Pr = cond['Pr']                    # bulk gas Prandtl number [-]

    # Step 5. Specific heat capacity from Prandtl number
    cp = Pr * k / mu                   # specific heat [J/kg/K]

    # Step 6. Mass flow rate per parallel path
    m_dot = m_dot_total / n_parallel   # [kg/s]

    # Step 7. Geometry calculations (ACHP PlainFins order)
    # 7.1 Fin pitch and fin collar diameter
    F_p = 1.0 / f_density              # fin pitch [m]
    D_c = D_o + 2.0 * t_fin            # fin collar outside diameter [m]

    # 7.2 Number of fins
    N_fin = L_tube * f_density         # number of fins [-]

    # 7.3 Fin sheet height (extends beyond outer tubes)
    height = P_t * (tubes_per_row + 1) # [m]

    # 7.4 Face area
    A_duct = height * L_tube           # face area [m2]

    # 7.5 Minimum free flow area
    A_c = (A_duct
           - t_fin * N_fin * (height - D_c * tubes_per_row)
           - tubes_per_row * D_c * L_tube)  # [m2]

    # 7.6 Maximum velocity and Reynolds number on collar diameter
    u_max = m_dot / (rho * A_c)        # [m/s]
    Re_Dc = rho * u_max * D_c / mu     # [-]

    # 7.7 Flow depth
    depth = (tube_rows + 1) * P_l      # [m]

    # 7.8 Area of one fin (both sides, minus tube cutouts)
    A_one_fin = 2.0 * (height * P_l * (tube_rows + 1)
                       - tubes_per_row * tube_rows * np.pi * D_c**2 / 4.0)

    # 7.9 Total fin area
    A_f = N_fin * A_one_fin            # [m2]

    # 7.10 Total outer area (fins + exposed tube)
    A_o = (A_f
           + tubes_per_row * tube_rows * np.pi * D_c * (L_tube - N_fin * t_fin))  # [m2]

    # 7.11 Hydraulic diameter
    D_h = 4.0 * A_c * depth / A_o      # [m]

    # Step 8. Colburn j factor (Wang et al. 2000, as in ACHP)
    if tube_rows == 1:
        # Single row correlation
        log_Re = np.log(Re_Dc)
        P1 = 1.9 - 0.23 * log_Re
        P2 = -0.236 + 0.126 * log_Re

        term1 = 0.108
        term2 = Re_Dc ** -0.29
        term3 = (P_t / P_l) ** P1
        term4 = (F_p / D_c) ** -1.084
        term5 = (F_p / D_h) ** -0.786
        term6 = (F_p / P_t) ** P2
        j = term1 * term2 * term3 * term4 * term5 * term6
    else:
        # Multi-row correlation (tube_rows >= 2)
        N = tube_rows
        log_Re = np.log(Re_Dc)

        P3 = -0.361 - 0.042 * N / log_Re + 0.158 * np.log(N * (F_p / D_c) ** 0.41)
        P4 = -1.224 - 0.076 * (P_l / D_h) ** 1.42 / log_Re
        P5 = -0.083 + 0.058 * N / log_Re
        P6 = -5.735 + 1.21 * np.log(Re_Dc / N)

        term1 = 0.086
        term2 = Re_Dc ** P3
        term3 = N ** P4
        term4 = (F_p / D_c) ** P5
        term5 = (F_p / D_h) ** P6
        term6 = (F_p / P_t) ** -0.93
        j = term1 * term2 * term3 * term4 * term5 * term6

    # Step 9. Convective heat transfer coefficient on total outer area A_o
    h = j * rho * u_max * cp / (Pr ** (2.0 / 3.0))  # [W/m2/K]

    # Step 10. Fin efficiency (Schmidt 1949, staggered layout, as in ACHP)
    # 10.1 Tube radius (uses tube OD, not collar diameter)
    r = D_o / 2.0

    # 10.2 Equivalent fin geometry
    X_D = np.sqrt(P_l**2 + P_t**2 / 4.0) / 2.0
    X_T = P_t / 2.0

    # 10.3 Equivalent fin radius ratio
    ratio = 1.27 * X_T / r * np.sqrt(X_D / X_T - 0.3)

    # 10.4 Fin parameter m
    m = np.sqrt(2.0 * h / (k_f * t_fin))

    # 10.5 Phi factor
    phi = (ratio - 1.0) * (1.0 + 0.35 * np.log(ratio))

    # 10.6 Fin efficiency
    m_r_phi = m * r * phi
    if m_r_phi == 0.0:
        eta_fin = 1.0
    else:
        eta_fin = np.tanh(m_r_phi) / m_r_phi

    # 10.7 Overall surface efficiency
    eta_o = 1.0 - (A_f / A_o) * (1.0 - eta_fin)

    # Step 11. Reduced heat transfer coefficient (referred to total outer area)
    h_reduced = eta_o * h

    return h_reduced
