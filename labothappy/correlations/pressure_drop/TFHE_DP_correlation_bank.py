# -*- coding: utf-8 -*-
"""
Finned-side pressure drop of tube & fins heat exchangers (TFHE), (AS, geom, cond) formalism.

Source : Handbook for transversely finned tube heat exchanger design,
         E. Pis'mennyi, G. Polupan, I. Carvajal-Mariscal, F. Sanchez-Silva, I. Pioro.
         Staggered bank of transversely finned tubes (annular or square fins: same formula).

Returns the pressure drop [Pa] of the WHOLE bank (all rows crossed by the gas), friction only.

geom : D_t_o, L_t, D_f_o, t_f, f_density, n_t, n_rows, n_passes, w (bundle width),
       P_V, P_H, n_parallel and n_series (optional, default 1)
cond : m_dot (WHOLE finned side), rho, mu (gas bulk properties)
"""
import warnings
import numpy as np
from labothappy.correlations.convection.correlation_warnings import ValidityRangeWarning

_NEEDED = ('D_t_o', 'L_t', 'D_f_o', 't_f', 'f_density', 'n_t', 'n_rows', 'n_passes', 'w',
           'P_V', 'P_H')


def TFHE_pismennyi_DP(AS, geom, cond):
    miss = [k for k in _NEEDED if k not in geom] + [k for k in ('m_dot', 'rho', 'mu') if k not in cond]
    if miss:
        raise KeyError(f"TFHE DP: missing {miss} in geom / cond.")
    D, t, L = geom['D_t_o'], geom['t_f'], geom['L_t']
    S_1, S_2 = geom['P_V'], geom['P_H']
    S_m = np.sqrt(S_1 * S_2)
    N_rows = geom['n_rows']
    rho, mu = cond['rho'], cond['mu']

    "Gas velocity - flow area"
    Fin_L = (geom['D_f_o'] - D) / 2
    n_f = L * geom['f_density'] - 1
    s = (L - n_f * t) / (n_f + 1)                                   # fin spacing
    n_tpr = geom['n_t'] / (N_rows * geom['n_passes'])                # tubes per row
    D_fin_c = D + 2 * Fin_L * t / s                                       # conventional fin diameter
    diag = np.sqrt(2) * S_m
    psi_c = (S_m - D_fin_c) / (diag - D_fin_c)
    S_free = L * geom['w'] - n_tpr * D_fin_c * L
    S_flow = S_free * (2 / psi_c) if psi_c > 2 else S_free
    u = cond['m_dot'] / geom.get('n_parallel', 1) / rho / S_flow

    "Resistance coefficient of one tube row"
    fact_1 = np.pi * (D * s + 2 * Fin_L * t + 2 * Fin_L * (Fin_L + D))
    fact_2 = S_m * s - (D * s + 2 * Fin_L * t)
    A_totF = fact_1 / fact_2                                              # reduced developed surface
    n = 0.17 * A_totF ** 0.25 * (S_1 / S_2) ** 0.57 * np.exp(-0.36 * (S_1 / S_2))
    C_r = 2.8 * A_totF ** 0.53 * (S_1 / S_2) ** 1.3 * np.exp(-0.9 * (S_1 / S_2))
    C_z = 1.0 if N_rows >= 6 else np.exp(0.1 * (6 / N_rows - 1))
    D_eq = 2 * (s * (S_m - D) - 2 * Fin_L * t) / (2 * Fin_L + s)
    xi_o = C_z * C_r * (u * D_eq / mu) ** (-n)

    C_op = 1.1       # real operating conditions of the heat-transfer surface (value of the handbook)
    return geom.get('n_series', 1) * C_op * xi_o * geom['n_passes'] * N_rows * rho * u ** 2 / 2




#%% CIRCULAR (ANNULAR) FINS: ESDU high fin

def TFHE_esdu_high_fin_DP(AS, geom, cond):
    """
    Gas-side pressure drop of a staggered bank of high-fin tubes with circular (annular) fins,
    using the ESDU method.

    DP = n_series (K_acc + n_rows_crossed K_f) rho V_max^2 / 2
    K_f = 4.567 Re^-0.242 A_increase^0.504 (P_normal / D_t_o)^-0.376 (P_parallel / D_t_o)^-0.546
    K_acc = 1 + (A_min / A_face)^2 (entrance + exit), A_face = w L_t (ht / fluids use
    P_normal L_t (tubes_per_row + 1), the same for a bundle without corbels).

    Definitions: Re = rho V_max D_t_o / mu, V_max = gas velocity in the minimum flow area A_min
    (min of the transverse and diagonal gaps, fins included, as ht / fluids); bulk gas
    properties. Returns the pressure drop [Pa] of the WHOLE bank (all rows) INCLUDING the
    entrance / exit term K_acc (the Pis'mennyi DP of this file is friction only).

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        Geometry dictionary with the following keys (SI units):
        - 'fin_type' : str, expected 'annular' (case-insensitive)
        - 'D_t_o' : float, tube outer (root) diameter [m]
        - 'D_f_o' : float, fin outer diameter [m]
        - 't_f' : float, fin thickness [m]
        - 'f_density' : float, number of fins per metre [1/m]
        - 'L_t' : float, tube length [m]
        - 'n_rows' : float, number of tube rows [-]
        - 'n_passes' : float, number of tube passes [-]
        - 'n_t' : float, total number of tubes [-]
        - 'P_V' : float, transverse pitch (90 deg to gas flow) [m]
        - 'P_H' : float, longitudinal pitch (along gas flow) [m]
        - 'w' : float, bundle width [m]
        - 'n_parallel' : float, optional, number of parallel banks [-]
        - 'n_series' : float, optional, number of banks in series [-]
        - 'theta' : float, optional, layout angle [deg]: 0/90 inline, 30/45/60 staggered (default 45)
    cond : dict
        Condition dictionary with the following keys (SI units):
        - 'm_dot' : float, total gas mass flow rate through the finned side [kg/s]
        - 'rho' : float, bulk gas density [kg/m3]
        - 'mu' : float, bulk gas dynamic viscosity [Pa s]

    Returns
    -------
    float
        Total pressure drop [Pa] of the whole bank (all rows) including the entrance/exit
        acceleration term.

    Reference
    ---------
    ESDU 86022 (1986) High-fin staggered tube banks: heat transfer and pressure drop for turbulent
    single phase gas flow; as given in Hewitt, G.F., Shires, G.L., Bott, T.R. (1994) Process Heat
    Transfer, CRC Press.
    Implementation follows ht.air_cooler.dP_ESDU_high_fin (C. Bell), which cites the reference
    above; areas and minimum flow area as fluids.geometry.AirCooledExchanger.

    Validity Range
    --------------
    - 157.48 <= fin_density <= 433.07 [1/m] (4 to 11 fins per inch)
    - 9.525 mm <= D <= 50.8 mm (3/8 to 2 inch)
    - 8.467 mm <= fin_height <= 15.875 mm (1/3 to 5/8 inch)
    - 1.2 <= D_fin / D <= 2.4
    - 5000 <= Re <= 50000
    - Staggered banks only.
    - Fluids: gases.
    """
    # Step 1. Check the fin shape
    fin_type = str(geom.get('fin_type', 'Annular')).strip().lower()
    if fin_type != 'annular':
        raise ValueError(
            f"TFHE DP: correlation is for circular (annular) fins only. "
            f"Received geom['fin_type'] = {geom.get('fin_type')!r}."
        )

    # Check tube arrangement validity
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
    if tube_arrang != 'staggered':
        warnings.warn(
            f"TFHE DP: ESDU 86022 correlation is for staggered tube banks only. "
            f"Received theta = {theta!r} deg.",
            ValidityRangeWarning
        )

    # Step 2. Inputs (SI units)
    D = geom['D_t_o']
    D_fin = geom['D_f_o']
    t_fin = geom['t_f']
    fin_density = geom['f_density']
    L_tube = geom['L_t']
    tube_rows = geom['n_rows'] * geom['n_passes']
    tubes_per_row = geom['n_t'] / tube_rows
    pitch_normal = geom['P_V']
    pitch_parallel = geom['P_H']
    n_parallel = geom.get('n_parallel', 1)
    m_dot = cond['m_dot'] / n_parallel
    rho = cond['rho']
    mu = cond['mu']

    # Step 3. Fin dimensions and areas of ONE tube
    fin_height = (D_fin - D) / 2
    fin_interval = 1 / fin_density
    A_bare_tube = np.pi * D * L_tube
    A_tube_showing = np.pi * D * L_tube * (1 - t_fin / fin_interval)
    A_per_fin = 2 * np.pi / 4 * (D_fin**2 - D**2) + np.pi * D_fin * t_fin
    A_fin = fin_density * L_tube * A_per_fin
    A_total = A_tube_showing + A_fin
    A_increase = A_total / A_bare_tube

    # Step 4. Minimum flow area, maximum velocity, and Reynolds number
    pitch_diagonal = np.sqrt(pitch_parallel**2 + (pitch_normal / 2)**2)
    blockage = D + 2 * fin_density * fin_height * t_fin
    A_diagonal = 2 * tubes_per_row * L_tube * (pitch_diagonal - blockage)
    A_normal = tubes_per_row * L_tube * (pitch_normal - blockage)
    A_min = min(A_diagonal, A_normal)
    V_max = m_dot / (A_min * rho)
    Re = rho * V_max * D / mu

    # Validity checks and warnings
    if not (157.48 <= fin_density <= 433.07):
        warnings.warn(
            f"TFHE DP: fin_density = {fin_density} [1/m] is outside the validity range "
            f"of [157.48, 433.07] [1/m].",
            ValidityRangeWarning
        )

    if not (0.009525 <= D <= 0.0508):
        warnings.warn(
            f"TFHE DP: tube outer diameter D = {D} [m] is outside the validity range "
            f"of [0.009525, 0.0508] [m].",
            ValidityRangeWarning
        )

    if not (0.008467 <= fin_height <= 0.015875):
        warnings.warn(
            f"TFHE DP: fin_height = {fin_height} [m] is outside the validity range "
            f"of [0.008467, 0.015875] [m].",
            ValidityRangeWarning
        )

    D_fin_ratio = D_fin / D
    if not (1.2 <= D_fin_ratio <= 2.4):
        warnings.warn(
            f"TFHE DP: fin diameter ratio D_fin / D = {D_fin_ratio} is outside the "
            f"validity range of [1.2, 2.4].",
            ValidityRangeWarning
        )

    if not (5000 <= Re <= 50000):
        warnings.warn(
            f"TFHE DP: Reynolds number Re = {Re} is outside the validity range "
            f"of [5000, 50000].",
            ValidityRangeWarning
        )

    # Step 1 (specific to DP): Face area of the bank
    # Bundle width * tube length (without corbels)
    w = geom['w']
    A_face = w * L_tube

    # Step 2 (specific to DP): Contraction ratio
    contraction_ratio = A_min / A_face

    # Step 3 (specific to DP): Friction coefficient per row
    term_re = 4.567 * (Re ** -0.242)
    term_area = A_increase ** 0.504
    term_transverse = (pitch_normal / D) ** -0.376
    term_parallel = (pitch_parallel / D) ** -0.546
    K_f = term_re * term_area * term_transverse * term_parallel

    # Step 4 (specific to DP): Entrance + exit (acceleration) coefficient
    K_acc = 1 + (contraction_ratio ** 2)

    # Step 5 (specific to DP): Total pressure drop of the bank
    dynamic_head = rho * (V_max ** 2) / 2
    DP_bank = (K_acc + tube_rows * K_f) * dynamic_head

    # Step 6 (specific to DP): Return pressure drop multiplied by banks in series
    n_series = geom.get('n_series', 1)
    return n_series * DP_bank


#%% CONTINUOUS PLAIN PLATE FINS (fin_type 'Square'): Wang, Chi, Chang (2000)

def TFHE_wang_plain_fin_DP(AS, geom, cond):
    """
    Gas-side frictional pressure drop of the whole bank (all rows) for a staggered bank of round
    tubes with continuous plain plate fins, using the Wang et al. (2000) correlation.

    Definitions: f is the FANNING friction factor; Re_Dc = rho u_max D_c / mu on the fin collar
    diameter D_c = D_t_o + 2 t_f, u_max = gas velocity in the minimum free flow area A_c; bulk
    gas properties. DP = f (A_o / A_c) G_c^2 / (2 rho), with A_o the total outer area (fins +
    exposed tube) and G_c = m_dot / A_c: frictional core pressure drop of the WHOLE bank (all rows),
    entrance / exit / acceleration NOT included.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Not used, kept for the common (AS, geom, cond) interface.
    geom : dict
        Geometry dictionary with the following keys (SI units):
        - 'fin_type' : str, expected 'Square' (case-insensitive)
        - 'theta' : float, tube layout angle [deg], expected 30, 45 or 60 (staggered)
        - 'D_t_o' : float, tube outer diameter [m]
        - 't_f' : float, fin thickness [m]
        - 'f_density' : float, number of fins per metre [1/m]
        - 'L_t' : float, tube length [m]
        - 'n_rows' : float, number of tube rows [-]
        - 'n_passes' : float, number of tube passes [-]
        - 'n_t' : float, total number of tubes [-]
        - 'P_V' : float, transverse pitch (normal to gas flow) [m]
        - 'P_H' : float, longitudinal pitch (along gas flow) [m]
        - 'n_parallel' : float, optional, number of parallel banks [-]
        - 'n_series' : float, optional, number of banks in series [-]
    cond : dict
        Condition dictionary with the following keys (SI units):
        - 'm_dot' : float, total gas mass flow rate through the finned side [kg/s]
        - 'rho' : float, bulk gas density [kg/m3]
        - 'mu' : float, bulk gas dynamic viscosity [Pa s]

    Returns
    -------
    float
        Frictional core pressure drop [Pa] of the whole bank (entrance / exit / acceleration
        not included).

    Reference
    ---------
    Wang, C.-C., Chi, K.-Y., Chang, C.-J. (2000) Heat transfer and friction characteristics of plain
    fin-and-tube heat exchangers, part II: Correlation. International Journal of Heat and Mass
    Transfer. TODO: check volume, pages, and DOI (paper not available).
    Equations as implemented in ACHP (Purdue), ACHP/FinCorrelations.py, function PlainFins.

    Validity Range
    --------------
    Validity: unknown (TODO, paper not available); air-side coils with small tubes
    (D_c about 6-13 mm). Staggered banks only, continuous plate fins.
    """
    # Check fin type geometry
    fin_type = str(geom.get('fin_type', 'Square')).strip().lower()
    if fin_type != 'square':
        warnings.warn(
            "Wang plain fin: correlation for continuous plate fins (fin_type 'Square').",
            ValidityRangeWarning
        )

    # Check tube arrangement
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
            "Wang plain fin: staggered banks only.",
            ValidityRangeWarning
        )

    # Read geometry parameters in SI units
    D_o = geom['D_t_o']
    t_fin = geom['t_f']
    f_density = geom['f_density']
    L_tube = geom['L_t']
    tube_rows = geom['n_rows'] * geom['n_passes']
    tubes_per_row = geom['n_t'] / tube_rows
    P_t = geom['P_V']
    P_l = geom['P_H']
    n_parallel = geom.get('n_parallel', 1)

    # Read bulk gas properties
    rho = cond['rho']
    mu = cond['mu']
    m_dot = cond['m_dot'] / n_parallel

    # Step 1. Fin pitch [m] and collar outside diameter [m]
    F_p = 1 / f_density
    D_c = D_o + (2 * t_fin)

    # Step 2. Number of fins [-]
    N_fin = L_tube * f_density

    # Step 3. Fin sheet height [m] (ACHP: the fin sheet extends beyond the outer tubes)
    height = P_t * (tubes_per_row + 1)

    # Step 4. Face area [m2]
    A_duct = height * L_tube

    # Step 5. Minimum free flow area [m2]
    term_fin_block = t_fin * N_fin * (height - (D_c * tubes_per_row))
    term_tube_block = tubes_per_row * D_c * L_tube
    A_c = A_duct - term_fin_block - term_tube_block

    # Step 6. Maximum gas velocity [m/s] and Reynolds number based on D_c [-]
    u_max = m_dot / (rho * A_c)
    Re_Dc = (rho * u_max * D_c) / mu

    # Step 7. Total surface area of one fin [m2]
    term_fin_geom = height * P_l * (tube_rows + 1)
    term_tube_geom = tubes_per_row * tube_rows * np.pi * (D_c ** 2) / 4
    A_one_fin = 2 * (term_fin_geom - term_tube_geom)

    # Step 8. Total fin area [m2]
    A_f = N_fin * A_one_fin

    # Step 9. Total outer heat transfer area [m2]
    term_bare_tube = tubes_per_row * tube_rows * np.pi * D_c * (L_tube - (N_fin * t_fin))
    A_o = A_f + term_bare_tube

    # Fanning friction factor coefficients (Wang et al. 2000, as in ACHP)
    N = tube_rows
    log_Re_Dc = np.log(Re_Dc)

    # Calculate intermediate friction factor parameters F1, F2, F3
    F1 = -0.764 + (0.739 * P_t / P_l) + (0.177 * F_p / D_c) - (0.00758 / N)
    F2 = -15.689 + (64.021 / log_Re_Dc)
    F3 = 1.696 - (15.695 / log_Re_Dc)

    # Fanning friction factor f [-]
    term_re_part = Re_Dc ** F1
    term_pitch_part = (P_t / P_l) ** F2
    term_fin_part = (F_p / D_c) ** F3
    f = 0.0267 * term_re_part * term_pitch_part * term_fin_part

    # Core mass flux [kg/m2/s]
    G_c = m_dot / A_c

    # Frictional core pressure drop of the whole bank [Pa]
    DP_bank = f * (A_o / A_c) * (G_c ** 2) / (2 * rho)

    # Return total pressure drop scaled by series banks
    n_series = geom.get('n_series', 1)
    return n_series * DP_bank
