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

geom : Tube_OD, Tube_t, Tube_L, Fin_OD, Fin_t, Fin_per_m, k_fin, n_tubes, n_rows, Tube_pass, w (bundle width),
       pitch_ratio (annular) or pitch_V, pitch_H (square), n_parallel (optional)
cond : m_dot (WHOLE finned side), rho, mu, k, Pr   (bulk properties of the gas)
"""
import numpy as np

_NEEDED = ('Tube_OD', 'Tube_t', 'Tube_L', 'Fin_OD', 'Fin_t', 'Fin_per_m', 'k_fin', 'n_tubes', 'n_rows',
           'Tube_pass', 'w')


def _need(geom, cond, extra=()):
    miss = [k for k in _NEEDED + tuple(extra) if k not in geom]
    if miss:
        raise KeyError(f"TFHE HTC: geom is missing {miss}.")
    miss = [k for k in ('m_dot', 'rho', 'mu', 'k', 'Pr') if k not in cond]
    if miss:
        raise KeyError(f"TFHE HTC: cond is missing {miss} (single-phase bulk properties of the gas).")


def _fin_efficiency_straight_rect(t_fin, L_fin, h_c, k_fin, D_fin, D_tube):
    """Efficiency eta_f of an annular fin treated as a straight rectangular fin of corrected length
    (Pis'mennyi), and the correction factor psi_eta. mL = sqrt(2 h/(k t)) L'."""
    if h_c < 0:
        raise ValueError(f"TFHE: negative heat transfer coefficient {h_c}.")
    m = np.sqrt(2 * h_c / (k_fin * t_fin))
    L_p = L_fin * (1 + (0.191 + 0.054 * (D_fin / D_tube)) * np.log(D_fin / D_tube))
    mL = m * L_p
    eta_f = 1.0 if mL == 0 else np.tanh(mL) / mL
    psi_eta = 1 - 0.016 * (D_fin / D_tube - 1) * (1 + np.tanh(2 * mL - 1))
    return eta_f, psi_eta


def _bank_h(geom, cond, area_fractions, Psi_f, sigma_1, sigma_2, pitch_for_diag):
    """Common part: gas velocity in the narrowest section, h_c, fin efficiency, reduced coefficient."""
    D, t = geom['Tube_OD'], geom['Fin_t']
    L, N_rows = geom['Tube_L'], geom['n_rows'] * geom['Tube_pass']
    N_fins = L * geom['Fin_per_m'] - 1
    s = (L - N_fins * t) / (N_fins + 1)                                   # fin spacing
    Fin_L = (geom['Fin_OD'] - D) / 2
    n_tpr = geom['n_tubes'] / N_rows                                      # tubes per row

    D_fin_c = D + 2 * Fin_L * t / s                                       # conventional fin diameter
    diag = np.sqrt(2) * pitch_for_diag                                    # staggered bank
    psi_c = (pitch_for_diag - D_fin_c) / (diag - D_fin_c)
    S_free = geom['Tube_L'] * geom['w'] - n_tpr * D_fin_c * L
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

    eta_f, psi_eta = _fin_efficiency_straight_rect(t, Fin_L, h_c, geom['k_fin'], geom['Fin_OD'], D)
    Af_A, At_A = area_fractions
    return (Af_A * eta_f * psi_eta + At_A) * h_c


def TFHE_pismennyi_annular_HTC(AS, geom, cond):
    """Annular (round) fins, staggered bank. sigma_1 = sigma_2 = geom['pitch_ratio'], pitch = pitch_ratio * Tube_OD."""
    _need(geom, cond, ('pitch_ratio',))
    D, t, L = geom['Tube_OD'], geom['Fin_t'], geom['Tube_L']
    N_fins = L * geom['Fin_per_m'] - 1
    s = (L - N_fins * t) / (N_fins + 1)
    # areas of ONE tube (the fractions do not depend on the number of tubes)
    A_tip = 2 * np.pi * (geom['Fin_OD'] / 2) * t * N_fins
    A_faces = 2 * np.pi * ((geom['Fin_OD'] / 2) ** 2 - (D / 2) ** 2) * N_fins
    A_tube = 2 * np.pi * (D / 2) * (L - N_fins * t)
    A_tot = A_tip + A_faces + A_tube
    Psi_f = (geom['Fin_OD'] ** 2 - D ** 2 + 2 * geom['Fin_OD'] * t) / (2 * D * s) + 1 - t / s
    sigma = geom['pitch_ratio']
    return _bank_h(geom, cond, ((A_tip + A_faces) / A_tot, A_tube / A_tot), Psi_f, sigma, sigma,
                   geom.get('pitch', sigma * D))


def TFHE_pismennyi_square_HTC(AS, geom, cond):
    """Square (plate) fins, staggered bank. sigma_1 = pitch_V / Tube_OD, sigma_2 = pitch_H / Tube_OD."""
    _need(geom, cond, ('pitch_V', 'pitch_H'))
    D, t, L = geom['Tube_OD'], geom['Fin_t'], geom['Tube_L']
    N_fins = L * geom['Fin_per_m'] - 1
    s = (L - N_fins * t) / (N_fins + 1)
    A_r = 2 * (geom['Fin_OD'] ** 2 - 0.785 * D ** 2 + 2 * geom['Fin_OD'] * t) * (L / s)   # fins of ONE tube
    A_t = np.pi * D * (L - N_fins * t)                                                    # bare tube of ONE tube
    A_tot = A_r + A_t
    Psi_f = 2 * (geom['Fin_OD'] ** 2 - 0.785 * D ** 2 + 2 * geom['Fin_OD'] * t) / (np.pi * D * s) + 1 - t / s
    return _bank_h(geom, cond, (A_r / A_tot, A_t / A_tot), Psi_f, geom['pitch_V'] / D, geom['pitch_H'] / D,
                   np.sqrt(geom['pitch_V'] * geom['pitch_H']))


def TFHE_pismennyi_HTC(AS, geom, cond):
    """Pis'mennyi correlation routed on geom['Fin_type'] : 'Annular' (round fins) or 'Square' (plate fins),
    case insensitive."""
    ft = str(geom.get('Fin_type', '')).strip().lower()
    if ft == 'annular':
        return TFHE_pismennyi_annular_HTC(AS, geom, cond)
    if ft == 'square':
        return TFHE_pismennyi_square_HTC(AS, geom, cond)
    raise ValueError(f"TFHE HTC: geom['Fin_type'] = {geom.get('Fin_type')!r}: expected 'Annular' or 'Square'.")