# -*- coding: utf-8 -*-
"""
Finned-side pressure drop of tube & fins heat exchangers (TFHE), (AS, geom, cond) formalism.

Source : Handbook for transversely finned tube heat exchanger design,
         E. Pis'mennyi, G. Polupan, I. Carvajal-Mariscal, F. Sanchez-Silva, I. Pioro.
         Staggered bank of transversely finned tubes (annular or square fins: same formula).

Returns the pressure drop [Pa] of the WHOLE bank (all rows crossed by the gas), friction only.

geom : Tube_OD, Tube_L, Fin_OD, Fin_t, Fin_per_m, n_tubes, n_rows, Tube_pass, w (bundle width),
       pitch_V, pitch_H, n_parallel and n_series (optional, default 1)
cond : m_dot (WHOLE finned side), rho, mu (gas bulk properties)
"""
import numpy as np

_NEEDED = ('Tube_OD', 'Tube_L', 'Fin_OD', 'Fin_t', 'Fin_per_m', 'n_tubes', 'n_rows', 'Tube_pass', 'w',
           'pitch_V', 'pitch_H')


def TFHE_pismennyi_DP(AS, geom, cond):
    miss = [k for k in _NEEDED if k not in geom] + [k for k in ('m_dot', 'rho', 'mu') if k not in cond]
    if miss:
        raise KeyError(f"TFHE DP: missing {miss} in geom / cond.")
    D, t, L = geom['Tube_OD'], geom['Fin_t'], geom['Tube_L']
    S_1, S_2 = geom['pitch_V'], geom['pitch_H']
    S_m = np.sqrt(S_1 * S_2)
    N_rows = geom['n_rows']
    rho, mu = cond['rho'], cond['mu']

    "Gas velocity - flow area"
    Fin_L = (geom['Fin_OD'] - D) / 2
    N_fins = L * geom['Fin_per_m'] - 1
    s = (L - N_fins * t) / (N_fins + 1)                                   # fin spacing
    n_tpr = geom['n_tubes'] / (N_rows * geom['Tube_pass'])                # tubes per row
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
    return geom.get('n_series', 1) * C_op * xi_o * geom['Tube_pass'] * N_rows * rho * u ** 2 / 2


