# -*- coding: utf-8 -*-
"""
Heat transfer coefficients of curved (helically coiled) pipes, (AS, geom, cond) formalism.
Used for the coil side of helical coil heat exchangers (HCHE); the single-phase curved-pipe function is also
the basis of the spiral plate bank (SPHE). Prefix HCHE_..._HTC; every function returns h [W/m2/K].

geom : Dh (or D, Tube_ID)   hydraulic / inner diameter [m]
       R_curv               radius of curvature of the flow path [m]
                            (or D_coil = coil diameter [m] and, optionally, coil_pitch [m]: the curvature
                            radius of a helix is R (1 + (pitch / (2 pi R))^2), R = D_coil / 2)
       L                    developed length of the flow path [m]  (or Tube_L)
cond : G [kg/m2/s]  mass flux in ONE tube
       single phase  : mu, k, Pr (+ mu_w, default mu)
       two phase     : p [Pa], x [-], q [W/m2] (boiling), optional mu_l, mu_v, k_l, k_v overrides
"""
import numpy as np
import CoolProp.CoolProp as CP

from labothappy.correlations.convection.pipe_HTC_correlation_bank import _gnielinski_core, _sat_props, _clamp_x

#%% HELPERS


def _D(geom):
    for key in ('Dh', 'D', 'Tube_ID'):
        if geom.get(key):
            return geom[key]
    if 'Tube_OD' in geom and 'Tube_t' in geom:
        return geom['Tube_OD'] - 2.0 * geom['Tube_t']
    raise KeyError("Curved pipe HTC: geom needs 'Dh' (or 'D', 'Tube_ID').")


def _R_curv(geom):
    if 'R_curv' in geom:
        return geom['R_curv']
    if 'D_coil' in geom:
        R = geom['D_coil'] / 2.0
        if 'coil_pitch' in geom:
            R *= 1.0 + (geom['coil_pitch'] / (2.0 * np.pi * R)) ** 2      # helix curvature radius
        return R
    raise KeyError("Curved pipe HTC: geom needs 'R_curv' (or 'D_coil' and optionally 'coil_pitch').")


def _L(geom):
    for key in ('L', 'Tube_L'):
        if geom.get(key):
            return geom[key]
    raise KeyError("Curved pipe HTC: geom needs the developed length 'L'.")


def _need_cond(cond, *keys):
    miss = [k for k in keys if k not in cond]
    if miss:
        raise KeyError(f"Curved pipe HTC: cond is missing {miss}.")


def _ito_f(Re, De):
    """Friction factor from the Ito correlation as used by Sigalotti et al. (Xhan correlation):
    f = 16 * 0.1033 sqrt(De) / (Re * [sqrt(1 + 1.729/De) - sqrt(1.729/De)]^3).
    NOT verified here (Fanning or Darcy convention to check against the source)."""
    num = 0.1033 * np.sqrt(De)
    den = (np.sqrt(1.0 + 1.729 / De) - np.sqrt(1.729 / De)) ** 3
    return 16.0 * num / (Re * den)


def _gnielinski_di_piazza(Re, Pr, D, L, f):
    return (((f / 8) * (Re - 1000) * Pr) / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))) * (1 + (D / L) ** (2 / 3))

#%% SINGLE PHASE


def _laminar(mu, Pr, mu_w, k, G, D, L, R):
    curv = D / (2 * R)
    Re = G * D / mu
    De = Re * np.sqrt(curv)
    if curv <= 0.01:
        return _gnielinski_core(mu, Pr, mu_w, k, G, D, L)
    if curv <= 0.1:                                    # Xhan correlation (Sigalotti et al.)
        Nu = (2.153 + 0.318 * De ** 0.643) * Pr ** 0.177
    else:                                              # Zhang et al.
        Nu = 0.8532 * De ** 0.4822 * curv ** (-0.03451)
    return Nu * k / D


def _transition(mu, Pr, mu_w, k, G, D, L, R):
    curv = D / (2 * R)
    Re = G * D / mu
    De = Re * np.sqrt(curv)
    if curv <= 0.01:
        return _gnielinski_core(mu, Pr, mu_w, k, G, D, L)
    f = _ito_f(Re, De)
    if curv <= 0.1:                                    # Xhan correlation (Sigalotti et al.)
        Nu = 0.04 * f * Re ** 1.13 * curv ** 0.077 * Pr ** 0.4
    else:                                              # Di Piazza & Ciofalo
        Nu = _gnielinski_di_piazza(Re, Pr, D, L, f)
    return Nu * k / D


def _turbulent(mu, Pr, mu_w, k, G, D, L, R):
    curv = D / (2 * R)
    Re = G * D / mu
    De = Re * np.sqrt(curv)
    if curv <= 0.01:
        return _gnielinski_core(mu, Pr, mu_w, k, G, D, L)
    if curv <= 0.2:                                    # Schmidt correction, Nellis & Klein
        return _gnielinski_core(mu, Pr, mu_w, k, G, D, L) * (1 + 3.6 * (1 - curv) * curv ** 0.8)
    f = _ito_f(Re, De)                                 # Di Piazza & Ciofalo
    return _gnielinski_di_piazza(Re, Pr, D, L, f) * k / D


def _curved_pipe_1P(mu, Pr, mu_w, k, G, D, L, R):
    """Laminar below Re_cr = 12730 (D/2R)^0.2 (Ito), turbulent above Re = 1e5, Xhan transition correlation between."""
    curv = D / (2 * R)
    Re = G * D / mu
    Re_cr = 12730 * curv ** 0.2
    if Re < Re_cr:
        return _laminar(mu, Pr, mu_w, k, G, D, L, R)
    h = _turbulent(mu, Pr, mu_w, k, G, D, L, R) if Re >= 1e5 else _transition(mu, Pr, mu_w, k, G, D, L, R)
    # Guard: curvature can only enhance heat transfer. The Ito friction factor (laminar-curved formula) used in the
    # transition / Di Piazza branches gives f below the straight-pipe value at high Re -> h below straight tube.
    return max(h, _gnielinski_core(mu, Pr, mu_w, k, G, D, L))


def curved_pipe_sigalotti_SP_HTC(AS, geom, cond):
    """Single phase in a curved (helically coiled) pipe: laminar / transition / turbulent regimes selected on Re and
    on the curvature ratio D/(2 R_curv); straight-pipe Gnielinski below 0.01.
    Reference: Sigalotti, Alvarado-Rodriguez, Rendon, Fluid flow in helically coiled pipes (2023), with
    Ito, Zhang et al., Nellis & Klein (Schmidt correction) and Di Piazza & Ciofalo.
    cond : G, mu, k, Pr (+ mu_w)       geom : Dh, R_curv (or D_coil), L"""
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    mu = cond['mu']
    return _curved_pipe_1P(mu, cond['Pr'], cond.get('mu_w', mu), cond['k'], cond['G'], _D(geom), _L(geom), _R_curv(geom))

#%% CONDENSING


def curved_pipe_zakeralhoseini_condensing_HTC(AS, geom, cond):
    """Condensation of R1234yf in helically coiled tubes (Zakeralhoseini et al.):
    h_tp = h_l * 0.27 (1 + 2.98 / X_tt^0.8) P_r^-0.75, h_l = curved-pipe single-phase h of the liquid flowing alone.
    cond : G, p, x      geom : Dh, R_curv (or D_coil), L"""
    _need_cond(cond, 'G', 'p', 'x')
    D, L, R, G, p = _D(geom), _L(geom), _R_curv(geom), cond['G'], cond['p']
    x = _clamp_x(cond['x'], 1e-4)
    p_crit = AS.p_critical()
    s = _sat_props(AS, p, cond)
    G_l = (1 - x) * G
    h_l = _curved_pipe_1P(s['mu_l'], s['Pr_l'], s['mu_l'], s['k_l'], G_l, D, L, R)
    inv_Xtt = (x / (1 - x)) ** 0.9 * (s['rho_l'] / s['rho_v']) ** 0.5 * (s['mu_v'] / s['mu_l']) ** 0.1
    return h_l * 0.27 * (1 + 2.98 / (1.0 / inv_Xtt) ** 0.8) * (p / p_crit) ** (-0.75)

#%% BOILING


def _boiling_inputs(AS, geom, cond):
    _need_cond(cond, 'G', 'p', 'x', 'q')
    s = _sat_props(AS, cond['p'], cond)
    x = _clamp_x(cond['x'], 1e-4)
    inv_Xtt = (x / (1 - x)) ** 0.9 * (s['rho_l'] / s['rho_v']) ** 0.5 * (s['mu_v'] / s['mu_l']) ** 0.1
    Bo = max(cond['q'], 1e-9) / (cond['G'] * s['h_lv'])
    return s, x, inv_Xtt, Bo


def curved_pipe_zhao_boiling_HTC(AS, geom, cond):
    """Convective boiling in small horizontal helically coiled tubes, once-through steam generator (Zhao et al. 2003):
    h_tp = (1.6 X_tt^-0.74 + 183000 Bo^1.46) h_liq, h_liq = Dittus-Boelter-type 0.023 Re_l^0.8 Pr_l^0.4 k_l/D, G_l = (1 - x) G.
    cond : G, p, x, q      geom : Dh"""
    s, x, inv_Xtt, Bo = _boiling_inputs(AS, geom, cond)
    D = _D(geom)
    Re_l = (1 - x) * cond['G'] * D / s['mu_l']
    h_liq = 0.023 * Re_l ** 0.8 * s['Pr_l'] ** 0.4 * s['k_l'] / D
    return (1.6 * inv_Xtt ** 0.74 + 183000 * Bo ** 1.46) * h_liq


def curved_pipe_nariai_boiling_HTC(AS, geom, cond):
    """Boiling in a helical coil at subcritical conditions, small-scale solar ORC (Lazova et al. 2017, Nariai-type):
    Nu = 170 Re_lo^0.8 Pr_l^(1/3) (Bo + 0.00015 X_tt^-0.67), Re_lo = G D / mu_l (all liquid).
    cond : G, p, x, q      geom : Dh"""
    s, x, inv_Xtt, Bo = _boiling_inputs(AS, geom, cond)
    D = _D(geom)
    Re_lo = cond['G'] * D / s['mu_l']
    return 170 * Re_lo ** 0.8 * s['Pr_l'] ** (1 / 3) * (Bo + 0.00015 * inv_Xtt ** 0.67) * s['k_l'] / D