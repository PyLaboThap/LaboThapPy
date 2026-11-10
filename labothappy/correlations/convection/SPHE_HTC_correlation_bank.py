# -*- coding: utf-8 -*-
"""
Heat transfer coefficients of spiral plate heat exchangers (SPHE), (AS, geom, cond) formalism.
Every function returns h [W/m2/K] for one channel (single phase).

geom : Dh               hydraulic diameter of the channel [m]   (Dh = 2 b w / (b + w) for a rectangular channel)
       R_curv           mean radius of curvature of the channel [m]
       ch_w             channel width (plate height, direction of the axis) [m]
       ch_b  (or H)     channel spacing (gap between plates) [m]
       L                developed length of the channel [m]
cond : G [kg/m2/s] mass flux of ONE channel, mu, k, Pr (+ mu_w, default mu)

Notation kept from the original file: in Picon-Nunez / Rajavel, 'ch_w' is the dimension used in the
curvature term (ch_w / 2R)^(1/6) and in H/ch_w. CHECK which dimension your geometry assigns to it.
"""
import numpy as np
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _gnielinski_core
from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import (
    _curved_pipe_1P, _R_curv, _L, _D, _need_cond)


def _AR(geom):
    b = geom.get('ch_b', geom.get('H'))
    w = geom.get('ch_w')
    if b is None or w is None:
        raise KeyError("SPHE rectangular duct: geom needs 'ch_w' and 'ch_b' (or 'H').")
    return min(b, w) / max(b, w)


def _Nu_lam_rect(AR):
    """Shah & London (1978), fully developed laminar, constant wall temperature."""
    return 7.541 * (1 - 2.61 * AR + 4.97 * AR ** 2 - 5.119 * AR ** 3 + 2.702 * AR ** 4 - 0.548 * AR ** 5)


def SPHE_rect_duct_HTC(AS, geom, cond):
    """Straight rectangular duct: Shah & London laminar (Re<2300), Gnielinski/Stephan-Preusser-Sieder-Tate
    (Re>10000), linear interpolation on h between. cond: G, mu, k, Pr (+mu_w)."""
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    D, L, G, mu, k, Pr = _D(geom), _L(geom), cond['G'], cond['mu'], cond['k'], cond['Pr']
    mu_w = cond.get('mu_w', mu)
    Re = G * D / mu
    h_lam = _Nu_lam_rect(_AR(geom)) * k / D
    if Re < 2300:
        return h_lam
    h_turb = _gnielinski_core(mu, Pr, mu_w, k, G, D, L)
    if Re > 10000:
        return h_turb
    return h_lam + (h_turb - h_lam) / (10000 - 2300) * (Re - 2300)


def SPHE_curved_rect_duct_HTC(AS, geom, cond):
    """Curved rectangular duct: Facao & Oliveira laminar (Re<2300), curved-pipe correlations otherwise.
    Bug fixed vs original: curvature ratio is Dh/(2R), not Dh/2*R."""
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    D, L, R, G, mu, k, Pr = _D(geom), _L(geom), _R_curv(geom), cond['G'], cond['mu'], cond['k'], cond['Pr']
    mu_w = cond.get('mu_w', mu)
    Re = G * D / mu
    if Re >= 2300:
        return _curved_pipe_1P(mu, Pr, mu_w, k, G, D, L, R)
    De = Re * np.sqrt(D / (2 * R))
    Nu_1000, Nu_1200 = 0.8508 * 1000 ** 0.5207, 38.5
    if De < 1000:
        Nu = 0.8508 * De ** 0.5207
    elif De > 1200:
        Nu = Nu_1200
    else:
        Nu = Nu_1000 + (Nu_1200 - Nu_1000) / 200 * (De - 1000)
    return Nu * k / D


def SPHE_picon_nunez_HTC(AS, geom, cond):
    """Shortcut design approach for spiral heat exchangers, Picon-Nunez et al. (2007), after Minton:
    Re <  Re_c : h = 1.86 cp G Re^(-2/3) Pr^(-2/3) (ch_w / 2R)^(1/6)                 (laminar)
    Re >= Re_c : h = (1 + 3.54 Dh/2R) 0.023 cp G Re^(-0.2) Pr^(-2/3)                 (turbulent)
    Re_c = 20000 (Dh/2R)^0.32.
    Regimes are SWAPPED vs the original file (laminar form was used above Re_cr, turbulent form below).
    Re_c: the 12730 (Dh/2R)^0.2 of Ito is a straight-pipe-like critical Re and is not what the source uses:
    the Minton value is taken (check against the paper)."""
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    D, R, G, mu, k, Pr = _D(geom), _R_curv(geom), cond['G'], cond['mu'], cond['k'], cond['Pr']
    curv = D / (2 * R)
    Re = G * D / mu
    cp = Pr * k / mu
    if Re < 20000 * curv ** 0.32:
        ch_w = geom.get('ch_w')
        if ch_w is None:
            raise KeyError("SPHE Picon-Nunez laminar branch: geom needs 'ch_w'.")
        return 1.86 * cp * G * Re ** (-2 / 3) * Pr ** (-2 / 3) * (ch_w / (2 * R)) ** (1 / 6)
    return (1 + 3.54 * curv) * 0.023 * cp * G * Re ** (-0.2) * Pr ** (-2 / 3)


def SPHE_rajavel_HTC(AS, geom, cond):
    """Rajavel (2014), experimental/numerical spiral plate HX: Nu = 0.02846 Re^0.828 Pr^0.336 (H/ch_w)^-0.0871.
    WATER ONLY, channels 4-6 mm, 50.6 < H/b < 76 (warning outside). The original returned nothing (missing return)."""
    import warnings
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    D, G, mu, k, Pr = _D(geom), cond['G'], cond['mu'], cond['k'], cond['Pr']
    H, ch_w = geom.get('H'), geom.get('ch_w')
    if H is None or ch_w is None:
        raise KeyError("SPHE Rajavel: geom needs 'H' and 'ch_w'.")
    if AS is not None and 'water' not in str(AS.name()).lower():
        warnings.warn("SPHE Rajavel: correlation established for water only.")
    Re = G * D / mu
    Nu = 0.02846 * Re ** 0.828 * Pr ** 0.336 * (H / ch_w) ** (-0.0871)
    return Nu * k / D
