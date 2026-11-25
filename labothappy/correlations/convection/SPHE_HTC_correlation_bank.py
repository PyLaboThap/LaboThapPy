# -*- coding: utf-8 -*-
"""
Heat transfer coefficients of spiral plate heat exchangers (SPHE), (AS, geom, cond) formalism.
Every function returns h [W/m2/K] for one channel (single phase).

geom : Dh               hydraulic diameter of the channel [m]   (Dh = 2 b H / (b + H) for a rectangular channel)
       R_curv           mean radius of curvature of the channel [m]
       b                channel gap = spacing between the plates [m] (set by the SPHE geometry)
       H                plate height = channel height, direction of the spiral axis [m]
       L                developed length of the channel [m]
cond : G [kg/m2/s] mass flux of ONE channel, mu, k, Pr (+ mu_w, default mu)

In Picon-Nunez / Rajavel, 'b' is the channel gap: curvature term (b / 2R)^(1/6) and ratio H/b.
"""
import numpy as np
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _gnielinski_core
from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import (
    _curved_pipe_1P, _R_curv, _L, _D, _need_cond)


def _AR(geom):
    """Aspect ratio of the rectangular channel, short side / long side [-] (<= 1).
    gap = geom['b'], height = geom['H'] (plate height)."""
    gap = geom.get('b')
    height = geom.get('H')
    if gap is None or height is None:
        raise KeyError("SPHE rectangular duct: geom needs the gap 'b' "
                       "and the plate height 'H'.")
    return min(gap, height) / max(gap, height)


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
    Re <  Re_c : h = 1.86 cp G Re^(-2/3) Pr^(-2/3) (b / 2R)^(1/6)                 (laminar)
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
        b = geom.get('b')
        if b is None:
            raise KeyError("SPHE Picon-Nunez laminar branch: geom needs 'b'.")
        return 1.86 * cp * G * Re ** (-2 / 3) * Pr ** (-2 / 3) * (b / (2 * R)) ** (1 / 6)
    return (1 + 3.54 * curv) * 0.023 * cp * G * Re ** (-0.2) * Pr ** (-2 / 3)


def SPHE_rajavel_HTC(AS, geom, cond):
    """Rajavel (2014), experimental/numerical spiral plate HX: Nu = 0.02846 Re^0.828 Pr^0.336 (H/b)^-0.0871.
    WATER ONLY, channels 4-6 mm, 50.6 < H/b < 76 (warning outside). The original returned nothing (missing return)."""
    import warnings
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    D, G, mu, k, Pr = _D(geom), cond['G'], cond['mu'], cond['k'], cond['Pr']
    H, b = geom.get('H'), geom.get('b')
    if H is None or b is None:
        raise KeyError("SPHE Rajavel: geom needs 'H' and 'b'.")
    if AS is not None and 'water' not in str(AS.name()).lower():
        warnings.warn("SPHE Rajavel: correlation established for water only.")
    Re = G * D / mu
    Nu = 0.02846 * Re ** 0.828 * Pr ** 0.336 * (H / b) ** (-0.0871)
    return Nu * k / D


def SPHE_morimoto_hotta_HTC(AS, geom, cond):
    """
    Morimoto and Hotta (1986) single-phase heat transfer coefficient for a spiral plate heat
    exchanger.

    Parameters:
    AS : object
        Fluid state object (unused in correlation, but part of the signature).
    geom : dict
        Geometry dictionary containing:
        - 'Dh' (or 'D_t_i'): hydraulic diameter of the channel, 2 H b / (H + b) [m]
        - 'R_curv': mean radius of curvature of the spiral channel [m]
    cond : dict
        Condition dictionary containing flow properties:
        - 'G': mass flux of ONE channel [kg/m2/s]
        - 'mu': dynamic viscosity [Pa*s]
        - 'k': thermal conductivity [W/m/K]
        - 'Pr': Prandtl number [-]

    Returns:
    float
        Heat transfer coefficient h [W/m2/K].

    Reference:
    Morimoto, E., Hotta, K. (1986) Study of geometric structure and heat transfer characteristics of
    spiral plate heat exchanger. Transactions of the Japan Society of Mechanical Engineers Series B
    52(474), 926-933 (in Japanese). doi:10.1299/kikaib.52.926

    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.conv_internal.Morimoto_Hotta, which cites the reference above (also Bidabadi, Sadaghiani,
    Azad 2013, Scientia Iranica 20(5), 1445-1454; Turgut & Coban 2016, Heat and Mass Transfer,
    doi:10.1007/s00231-016-1861-y).

    Validity range:
    Not given by the reference library (paper in Japanese). A warning is raised if Re < 2000: this
    lower limit is a labothappy safeguard (turbulent-type correlation), not from the paper.

    Definitions:
    Re and Nu based on hydraulic diameter (Dh) with bulk properties. R_m is the mean spiral radius.
    """
    import warnings

    # 1. Check required conditions and extract properties
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Extract hydraulic diameter and mean radius of curvature
    Dh = _D(geom)
    R_m = _R_curv(geom)

    # Calculate Reynolds number based on hydraulic diameter and bulk properties
    Re = G * Dh / mu

    # Emit warning if Re is below the turbulent safeguard limit
    if Re < 2000:
        warnings.warn(
            "SPHE Morimoto & Hotta: Re < 2000. This lower limit is a labothappy "
            "safeguard (turbulent-type correlation), not from the paper."
        )

    # Calculate Nusselt number using the Morimoto & Hotta (1986) correlation
    # Nu = 0.0239 * (1 + 5.54 * Dh / R_m) * Re**0.806 * Pr**0.268
    curvature_term = 1.0 + 5.54 * Dh / R_m
    reynolds_term = Re ** 0.806
    prandtl_term = Pr ** 0.268
    Nu = 0.0239 * curvature_term * reynolds_term * prandtl_term

    # Calculate heat transfer coefficient from Nusselt number
    h = Nu * k / Dh

    return h
