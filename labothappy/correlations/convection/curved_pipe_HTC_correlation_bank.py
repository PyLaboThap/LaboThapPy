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
import warnings
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


def curved_pipe_schmidt_SP_HTC(AS, geom, cond):
    """Single-phase turbulent heat transfer inside a helically coiled tube, Schmidt (1967).

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Fluid state object (unused, kept for API consistency).
    geom : dict
        Geometry dictionary with keys:
        - 'Dh' (or 'D', 'Tube_ID') : inner diameter [m]
        - 'R_curv' (or 'D_coil' with optional 'coil_pitch') : radius of curvature [m]
    cond : dict
        Flow conditions dictionary with keys:
        - 'G' : mass flux in ONE tube [kg/m2/s]
        - 'mu' : dynamic viscosity [Pa s]
        - 'k' : thermal conductivity [W/m/K]
        - 'Pr' : Prandtl number [-] (bulk properties)

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K] based on inner diameter.

    Reference
    ---------
    Schmidt, E.F. (1967) Waermeuebergang und Druckverlust in Rohrschlangen.
    Chemie Ingenieur Technik 39(13), 781-789. doi:10.1002/cite.330391302

    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.conv_internal.helical_turbulent_Nu_Schmidt, which cites the
    reference above (also Hardik, Baburajan, Prabhu 2015, IJHMT 89, 522-538,
    and Rohsenow, Handbook of Heat Transfer 3E). See also the review by
    El-Genk & Schriener (2016), Heat Transfer Engineering, doi:10.1080/01457632.2016.1194693.

    Validity
    --------
    Re_crit < Re < 1.5e5 (turbulent form); curv < 0.14 (range of Schmidt's Re_crit).

    Definitions
    -----------
    Re = G * D / mu (inner diameter D). Properties at bulk temperature.
    Single-phase only (no quality, no heat flux).
    """
    # --- read inputs ---------------------------------------------------------
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    G = cond['G']          # [kg/m2/s]
    mu = cond['mu']        # [Pa s]
    k = cond['k']          # [W/m/K]
    Pr = cond['Pr']        # [-]

    D = _D(geom)           # [m] inner diameter
    R = _R_curv(geom)      # [m] radius of curvature of the coil

    # --- geometry and flow parameters ----------------------------------------
    curv = D / (2.0 * R)   # [-] curvature ratio (Di / Dc)
    Re = G * D / mu        # [-] Reynolds number on inner diameter

    # --- Schmidt critical Reynolds number (transition to turbulence) ---------
    # Re_crit = 2300 * (1 + 8.6 * curv**0.45)  (from fluids library)
    Re_crit = 2300.0 * (1.0 + 8.6 * curv ** 0.45)

    # --- validity warnings ---------------------------------------------------
    if Re <= Re_crit:
        warnings.warn(
            f"curved_pipe_schmidt_SP_HTC: Re={Re:.0f} <= Re_crit={Re_crit:.0f} "
            "(below turbulent regime of the correlation)."
        )
    if Re >= 1.5e5:
        warnings.warn(
            f"curved_pipe_schmidt_SP_HTC: Re={Re:.0f} >= 1.5e5 "
            "(above stated validity range of the turbulent correlation)."
        )
    if curv >= 0.14:
        warnings.warn(
            f"curved_pipe_schmidt_SP_HTC: curv={curv:.3f} >= 0.14 "
            "(outside Schmidt's curvature range for Re_crit)."
        )

    # --- Nusselt number correlation ------------------------------------------
    if Re <= 2.2e4:
        # Low-Re turbulent branch
        coefficient = 1.0 + 14.8 * (1.0 + curv) * curv ** (1.0 / 3.0)
        exponent = 0.8 - 0.22 * curv ** 0.1
        Nu = 0.023 * coefficient * Re ** exponent * Pr ** (1.0 / 3.0)
    else:
        # High-Re turbulent branch
        coefficient = 1.0 + 3.6 * (1.0 - curv) * curv ** 0.8
        Nu = 0.023 * coefficient * Re ** 0.8 * Pr ** (1.0 / 3.0)

    # --- heat transfer coefficient -------------------------------------------
    h = Nu * k / D         # [W/m2/K]

    return h


def curved_pipe_mori_nakayama_SP_HTC(AS, geom, cond):
    """
    Calculate the single-phase turbulent heat transfer coefficient inside a curved pipe
    (helical coil) using the correlation by Mori and Nakayama (1967).

    Parameters
    ----------
    AS : object
        CoolProp abstract state or similar fluid property handler (not explicitly used
        for single-phase properties here, but kept for interface consistency).
    geom : dict
        Geometry dictionary containing keys for diameter ('Dh', 'D', or 'Tube_ID') [m],
        and curvature radius ('R_curv' or 'D_coil', optionally 'coil_pitch') [m].
    cond : dict
        Conditions dictionary containing:
        - G : mass flux in one tube [kg/m^2/s]
        - mu : dynamic viscosity [Pa s]
        - k : thermal conductivity [W/m/K]
        - Pr : Prandtl number [-]

    Returns
    -------
    h : float
        Single-phase heat transfer coefficient on the inner diameter [W/m^2/K].

    Reference
    ---------
    Mori, Y., Nakayama, W. (1967) Study on forced convective heat transfer in curved pipes.
    International Journal of Heat and Mass Transfer 10(5), 681-695.
    doi:10.1016/0017-9310(67)90113-5

    Review reference:
    El-Genk, M.S., Schriener, T.M. (2016) A review and correlations for convection heat transfer
    and pressure losses in toroidal and helically coiled tubes. Heat Transfer Engineering
    (online 2016), doi:10.1080/01457632.2016.1194693 (as cited by ht).

    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.conv_internal.helical_turbulent_Nu_Mori_Nakayama (as coded), which cites the reference above.

    Validity range
    --------------
    - Turbulent flow regime.
    - Reference library applicability limit: Re * curv**2 > 0.1. A warning is emitted
      if this condition is violated.

    Definitions
    -----------
    - Diameter used in Reynolds number: inner diameter D [m].
    - Properties used: bulk fluid properties (mu, k, Pr).
    """
    # 1. Read required conditions and geometry parameters
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    D = _D(geom)
    R = _R_curv(geom)

    # 2. Compute curvature ratio curv = D / (2 * R)
    curv = D / (2.0 * R)

    # 3. Compute Reynolds number Re = G * D / mu
    Re = G * D / mu

    # Check applicability limit given by the reference library
    applicability_parameter = Re * (curv ** 2.0)
    if applicability_parameter <= 0.1:
        warnings.warn(
            f"Re * curv**2 = {applicability_parameter} is <= 0.1. "
            "This is outside the validity range specified for the Mori & Nakayama correlation.",
            UserWarning
        )

    # 4. Calculate Nusselt number branches based on Prandtl number
    if Pr < 1.0:
        # Equation for Prandtl numbers less than 1
        term1 = Pr / (26.2 * (Pr ** (2.0 / 3.0) - 0.074)) * (Re ** 0.8) * (curv ** 0.1)
        term2 = 1.0 + 0.098 * ((Re * (curv ** 2.0)) ** -0.2)
    else:
        # Equation for Prandtl numbers greater than or equal to 1
        term1 = (Pr ** 0.4) / 41.0 * (Re ** (5.0 / 6.0)) * (curv ** (1.0 / 12.0))
        term2 = 1.0 + 0.061 / ((Re * (curv ** 2.5)) ** (1.0 / 6.0))

    # 5. Calculate final Nusselt number and heat transfer coefficient
    Nu = term1 * term2
    h = Nu * k / D

    return h


def curved_pipe_xin_ebadian_SP_HTC(AS, geom, cond):
    """
    Calculate the single-phase turbulent convective heat transfer coefficient inside
    a helical pipe using the correlation by Xin and Ebadian (1997).

    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.conv_internal.helical_turbulent_Nu_Xin_Ebadian, which cites the reference below.
    This correlation is also reviewed by El-Genk & Schriener (2016).

    Parameters
    ----------
    AS : object
        CoolProp abstract state object (or similar fluid property interface).
    geom : dict
        Geometry dictionary containing:
        - 'Dh' (or 'D', 'Tube_ID') : inner diameter of the pipe [m]
        - 'R_curv' (or 'D_coil' and optionally 'coil_pitch') : radius of curvature [m]
    cond : dict
        Conditions dictionary containing:
        - 'G' : mass flux in ONE tube [kg/m2/s]
        - 'mu' : dynamic viscosity [Pa s]
        - 'k' : thermal conductivity [W/m/K]
        - 'Pr' : Prandtl number [-]

    Returns
    -------
    float
        Convective heat transfer coefficient, h [W/m2/K], evaluated on the inner diameter.

    Reference
    ---------
    Xin, R.C., Ebadian, M.A. (1997) The effects of Prandtl numbers on local and average
    convective heat transfer characteristics in helical pipes. Journal of Heat Transfer 119(3),
    467-473.
    doi:10.1115/1.2824120

    Review Reference:
    El-Genk, M.S., Schriener, T.M. (2016) A review and correlations for convection heat transfer
    and pressure losses in toroidal and helically coiled tubes. Heat Transfer Engineering
    (online 2016), doi:10.1080/01457632.2016.1194693 (as cited by ht).

    Validity Range
    --------------
    - Reynolds number (Re) : 5e3 < Re < 1e5 (Xin & Ebadian 1997, Eq. 6)
    - Prandtl number (Pr) : 0.7 < Pr < 5
    - Curvature ratio (curv) : 0.0267 < curv < 0.0884
    """
    # Step 1: Ensure required conditions are present and extract properties
    _need_cond(cond, 'G', 'mu', 'k', 'Pr')
    G = cond['G']
    mu = cond['mu']
    k = cond['k']
    Pr = cond['Pr']

    # Extract inner diameter and curvature radius from geometry
    D = _D(geom)
    R = _R_curv(geom)

    # Step 2: Calculate the curvature ratio (curv = D / (2 * R))
    curv = D / (2.0 * R)

    # Step 3: Calculate the Reynolds number (Re = G * D / mu)
    Re = G * D / mu

    # Check validity ranges and warn if outside
    if not (5e3 < Re < 1e5):
        warnings.warn(
            f"Xin & Ebadian SP HTC: Re = {Re:.2e} is outside the validity range (5e3 < Re < 1e5, "
            "Eq. 6)."
        )
    if not (0.7 < Pr < 5.0):
        warnings.warn(
            f"Xin & Ebadian SP HTC: Pr = {Pr:.2f} is outside the validity range (0.7 < Pr < 5)."
        )
    if not (0.0267 < curv < 0.0884):
        warnings.warn(
            f"Xin & Ebadian SP HTC: curv = {curv:.4f} is outside the validity range "
            "(0.0267 < curv < 0.0884)."
        )

    # Step 4: Calculate Nusselt number using Xin & Ebadian (1997) equation
    # Nu = 0.00619 * Re^0.92 * Pr^0.4 * (1 + 3.455 * curv)
    term_re = Re ** 0.92
    term_pr = Pr ** 0.4
    term_curv = 1.0 + 3.455 * curv
    Nu = 0.00619 * term_re * term_pr * term_curv

    # Calculate the heat transfer coefficient h = Nu * k / D
    h = Nu * k / D

    return h

#%% CONDENSING


def curved_pipe_zakeralhoseini_condensing_HTC(AS, geom, cond):
    """Condensation of R1234yf in helically coiled tubes (Zakeralhoseini et al. 2020), Eq. (25):
    h_tp = h_l * 0.27 (1 + 2.98 / X_tt^0.8) p_r^-0.75, with h_l = Nu_l k_l / d and Nu_l from
    Mori & Nakayama, Eq. (22): Nu_l = (1/41) Re_l^(5/6) Pr_l^0.4 (d/D)^(1/12)
    {1 + 0.061 / [Re_l (d/D)^2.5]^(1/6)}, Re_l = G (1 - x) d / mu_l (Eq. 23), D = coil diameter
    (the printed "p_r^0.4" of Eq. 22 is read as Pr_l^0.4, as in Mori-Nakayama).
    cond : G, p, x      geom : Dh, R_curv (or D_coil), L
    Reference: Zakeralhoseini, S., Sajadi, B., et al. (2020) Experimental investigation of the heat
    transfer coefficient and pressure drop of R1234yf during flow condensation in helically coiled
    tubes, Int. J. Thermal Sciences 157, 106516, Eqs. (21)-(25).
    Validity: R1234yf, d 8.2 mm, coil 95.3 mm, pitch 16.7 mm, Tsat 20 / 30 C, G 95-285 kg/m2/s,
    x 0.1-0.8; 95 % of data within +-15 %. Eq. (22): Pr_l > 1 and Re_l (d/D)^2.5 > 0.4.
    (2026-11-21: h_l was the 'Sigalotti' single-phase curved h; corrected to Eq. (22).)"""
    _need_cond(cond, 'G', 'p', 'x')
    D, L, R, G, p = _D(geom), _L(geom), _R_curv(geom), cond['G'], cond['p']
    x = _clamp_x(cond['x'], 1e-2)
    p_crit = AS.p_critical()
    s = _sat_props(AS, p, cond)
    # Validity warnings (paper range)
    if not 95 <= G <= 285:
        warnings.warn(f"Zakeralhoseini 2020: G = {G:.0f} kg/m2/s outside 95-285.")
    if not 0.1 <= x <= 0.8:
        warnings.warn(f"Zakeralhoseini 2020: x = {x:.2f} outside 0.1-0.8.")
    # Liquid Reynolds number, Eq. (23), and tube / coil diameter ratio
    Re_l = G * (1 - x) * D / s['mu_l']
    d_over_Dc = D / (2.0 * R)
    Pr_l = s['Pr_l']
    if Pr_l <= 1 or Re_l * d_over_Dc ** 2.5 <= 0.4:
        warnings.warn("Zakeralhoseini 2020: Mori-Nakayama Eq. (22) used outside Pr_l > 1, "
                      "Re_l (d/D)^2.5 > 0.4.")
    # Liquid Nusselt number, Mori & Nakayama, Eq. (22)
    term_main = (1.0 / 41.0) * Re_l ** (5.0 / 6.0) * Pr_l ** 0.4 * d_over_Dc ** (1.0 / 12.0)
    term_corr = 1.0 + 0.061 / (Re_l * d_over_Dc ** 2.5) ** (1.0 / 6.0)
    Nu_l = term_main * term_corr
    h_l = Nu_l * s['k_l'] / D
    # Martinelli parameter and Eq. (25)
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
    h_tp = (1.6 X_tt^-0.74 + 183000 Bo^1.46) h_lo, h_lo = Dittus-Boelter 0.023 Re_lo^0.8 Pr_l^0.4 k_l/D,
    Re_lo = G D / mu_l: TOTAL flow assumed liquid (text after Eq. 10 of the paper).
    cond : G, p, x, q      geom : Dh
    Reference: Zhao, L., Guo, L., Bai, B., Hou, Y., Zhang, X. (2003) Convective boiling heat transfer and
    two-phase flow characteristics inside a small horizontal helically coiled tubing once-through steam
    generator, Int. J. Heat Mass Transfer 46, 4779-4788, Eqs. (10), (21), (22).
    Validity: steam-water, d 9 mm, coil D 292 mm, pitch 30 mm, p 0.5-3.5 MPa, G 236-943 kg/m2/s,
    q 0-900 kW/m2, x up to 0.95; mean deviation 12 %.
    (2026-11-21: h_lo was computed on (1 - x) G; corrected to the total flow, paper definition.)"""
    s, x, inv_Xtt, Bo = _boiling_inputs(AS, geom, cond)
    D = _D(geom)
    # Validity warnings (paper range)
    if not 0.5e6 <= cond['p'] <= 3.5e6:
        warnings.warn(f"Zhao 2003: p = {cond['p'] / 1e6:.2f} MPa outside 0.5-3.5 MPa (steam-water data).")
    if not 236 <= cond['G'] <= 943:
        warnings.warn(f"Zhao 2003: G = {cond['G']:.0f} kg/m2/s outside 236-943.")
    # Zhao et al. (2003): h_lo for the total flow assumed liquid (Re = G d / mu_l, nomenclature)
    Re_lo = cond['G'] * D / s['mu_l']
    h_lo = 0.023 * Re_lo ** 0.8 * s['Pr_l'] ** 0.4 * s['k_l'] / D
    # Eq. (21)
    return (1.6 * inv_Xtt ** 0.74 + 183000 * Bo ** 1.46) * h_lo


def curved_pipe_nariai_boiling_HTC(AS, geom, cond):
    """Boiling in a helical coil at subcritical conditions, small-scale solar ORC (Lazova et al. 2017, Nariai-type):
    Nu = 170 Re_lo^0.8 Pr_l^(1/3) (Bo + 0.00015 X_tt^-0.67), Re_lo = G D / mu_l (all liquid).
    cond : G, p, x, q      geom : Dh"""
    s, x, inv_Xtt, Bo = _boiling_inputs(AS, geom, cond)
    D = _D(geom)
    Re_lo = cond['G'] * D / s['mu_l']
    return 170 * Re_lo ** 0.8 * s['Pr_l'] ** (1 / 3) * (Bo + 0.00015 * inv_Xtt ** 0.67) * s['k_l'] / D