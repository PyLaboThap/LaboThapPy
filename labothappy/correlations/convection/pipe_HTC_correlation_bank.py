# -*- coding: utf-8 -*-
"""
Tube-side (in-tube) heat transfer coefficients, in the (AS, geom, cond) formalism of the
BPHX and shell-side STHE correlations.  Every correlation returns h [W/(m^2 K)].

AS   : CoolProp AbstractState of the fluid (it is moved by the saturation / wall property calls:
       build the bulk properties of the cell BEFORE calling).
geom : Tube_ID [m]   (or Tube_OD and Tube_t [m] -> Tube_ID = Tube_OD - 2 Tube_t)
       Tube_L [m]    straight tube length (Gnielinski / Steiner-Taborek entrance effect)
       n_tubes, n_passes (only to build G from m_dot when cond has no 'G')
       Ra [m]        wall roughness (VDI boiling, optional, default 3e-6)
cond : G [kg/m2/s]   mass flux inside ONE tube (or m_dot [kg/s] of the whole side + geom n_tubes, n_passes)
       p [Pa], x [-], q [W/m2] (local heat flux), T_w [K], mu, mu_w, k, Pr, cp, rho (single phase)

Naming: pipe_<author>_<phase>_HTC.  Pool boiling / film boiling / boiling curve (outside
tubes, T-based signature) are kept unchanged at the end of the file.
"""

from math import inf
import numpy as np
import warnings
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP

#%% HELPERS

def _tube_ID(geom):
    if 'Tube_ID' in geom:
        return geom['Tube_ID']
    if 'Tube_OD' in geom and 'Tube_t' in geom:
        return geom['Tube_OD'] - 2.0 * geom['Tube_t']
    raise KeyError("Tube-side HTC: geom needs 'Tube_ID' (or 'Tube_OD' and 'Tube_t').")


def _tube_G(geom, cond):
    """Mass flux inside one tube [kg/m2/s]: cond['G'], else m_dot / (tubes per pass * tube section)."""
    if 'G' in cond:
        return cond['G']
    D = _tube_ID(geom)
    n_tubes_pass = geom['n_tubes'] / geom.get('n_passes', 1)
    return cond['m_dot'] / (n_tubes_pass * np.pi * D**2 / 4.0)


def _sigma(AS, p):
    """Surface tension at saturation pressure p [N/m] (AS first, PropsSI as a fallback)."""
    AS.update(CP.PQ_INPUTS, p, 0.5)
    try:
        return AS.surface_tension()
    except ValueError:
        # Tabular backends (BICUBIC&HEOS, used by hex_MB) can refuse surface_tension():
        # same fallback as BPHE_HTC_correlation_bank._compute_sat_props
        return PropsSI('surface_tension', 'P', p, 'Q', 0, AS.fluid_names()[0])

def _sat_props(AS, p, cond=None):
    """Saturated liquid / vapor properties at p (moves AS).
    cond['mu_l'], cond['mu_v'], cond['k_l'], cond['k_v'] override the CoolProp transport properties
    (some fluids, e.g. R1233zd(E) in some CoolProp versions, have no transport model: give them then).
    Pr = cp mu / k."""
    cond = cond or {}
    def transport(key, fn):
        if key in cond:
            return cond[key]
        try:
            return fn()
        except ValueError:
            raise ValueError(f"CoolProp has no transport model for {AS.fluid_names()[0]}: give cond[{key!r}].")
    AS.update(CP.PQ_INPUTS, p, 0.0)
    liq = dict(rho_l=AS.rhomass(), mu_l=transport('mu_l', AS.viscosity), k_l=transport('k_l', AS.conductivity),
               cp_l=AS.cpmass(), h_l=AS.hmass(), T_bub=AS.T())
    liq['Pr_l'] = liq['cp_l'] * liq['mu_l'] / liq['k_l']
    AS.update(CP.PQ_INPUTS, p, 1.0)
    vap = dict(rho_v=AS.rhomass(), mu_v=transport('mu_v', AS.viscosity), k_v=transport('k_v', AS.conductivity),
               h_v=AS.hmass(), T_dew=AS.T())
    vap['Pr_v'] = AS.cpmass() * vap['mu_v'] / vap['k_v']
    s = {**liq, **vap}
    s['h_lv'] = s['h_v'] - s['h_l']
    s['T_sat'] = 0.5 * (s['T_bub'] + s['T_dew'])
    try:
        s['sigma'] = cond['sigma'] if 'sigma' in cond else _sigma(AS, p)
    except ValueError:
        s['sigma'] = None                       # checked by the correlations that need it (_need_sigma)
    return s


def _need_sigma(AS, s):
    if s['sigma'] is None:
        raise ValueError(f"CoolProp has no surface tension for {AS.fluid_names()[0]}: give cond['sigma'].")


def _clamp_x(x, eps=1e-6):
    """Quality kept strictly inside ]0, 1[ for the correlations built on X_tt / Co."""
    return min(max(x, eps), 1.0 - eps)

#%% SINGLE PHASE

def _gnielinski_core(mu, Pr, mu_w, k, G, Dh, L):
    """Laminar Stephan-Preusser / transition Gnielinski / turbulent Sieder-Tate."""
    def stephan_preusser(Re):
        Nu_1 = 0.0677 * (Re * Pr * (Dh / L))**1.33
        Nu_2 = 1 + 0.1 * Pr * (Re * Dh / L)**0.3
        return 3.657 + Nu_1 / Nu_2

    def gnielinski_transition(Re):
        f = (1.82 * np.log10(Re) - 1.64)**(-2)
        return (((f / 8) * (Re - 1000) * Pr) / (1 + 12.7 * (f / 8)**0.5 * (Pr**(2 / 3) - 1))) * (1 + (Dh / L)**(2 / 3))

    def sieder_tate(Re):
        return 0.027 * Re**0.8 * Pr**(1 / 3) * (mu / mu_w)**0.14

    Re = G * Dh / mu
    if Re > 10000:                      # fully turbulent
        Pr_min, Pr_max = 0.1, 1000
        Nu = sieder_tate(Re)
    elif Re < 2300:                     # fully laminar
        Pr_min, Pr_max = 0.6, inf
        Nu = stephan_preusser(Re)
    else:                               # transition zone
        Pr_min, Pr_max = 0.1, 1000
        Nu = gnielinski_transition(Re)

    if Re >= 1e6 or Re <= 0:
        warnings.warn('Gnielinski single-phase: Reynolds out of validity range !!!')
    if Pr >= Pr_max or Pr <= Pr_min:
        warnings.warn('Gnielinski single-phase: Prandtl out of validity range !!!')
    return Nu * k / Dh


def pipe_gnielinski_SP_HTC(AS, geom, cond):
    """Gnielinski (transition) / Stephan-Preusser (laminar) / Sieder-Tate (turbulent).
    cond : G (or m_dot), mu, Pr, k, mu_w (optional, mu_w = mu if absent)
    geom : Tube_ID, Tube_L
    References: Ammar & Park, Validation of the Gnielinski correlation for enhanced tubes;
                Caputo et al., Heat exchanger design based on economic optimisation."""
    mu = cond['mu']
    return _gnielinski_core(mu, cond['Pr'], cond.get('mu_w', mu), cond['k'],
                            _tube_G(geom, cond), _tube_ID(geom), geom['Tube_L'])


def pipe_dittus_boelter_SP_HTC(AS, geom, cond):
    """Dittus-Boelter, Nu = 0.023 Re^0.8 Pr^n, n = 0.4 (fluid heated) or 0.3 (fluid cooled).
    cond : G (or m_dot), mu, Pr, k, heating (bool, default True: tube fluid is heated)
    geom : Tube_ID
    Reference: Incropera, Fundamentals of heat transfer."""
    D = _tube_ID(geom)
    Re = _tube_G(geom, cond) * D / cond['mu']
    n = 0.4 if cond.get('heating', True) else 0.3
    return 0.023 * Re**0.8 * cond['Pr']**n * cond['k'] / D


def pipe_meshram_HTC(AS, geom, cond):
    """Meshram et al., sCO2 in a printed circuit heat exchanger: Nu = 0.0177 Re^0.871 Pr^0.7.
    (Printed circuit heat exchanger correlation, kept here for convenience.)
    cond : G, mu, Pr, k      geom : Dh (or Tube_ID)
    Reference: Meshram et al., Modeling and analysis of a printed circuit heat exchanger for
    supercritical CO2 power cycle applications."""
    Dh = geom.get('Dh', None) or _tube_ID(geom)
    Re = cond['G'] * Dh / cond['mu']
    return 0.0177 * Re**0.871 * cond['Pr']**0.7 * cond['k'] / Dh

#%% TRANSCRITICAL (sCO2 cooling)

def pipe_liu_sCO2_HTC(AS, geom, cond):
    """Liu et al. (2014), sCO2 cooled in a large tube.
    Validity (Table 2 of the paper): horizontal tubes, CO2 cooling, D 10.7 mm (p 7.5 / 8.0 / 8.5 MPa,
    0.4-0.8 kg/min, i.e. G ~74-148 kg/m2/s, T_in 25-67 C), D 6 mm (8.0 MPa, 0.35 kg/min, T_in 25-58 C)
    and D 4 mm (8.0 MPa, 0.6 kg/min, T_in 32-53 C); Eq. (29) is presented for the 10.7 mm tube
    (most data within +-15 %).
    cond : G, p, T_w, k, rho, mu, cp  (bulk values); T (bulk, optional: cooling check)
    geom : Tube_ID
    Definitions: Nu_w = h d / lambda (Eq. 8: lambda of CO2, wall or bulk not stated; bulk k used here),
    Re_w = G d / mu_w, Pr_w, rho_w, cp_w at the wall temperature, rho_b, cp_b bulk.
    Reference: Liu, Z.-B., He, Y.-L., Yang, Y.-F., Fei, J.-Y. (2014) Experimental study on heat transfer
    and pressure drop of supercritical CO2 cooled in a large tube, Appl. Therm. Eng. 70, 307-315,
    Eq. (29)."""
    D = _tube_ID(geom)
    G, rho, cp = _tube_G(geom, cond), cond['rho'], cond['cp']
    # Validity warnings (test matrix of Table 2)
    if not 7.5e6 <= cond['p'] <= 8.5e6:
        warnings.warn(f"Liu 2014 sCO2: p = {cond['p'] / 1e6:.2f} MPa outside 7.5-8.5 MPa (Table 2).")
    if not 4e-3 <= D <= 10.7e-3:
        warnings.warn(f"Liu 2014 sCO2: D = {D * 1e3:.2f} mm outside 4-10.7 mm (Table 2).")
    if cond.get('T') is not None and cond['T_w'] > cond['T']:
        warnings.warn("Liu 2014 sCO2: correlation for CO2 COOLING (wall colder than the bulk).")
    AS.update(CP.PT_INPUTS, cond['p'], cond['T_w'])           # wall state (moves AS)
    rho_w, cp_w, mu_w, Pr_w = AS.rhomass(), AS.cpmass(), AS.viscosity(), AS.Prandtl()
    Re_w = G * D / mu_w
    Nu = 0.01 * Re_w**0.9 * Pr_w**0.5 * (rho_w / rho)**0.906 * (cp_w / cp)**(-0.585)
    return Nu * cond['k'] / D


def pipe_cheng_sCO2_HTC(AS, geom, cond):
    """Cheng et al. (2024), sCO2 in a horizontal tube, Froude number analysis.
    Validity: D = 10 mm, G 496.7-1346.2 kg/m2/s, q 97.4-400.3 kW/m2, P 7.53-23.51 MPa, CO2.
    cond : G, q, T_w, p, mu, k, h_in, h_out (cell enthalpies; h if h_in/h_out absent)
    geom : Tube_ID
    References: Cheng, Xu, Cao, Zhou, Liu, Supercritical carbon dioxide heat transfer in horizontal
    tube based on the Froude number analysis (2024), Energy 294, 130980, Eqs. (28)-(33), (38)-(40);
    Banuti (2015, 2019) for T_minus / T_plus; Yang & Yang (2012) for the pseudo-critical temperature.
    Definitions: Re_b = G d / mu_b; Pr_b,ave = mu_b cp_ave / k_b with cp_ave = (i_w - i_b)/(T_w - T_b),
    i_b = mean of the cell inlet / outlet enthalpies; x = pseudo-vapor quality (cell average);
    Fr_VL,w = G^2 x^2 / (rho_w^2 g d) (Eq. 33); SBO_w = q / (G i_w) (Eq. 32) uses an ABSOLUTE enthalpy:
    CoolProp default reference state (the paper used NIST data). The paper gives Nu_top and Nu_bot
    (Eqs. 39-40); labothappy returns their mean (perimeter-averaged estimate, labothappy choice).
    Paper ranges: Re_b 6.18e4-5.35e5, Fr_VL,w 7.58e-5-1834 (warnings).
    (2026-11-21: Fr_VL,w used x_in not squared and cp_ave the inlet state: corrected.)"""
    fluid = AS.fluid_names()[0]
    D_in = _tube_ID(geom)
    G, q, T_w, P = _tube_G(geom, cond), cond['q'], cond['T_w'], cond['p']
    mu, k = cond['mu'], cond['k']
    h_in = cond.get('h_in', cond.get('h'))
    h_out = cond.get('h_out', h_in)
    if h_in is None:
        raise KeyError("Cheng sCO2: cond needs 'h_in' and 'h_out' (or 'h').")

    def find_pseudo_crit_T(P):                       # P in [7.5, 14] MPa (Yang & Yang, 2012)
        p = P * 1e-6
        T_pc = -31.40 + 12.15 * p - 0.6927 * p**2 + 0.03160 * p**3 - 0.0007521 * p**4
        return T_pc + 273.15

    def find_h_cp_L(fluid):                          # Banuti (2019)
        P_L = PropsSI('PCRIT', fluid) * 0.1
        if P_L < PropsSI('PTRIPLE', fluid):
            P_L = PropsSI('PTRIPLE', fluid) * 1.05
        T_L = PropsSI('T', 'P', P_L, 'Q', 0, fluid)
        SC = 0.01
        cp_L = PropsSI('CPMASS', 'P', P_L, 'T', T_L - SC, fluid)
        h_L = PropsSI('H', 'P', P_L, 'T', T_L - SC, fluid)
        return cp_L, h_L - cp_L * T_L

    def T_plus_minus(P, fluid):                      # Banuti (2015)
        T_pc = find_pseudo_crit_T(P)
        cp_pc = PropsSI('CPMASS', 'T', T_pc, 'P', P, fluid)
        h_0_pc = PropsSI('H', 'T', T_pc, 'P', P, fluid)
        cp_L, h_L_0 = find_h_cp_L(fluid)
        T_max = PropsSI('TMAX', fluid)
        cp_IG = PropsSI('CPMASS', 'T', T_max, 'P', P, fluid)
        T_plus = (h_0_pc - cp_pc * T_pc) / (cp_IG - cp_pc)
        T_minus = (h_L_0 - h_0_pc + cp_pc * T_pc) / (cp_pc - cp_L)
        return T_pc, T_minus, T_plus

    T_pc, T_minus, T_plus = T_plus_minus(P, fluid)

    h_w = PropsSI('H', 'T', T_w, 'P', P, fluid)
    h_pc = PropsSI('H', 'T', T_pc, 'P', P, fluid)
    rho_w = PropsSI('D', 'T', T_w, 'P', P, fluid)

    SBO_w = q / (G * h_w)

    h_LL = PropsSI('H', 'T', T_minus, 'P', P, fluid)
    h_LV = PropsSI('H', 'T', T_plus, 'P', P, fluid)
    x_pb_in = (h_in - h_LL) / (h_LV - h_LL)          # pseudo-vapor quality
    x_pb_out = (h_out - h_LL) / (h_LV - h_LL)
    x_ave = 0.5 * (x_pb_in + x_pb_out)

    Re = G * D_in / mu
    if not 6.18e4 <= Re <= 5.35e5:
        warnings.warn(f"Cheng 2024 sCO2: Re_b = {Re:.3g} outside 6.18e4-5.35e5 (paper range).")
    # Bulk enthalpy and temperature of the cell (mean of inlet and outlet enthalpies)
    h_b = 0.5 * (h_in + h_out)
    T_b = PropsSI('T', 'H', h_b, 'P', P, fluid)
    dT = T_w - T_b
    if abs(dT) < 1e-6:
        dT = 1e-6
    cp_ave = (h_w - h_b) / dT                        # Eq. after (40): cp = (i_w - i_b) / (T_w - T_b)
    Pr_ave = mu * cp_ave / k

    if x_ave < 0 or x_ave > 1:                       # vapor-like / liquid-like
        return 0.019 * (Re**0.79) * (Pr_ave**0.51) * k / D_in

    g = 9.81
    Fr_VL_w = G**2 * x_ave**2 / (rho_w**2 * g * D_in)            # two-phase-like, Eq. (33)
    if not 7.58e-5 <= Fr_VL_w <= 1834:
        warnings.warn(f"Cheng 2024 sCO2: Fr_VL,w = {Fr_VL_w:.3g} outside 7.58e-5-1834 (paper range).")
    Nu_top = 0.127 * (Re**1.065) * (Pr_ave**1.053) * (SBO_w**0.66) * (Fr_VL_w**(-0.012))
    Nu_bot = 1.512 * (Re**0.878) * (Pr_ave**0.685) * (SBO_w**0.69) * (Fr_VL_w**(-0.0213))
    return 0.5 * (Nu_top + Nu_bot) * k / D_in


def pipe_petukhov_supercritical_HTC(AS, geom, cond):
    """
    Petukhov et al. (1983) correlation for single-phase turbulent flow in a tube at
    supercritical / near-critical pressure, with wall-property corrections of the
    Darcy friction factor.

    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.conv_supercritical.Nu_Petukhov, which cites the reference below
    (assessed in Chen, Fang, Xu, Su 2015, Annals of Nuclear Energy 76, 451-460,
    doi:10.1016/j.anucene.2014.10.027: second best in the normal heat transfer
    category, MAD 12.0 %).

    Arguments
    ---------
    AS : CoolProp AbstractState
        Fluid state object (will be moved to wall state (p, T_w)).
    geom : dict
        Geometry dictionary with keys:
        - 'Tube_ID' [m]  (or 'Tube_OD' [m] and 'Tube_t' [m])
    cond : dict
        Operating conditions dictionary with keys:
        - 'G' [kg/m2/s]  mass flux in ONE tube (or 'm_dot' [kg/s] + geom n_tubes, n_passes)
        - 'p' [Pa]       bulk pressure
        - 'T_w' [K]      wall temperature
        - 'rho' [kg/m3]  bulk density
        - 'mu' [Pa s]    bulk dynamic viscosity
        - 'Pr' [-]       bulk Prandtl number (uses bulk cp of the cell, not integral-averaged)
        - 'k' [W/m/K]    bulk thermal conductivity

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K].

    Reference
    ---------
    Petukhov, B.S., Kurganov, V.A., Ankudinov, V.B. (1983) Heat transfer and flow
    resistance in the turbulent pipe flow of a fluid with near-critical state
    parameters. High Temperature 21(1), 81-89.

    Validity
    --------
    Turbulent correlation. Warns if Re_b < 1e4 (labothappy safeguard, not from the paper).
    Pr_b uses the bulk cp of the cell (as in ht library), not an integral-averaged cp.

    Definitions
    -----------
    G : mass flux in one tube [kg/m2/s]
    Re_b, Pr_b, k_b : bulk properties of the cell (from cond)
    rho_w, mu_w : wall properties at (p, T_w)
    f : Darcy friction factor
    Nu : Nusselt number based on tube inner diameter D
    AS : moved to wall state (p, T_w) inside the function
    """
    # 1. Geometry and bulk conditions
    G = _tube_G(geom, cond)                     # mass flux in one tube [kg/m2/s]
    D = _tube_ID(geom)                          # tube inner diameter [m]
    p = cond['p']                               # bulk pressure [Pa]
    T_w = cond['T_w']                           # wall temperature [K]

    # Bulk properties of the cell (from cond dictionary)
    rho_b = cond['rho']                         # bulk density [kg/m3]
    mu_b = cond['mu']                           # bulk dynamic viscosity [Pa s]
    Pr_b = cond['Pr']                           # bulk Prandtl number [-]
    k_b = cond['k']                             # bulk thermal conductivity [W/m/K]

    # 2. Wall properties (move AS to wall state)
    AS.update(CP.PT_INPUTS, p, T_w)             # set AS to (p, T_w)
    rho_w = AS.rhomass()                        # wall density [kg/m3]
    mu_w = AS.viscosity()                       # wall dynamic viscosity [Pa s]

    # 3. Bulk Reynolds number
    Re_b = G * D / mu_b                         # bulk Reynolds number [-]

    # 4. Isothermal friction factor (Filonenko, Darcy)
    f_0 = (1.82 * np.log10(Re_b) - 1.64) ** -2  # Darcy friction factor [-]

    # 5. Friction factor with property corrections
    f = f_0 * (rho_w / rho_b) ** 0.4 * (mu_w / mu_b) ** 0.2  # corrected Darcy friction factor [-]

    # 6. Nusselt number (Petukhov equation)
    numerator = (f / 8.0) * Re_b * Pr_b
    denominator = 1.0 + 900.0 / Re_b + 12.7 * (f / 8.0) ** 0.5 * (Pr_b ** (2.0 / 3.0) - 1.0)
    Nu = numerator / denominator                # Nusselt number [-]

    # 7. Heat transfer coefficient
    h = Nu * k_b / D                            # HTC [W/m2/K]

    # Validity warning (turbulent correlation safeguard)
    if Re_b < 1e4:
        warnings.warn(
            'Petukhov supercritical HTC: Re_b = {:.0f} < 1e4. '
            'Correlation is for turbulent flow; this limit is a labothappy safeguard.'.format(Re_b)
        )

    return h

#%% CONDENSATION INSIDE TUBES

def pipe_el_hajal_flow_pattern(AS, geom, cond):
    """
    Flow pattern map for condensation inside horizontal plain tubes.
    Returns a dictionary containing the local flow regime and the quantities
    needed by the Thome condensation model. Bubbly flow is not identified
    (no heat transfer model for it in Part 2).

    Arguments
    ---------
    AS : CoolProp AbstractState
        Fluid state object (moved by saturation property calls).
    geom : dict
        Geometry dictionary with keys:
        - 'Tube_ID' [m] (or 'Tube_OD' [m] and 'Tube_t' [m])
        - 'n_tubes' [-], 'n_passes' [-] (optional, used if 'G' is absent in cond)
    cond : dict
        Operating conditions dictionary with keys:
        - 'p' [Pa]       saturation pressure
        - 'x' [-]        vapour quality (clamped to ]0, 1[)
        - 'G' [kg/m2/s]  mass flux in ONE tube (optional if 'm_dot' [kg/s] and geom are provided)
        - 'sigma' [N/m]  surface tension (optional if CoolProp provides it)

    Returns
    -------
    out : dict
        Dictionary with keys:
        - 'regime'      [-]  local flow regime ('S', 'SW', 'MF', 'A', 'I')
        - 'eps'         [-]  local void fraction (log-mean epsilon)
        - 'theta_strat' [rad] local stratified angle (Biberg expression)
        - 'G_wavy'      [kg/m2/s] local wavy mass flux transition threshold
        - 'G_strat'     [kg/m2/s] local stratified mass flux transition threshold
        - 'G_mist'      [kg/m2/s] local mist mass flux transition threshold
        - 'x_IA'        [-]  annular-intermittent transition quality

    Reference
    ---------
    El Hajal, J., Thome, J.R., Cavallini, A. (2003) Condensation in horizontal tubes,
    part 1: two-phase flow pattern map. International Journal of Heat and Mass Transfer
    46, 3349-3363.
    Equation numbers refer to this paper.

    Validity range and fluids tested
    --------------------------------
    16 < G < 1532 kg/m2/s; 3.14 mm < D < 21.4 mm; 0.02 < p / p_crit < 0.8;
    76 < (We/Fr)_L < 884. Twenty fluids including ammonia and refrigerants (paper Section 9).

    Definitions
    -----------
    - G : mass flux in one tube [kg/m2/s]
    - D : tube inner diameter [m]
    - p : pressure [Pa]
    - x : vapor quality [-]
    - rho_L, rho_V : liquid and vapor densities [kg/m3]
    - mu_L, mu_V : liquid and vapor dynamic viscosities [Pa s]
    - sigma : surface tension [N/m]
    - g : gravitational acceleration [m/s2] = 9.81 m/s2
    """
    # Void fraction models of labothappy (void_fraction.py): the functions behind
    # compute_void_fraction (called directly: it rebuilds two AbstractStates at each call)
    from labothappy.correlations.void_fraction.void_fraction import (
        void_fraction_homogeneous, void_fraction_rouhani_axelsson_steiner)

    # 1. Inputs and properties
    G = _tube_G(geom, cond)                     # mass flux inside one tube [kg/m2/s]
    D = _tube_ID(geom)                          # tube inner diameter [m]
    p = cond['p']                               # saturation pressure [Pa]
    x_local = _clamp_x(cond['x'])               # local quality [-]

    s = _sat_props(AS, p, cond)                 # saturated properties at p
    _need_sigma(AS, s)                          # ensure surface tension is available
    rho_L = s['rho_l']                          # liquid density [kg/m3]
    rho_V = s['rho_v']                          # vapor density [kg/m3]
    mu_L = s['mu_l']                            # liquid viscosity [Pa s]
    mu_V = s['mu_v']                            # vapor viscosity [Pa s]
    sigma = s['sigma']                          # surface tension [N/m]
    g = 9.81                                    # gravitational acceleration [m/s2]

    # 2. Quality grid (Step 1)
    x_grid = np.linspace(0.001, 0.999, 999)     # grid of qualities [-]
    x_all = np.append(x_grid, x_local)          # grid + local quality (local is last element) [-]

    # 3. Void fraction, LM epsilon (Eqs. 1, 2, 6) and geometry (Eqs. 7-12, Part 2 Eq. 19) (Step 2)
    # The labothappy void fraction functions take one quality at a time: loop over x_all
    eps_h = np.zeros(len(x_all))
    eps_ra = np.zeros(len(x_all))
    for i in range(len(x_all)):
        eps_h[i] = void_fraction_homogeneous(rho_L, rho_V, x_all[i])             # Eq. 1
        # Eq. 2, Steiner horizontal Rouhani-Axelsson
        eps_ra[i] = void_fraction_rouhani_axelsson_steiner(rho_L, rho_V, x_all[i], G, sigma)
    eps = (eps_h - eps_ra) / np.log(eps_h / eps_ra)                             # Eq. 6 (LM epsilon)

    A_Ld = np.pi * (1 - eps) / 4  # A_L / D^2 (Eqs. 7-9)
    A_Vd = np.pi * eps / 4                                                      # A_V / D^2

    term_biberg_1 = np.pi * (1 - eps)
    bracket_2 = 1 - 2 * (1 - eps) + (1 - eps)**(1 / 3) - eps**(1 / 3)
    term_biberg_2 = (3 * np.pi / 2)**(1 / 3) * bracket_2
    bracket_3 = (1 - 2 * (1 - eps)) * (1 + 4 * ((1 - eps)**2 + eps**2))
    term_biberg_3 = (1 / 200) * (1 - eps) * eps * bracket_3
    # Stratified angle, explicit expression of Biberg (Part 2 Eq. 19)
    theta_strat = 2 * np.pi - 2 * (term_biberg_1 + term_biberg_2 - term_biberg_3)

    h_Ld = 0.5 * (1 - np.cos((2 * np.pi - theta_strat) / 2))                    # Eq. 11

    # 4. Transition curves (Step 3)
    We_Fr_L = g * D**2 * rho_L / sigma                                          # Eq. 20 (scalar)

    # Stratified-wavy to intermittent / annular transition, condensation (Eq. 16)
    height_term = (1 - (2 * h_Ld - 1)**2)**0.5
    term_1 = 16 * A_Vd**3 * g * D * rho_L * rho_V / (x_all**2 * np.pi**2 * height_term)
    term_2 = np.pi**2 / (25 * h_Ld**2) * We_Fr_L**(-1.023) + 1
    exponent = -(x_all**2 - 0.97)**2 / (x_all * (1 - x_all))
    G_wavy = (term_1 * term_2)**0.5 + 50 - 75 * np.exp(exponent)

    term_strat_num = 226.3**2 * A_Ld * A_Vd**2 * rho_V * (rho_L - rho_V) * mu_L * g
    term_strat_den = x_all**2 * (1 - x_all) * np.pi**3
    G_strat = (term_strat_num / term_strat_den)**(1 / 3) + 20 * x_all           # Eq. 17

    x_IA_denom = 0.2914 * (rho_V / rho_L)**(-1 / 1.75) * (mu_L / mu_V)**(-1 / 7) + 1
    x_IA = 1 / x_IA_denom                                                       # Eq. 18 (scalar)

    xi = (1.138 + 2 * np.log10(np.pi / (1.5 * A_Ld)))**(-2)                     # Eq. 21
    G_mist_num = 7680 * A_Vd**2 * g * D * rho_L * rho_V
    G_mist_den = x_all**2 * np.pi**2 * xi * We_Fr_L
    G_mist = (G_mist_num / G_mist_den)**0.5                                     # Eq. 19

    # 5. Minima held for x > x_min (paper p. 3356) (Step 4)
    grid_len = len(x_grid)

    idx_min_wavy = np.argmin(G_wavy[:grid_len])
    x_min_wavy = x_grid[idx_min_wavy]
    G_wavy_min = G_wavy[idx_min_wavy]
    if x_local > x_min_wavy:
        G_wavy_local = G_wavy_min
    else:
        G_wavy_local = G_wavy[-1]

    idx_min_mist = np.argmin(G_mist[:grid_len])
    x_min_mist = x_grid[idx_min_mist]
    G_mist_min = G_mist[idx_min_mist]
    if x_local > x_min_mist:
        G_mist_local = G_mist_min
    else:
        G_mist_local = G_mist[-1]

    G_strat_local = G_strat[-1]
    eps_local = eps[-1]
    theta_strat_local = theta_strat[-1]

    # 6. Regime determination (paper p. 3356, logic in this order) (Step 5)
    if G < G_strat_local:
        regime = 'S'        # fully stratified
    elif G < G_wavy_local:
        regime = 'SW'       # stratified-wavy
    elif G > G_mist_local:
        regime = 'MF'       # mist
    elif x_local >= x_IA:
        regime = 'A'        # annular
    else:
        regime = 'I'        # intermittent

    # 7. Validity checks and warnings (paper Section 9)
    p_crit = AS.p_critical()
    p_red = p / p_crit

    outside_validity = False
    warning_messages = []

    if G <= 16 or G >= 1532:
        outside_validity = True
        warning_messages.append(f"Mass flux G = {G:.1f} kg/m2/s outside [16, 1532]")

    if D <= 3.14e-3 or D >= 21.4e-3:
        outside_validity = True
        warning_messages.append(f"Diameter D = {D * 1e3:.2f} mm outside [3.14, 21.4] mm")

    if p_red <= 0.02 or p_red >= 0.8:
        outside_validity = True
        warning_messages.append(f"Reduced pressure p_red = {p_red:.3f} outside [0.02, 0.8]")

    if We_Fr_L <= 76 or We_Fr_L >= 884:
        outside_validity = True
        warning_messages.append(f"(We/Fr)_L = {We_Fr_L:.1f} outside [76, 884]")

    if outside_validity:
        full_message = ("El Hajal flow pattern map: " + ", ".join(warning_messages)
                        + " (validity range exceeded): extrapolated.")
        warnings.warn(full_message)

    # 8. Return dictionary
    return {
        'regime': regime,
        'eps': eps_local,
        'theta_strat': theta_strat_local,
        'G_wavy': G_wavy_local,
        'G_strat': G_strat_local,
        'G_mist': G_mist_local,
        'x_IA': x_IA
    }


def pipe_thome_condensing_HTC(AS, geom, cond):
    """
    In-tube condensation heat transfer coefficient, Thome, El Hajal & Cavallini (2003), Part 2,
    with the flow pattern map of Part 1.

    Parameters
    ----------
    AS : CoolProp AbstractState
        Fluid state object (will be moved by saturation and wall property calls).
    geom : dict
        Geometry dictionary containing:
        - 'Tube_ID' [m] (or 'Tube_OD' [m] and 'Tube_t' [m])
        - 'n_tubes' [-], 'n_passes' [-] (optional, used if 'G' is absent in cond)
    cond : dict
        Conditions dictionary containing:
        - 'p' [Pa]       saturation pressure
        - 'x' [-]        thermodynamic vapor quality
        - 'T_w' [K]      inner wall temperature (required, raises KeyError if absent)
        - 'G' [kg/m2/s]  mass flux inside ONE tube (or 'm_dot' [kg/s] + geom n_tubes, n_passes)
        - 'm_dot' [kg/s] total mass flow rate (optional, used if 'G' is absent)

    Returns
    -------
    h : float
        Two-phase heat transfer coefficient [W/m2/K] on the inner surface.

    References
    ----------
    Thome, J.R., El Hajal, J., Cavallini, A. (2003) Condensation in horizontal tubes, part 2:
    new heat transfer model based on flow regimes. International Journal of Heat and Mass
    Transfer 46, 3365-3387.
    Flow pattern map: El Hajal, J., Thome, J.R., Cavallini, A. (2003) part 1, IJHMT 46, 3349-3363.

    Validity range and fluids tested
    --------------------------------
    24 <= G <= 1022 kg/m2/s; 3.1 <= D <= 21.4 mm; 0.02 <= p / p_crit <= 0.80.
    15 fluids tested: refrigerants (R-11, R-12, R-22, R-32, R-113, R-125, R-134a, R-236ea,
    R-32/R-125, R-404A, R-410A) and hydrocarbons (propane, n-butane, iso-butane, propylene).
    85% of refrigerant data within +-20%.

    Definitions and implementation details
    --------------------------------------
    - Void fraction: LM epsilon of Part 1 (log-mean of the homogeneous and of the Steiner
      horizontal Rouhani-Axelsson void fractions), from the labothappy void_fraction models,
      through pipe_el_hajal_flow_pattern.
    - Regimes and falling film angle theta:
      * 'A', 'I', 'MF': theta = 0.0 (fully annular, intermittent, mist)
      * 'SW': stratified-wavy, theta = theta_strat * ((G_wavy - G) / (G_wavy - G_strat))**0.5
      * 'S': fully stratified, theta = theta_strat; f_i multiplied by G / G_strat (Eq. 16)
    - delta = D/2 rule: whenever delta > d/2, delta is set equal to d/2; negative argument
      indicates liquid does not fit in the wetted sector even with delta = d/2.
    - Quality limits (Part 2 p. 3373): x > 0.99 is reset to 0.99; x < 0.01 (lower limit of
      applicability) emits a warning, x kept >= 1e-6. Data of the paper: 0.03 < x < 0.97.
    - Wall temperature: dT = max(T_sat - T_w, 1e-3) (labothappy choice to avoid a division by zero
      if T_w >= T_sat).
    """
    # 1. Inputs: G, D, p, T_w, g
    G = _tube_G(geom, cond)
    D = _tube_ID(geom)
    p = cond['p']
    if 'T_w' not in cond:
        raise KeyError("Thome condensing HTC: cond needs 'T_w' [K].")
    T_w = cond['T_w']
    g = 9.81

    # 2. Quality limits (Part 2 p. 3373)
    x = cond['x']
    if x > 0.99:
        x = 0.99
    elif x < 0.01:
        warnings.warn(f"Thome condensing HTC: quality x = {x} below the lower limit of "
                      "applicability x = 0.01 of Thome et al. (2003).")
        x = max(x, 1e-6)

    # 3. Saturated properties
    s = _sat_props(AS, p, cond)
    _need_sigma(AS, s)
    rho_L = s['rho_l']
    rho_V = s['rho_v']
    mu_L = s['mu_l']
    k_L = s['k_l']
    Pr_L = s['Pr_l']
    h_LV = s['h_lv']
    sigma = s['sigma']
    T_sat = s['T_sat']

    # 4. Flow pattern and void fraction (calls pipe_el_hajal_flow_pattern from the same module)
    fp = pipe_el_hajal_flow_pattern(AS, geom, {**cond, 'x': x})
    eps = fp['eps']
    theta_strat = fp['theta_strat']
    regime = fp['regime']

    # 5. Falling film angle theta (Eq. 6 and implementation steps 4-6)
    if regime in ('A', 'I', 'MF'):
        theta = 0.0
    elif regime == 'SW':
        numerator_SW = fp['G_wavy'] - G
        denominator_SW = fp['G_wavy'] - fp['G_strat']
        theta = theta_strat * (numerator_SW / denominator_SW)**0.5
    else:
        theta = theta_strat             # regime 'S', fully stratified

    # 6. Liquid area and film thickness (Eqs. 4, 10)
    A = np.pi * D**2 / 4.0
    A_L = (1.0 - eps) * A
    argument = D**2 - 8.0 * A_L / (2.0 * np.pi - theta)
    if argument < 0:
        delta = D / 2.0
    else:
        delta = (D - np.sqrt(argument)) / 2.0

    if delta > D / 2.0:
        delta = D / 2.0

    # 7. Velocities (Eqs. 11, 12)
    u_L = G * (1.0 - x) / (rho_L * (1.0 - eps))
    u_V = G * x / (rho_V * eps)

    # 8. Liquid film Reynolds number (Eq. 8)
    Re_L = 4.0 * G * (1.0 - x) * delta / ((1.0 - eps) * mu_L)

    # 9. Interfacial roughness factor (Eq. 15, Eq. 16 for regime 'S')
    roughness_term = (u_V / u_L)**0.5 * ((rho_L - rho_V) * g * delta**2 / sigma)**0.25
    if regime == 'S':
        f_i = 1.0 + roughness_term * (G / fp['G_strat'])
    else:
        f_i = 1.0 + roughness_term

    # 10. Convective condensation (Eq. 7)
    c = 0.003
    n = 0.74
    m = 0.5
    alpha_c = c * Re_L**n * Pr_L**m * k_L / delta * f_i

    # 11. Falling film condensation, Nusselt (Eq. 17)
    dT = max(T_sat - T_w, 1e-3)
    alpha_f = 0.728 * (rho_L * (rho_L - rho_V) * g * h_LV * k_L**3 / (mu_L * D * dT))**0.25

    # 12. Local coefficient (Eq. 2)
    r = D / 2.0
    alpha_tp = (alpha_f * r * theta + (2.0 * np.pi - theta) * r * alpha_c) / (2.0 * np.pi * r)

    # Validity range checks (warnings emitted, calculation completed anyway)
    p_crit = AS.p_critical()
    p_red = p / p_crit
    if G < 24 or G > 1022:
        warnings.warn(f"Thome condensing HTC: mass flux G = {G:.1f} kg/m2/s outside the "
                      "tested range [24, 1022].")
    if D < 0.0031 or D > 0.0214:
        warnings.warn(f"Thome condensing HTC: tube diameter D = {D * 1e3:.2f} mm outside the "
                      "tested range [3.1, 21.4] mm.")
    if p_red < 0.02 or p_red > 0.80:
        warnings.warn(f"Thome condensing HTC: reduced pressure p_red = {p_red:.3f} outside the "
                      "tested range [0.02, 0.80].")

    return alpha_tp


def pipe_dobson_chato_condensing_HTC(AS, geom, cond):
    """Dobson & Chato: gravity-driven film (Nu = 0.555 [...]^0.25) for Re_v <= 35000, annular
    Nu = 0.023 Re_l^0.8 Pr_l^0.4 (1 + 2.22 / X_tt^0.89) above.
    cond : G, p, x, T_w      geom : Tube_ID
    Reference: Incropera, Principles of heat transfer."""
    g = 9.81
    D_in = _tube_ID(geom)
    G, P_sat, x = _tube_G(geom, cond), cond['p'], min(max(cond['x'], 0.0), 1.0)
    s = _sat_props(AS, P_sat, cond)
    rho_l, rho_v, mu_l, mu_v = s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v']
    k_l, cp_l, Pr_l, h_fg, T_sat = s['k_l'], s['cp_l'], s['Pr_l'], s['h_lv'], s['T_sat']
    dT = max(T_sat - cond['T_w'], 1e-3)

    Re_v = G * D_in / mu_v                                    # all the flow as vapor
    def dobson_chato():
        h_2_fg = h_fg + 0.375 * cp_l * dT
        return 0.555 * ((rho_l * g * (rho_l - rho_v) * h_2_fg * D_in**3) / (mu_l * k_l * dT))**0.25

    if Re_v <= 35000:
        Nu = dobson_chato()
    else:
        Re_Dl = G * (1 - x) * D_in / mu_l
        if x > 0:
            x_tt = ((1 - x) / x)**0.9 * (rho_v / rho_l)**0.5 * (mu_l / mu_v)**0.1
            Nu = 0.023 * Re_Dl**0.8 * Pr_l**0.4 * (1 + 2.22 / x_tt**0.89) if x_tt > 0 else dobson_chato()
        else:
            Nu = dobson_chato()
    return Nu * k_l / D_in


def pipe_shah_condensing_HTC(AS, geom, cond):
    """
    Shah (1979, 1981) in-tube film condensation heat transfer coefficient.

    Implementation follows the independent implementation in the ht library
    (C. Bell), ht.condensation.Shah, which cites the references below (also given
    in Kakac, Boilers, Evaporators and Condensers, Wiley 1991).

    Parameters
    ----------
    AS : CoolProp AbstractState
        Fluid state object.
    geom : dict
        Geometry dictionary containing 'Tube_ID' [m] (or 'Tube_OD' [m] and 'Tube_t' [m]),
        straight tube length 'Tube_L' [m] (if needed), and optionally 'n_tubes' [-]
        and 'n_passes' [-].
    cond : dict
        Conditions dictionary containing:
        - 'p' : saturation pressure [Pa]
        - 'x' : vapor quality [-]
        - 'G' : mass flux inside ONE tube [kg/(m^2 s)] (optional if 'm_dot' [kg/s]
          and geom are provided)

    Returns
    -------
    float
        Heat transfer coefficient h [W/(m^2 K)] on the inner surface.

    References
    ----------
    Shah, M.M. (1979) A general correlation for heat transfer during film
    condensation inside pipes. International Journal of Heat and Mass Transfer
    22(4), 547-556. doi:10.1016/0017-9310(79)90058-9.
    Shah, M.M. (1981) ASHRAE Transactions 87(3), 1086-1100.

    Validity range
    --------------
    0.002 <= p_red <= 0.44; Re_lo >= 350.
    Fluids tested in the original 1979 data: water, R11, R12, R22, R113, methanol,
    ethanol, benzene, toluene, trichloroethylene.
    Diameter D: 7-40 mm.
    Mass flux G: 11-211 kg/(m^2 s) (ranges of the 1979 data "as usually reported").

    Definitions
    -----------
    - G of one tube [kg/(m^2 s)]
    - Re_lo with the total flow as liquid [-]
    - Saturated-liquid properties (mu_l, k_l, Pr_l)
    - Reduced pressure p_red = p / p_crit [-]
    - No vapour property used
    - Returns h [W/(m^2 K)] on the inner surface
    """
    # 1. Mass flux in one tube, internal diameter, pressure, and clamped quality
    G = _tube_G(geom, cond)
    D = _tube_ID(geom)
    p = cond['p']
    x = _clamp_x(cond['x'])

    # 2. Saturated-liquid properties at p
    s = _sat_props(AS, p, cond)
    mu_l = s['mu_l']
    k_l = s['k_l']
    Pr_l = s['Pr_l']

    # 3. Reduced pressure
    p_crit = AS.p_critical()
    p_red = p / p_crit

    # 4. Liquid-only Reynolds number (all the flow as liquid)
    Re_lo = G * D / mu_l

    # 5. Liquid-only HTC, Dittus-Boelter correlation
    h_lo = 0.023 * Re_lo**0.8 * Pr_l**0.4 * k_l / D

    # 6. Two-phase multiplier terms
    term_liquid = (1 - x)**0.8
    term_vapour = 3.8 * x**0.76 * (1 - x)**0.04 / p_red**0.38

    # 7. Two-phase heat transfer coefficient
    h = h_lo * (term_liquid + term_vapour)

    # Validity checks and warnings
    outside_validity = False
    warning_messages = []

    if p_red < 0.002 or p_red > 0.44:
        outside_validity = True
        warning_messages.append(f"Reduced pressure p_red = {p_red:.4f} outside [0.002, 0.44]")

    if Re_lo < 350:
        outside_validity = True
        warning_messages.append(f"Liquid-only Reynolds Re_lo = {Re_lo:.1f} below 350")

    if outside_validity:
        full_message = ("Shah condensing HTC: " + ", ".join(warning_messages)
                        + " (validity range exceeded): extrapolated.")
        warnings.warn(full_message)

    return h


def pipe_cavallini_smith_zecchin_condensing_HTC(AS, geom, cond):
    """Cavallini, Smith & Zecchin (1974) forced-convection condensation inside tubes.
    Single correlation for annular flow regime using an equivalent Reynolds number.
    cond : G [kg/m2/s] (or m_dot [kg/s] + geom n_tubes, n_passes), p [Pa], x [-]
    geom : Tube_ID [m] (or Tube_OD [m] and Tube_t [m])
    Returns h [W/m2/K].
    Reference: Cavallini, A., Smith, J.R., Zecchin, R. (1974) A dimensionless correlation for
    heat transfer in forced convection condensation. 6th International Heat Transfer Conference,
    Tokyo, 309-313.
    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.condensation.Cavallini_Smith_Zecchin, which cites the reference above (also given in
    Kakac, Boilers, Evaporators and Condensers, Wiley 1991).
    Validity range not given in the reference library; forced-convection (annular) condensation.
    Labothappy safeguard: warns if Re_eq < 1000 (turbulent Dittus-Boelter-type form),
    not from the paper.
    Definitions: G per tube; Re_l, Re_v with each phase flowing alone; saturated properties at p;
    Nu based on tube ID with liquid thermal conductivity; quality clamped to [1e-6, 1-1e-6]."""
    # Step 1: mass flux, diameter, pressure, quality
    G = _tube_G(geom, cond)
    D = _tube_ID(geom)
    p = cond['p']
    x = _clamp_x(cond['x'])

    # Step 2: saturated properties at pressure p
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    mu_v = s['mu_v']
    k_l = s['k_l']
    Pr_l = s['Pr_l']

    # Step 3: liquid-phase Reynolds number (liquid flowing alone)
    Re_l = G * (1.0 - x) * D / mu_l

    # Step 4: vapour-phase Reynolds number (vapour flowing alone)
    Re_v = G * x * D / mu_v

    # Step 5: equivalent Reynolds number (Cavallini, Smith, Zecchin 1974, as implemented in ht)
    Re_eq = Re_v * (mu_v / mu_l) * (rho_l / rho_v)**0.5 + Re_l

    # Validity safeguard: turbulent correlation form
    if Re_eq < 1000:
        warnings.warn('Cavallini-Smith-Zecchin: Re_eq < 1000, correlation is for turbulent '
                      'annular flow.')

    # Step 6: Nusselt number (liquid properties, Dittus-Boelter type)
    Nu = 0.05 * Re_eq**0.8 * Pr_l**0.33

    # Step 7: heat transfer coefficient
    h = Nu * k_l / D
    return h

#%% BOILING INSIDE TUBES

# Reference nucleate boiling HTC a_0 [W/m2/K] at q_0, P* = 0.1, D = 10 mm, Ra = 1 um (VDI Heat Atlas)
_A0_VDI = {
    "Methane": 8060, "CH4": 8060, "Ethane": 5210, "C2H6": 5210, "Propane": 4000, "C3H8": 4000,
    "n-Butane": 3300, "C4H10": 3300, "n-Pentane": 3070, "C5H12": 3070, "Cyclopentane": 3070,
    "Isopentane": 2940, "n-Hexane": 2840, "C6H14": 2840, "n-Heptane": 2420, "C7H16": 2420,
    "Cyclohexane": 2420, "C6H12": 2420, "Benzene": 2730, "C6H6": 2730, "Toluene": 2910, "C7H8": 2910,
    "Diphenyl": 2030, "C12H10": 2030, "Methanol": 2770, "CH4O": 2770, "Ethanol": 3690, "C2H6O": 3690,
    "n-Propanol": 3170, "C3H8O": 3170, "Isopropanol": 2920, "n-Butanol": 2750, "C4H10O": 2750,
    "Isobutanol": 2940, "Acetone": 3270, "C3H6O": 3270, "R11": 2690, "CCl3F": 2690,
    "R12": 3290, "CCl2F2": 3290, "R13": 3910, "CClF3": 3910, "CBrF3": 3380, "R22": 3930, "CHClF2": 3930,
    "R23": 4870, "CHF3": 4870, "R113": 2180, "C2Cl3F3": 2180, "R114": 2460, "C2Cl2F4": 2460,
    "R115": 2890, "C2ClF5": 2890, "R123": 2600, "C2H2Cl2F3": 2600, "R134a": 3500, "C2H2F4": 3500,
    "R152a": 4000, "C2H4F2": 4000, "R226": 3700, "C3HClF6": 3700, "R227": 3800, "C3HF7": 3800,
    "RC318": 2710, "C4F8": 2710, "R502": 2900, "C2F5Cl": 2900, "Chloromethane": 4790, "CH3Cl": 4790,
    "Tetrachloromethane": 2320, "CCl4": 2320, "Tetrafluoromethane": 4500, "CF4": 4500,
    "Helium": 1990, "He": 1990, "Hydrogen": 12220, "H2": 12220, "Neon": 8920, "Ne": 8920,
    "Nitrogen": 4380, "N2": 4380, "Argon": 3870, "Ar": 3870, "Oxygen": 4120, "O2": 4120,
    "CarbonDioxide": 18890, "CO2": 18890,
    "n-Propane": 4000, "Water": 25580, "H2O": 25580, "Ammonia": 36640, "NH3": 36640,
    "SulfurHexafluoride": 12230, "SF6": 12230, "R13B1": 3380,
}  # fluids missing here (R1233zd(E), ...): give cond['a_0']
# 2026-11-21: values checked against Steiner & Taborek (1992) Table 2 (from VDI-Waermeatlas 1991):
# NH3 was 7000 and CO2 5500 (wrong); Table 2 gives 36,640 and 18,890 at q_o = 150 kW/m2.
# The a_0 of water, NH3, CO2, SF6 are given at q_o = 150 kW/m2, cryogenics at 10 kW/m2, helium at
# 1 kW/m2, R13B1 at 30 kW/m2, all others at 20 kW/m2 (Table 1 / Table 2).


def _a0(fluid, cond):
    if 'a_0' in cond:
        return cond['a_0']
    if fluid in _A0_VDI:
        return _A0_VDI[fluid]
    raise KeyError(f"No reference nucleate boiling HTC a_0 for {fluid!r}: give cond['a_0'].")


def pipe_steiner_taborek_boiling_HTC(AS, geom, cond):
    """Steiner-Taborek asymptotic model (Wolverine Engineering Data Book, Boiling inside plain tubes).
    cond : G, p, x, q (local heat flux); h_nb_o (optional, reference nucleate boiling HTC at the
           normalized heat flux q_o of the fluid group below and p_r = 0.1;
           default: VDI a_0 of the fluid, or the cyclopentane estimate of the original code)
    geom : Tube_ID, Tube_L
    Convective part: liquid-only (total flow) Gnielinski/Sieder-Tate HTC, no wall viscosity correction.
    Reference: Steiner, D., Taborek, J. (1992) Flow boiling heat transfer in vertical tubes correlated
    by an asymptotic model, Heat Transfer Engineering 13(2), 43-69, Eqs. (3)-(15), Tables 1-3.
    Normalized heat flux q_o (Table 1): 150 kW/m2 inorganic fluids (water, NH3, CO2, SF6), 20 kW/m2
    hydrocarbons / refrigerants / organics (30 kW/m2 R13B1, Table 2), 10 kW/m2 cryogenics (H2 to O2),
    1 kW/m2 helium. Validity: vertical tubes, rho_L / rho_G 3.75-5000; Eq. (5) F_tp for x <= 0.6.
    (2026-11-21 corrections vs the paper: q_o by fluid group (was 20 kW/m2 for all); convective-only
    branch Eq. (5a) (alpha_LO multiplier, exponent 0.01, density ratio placement); F(M) <= 2.5 and
    Table 3 values for H2, He, NH3; nf Eq. (10b) for cryogenics.)"""
    fluid = AS.fluid_names()[0]
    D_in, L = _tube_ID(geom), geom['Tube_L']
    G, P_sat, x, q_act = _tube_G(geom, cond), cond['p'], _clamp_x(cond['x']), cond['q']

    s = _sat_props(AS, P_sat, cond)
    _need_sigma(AS, s)
    rho_l, rho_v, mu_l, k_l, Pr_l = s['rho_l'], s['rho_v'], s['mu_l'], s['k_l'], s['Pr_l']
    mu_v, k_v, Pr_v = s['mu_v'], s['k_v'], s['Pr_v']
    sigma, T_sat, h_lv = s['sigma'], s['T_bub'], s['h_lv']
    P_r = P_sat / AS.p_critical()

    h_L = _gnielinski_core(mu_l, Pr_l, mu_l, k_l, G, D_in, L)

    r_o = 0.3e-6
    q_ONB = (2 * sigma * T_sat * h_L) / (r_o * rho_v * h_lv)         # onset of nucleate boiling

    if q_act < q_ONB:                                                # convective boiling only
        # Eq. (5a): alpha_cb / alpha_LO = F_tp, h_v = alpha_GO (total flow as gas)
        h_v = _gnielinski_core(mu_v, Pr_v, mu_v, k_v, G, D_in, L)
        a = (1 - x)**1.5 + 1.9 * x**0.6 * (1 - x)**0.01 * (rho_l / rho_v)**0.35
        b = (h_v / h_L) * x**0.01 * (1 + 8 * (1 - x)**0.7 * (rho_l / rho_v)**0.67)
        F_tp = 1 if (a == 0 or b == 0) else (a**(-2.2) + b**(-2))**(-0.5)
        return F_tp * h_L                                            # Eq. (4)

    # Nucleate + convective (asymptotic, exponent 3)
    M = AS.molar_mass() * 1e3                                        # g/mol

    # Fluid groups of Table 1 (normalized heat flux q_o) and Eq. (10b) (cryogenics)
    inorganic = ('Water', 'H2O', 'Ammonia', 'NH3', 'CarbonDioxide', 'CO2', 'SulfurHexafluoride', 'SF6')
    cryogenic = ('Hydrogen', 'H2', 'ParaHydrogen', 'Neon', 'Ne', 'Nitrogen', 'N2', 'Argon', 'Ar',
                 'Oxygen', 'O2')
    helium = ('Helium', 'He')
    if fluid in helium:
        q_o = 1000
    elif fluid in cryogenic:
        q_o = 10000
    elif fluid in inorganic:
        q_o = 150000
    elif fluid in ('R13B1', 'CBrF3'):
        q_o = 30000
    else:
        q_o = 20000

    # Eq. (5) F_tp is given for x <= 0.6
    if x > 0.6:
        warnings.warn(f"Steiner-Taborek: x = {x:.2f} > 0.6, Eq. (5) F_tp used outside its stated "
                      "range (Steiner & Taborek 1992).")
    if 'h_nb_o' in cond:
        h_nb_o = cond['h_nb_o']
    elif fluid == 'Cyclopentane':
        h_nb_o = 3070 * (2420 / 2840)    # assumption: h_nb,o(cyclopentane) = h_nb,o(cyclohexane) * (n-pentane / n-hexane)
    else:
        h_nb_o = _a0(fluid, cond)
    D_in_o = 0.01

    F_tp = ((1 - x)**1.5 + 1.9 * x**0.6 * (rho_l / rho_v)**0.35)**1.1
    F_pf = 2.816 * P_r**0.45 + (3.4 + (1.7 / (1 - P_r**7))) * P_r**3.7
    if fluid in cryogenic or fluid in helium:
        nf = 0.7 - 0.13 * np.exp(1.105 * P_r)                        # Eq. (10b), cryogenics
    else:
        nf = 0.8 - 0.1 * np.e**(1.75 * P_r)                          # Eq. (10a)

    # Residual factor F(M): Table 3 for H2, He, NH3; otherwise Eq. (14), limited to 2.5
    if fluid in ('Hydrogen', 'H2', 'ParaHydrogen'):
        F_M = 0.35
    elif fluid in helium:
        F_M = 0.86
    elif fluid in ('Ammonia', 'NH3'):
        F_M = 1.24
    else:
        F_M = min(0.377 + 0.199 * np.log(M) + 0.000028427 * M**2, 2.5)
    R_p, R_p_o = 1e-6, 1e-6
    F_nb = F_pf * (q_act / q_o)**nf * (R_p / R_p_o)**0.133 * (D_in / D_in_o)**(-0.4) * F_M
    return ((h_nb_o * F_nb)**3 + (h_L * F_tp)**3)**(1 / 3)


def pipe_vdi_boiling_HTC(AS, geom, cond):
    """VDI Heat Atlas (2010) flow boiling in horizontal tubes (Gorenflo-type nucleate form).
    cond : G, p, x, q; a_0 (optional if the fluid is not in the VDI table)
    geom : Tube_ID; Ra (wall roughness, default 3e-6 m)"""
    fluid = AS.fluid_names()[0]
    D_in = _tube_ID(geom)
    G, P_sat, q = _tube_G(geom, cond), cond['p'], cond['q']
    x = min(max(cond['x'], 0.0), 1.0)

    if fluid in ('Helium', 'He'):
        q_0 = 1000
    elif fluid in ('H2', 'N2', 'Neon', 'Air', 'O2', 'Ar', 'Hydrogen', 'Nitrogen', 'Oxygen', 'Argon'):
        q_0 = 10000
    elif fluid in ('H2O', 'NH3', 'CO2', 'SF6', 'Water', 'Ammonia', 'CarbonDioxide',
                   'SulfurHexafluoride'):     # CoolProp names added 2026-11-21
        q_0 = 150000
    else:
        q_0 = 20000
    a_0 = _a0(fluid, cond)

    P_crit = AS.p_critical()
    P_01 = 0.1 * P_crit

    AS.update(CP.PQ_INPUTS, P_01, 0)
    MM_fluid, rho_l_01, h_l_01 = AS.molar_mass(), AS.rhomass(), AS.hmass()
    AS.update(CP.PQ_INPUTS, P_01, 1)
    h_v_01, rho_v_01 = AS.hmass(), AS.rhomass()
    sigma_01 = 0.010449499360652493 if fluid == 'R134a' else _sigma(AS, P_01)

    P_star = P_sat / P_crit
    n = 0.9 - 0.36 * P_star**0.13
    MM_H2 = 0.00201588
    C_F = 0.789 * (MM_fluid / MM_H2)**0.11
    F_p = 2.692 * P_star**0.43 + (1.65 * P_star**6.5 / (1 - P_star**4.4))
    F_d = (1e-2 / D_in)**0.5
    F_W = (geom.get('Ra', 3e-6) / 1e-6)**0.133
    F_G = (G / 100.0)**0.25

    g = 9.81
    q_cr_01 = 0.13 * (h_v_01 - h_l_01) * rho_v_01**0.5 * (sigma_01 * g * (rho_l_01 - rho_v_01))**0.25
    q_cr_PB = 2.79 * q_cr_01 * P_star**0.4 * (1 - P_star)
    F_x = 1 - P_star**0.1 * (q / q_cr_PB)**0.3 * x

    return a_0 * C_F * (q / q_0)**n * F_p * F_d * F_W * F_G * F_x


def pipe_gungor_winterton_boiling_HTC(AS, geom, cond):
    """Gungor & Winterton (1986): h = E h_l + S h_nb, h_nb = Cooper pool boiling.
    cond : G, p, x, q      geom : Tube_ID
    Reference: Gungor, K.E., Winterton, R.H.S. (1986) A general correlation for flow boiling in tubes
    and annuli, Int. J. Heat Mass Transfer 29(3), 351-358, Eqs. (2)-(12).
    Definitions: Re_l = G (1 - x) d / mu_l (liquid fraction) in h_l AND in the suppression factor S
    (Eq. 10); Fr = G^2 / (rho_l^2 g d); properties at saturation; returns h [W/m2/K].
    (2026-11-21: S used a two-phase Re_tp = G d / mu_tp before; corrected to Re_l, Eq. 10.)"""
    D_in = _tube_ID(geom)
    G, P_sat, q = _tube_G(geom, cond), cond['p'], cond['q']
    x = _clamp_x(cond['x'])

    s = _sat_props(AS, P_sat, cond)
    rho_l, rho_v, mu_l, mu_v, k_l, Pr_l, h_lv = s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v'], s['k_l'], s['Pr_l'], s['h_lv']
    P_crit, MM = AS.p_critical(), AS.molar_mass()
    P_r = P_sat / P_crit

    # liquid-fraction convective HTC: Dittus-Boelter with Re_l = G (1 - x) D / mu_l
    Re_l = G * (1 - x) * D_in / mu_l
    h_f = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D_in

    g = 9.81
    Fr_l = G**2 / (rho_l**2 * g * D_in)                       # liquid Froude number (threshold 0.05)
    X_tt = ((1 - x) / x)**0.9 * (rho_v / rho_l)**0.5 * (mu_l / mu_v)**0.1
    Bo = q / (G * h_lv)

    if Fr_l <= 0.05:
        E_2, S_2 = Fr_l**(0.1 - 2 * Fr_l), Fr_l**0.5
    else:
        E_2, S_2 = 1.0, 1.0

    E = E_2 * (1 + 24000 * Bo**1.16 + 1.37 * (1 / X_tt)**0.86)
    # Gungor & Winterton (1986) Eq. (10): S uses the liquid Reynolds number Re_l = G (1 - x) d / mu_l
    S = S_2 / (1 + 1.15e-6 * E**2 * Re_l**1.17)

    h_nb = 55 * P_r**0.12 * (-0.4343 * np.log(P_r))**(-0.55) * (MM * 1e3)**(-0.5) * q**0.67   # Cooper
    return E * h_f + S * h_nb


def pipe_choi_boiling_HTC(AS, geom, cond):
    """Choi-type flow boiling: h = E h_l + S F_m h_nb (Lockhart-Martinelli E, convection number S).
    cond : G, p, x, q      geom : Tube_ID"""
    D_i = _tube_ID(geom)
    G, P_sat, q = _tube_G(geom, cond), cond['p'], cond['q']
    x = _clamp_x(cond['x'])

    s = _sat_props(AS, P_sat, cond)
    _need_sigma(AS, s)
    rho_l, rho_v, mu_l, mu_v, k_l, Pr_l, cp_l, h_lv = (s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v'],
                                                      s['k_l'], s['Pr_l'], s['cp_l'], s['h_lv'])
    T_sat_bubble, T_sat_dew, T_sat, sigma = s['T_bub'], s['T_dew'], s['T_sat'], s['sigma']

    Re_l = G * (1 - x) * D_i / mu_l
    h_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * (k_l / D_i)

    Bo = q / (G * h_lv)
    X_tt = ((1 - x) / x)**0.9 * (rho_v / rho_l)**0.5 * (mu_l / mu_v)**0.1
    E = 49.971 * Bo**0.383 * X_tt**(-0.758)

    Co = ((1 - x) / x)**0.9 * (rho_v / rho_l)**0.5
    S = 0.909 * Co**0.301

    g, beta = 9.81, 35                                         # beta: contact angle [deg]
    bd = 0.0146 * beta * (2 * sigma * (rho_l - rho_v) / g)**0.5
    h_nb = 207 * k_l / bd * ((q * bd) / (k_l * T_sat))**0.674 * (rho_v / rho_l)**0.581 * Pr_l**0.533

    F_m = (1 + 0.039 * cp_l * (T_sat_dew - T_sat_bubble) / h_lv)**(-1)
    return E * h_l + S * F_m * h_nb


def pipe_chen_boiling_HTC(AS, geom, cond):
    """
    Chen (1966) saturated flow boiling in tubes with Edelstein fit of the suppression factor S.
    Superposition model: h = S * h_nb + F * h_l.

    Arguments
    ---------
    AS   : CoolProp AbstractState (will be moved by saturation / wall property calls).
    geom : dict with keys
           'Tube_ID' [m]  (or 'Tube_OD' [m] and 'Tube_t' [m]).
    cond : dict with keys
           'G'       [kg/m2/s]  mass flux in ONE tube (or 'm_dot' [kg/s] + geom n_tubes, n_passes),
           'p'       [Pa]       saturation pressure,
           'x'       [-]        vapour quality (clamped to ]0, 1[),
           'T_w'     [K]        wall temperature (required),
           'q'       [W/m2]     NOT used by this correlation (present for API consistency).

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m2/K].

    Reference
    ---------
    Chen, J.C. (1966) Correlation for boiling heat transfer to saturated fluids in convective flow.
    Industrial & Engineering Chemistry Process Design and Development 5(3), 322-329.
    doi:10.1021/i260019a023.
    S factor fit: Edelstein, S., Perez, A.J., Chen, J.C. (1984) AIChE Journal 30(5), 840-841,
    doi:10.1002/aic.690300528.
    Nucleate term: Forster & Zuber (1955).
    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.boiling_flow.Chen_Edelstein, which cites the references above (model assembled as in
    Bertsch, Groll, Garimella 2008, Nanoscale Microscale Thermophys. Eng. 12(3), 187-227).

    Validity
    --------
    Quality range of the original water / hydrocarbon data as usually reported (0.01-0.71);
    not given by ht. Warning emitted outside this range.

    Definitions
    -----------
    G          : mass flux in one tube [kg/m2/s].
    Re_l       : liquid-only Reynolds number (liquid flowing alone in the tube).
    X_tt       : Martinelli parameter (turbulent-turbulent).
    Properties : saturated at pressure p (T_sat = mean of bubble and dew temperatures for mixtures).
    dT_sat     : T_w - T_sat. If dT_sat <= 0: no nucleate term (labothappy choice).
    dp_sat     : p_sat(T_w) - p.
    AS         : is moved by the saturation and wall property calls.
    """
    # 1. Geometry, flow, pressure, quality, wall temperature
    G = _tube_G(geom, cond)                     # [kg/m2/s]
    D = _tube_ID(geom)                          # [m]
    p = cond['p']                               # [Pa]
    x = _clamp_x(cond['x'])                     # [-]
    if 'T_w' not in cond:
        raise KeyError("Chen boiling: cond must contain 'T_w' [K] (wall temperature).")
    T_w = cond['T_w']                           # [K]

    # 2. Saturated properties at pressure p
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']                          # [kg/m3]
    rho_v = s['rho_v']                          # [kg/m3]
    mu_l = s['mu_l']                            # [Pa s]
    mu_v = s['mu_v']                            # [Pa s]
    k_l = s['k_l']                              # [W/m/K]
    cp_l = s['cp_l']                            # [J/kg/K]
    Pr_l = s['Pr_l']                            # [-]
    h_lv = s['h_lv']                            # [J/kg]
    sigma = s['sigma']                          # [N/m]
    T_sat = s['T_sat']                          # [K]
    _need_sigma(AS, s)                          # ensure sigma is available

    # 3. Liquid-phase Reynolds number (liquid flowing alone)
    Re_l = D * G * (1.0 - x) / mu_l             # [-]

    # 4. Liquid-alone HTC (Dittus-Boelter, Eq. in Chen 1966)
    h_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D   # [W/m2/K]

    # 5. Martinelli parameter (turbulent-turbulent)
    X_tt = ((1.0 - x) / x)**0.9 * (rho_v / rho_l)**0.5 * (mu_l / mu_v)**0.1   # [-]

    # 6. Enhancement factor F (Chen 1966)
    F = (1.0 + X_tt**-0.5)**1.78                # [-]

    # 7. Two-phase Reynolds number
    Re_tp = Re_l * F**1.25                      # [-]

    # 8. Suppression factor S (Edelstein fit, 1984)
    S = 0.9622 - 0.5822 * np.arctan(Re_tp / 6.18e4)   # [-]

    # 9. Wall superheat
    dT_sat = T_w - T_sat                        # [K]

    # 10. Nucleate boiling HTC (Forster-Zuber, 1955)
    if dT_sat <= 0.0:
        # Labothappy choice: no nucleate boiling without positive superheat
        # (Forster-Zuber term needs dT_sat > 0 and dp_sat > 0).
        h_nb = 0.0
    else:
        # Saturation pressure at wall temperature
        AS.update(CP.QT_INPUTS, 0, T_w)         # move AS to saturated liquid at T_w
        p_sat_w = AS.p()                        # [Pa]
        dp_sat = p_sat_w - p                    # [Pa]

        # Forster-Zuber written on several lines for readability
        numerator = k_l**0.79 * cp_l**0.45 * rho_l**0.49
        denominator = sigma**0.5 * mu_l**0.29 * h_lv**0.24 * rho_v**0.24
        h_nb = 0.00122 * numerator / denominator * dT_sat**0.24 * dp_sat**0.75   # [W/m2/K]

    # 11. Total HTC (Chen superposition)
    h = S * h_nb + F * h_l                      # [W/m2/K]

    # Validity warning (quality range from original data, not given by ht)
    if x < 0.01 or x > 0.71:
        warnings.warn(
            f"Chen boiling: quality x = {x:.3f} outside usual range 0.01-0.71 "
            "(original water/hydrocarbon data range as usually reported; not given by ht)."
        )

    return h


def pipe_liu_winterton_boiling_HTC(AS, geom, cond):
    """
    Saturated flow boiling in tubes using the Liu and Winterton (1991) asymptotic superposition
    model with the Cooper pool boiling correlation.

    Arguments:
    AS   : CoolProp AbstractState of the fluid (modified in-place for saturation properties)
    geom : dict containing tube geometry with SI units:
           'Tube_ID' [m] (or 'Tube_OD' [m] and 'Tube_t' [m])
           'n_tubes' [-], 'n_passes' [-] (optional, used if 'G' is absent in cond)
    cond : dict containing flow conditions with SI units:
           'p' [Pa] - saturation pressure
           'x' [-] - thermodynamic vapor quality
           'T_w' [K] - inner wall temperature
           'G' [kg/(m^2 s)] - mass flux inside one tube (optional if 'm_dot' and geom are provided)
           'm_dot' [kg/s] - total mass flow rate (optional, used if 'G' is absent)

    Returns:
    h : float - two-phase heat transfer coefficient [W/(m^2 K)]

    Reference:
    Liu, Z., Winterton, R.H.S. (1991) A general correlation for saturated and subcooled flow
    boiling in tubes and annuli, based on a nucleate pool boiling equation. International
    Journal of Heat and Mass Transfer 34(11), 2759-2766. doi:10.1016/0017-9310(91)90234-6.

    Validity range and fluids tested:
    About 5000 data points (water, refrigerants, ethylene glycol; saturated and subcooled boiling)
    according to the paper title and abstract; range not checked.

    Definitions:
    G of one tube; Re_lo with the total flow as liquid; saturated properties; dT_sat = T_w - T_sat
    (if dT_sat <= 0: no nucleate term, labothappy choice); Cooper with Rp = 1 um and M in g/mol;
    returns h [W/m^2/K].

    Implementation follows the independent implementation in the ht library (C. Bell),
    ht.boiling_flow.Liu_Winterton (wall-superheat form), which cites the reference above
    (reviewed in Fang, Zhou, Li 2013, IJR 36(8), 2017-2039). The low-Froude correction for
    horizontal tubes of the paper is not implemented (as in ht).
    """
    # 1. Read geometry, mass flux, pressure, quality, and wall temperature
    G = _tube_G(geom, cond)
    D = _tube_ID(geom)
    p = cond['p']

    if 'T_w' not in cond:
        raise KeyError("Liu-Winterton boiling HTC: cond needs 'T_w' [K].")
    T_w = cond['T_w']

    x = _clamp_x(cond['x'])

    # 2. Get saturated liquid properties
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']
    rho_v = s['rho_v']
    mu_l = s['mu_l']
    k_l = s['k_l']
    Pr_l = s['Pr_l']
    T_sat = s['T_sat']

    # 3. Critical pressure and molar mass
    p_crit = AS.p_critical()
    p_red = p / p_crit
    if p_red >= 0.95:
        warnings.warn("Liu-Winterton boiling HTC: reduced pressure p_red >= 0.95 (validity limit).")

    M = AS.molar_mass() * 1000.0

    # 4. Liquid-only Reynolds number
    Re_lo = G * D / mu_l

    # 5. Liquid-only heat transfer coefficient (Dittus-Boelter equation)
    h_lo = 0.023 * Re_lo**0.8 * Pr_l**0.4 * k_l / D

    # 6. Enhancement factor
    F = (1.0 + x * Pr_l * (rho_l / rho_v - 1.0))**0.35

    # 7. Suppression factor
    S = 1.0 / (1.0 + 0.055 * F**0.1 * Re_lo**0.16)

    # 8. Wall superheat and nucleate boiling contribution (Cooper pool boiling)
    dT_sat = T_w - T_sat
    if dT_sat <= 0.0:
        # Labothappy choice: no nucleate boiling without positive superheat
        h_nb = 0.0
    else:
        # Cooper pool boiling, wall-superheat form, roughness Rp = 1 um (log10(Rp) = 0)
        inner = 55.0 * dT_sat**0.67 * p_red**0.12 * (-np.log10(p_red))**-0.55 * M**-0.5
        h_nb = inner**(1.0 / 0.33)

    # 9. Asymptotic superposition combining convective and nucleate boiling terms
    h = ((F * h_lo)**2.0 + (S * h_nb)**2.0)**0.5

    return h

def pipe_shah_1982_boiling_HTC(AS, geom, cond):
    """
    Saturated flow boiling heat transfer coefficient inside a tube using the Shah (1982)
    chart correlation in equation form.

    Equations numbered as in the independent implementation of the ThermoCycle library
    (Shah1982.mo, Eqs. 1-14). ThermoCycle smooths the switches between regimes with
    transition functions; they are NOT reproduced here: sharp switches as in Shah (1982)
    and in the ACHP implementation (Correlations.ShahEvaporation_Average).

    Arguments
    ---------
    AS : CoolProp AbstractState
        Fluid state object (will be moved by saturation property calls).
    geom : dict
        Geometry dictionary with keys:
        - 'Tube_ID' [m] (or 'Tube_OD' [m] and 'Tube_t' [m])
        - 'theta' [deg] (optional, tube inclination from the horizontal, default 0.0)
    cond : dict
        Operating conditions dictionary with keys:
        - 'p' [Pa]       saturation pressure
        - 'x' [-]        vapour quality (clamped to ]0, 1[)
        - 'q' [W/m2]     local wall heat flux (required, must be >= 0)
        - 'G' [kg/m2/s]  mass flux in ONE tube (optional if 'm_dot' [kg/s] and geom are provided)

    Returns
    -------
    h : float
        Saturated flow boiling heat transfer coefficient [W/m2/K].

    Reference
    ---------
    Shah, M.M. (1982) Chart correlation for saturated boiling heat transfer: equations and further
    study. ASHRAE Transactions 88(1), 185-196.

    Validity range
    --------------
    Validity range: not given by the reference implementations (TODO from Shah 1982);
    saturated boiling, water and refrigerants (Shah's database).
    Inclined tubes (theta not 0 or 90 deg): not covered by Shah, treated as horizontal (warning).

    Definitions
    -----------
    - Re_l: liquid alone, G * (1 - x) * D / mu_l [-], with properties saturated at p
      (T_sat = mean of bubble and dew temperatures for mixtures, as _sat_props).
    - Bo: boiling number with the wall heat flux q and h_lv [-].
    - Fr_l: liquid Froude number of the whole flow as liquid [-].
    - q: wall heat flux [W/m2].
    - Returns the local h [W/m2/K] at quality x.
    """
    # Check required inputs
    if 'q' not in cond:
        raise KeyError("Shah 1982 boiling HTC: cond needs 'q' [W/m2] (local wall heat flux).")

    q = cond['q']
    if q < 0:
        raise ValueError("Shah 1982 is a boiling correlation: q must be >= 0.")

    # 1. Geometry, mass flux, pressure, and quality
    G = _tube_G(geom, cond)                     # mass flux in one tube [kg/m2/s]
    D = _tube_ID(geom)                          # tube inner diameter [m]
    p = cond['p']                               # saturation pressure [Pa]
    x = _clamp_x(cond['x'])                     # vapor quality [-]
    theta = geom.get('theta', 0.0)              # tube inclination [deg]

    # 2. Saturated properties at pressure p
    s = _sat_props(AS, p, cond)
    rho_l = s['rho_l']                          # saturated liquid density [kg/m3]
    rho_v = s['rho_v']                          # saturated vapor density [kg/m3]
    mu_l = s['mu_l']                            # saturated liquid dynamic viscosity [Pa s]
    k_l = s['k_l']                              # saturated liquid thermal conductivity [W/m/K]
    Pr_l = s['Pr_l']                            # saturated liquid Prandtl number [-]
    h_lv = s['h_lv']                            # latent heat of vaporization [J/kg]

    # 3. Liquid-alone heat transfer coefficient, Dittus-Boelter (Eq. 5)
    Re_l = G * (1.0 - x) * D / mu_l             # liquid-alone Reynolds number [-]
    h_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D   # liquid-alone HTC [W/m2/K]

    # 4. Convection number (Eq. 2)
    Co = ((1.0 - x) / x)**0.8 * (rho_v / rho_l)**0.5   # convection number [-]

    # 5. Boiling number (Eq. 3)
    Bo = q / (G * h_lv)                         # boiling number [-]

    # 6. Liquid Froude number (Eq. 4)
    g = 9.81                                    # gravitational acceleration [m/s2]
    Fr_l = G**2 / (rho_l**2 * g * D)            # liquid Froude number [-]

    # 7. Parameter N, depends on the tube orientation
    if abs(theta) != 0.0 and abs(theta) != 90.0:
        warnings.warn(
            f"Shah (1982) gives horizontal and vertical tubes only: inclined tube "
            f"(theta = {theta} deg) treated as horizontal."
        )
    if abs(theta) == 90.0:
        N = Co                                  # vertical tube, all Fr_l (Eq. 13)
    elif Fr_l >= 0.04:
        N = Co                                  # horizontal tube, Fr_l >= 0.04 (Eq. 13)
    else:
        N = 0.38 * Fr_l**-0.3 * Co              # horizontal tube, Fr_l < 0.04 (Eq. 14)

    # 8. Constant F (Eqs. 11-12)
    if Bo >= 11e-4:
        F = 14.70                               # Eq. 11
    else:
        F = 15.43                               # Eq. 12

    # 9. Convective boiling enhancement factor psi_cb (Eq. 8)
    psi_cb = 1.8 / N**0.8                       # [-]

    # 10. Flow boiling regimes and two-phase multiplier psi
    if N > 1.0:
        # Nucleate boiling regime
        if Bo > 0.3e-4:
            psi_nb = 230.0 * Bo**0.5            # Eq. 6
        else:
            psi_nb = 1.0 + 46.0 * Bo**0.5       # Eq. 7
        psi = max(psi_nb, psi_cb)
    elif 0.1 < N <= 1.0:
        # Bubble suppression regime
        psi_bs = F * Bo**0.5 * np.exp(2.74 * N**-0.1)   # Eq. 9
        psi = max(psi_bs, psi_cb)
    else:
        # Convective boiling / low convection number regime (N <= 0.1)
        psi_bs = F * Bo**0.5 * np.exp(2.47 * N**-0.15)  # Eq. 10
        psi = max(psi_bs, psi_cb)

    # 11. Two-phase heat transfer coefficient (Eq. 1: psi = h_TP / h_l)
    h = psi * h_l                               # [W/m2/K]

    # TODO: validity range from Shah (1982).

    return h


#%% EXTERNAL BOILING (outside the tubes): unchanged, T-based signatures

def pool_boiling(fluid, T_sat, T_tube):
    """
    fluid : fluid name, T_sat : external fluid saturation temperature [K], T_tube : tube temperature [K]
    Returns the pool boiling HEAT FLUX q [W/m^2] (not an HTC: h = q / (T_tube - T_sat)).
    Reference: Van Long Le - Heat Pipe Code; VDI Heat Atlas.
    """
    g = 9.81
    rho_l = PropsSI('D', 'T', T_sat, 'Q', 0, fluid)
    rho_v = PropsSI('D', 'T', T_sat, 'Q', 1, fluid)
    mu_l = PropsSI('V', 'T', T_sat, 'Q', 0, fluid)
    cp_l = PropsSI('C', 'T', T_sat, 'Q', 0, fluid)
    sigma = PropsSI('I', 'T', T_sat, 'Q', 0.5, fluid)
    Dh_evap = PropsSI('H', 'T', T_sat, 'Q', 1, fluid) - PropsSI('H', 'T', T_sat, 'Q', 0, fluid)

    C_nb = 0.013       # pentane on polished nickel
    n = 1.7            # Prandtl exponent, fluids other than water
    DT_e = T_tube - T_sat
    Pr_l = PropsSI('Prandtl', 'T', T_sat, 'Q', 0, fluid)

    coef_1 = mu_l * Dh_evap
    coef_2 = g * (rho_l - rho_v) / sigma
    coef_3 = (cp_l * DT_e) / (C_nb * Dh_evap * Pr_l**n)
    return coef_1 * coef_2**(1 / 2) * coef_3**3


def film_boiling(D_out, fluid, T_sat, T_tube, P_sat):
    """Film boiling on a horizontal cylinder [W/(m^2 K)]. Reference: Incropera."""
    g = 9.81
    e = 0.9
    sigma = 5.67e-8
    T_f = (T_sat + T_tube) / 2
    (k_v, rho_v, nu_v, cp_v) = PropsSI(('L', 'D', 'V', 'C'), 'T', T_f, 'P', P_sat, fluid)
    h_v = PropsSI('H', 'T', T_sat, 'Q', 1, fluid)
    (rho_l, h_l) = PropsSI(('D', 'H'), 'T', T_sat, 'Q', 0, fluid)
    C = 0.62

    Dh_evap = h_v - h_l
    Dh_evap_corrected = Dh_evap + 0.8 * cp_v * (T_tube - T_sat)

    coef = (g * (rho_l - rho_v) * Dh_evap_corrected * D_out**3) / (nu_v * k_v * (T_tube - T_sat))
    h_conv = (D_out / k_v) * C * coef**(1 / 4)
    h_rad = e * sigma * (T_tube**4 - T_sat**4) / (T_tube - T_sat)

    def equation(x):
        return x**(4 / 3) - (h_conv**(4 / 3) + h_rad * x**(1 / 3))

    return fsolve(equation, x0=[100])


def boiling_curve(D_out, fluid, T_sat, P_sat):
    """Boiling curve of a horizontal cylinder, wall superheat 0.01-1000 K. Reference: Incropera.
    Returns (h_final, DT) [W/(m^2 K)], [K]."""
    DT = np.linspace(1, 1000, 1000)
    q_pool = np.zeros(len(DT))
    q_film = np.zeros(len(DT))
    q_final = np.zeros(len(DT))

    g = 9.81
    C_crit = 0.15
    C_mhf = 0.09

    rho_l = PropsSI('D', 'T', T_sat, 'Q', 0, fluid)
    rho_v = PropsSI('D', 'T', T_sat, 'Q', 1, fluid)
    sigma = PropsSI('I', 'T', T_sat, 'Q', 0, fluid)
    Dh_evap = PropsSI('H', 'T', T_sat, 'Q', 1, fluid) - PropsSI('H', 'T', T_sat, 'Q', 0, fluid)

    q_crit = C_crit * Dh_evap * rho_v * ((sigma * g * (rho_l - rho_v)) / rho_v**2)**(1 / 4)
    q_min = C_mhf * rho_v * Dh_evap * ((g * sigma * (rho_l - rho_v)) / (rho_l + rho_v)**2)**(1 / 4)

    for i in range(len(DT)):
        T_surf = DT[i] + T_sat
        q_pool[i] = pool_boiling(fluid, T_sat, T_surf)
        h_film = film_boiling(D_out, fluid, T_sat, T_surf, P_sat)
        q_film[i] = h_film * DT[i]

    DT_crit = np.argmin(np.abs(q_pool - q_crit))
    DT_min = np.argmin(np.abs(q_film - q_min))
    Delta = (q_min - q_crit) / (DT_min - DT_crit)

    for i in range(len(DT)):
        if DT[i] <= DT_crit:
            q_final[i] = q_pool[i]
        elif DT[i] > DT_min:
            q_final[i] = q_film[i]
        else:
            q_final[i] = q_crit + Delta * (DT[i] - DT_crit)

    h_final = q_final / DT
    DT = np.concatenate((np.array([0]), DT), axis=0)
    h_final = np.concatenate((np.array([0]), h_final), axis=0)
    return h_final, DT