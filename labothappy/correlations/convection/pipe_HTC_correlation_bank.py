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
    try:
        AS.update(CP.PQ_INPUTS, p, 0.5)
        return AS.surface_tension()
    except Exception:
        return PropsSI('I', 'P', p, 'Q', 0.5, AS.fluid_names()[0])


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
        Nu_1 = 0.0677 * (Re * Pr * (Dh / L)**1.33)
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
    Validity: D 4-10.7 mm, T_in 25-67 C, P 7.5-8.5 MPa, G 74-200 kg/m2/s, CO2.
    cond : G, p, T_w, k, rho, mu, cp  (bulk values)      geom : Tube_ID
    Reference: Liu, He, Yang, Fei, Experimental study on heat transfer and pressure drop of
    supercritical CO2 cooled in a large tube, 2014."""
    D = _tube_ID(geom)
    G, rho, cp = _tube_G(geom, cond), cond['rho'], cond['cp']
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
    tube based on the Froude number analysis (2024); Banuti (2015, 2019) for T_minus / T_plus;
    Yang & Yang (2012) for the pseudo-critical temperature."""
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
    T_in = PropsSI('T', 'H', h_in, 'P', P, fluid)
    dT = T_w - T_in
    if abs(dT) < 1e-6:
        dT = 1e-6
    cp_ave = (h_w - h_in) / dT
    Pr_ave = mu * cp_ave / k

    if x_ave < 0 or x_ave > 1:                       # vapor-like / liquid-like
        return 0.019 * (Re**0.79) * (Pr_ave**0.51) * k / D_in

    g = 9.81
    Fr_VL_w = G**2 * x_pb_in / (rho_w**2 * g * D_in)             # two-phase-like
    Nu_top = 0.127 * (Re**1.065) * (Pr_ave**1.053) * (SBO_w**0.66) * (Fr_VL_w**(-0.012))
    Nu_bot = 1.512 * (Re**0.878) * (Pr_ave**0.685) * (SBO_w**0.69) * (Fr_VL_w**(-0.0213))
    return 0.5 * (Nu_top + Nu_bot) * k / D_in

#%% CONDENSATION INSIDE TUBES

def pipe_thome_condensing_HTC(AS, geom, cond):
    """Thome / El Hajal (stratified-annular, Biberg void-fraction geometry) in-tube condensation.
    cond : G, p, x, T_w (wall temperature [K])      geom : Tube_ID
    Void fraction: labothappy.correlations.void_fraction.void_fraction_old.void_fraction.
    Reference: Thome, El Hajal, Cavallini, Condensation in horizontal tubes (2003)."""
    from labothappy.correlations.void_fraction.void_fraction_old import void_fraction
    g = 9.81
    D_i = _tube_ID(geom)
    G, P_sat, x = _tube_G(geom, cond), cond['p'], min(max(cond['x'], 0.0), 1.0)
    if 'T_w' not in cond:
        raise KeyError("Thome condensation: cond needs 'T_w'.")

    s = _sat_props(AS, P_sat, cond)
    _need_sigma(AS, s)
    rho_L, rho_G, mu_L, lambda_L, Pr_L = s['rho_l'], s['rho_v'], s['mu_l'], s['k_l'], s['Pr_l']
    sigma, h_LG, T_sat = s['sigma'], s['h_lv'], s['T_sat']
    dT = max(T_sat - cond['T_w'], 1e-3)

    r = D_i / 2
    A = np.pi * r**2
    eps = void_fraction(x, rho_G, rho_L)[0]
    A_L = (1 - eps) * A

    u_L = G * (1 - x) / (rho_L * (1 - eps)) if eps != 1 else G / rho_L
    u_G = G * x / (rho_G * eps) if eps != 0 else G / rho_G

    # Biberg expression for the stratified (dry) angle
    theta = 2 * np.pi - 2 * (np.pi * (1 - eps) + (1.5 * np.pi)**(1 / 3) * (1 - 2 * (1 - eps) + (1 - eps)**(1 / 3) - eps**(1 / 3))
                             - 0.005 * (1 - eps) * eps * (1 - 2 * (1 - eps)) * (1 + 4 * ((1 - eps)**2 + eps**2)))

    # Liquid film thickness
    if (2 * np.pi - theta) / (2 * np.pi) >= 1e-3:
        delta = (D_i - np.sqrt(D_i**2 - A_L / ((2 * np.pi - theta) / 8))) / 2
    else:
        delta = D_i

    c, m, n = 0.003, 0.5, 0.74
    Re_L = (4 * G * (1 - x) * delta) / ((1 - eps) * mu_L) if eps != 1 else 0.0
    f_i = 1 + (u_G / u_L)**0.5 * ((rho_L - rho_G) * g * delta**2 / sigma)**0.25
    h_c = c * (Re_L**n) * (Pr_L**m) * (lambda_L / delta) * f_i if delta != 0 else 0.0

    h_f = 0.728 * ((rho_L * (rho_L - rho_G) * g * h_LG * lambda_L**3) / (mu_L * D_i * dT))**0.25

    return (h_f * r * theta + (2 * np.pi - theta) * r * h_c) / (2 * np.pi * r)


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
    "CarbonDioxide": 5500, "CO2": 5500,
    "n-Propane": 4000, "Water": 25580, "H2O": 25580, "Ammonia": 7000, "NH3": 7000,   # added (from memory of the VDI table): to check
}  # fluids missing here (R1233zd(E), ...): give cond['a_0']


def _a0(fluid, cond):
    if 'a_0' in cond:
        return cond['a_0']
    if fluid in _A0_VDI:
        return _A0_VDI[fluid]
    raise KeyError(f"No reference nucleate boiling HTC a_0 for {fluid!r}: give cond['a_0'].")


def pipe_steiner_taborek_boiling_HTC(AS, geom, cond):
    """Steiner-Taborek asymptotic model (Wolverine Engineering Data Book, Boiling inside plain tubes).
    cond : G, p, x, q (local heat flux); h_nb_o (optional, reference nucleate boiling HTC at q_o = 20 kW/m2;
           default: VDI a_0 of the fluid, or the cyclopentane estimate of the original code)
    geom : Tube_ID, Tube_L
    Convective part: liquid-only (total flow) Gnielinski/Sieder-Tate HTC, no wall viscosity correction."""
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
        h_v = _gnielinski_core(mu_v, Pr_v, mu_v, k_v, G, D_in, L)
        a = (1 - x)**1.5 + 1.9 * x**0.6 * (1 - x)**0.1 * (rho_l / rho_v)**0.35
        b = (h_v / h_L) * x**0.01 * (1 + 8 * (1 - x)**0.7) * (rho_l / rho_v)**0.67
        F_tp = 1 if (a == 0 or b == 0) else (a**(-2.2) + b**(-2))**(-0.5)
        return F_tp * h_v

    # Nucleate + convective (asymptotic, exponent 3)
    M = AS.molar_mass() * 1e3                                        # g/mol
    q_o = 20000
    if 'h_nb_o' in cond:
        h_nb_o = cond['h_nb_o']
    elif fluid == 'Cyclopentane':
        h_nb_o = 3070 * (2420 / 2840)    # assumption: h_nb,o(cyclopentane) = h_nb,o(cyclohexane) * (n-pentane / n-hexane)
    else:
        h_nb_o = _a0(fluid, cond)
    D_in_o = 0.01

    F_tp = ((1 - x)**1.5 + 1.9 * x**0.6 * (rho_l / rho_v)**0.35)**1.1
    F_pf = 2.816 * P_r**0.45 + (3.4 + (1.7 / (1 - P_r**7))) * P_r**3.7
    nf = 0.8 - 0.1 * np.e**(1.75 * P_r)
    F_M = 0.377 + 0.199 * np.log(M) + 0.000028427 * M**2
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
    elif fluid in ('H2O', 'NH3', 'CO2', 'SF6', 'Water', 'Ammonia'):
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
    Reference: Gungor & Winterton, A general correlation for flow boiling in tubes and annuli (1986)."""
    D_in = _tube_ID(geom)
    G, P_sat, q = _tube_G(geom, cond), cond['p'], cond['q']
    x = _clamp_x(cond['x'])

    s = _sat_props(AS, P_sat, cond)
    rho_l, rho_v, mu_l, mu_v, k_l, Pr_l, h_lv = s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v'], s['k_l'], s['Pr_l'], s['h_lv']
    P_crit, MM = AS.p_critical(), AS.molar_mass()
    P_r = P_sat / P_crit

    # two-phase viscosity at (p, x): CoolProp when defined, McAdams otherwise
    try:
        AS.update(CP.PQ_INPUTS, P_sat, x)
        mu_tp = AS.viscosity()
    except Exception:
        mu_tp = 1.0 / (x / mu_v + (1 - x) / mu_l)

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
    Re_tp = G * D_in / mu_tp                                  # (alternative: Re_l with 2.53e-6)
    S = S_2 / (1 + 0.00000115 * E**2 * Re_tp**1.17)

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