
import CoolProp.CoolProp as CP
import numpy as np

def pseudo_critical_temperature(P, fluid):
    """
    Fluid-general pseudo-critical temperature: the temperature of maximum
    isobaric specific heat (cp) along the isobar at P, which is the
    defining property of the Widom line / pseudo-critical point.

    Self-contained: creates its own CoolProp AbstractState internally, so
    no AbstractState is needed from the caller.

    Inputs
    ------
    P     : Pressure [Pa]
    fluid : CoolProp fluid name (e.g. 'CO2')

    Outputs
    -------
    T_pc : Pseudo-critical temperature [K]
    """
    AS = CP.AbstractState('HEOS', fluid)

    T_crit = AS.keyed_output(CP.iT_critical)
    T_search = np.linspace(T_crit * 0.9, T_crit * 1.3, 500)

    cp_search = []
    for T in T_search:
        AS.update(CP.PT_INPUTS, P, T)
        cp_search.append(AS.cpmass())

    T_pc = T_search[np.argmax(cp_search)]

    return T_pc


def pseudo_boiling_bounds(P, fluid, T_pc_func=None):
    """
    Pseudo-critical temperature and the pseudo-boiling transition bounds
    (T_minus, T_plus) for a supercritical fluid, based on Banuti's linear
    pseudo-boiling model.

    Self-contained: creates its own CoolProp AbstractState internally, so
    no AbstractState is needed from the caller.

    Inputs
    ------
    P         : Pressure [Pa]
    fluid     : CoolProp fluid name (e.g. 'CO2')
    T_pc_func : Optional. A function(P) -> T_pc [K] to use instead of the
                default fluid-general cp-maximum search. Pass a
                fluid-specific empirical correlation here (e.g.
                pseudo_critical_temperature_CO2) for speed or accuracy
                on a fluid it was fitted to.

    Outputs
    -------
    T_pc    : Pseudo-critical temperature [K]
    T_minus : Lower (liquid-like) transition bound [K]
    T_plus  : Upper (vapor-like) transition bound [K]

    Reference
    ---------
    Banuti, D. (2019). "The Latent Heat of Supercritical Fluids."
    Periodica Polytechnica Chemical Engineering, 63(2), 270-275.
    """
    AS = CP.AbstractState('HEOS', fluid)

    if T_pc_func is None:
        T_pc = pseudo_critical_temperature(P, fluid)
    else:
        T_pc = T_pc_func(P)

    AS.update(CP.PT_INPUTS, P, T_pc)
    cp_pc = AS.cpmass()
    h_0_pc = AS.hmass()

    # --- Liquid-like reference state (low-temperature anchor) ---
    P_crit = AS.keyed_output(CP.iP_critical)
    P_triple = AS.keyed_output(CP.iP_triple)

    P_L = P_crit * 0.1
    if P_L < P_triple:
        P_L = P_triple * 1.05

    AS.update(CP.PQ_INPUTS, P_L, 0)
    T_L = AS.T()

    SC = 0.01  # K, subcooling offset to stay on the liquid branch

    AS.update(CP.PT_INPUTS, P_L, T_L - SC)
    cp_L = AS.cpmass()
    h_L = AS.hmass()

    h_L_0 = h_L - cp_L * T_L
    # --------------------------------------------------------------

    # Ideal gas cp approximation (vapor-like reference)
    T_max = AS.keyed_output(CP.iT_max)
    AS.update(CP.PT_INPUTS, P, T_max)
    cp_IG = AS.cpmass()

    T_plus = (h_0_pc - cp_pc * T_pc) / (cp_IG - cp_pc)
    T_minus = (h_L_0 - h_0_pc + cp_pc * T_pc) / (cp_pc - cp_L)

    return T_pc, T_minus, T_plus

def pseudo_critical_temperature_CO2(P):
    """
    Pseudo-critical temperature of CO2 as a function of pressure.

    Assumptions
    -----------
    P in [7.5,14] MPa

    Inputs
    ------
    P : Pressure [Pa]

    Outputs
    -------
    T_pc : Pseudo-critical temperature [K]
    
    Reference
    ---------
    Investigation on the Properties of Supercritical CO2 Fluid and its Heat
    Transfer Characteristics (2012)
    
    Z. Yang & J. Yang
    """
    
    p = P*1e-6 # MPa
    T_pc = - 31.40 + 12.15*p - 0.6927*p**2 + 0.03160 * p**3 - 0.0007521 * p**4
    
    return T_pc + 273.15