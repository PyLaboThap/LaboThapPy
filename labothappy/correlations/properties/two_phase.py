
import CoolProp.CoolProp as CP

EPS = 1e-12


#=======================================================================
# SHARED TWO-PHASE PROPERTY EXTRACTION
#=======================================================================

def get_saturated_phase_properties(AS):
    """
    Extract saturated liquid/vapor properties and quality from a two-phase
    CoolProp AbstractState.

    Centralizes the AS -> (x, rho_l, rho_v, mu_l, mu_v, sigma) derivation so
    that callers outside this module (e.g. pressure_drop correlations) stay
    in sync with what `compute_void_fraction` uses internally.

    Parameters
    ----------
    AS : CoolProp.AbstractState
        Two-phase fluid state (quality and pressure already set).

    Returns
    -------
    dict with keys: x, P, rho_l, rho_v, mu_l, mu_v, sigma
    """
    x = AS.Q()
    x = max(EPS, min(1.0 - EPS, x))

    P = AS.p()

    AS_l = CP.AbstractState(AS.backend_name(), AS.fluid_names()[0])
    AS_l.update(CP.PQ_INPUTS, P, 0)

    AS_v = CP.AbstractState(AS.backend_name(), AS.fluid_names()[0])
    AS_v.update(CP.PQ_INPUTS, P, 1)

    return {
        "x": x,
        "P": P,
        "rho_l": AS_l.rhomass(),
        "rho_v": AS_v.rhomass(),
        "mu_l": AS_l.viscosity(),
        "mu_v": AS_v.viscosity(),
        "sigma": AS.surface_tension(),
    }

def compute_two_phase_density(x, rho_l, rho_g, alpha=None):
    """
    Compute two-phase density.
    
    If alpha is provided: ρ_mix = 1 / (α/ρ_g + (1-α)/ρ_l)
    Otherwise: use homogeneous model
    
    Parameters
    ----------
    x : float
        Vapor quality [-]
    rho_l : float
        Liquid density [kg/m³]
    rho_g : float
        Vapor density [kg/m³]
    alpha : float, optional
        Void fraction [-]. If None, computed from homogeneous model.
    
    Returns
    -------
    float
        Mixture density [kg/m³]
    """
    x = max(EPS, min(1.0 - EPS, x))

    if alpha is None:
        # Deferred import: void_fraction.py imports get_saturated_phase_properties
        # from this module, so importing at module level here would be circular.
        from labothappy.correlations.void_fraction.void_fraction import void_fraction_homogeneous
        alpha = void_fraction_homogeneous(rho_l, rho_g, x)
    
    # ρ_mix = 1 / (α/ρ_g + (1-α)/ρ_l)
    return 1.0 / (alpha / rho_g + (1.0 - alpha) / rho_l)