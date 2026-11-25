import numpy as np
import CoolProp.CoolProp as CP

EPS = 1e-12


#=======================================================================
# SHARED TWO-PHASE PROPERTY EXTRACTION
#=======================================================================

# HEOS states used for the surface tension: one per fluid, created once and kept
_HEOS_FOR_SIGMA = {}


def _get_surface_tension(AS, P, x):
    """Surface tension [N/m] at saturation pressure P, always from a HEOS state of the fluid
    (AS itself if its backend is HEOS). The tabular backends (BICUBIC&HEOS) are not used: with
    CoolProp 7.2 their surface_tension() raises or returns the value of the first state for every
    later state. nan if CoolProp has no surface tension for the fluid."""
    try:
        if AS.backend_name() == 'HelmholtzEOSBackend':
            return AS.surface_tension()
        fluid = AS.fluid_names()[0]
        if fluid not in _HEOS_FOR_SIGMA:
            _HEOS_FOR_SIGMA[fluid] = CP.AbstractState("HEOS", fluid)
        AS_s = _HEOS_FOR_SIGMA[fluid]
        AS_s.update(CP.PQ_INPUTS, P, x)
        return AS_s.surface_tension()
    except ValueError:
        return np.nan   # not available in CoolProp for this fluid

def _get_viscosity(AS, P, x):
    try:
        return AS.viscosity()
    except ValueError:
        pass
    try:
        AS_s = CP.AbstractState("HEOS", AS.fluid_names()[0])
        AS_s.update(CP.PQ_INPUTS, P, x)
        return AS_s.viscosity()
    except ValueError:
        return np.nan   # not available in CoolProp for this fluid

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
    sigma = _get_surface_tension(AS, P, x)
    mu_l = _get_viscosity(AS_l, P, x)
    mu_v = _get_viscosity(AS_v, P, x)

    return {
        "x": x,
        "P": P,
        "rho_l": AS_l.rhomass(),
        "rho_v": AS_v.rhomass(),
        "mu_l": mu_l,
        "mu_v": mu_v,
        "sigma": sigma,
    }

def compute_two_phase_density(x, rho_l, rho_g, alpha=None):
    """
    Compute two-phase density.
    
    In-situ (void-fraction weighted) density: ρ_mix = α·ρ_g + (1-α)·ρ_l
    If alpha is not provided, the homogeneous void fraction is used, which
    gives ρ_mix = 1 / (x/ρ_g + (1-x)/ρ_l).
    
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
    
    # ρ_mix = α·ρ_g + (1-α)·ρ_l  (volume fractions weight the densities)
    return alpha * rho_g + (1.0 - alpha) * rho_l