"""
Single entry point for all pressure drop correlations.

    dp = pressure_drop(AS, geom, G, correlation, **extra)

AS          : CoolProp AbstractState, already set at the flow state
geom        : equivalent channel geometry, built by the component model
            'D_h' hydraulic diameter [m], 'L' flow length [m], 'K', roughness [m] (optional)
            + correlation-specific keys (chevron-angle, ...)
G           : Mass flux [kg/(m^2 s)]
correlation : name, see available_correlations()
**extra     : only for correlations that need more arguments

Returns the frictional pressure drop [Pa] over geom['L].
"""


from labothappy.correlations.pressure_drop.pipe_DP import (
    friction_factor_churchill, friction_factor_swamee_jain, friction_factor_haaland,
    friction_factor_konakov, friction_factor_petukhov, friction_factor_cheng_CO2,
    pressure_drop_friedel, pressure_drop_muller_steinhagen_heck
)
from labothappy.correlations.pressure_drop.plate_DP import han_BPHEX_DP
from labothappy.correlations.pressure_drop.shell_and_tube_DP import (
    shell_DP_kern, shell_bell_delaware_DP)
from labothappy.correlations.properties.two_phase import get_saturated_phase_properties

# =============================================================================
# 1) SINGLE PHASE - friction factor correlations (all share Darcy-Weisbach)
# =============================================================================

def _darcy_weisbach(friction_factor):
    """Turn a friction factor f(AS, geom, G, Re) into a pressure drop function."""
    def dp(AS, geom, G, **extra):
        rho, mu = AS.rhomass(), AS.viscosity()
        D, L = geom['D'], geom['L']
        Re = G * D / mu
        f = friction_factor(AS, geom, G, Re)
        return f * (L / D) * G**2 / (2 * rho)
    return dp

FRICTION_FACTORS = {
    "Churchill":   lambda geom, Re: friction_factor_churchill(geom.get('K', 0.0), geom['D'], Re),
    "Swamee-Jain": lambda geom, Re: friction_factor_swamee_jain(geom.get('K', 0.0), geom['D'], Re),
    "Haaland":     lambda geom, Re: friction_factor_haaland(geom.get('K', 0.0), geom['D'], Re),
    "Konakov":     lambda Re: friction_factor_konakov(Re),
    "Petukhov":    lambda Re: friction_factor_petukhov(Re),
    "Cheng-CO2":   lambda AS, geom, G: friction_factor_cheng_CO2(G, geom['D'], AS.p(), AS.hmass(), AS.viscosity()),
}


# =============================================================================
# 2) SINGLE PHASE - correlations giving the pressure drop directly
# =============================================================================

def _dp_shell_kern(AS, geom, G, m_dot, T_wall, **extra):
    return shell_DP_kern(m_dot, T_wall, AS.hmass(), AS.p(), AS, geom)

def _dp_shell_bell_delaware(AS, geom, G, m_dot, **extra):
    return shell_bell_delaware_DP(m_dot, AS.hmass(), AS.p(), AS, geom)


# =============================================================================
# 3) TWO PHASE
# =============================================================================

def _dp_friedel(AS, geom, G, **extra):
    s = get_saturated_phase_properties(AS)
    return pressure_drop_friedel(G, s['x'], s['rho_l'], s['rho_v'], s['mu_l'], s['mu_v'],
                                 s['sigma'], geom['D'], geom['L'], K=geom.get('K', 0.0))

def _dp_msh(AS, geom, G, **extra):
    s = get_saturated_phase_properties(AS)
    return pressure_drop_muller_steinhagen_heck(G, s['x'], s['rho_l'], s['rho_v'], s['mu_l'],
                                                s['mu_v'], geom['D'], geom['L'], K=geom.get('K', 0.0))

def _dp_han_bphex(AS, geom, G, m_dot, **extra):
    s = get_saturated_phase_properties(AS)
    return han_BPHEX_DP(s['mu_l'], G, geom['D'], geom['chevron_angle'], geom['pitch_co'],
                        s['rho_v'], s['rho_l'], geom['L'], geom['n_canals'], m_dot, geom['D_port'])

# _dp_choi(AS, geom, G, AS_ex, **extra): same idea, AS_ex passed as an extra


# =============================================================================
# REGISTRY + DISTRIBUTOR
# =============================================================================

CORRELATIONS = {
    # name                 : (function,                      phase)
    **{name: (_darcy_weisbach(f), "single-phase") for name, f in FRICTION_FACTORS.items()},
    "Shell_Kern":           (_dp_shell_kern,                 "single-phase"),
    "Shell_Bell_Delaware":  (_dp_shell_bell_delaware,        "single-phase"),
    "Friedel":              (_dp_friedel,                    "two-phase"),
    "MSH":                  (_dp_msh,                        "two-phase"),
    "Han_BPHEX":            (_dp_han_bphex,                  "two-phase"),
}


def available_correlations(phase=None):
    """Names of the available correlations, optionally filtered by phase ('single-phase' or 'two-phase')."""
    return [n for n, (_, ph) in CORRELATIONS.items() if phase is None or ph == phase]


def pressure_drop(AS, geom, G, correlation, **extra):
    """Frictional pressure drop [Pa] over geom['L']."""
    try:
        func, _ = CORRELATIONS[correlation]
    except KeyError:
        raise ValueError(f"Unknown pressure drop correlation {correlation!r}. "
                         f"Available: {available_correlations()}") from None
    return func(AS, geom, G, **extra)
