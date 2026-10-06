"""
Uniform calling interface for the heat transfer coefficient (HTC) and pressure
drop (DP) correlations used by the moving-boundary heat exchanger model
(``hex_MB_charge_sensitive_bis``).

Every correlation is called the same way as the pipe correlations
(``pressure_drop_pipe_single_phase(AS, pipe_geom, G, correlation=...)``):

    alpha = heat_transfer_coefficient(AS, geom, G, correlation="Gnielinski", cond=cond)
    dP    = pressure_drop(AS, geom, G, correlation="Churchill", cond=cond)

Arguments
---------
AS : CoolProp.AbstractState
    Fluid state of the cell, already set by the caller at the cell mean
    state (h, p). Adapters may leave AS in another state on return, so the
    caller must re-set it before reusing it.
geom : dict
    Geometry of ONE side of the heat exchanger. It always contains the
    side-neutral keys below, plus every heat exchanger parameter:
        'D'        hydraulic diameter of the side [m]
        'L'        total flow length of the side [m]
        'K'        absolute roughness [m]
        'theta'    inclination [deg]
        'A'        heat transfer area of the side [m^2]
        'n_canals' number of channels (plates only)
        'canal_t'  channel thickness (plates only)
G : float
    Mass flux of the side [kg/(m^2 s)].
correlation : str
    Name of the correlation (see HTC_CORRELATIONS / DP_CORRELATIONS).
cond : CellConditions
    Operating data that AS alone cannot provide (wall temperature, quality,
    heat flux, flow rate, ...).

Return
------
alpha [W/(m^2 K)] for heat transfer, dP [Pa] over the FULL flow length of the
side for pressure drop (the heat exchanger model scales it to each cell).

Transition note
---------------
The existing correlations do not have this signature yet. Each entry of the
registries below is a thin adapter translating the uniform call into the
existing signature. When a correlation is rewritten with the uniform
signature, its adapter can simply be replaced by the function itself.
The existing correlation files are NOT modified by this module.

author: moving-boundary model cleaning (hex_MB_charge_sensitive_bis)
"""

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d

# Existing correlations (called unchanged through the adapters below)
from labothappy.correlations.convection.pipe_htc import (
    gnielinski_pipe_htc, boiling_curve, horizontal_flow_boiling,
    flow_boiling_gungor_winterton, Liu_sCO2, Cheng_sCO2, thome_condensation,
    choi_boiling, Meshram)
from labothappy.correlations.convection.plate_htc import (
    water_plate_HTC, martin_BPHEX_HTC, muley_manglik_BPHEX_HTC,
    han_boiling_BPHEX_HTC, han_cond_BPHEX_HTC, thonon_plate_HTC,
    kumar_plate_HTC, martin_holger_plate_HTC, amalfi_plate_HTC,
    shah_condensation_plate_HTC)
from labothappy.correlations.convection.shell_and_tube_htc import (
    shell_bell_delaware_htc, shell_htc_kern)
from labothappy.correlations.convection.tube_bank_htc import ext_tube_film_condens
from labothappy.correlations.convection.fins_htc import htc_tube_and_fins
from labothappy.correlations.convection.printed_circuit_htc import PCHE_Lee, PCHE_conv
from labothappy.correlations.pressure_drop.pipe_DP import (
    pressure_drop_pipe_single_phase, pressure_drop_pipe_frictional_two_phase)
from labothappy.correlations.pressure_drop.shell_and_tube_DP import (
    shell_DP_kern, shell_DP_donohue, shell_bell_delaware_DP)
from labothappy.correlations.pressure_drop.fins_DP import DP_tube_and_fins
from labothappy.correlations.properties.thermal_conductivity import conducticity_R1233zd

X_MIN = 1e-4   # quality clipping used for two-phase correlations
X_MAX = 1.0 - 1e-4


# =============================================================================
# CELL OPERATING CONDITIONS
# =============================================================================

@dataclass
class CellConditions:
    """Operating data of one side of one cell (everything AS does not hold)."""
    fluid: str
    m_dot: float            # mass flow rate of the side, whole heat exchanger [kg/s]
    p: float                # mean pressure of the cell [Pa]
    h: float                # mean specific enthalpy of the cell [J/kg]
    T: float                # mean bulk temperature of the cell [K]
    T_wall: float           # wall temperature estimate [K]
    x: float = float("nan")         # mean vapour quality (two-phase cells) [-]
    T_sat: float = float("nan")     # mean saturation temperature [K]
    h_in: float = float("nan")      # enthalpy at the cell boundary of lower enthalpy [J/kg]
    h_out: float = float("nan")     # enthalpy at the cell boundary of higher enthalpy [J/kg]
    p_su: float = float("nan")      # side inlet pressure [Pa]
    q: float = float("nan")         # heat flux estimate [W/m^2]
    Q: float = float("nan")         # cell heat rate [W]
    DT_lm: float = float("nan")     # F * LMTD of the cell [K]
    alpha_other: float = float("nan")  # HTC of the other side [W/(m^2 K)]


# =============================================================================
# PROPERTY HELPERS (shared by the adapters)
# =============================================================================

def _conductivity(AS, cond, T=None, p=None):
    """Thermal conductivity, with the R1233zd(E) fit CoolProp lacks."""
    if cond.fluid == 'R1233zd(E)':
        return conducticity_R1233zd(cond.T if T is None else T, cond.p if p is None else p)
    return AS.conductivity()


def bulk_properties(AS, cond):
    """mu, k, Pr, cp, rho at the bulk state (AS already set at (h, p))."""
    try:
        mu, cp, rho = AS.viscosity(), AS.cpmass(), AS.rhomass()
        k = _conductivity(AS, cond)
    except ValueError:
        # Tabular backends sometimes fail on transport properties: fall back to (p, T)
        name = AS.fluid_names()[0]
        mu, cp, rho = PropsSI(("V", "CPMASS", "D"), "P", cond.p, "T", cond.T, name)
        k = (conducticity_R1233zd(cond.T, cond.p) if cond.fluid == 'R1233zd(E)'
             else PropsSI("L", "P", cond.p, "T", cond.T, name))
    return {"mu": mu, "k": k, "Pr": mu * cp / k, "cp": cp, "rho": rho}


def wall_viscosity(AS, cond, mu_bulk):
    """Viscosity at wall temperature. Falls back to the bulk value if the
    wall state cannot be evaluated (the viscosity ratio then equals 1)."""
    for dT in (0.0, -1.0, 1.0):
        try:
            AS.update(CP.PT_INPUTS, cond.p, cond.T_wall + dT)
            return AS.viscosity()
        except ValueError:
            continue
    return mu_bulk


_SAT_CACHE = {}


def saturated_properties(AS, cond):
    """Saturated liquid / vapour properties at the cell pressure.

    They depend only on (fluid, backend, p): results are cached with p rounded
    to 1 Pa, which avoids recomputing them at every residual evaluation (this
    matters for fluids with slow property fallbacks such as R1233zd(E))."""
    key = (cond.fluid, AS.backend_name(), round(cond.p))
    props = _SAT_CACHE.get(key)
    if props is None:
        if len(_SAT_CACHE) > 20000:
            _SAT_CACHE.clear()
        props = _SAT_CACHE[key] = _saturated_properties(AS, cond, float(round(cond.p)))
    return props


def _saturated_properties(AS, cond, p):
    AS.update(CP.PQ_INPUTS, p, 0)
    try:
        mu_l = AS.viscosity()
    except ValueError:
        mu_l = PropsSI('V', 'P', p, 'Q', 0, cond.fluid)
    rho_l, h_l, cp_l = AS.rhomass(), AS.hmass(), AS.cpmass()
    k_l = _conductivity(AS, cond, T=AS.T(), p=p)
    AS.update(CP.PQ_INPUTS, p, 1)
    rho_v, h_v = AS.rhomass(), AS.hmass()
    return {"mu_l": mu_l, "rho_l": rho_l, "k_l": k_l, "cp_l": cp_l,
            "Pr_l": mu_l * cp_l / k_l, "rho_v": rho_v, "i_fg": h_v - h_l}


def _x(cond):
    return min(X_MAX, max(X_MIN, cond.x))


def _m_dot_per_shell(geom, cond):
    return cond.m_dot / geom.get('n_parallel', 1)


# =============================================================================
# HEAT TRANSFER COEFFICIENT ADAPTERS  -- signature (AS, geom, G, cond) -> alpha
# =============================================================================

# ---- single phase / supercritical, channel flow ----------------------------

def _htc_gnielinski(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    mu_w = wall_viscosity(AS, cond, b["mu"])
    return gnielinski_pipe_htc(b["mu"], b["Pr"], mu_w, b["k"], G, geom['D'], geom['L'])[0]


def _htc_meshram(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return Meshram(geom['D'], G, b["k"], b["mu"], b["Pr"])


def _htc_liu(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return Liu_sCO2(G, cond.p, cond.T_wall, b["k"], b["rho"], b["mu"], b["cp"], geom['D'], cond.fluid)


def _htc_cheng_sco2(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return Cheng_sCO2(G, cond.q, cond.T_wall, cond.p, cond.h_in, cond.h_out,
                      b["mu"], b["k"], geom['D'], cond.fluid)


def _htc_pche_lee(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return PCHE_Lee(geom['alpha'], geom['D_c'], G, b["k"], geom['L'], b["mu"], b["Pr"], b["rho"])


def _htc_pche_conv(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    mu_w = wall_viscosity(AS, cond, b["mu"])
    return PCHE_conv(geom['alpha'], geom['D_c'], G, b["k"], geom['L'], b["mu"], mu_w,
                     b["Pr"], cond.T, geom['type_channel'])


# ---- single phase, plates ----------------------------------------------------

def _htc_water_plate(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return water_plate_HTC(b["mu"], b["Pr"], b["k"], G, geom['D'])


def _htc_martin_holger_plate(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return martin_holger_plate_HTC(b["mu"], b["Pr"], b["k"], cond.m_dot, geom['n_canals'], cond.T,
                                   cond.p, cond.fluid, geom['D'], geom['l'], geom['w'],
                                   geom['amplitude'], geom['chevron_angle'])


def _htc_martin_bphex(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    mu_w = wall_viscosity(AS, cond, b["mu"])
    return martin_BPHEX_HTC(b["mu"], mu_w, b["Pr"], b["k"], G, geom['D'], geom['chevron_angle'])


def _htc_muley_manglik(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    mu_w = wall_viscosity(AS, cond, b["mu"])
    return muley_manglik_BPHEX_HTC(b["mu"], mu_w, b["Pr"], b["k"], G, geom['D'], geom['chevron_angle'])


def _htc_thonon_plate(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    return thonon_plate_HTC(b["mu"], b["Pr"], b["k"], G, geom['D'], geom['chevron_angle'])


def _htc_kumar_plate(AS, geom, G, cond):
    b = bulk_properties(AS, cond)
    mu_w = wall_viscosity(AS, cond, b["mu"])
    return kumar_plate_HTC(b["mu"], b["Pr"], b["k"], G, geom['D'], mu_w, geom['chevron_angle'])


# ---- shell side / finned side ---------------------------------------------------

def _htc_shell_kern(AS, geom, G, cond):
    return shell_htc_kern(_m_dot_per_shell(geom, cond), cond.T_wall, cond.T, cond.p, AS, geom)[0]


def _htc_shell_bell_delaware(AS, geom, G, cond):
    return shell_bell_delaware_htc(_m_dot_per_shell(geom, cond), cond.T, cond.T_wall, cond.p,
                                   cond.fluid, geom)


def _htc_tube_and_fins(AS, geom, G, cond):
    return htc_tube_and_fins(cond.fluid, geom, cond.p, cond.h, _m_dot_per_shell(geom, cond),
                             geom['Fin_type'])[0]


# ---- two phase -------------------------------------------------------------------

def _htc_han_cond_bphex(AS, geom, G, cond):
    s = saturated_properties(AS, cond)
    return han_cond_BPHEX_HTC(_x(cond), s["mu_l"], s["k_l"], s["Pr_l"], s["rho_l"], s["rho_v"], G,
                              geom['D'], geom['plate_pitch_co'], geom['chevron_angle'], geom['l_v'],
                              geom['n_canals'], cond.m_dot, geom['canal_t'])


def _htc_han_boiling_bphex(AS, geom, G, cond):
    s = saturated_properties(AS, cond)
    return float(np.squeeze(han_boiling_BPHEX_HTC(
        _x(cond), s["mu_l"], s["k_l"], s["Pr_l"], s["rho_l"], s["rho_v"], s["i_fg"], G,
        cond.DT_lm, cond.Q, cond.alpha_other, geom['D'], geom['chevron_angle'], geom['plate_pitch_co'])))


def _htc_gnielinski_liquid(AS, geom, G, cond):
    """Gnielinski evaluated with saturated-liquid properties (two-phase use)."""
    s = saturated_properties(AS, cond)
    return gnielinski_pipe_htc(s["mu_l"], s["Pr_l"], s["mu_l"], s["k_l"], G, geom['D'], geom['L'])[0]


def _htc_amalfi_plate(AS, geom, G, cond):
    return amalfi_plate_HTC(geom['D'], geom['l'], geom['w'], geom['amplitude'], geom['chevron_angle'],
                            geom['n_canals'], geom['A'], cond.m_dot, cond.p, AS)


def _htc_shah_condensation_plate(AS, geom, G, cond):
    return shah_condensation_plate_HTC(geom['D'], geom['l_v'], geom['w_v'], geom['amplitude'],
                                       geom['phi'], cond.m_dot, cond.p, geom['n_canals'], AS)


def _htc_flow_boiling(AS, geom, G, cond):
    return horizontal_flow_boiling(AS, G, cond.p, _x(cond), geom['D'], cond.q)


def _htc_gungor_winterton(AS, geom, G, cond):
    s = saturated_properties(AS, cond)
    return flow_boiling_gungor_winterton(cond.fluid, G, cond.p, _x(cond), geom['D'], cond.q,
                                         s["mu_l"], s["Pr_l"], s["k_l"])


def _htc_choi_boiling(AS, geom, G, cond):
    return choi_boiling(AS, cond.p, _x(cond), geom['D'], G, cond.q)


def _htc_thome_condensation(AS, geom, G, cond):
    return thome_condensation(AS, geom['D'], G, cond.p, cond.T_sat, cond.T_wall, _x(cond))


def _htc_ext_tube_film_condens(AS, geom, G, cond):
    rho = AS.rhomass()
    return ext_tube_film_condens(geom['Tube_OD'], cond.fluid, cond.T, cond.T_wall, G / rho)


def _htc_horizontal_tube_internal_condensation(AS, geom, G, cond):
    # Placeholder kept from the original model: the correlation call was
    # commented out there and a constant value returned instead.
    return 20000.0


@lru_cache(maxsize=64)
def _boiling_curve_interp(fluid, p_rounded, D_out):
    AS = CP.AbstractState("HEOS", fluid)
    AS.update(CP.PQ_INPUTS, p_rounded, 0)
    h_boil, DT = boiling_curve(D_out, fluid, AS.T(), p_rounded)
    return interp1d(DT, h_boil, bounds_error=False, fill_value=(h_boil[0], h_boil[-1]))


def _htc_boiling_curve(AS, geom, G, cond):
    """Pool boiling curve computed once at the side inlet pressure (cached)."""
    D_out = geom.get('Tube_OD')
    if D_out is None:
        warnings.warn("'Boiling_curve' needs 'Tube_OD'; using the original constant 20000 W/m2K.",
                      stacklevel=3)
        return 20000.0
    try:
        f = _boiling_curve_interp(cond.fluid, round(cond.p_su, -2), D_out)
    except Exception:
        return 20000.0   # same fallback as the original model
    return float(f(abs(cond.T_wall - cond.T)))


HTC_CORRELATIONS = {
    # single phase / supercritical
    "Gnielinski": _htc_gnielinski,
    "Meshram": _htc_meshram,
    "Liu": _htc_liu,
    "Cheng_sCO2": _htc_cheng_sco2,
    "Lee": _htc_pche_lee,
    "PCHE_conv": _htc_pche_conv,
    "water_plate_HTC": _htc_water_plate,
    "martin_holger_plate_HTC": _htc_martin_holger_plate,
    "martin_BPHEX_HTC": _htc_martin_bphex,
    "muley_manglik_BPHEX_HTC": _htc_muley_manglik,
    "thonon_plate_HTC": _htc_thonon_plate,
    "kumar_plate_HTC": _htc_kumar_plate,
    "Shell_Kern_HTC": _htc_shell_kern,
    "Shell_Bell_Delaware_HTC": _htc_shell_bell_delaware,
    "Tube_And_Fins": _htc_tube_and_fins,
    # two phase
    "Han_cond_BPHEX": _htc_han_cond_bphex,
    "Han_Boiling_BPHEX_HTC": _htc_han_boiling_bphex,
    "Boiling_curve": _htc_boiling_curve,
    "gnielinski_pipe_htc": _htc_gnielinski_liquid,
    "amalfi_plate_HTC": _htc_amalfi_plate,
    "shah_condensation_plate_HTC": _htc_shah_condensation_plate,
    "Flow_boiling": _htc_flow_boiling,
    "Flow_boiling_gungor_winterton": _htc_gungor_winterton,
    "choi_boiling": _htc_choi_boiling,
    "Thome_Condensation": _htc_thome_condensation,
    "ext_tube_film_condens": _htc_ext_tube_film_condens,
    "Horizontal_Tube_Internal_Condensation": _htc_horizontal_tube_internal_condensation,
}

# Alternative names accepted for backward compatibility with existing scripts
HTC_ALIASES = {
    "Liu_sCO2": "Liu",
    "PCHE_Lee": "PCHE_conv",   # the original model mapped the name 'PCHE_Lee' to PCHE_conv
}


# =============================================================================
# PRESSURE DROP ADAPTERS  -- signature (AS, geom, G, cond) -> dP over full length
# =============================================================================

_PIPE_1P = ('Churchill', 'Swamee-Jain', 'Haaland', 'Konakov', 'Petukhov', 'Cheng-CO2')
_PIPE_2P = ('Friedel', 'MSH')


def _make_pipe_1p(name):
    def _dp(AS, geom, G, cond):
        return pressure_drop_pipe_single_phase(AS, geom, G, correlation=name)
    _dp.__name__ = f"_dp_pipe_{name}"
    return _dp


def _make_pipe_2p(name):
    def _dp(AS, geom, G, cond):
        return pressure_drop_pipe_frictional_two_phase(AS, geom, G, correlation=name)
    _dp.__name__ = f"_dp_pipe_{name}"
    return _dp


def _dp_pipe_choi(AS, geom, G, cond):
    # Choi needs the state at the end of the flow path; the cell outlet is used.
    AS_ex = CP.AbstractState(AS.backend_name(), AS.fluid_names()[0])
    AS_ex.update(CP.HmassP_INPUTS, cond.h_out, cond.p)
    return pressure_drop_pipe_frictional_two_phase(AS, geom, G, correlation='Choi', AS_ex=AS_ex)


def _dp_shell_kern(AS, geom, G, cond):
    return shell_DP_kern(_m_dot_per_shell(geom, cond), cond.T_wall, cond.h, cond.p, AS, geom) \
        * geom.get('n_series', 1)


def _dp_shell_donohue(AS, geom, G, cond):
    return shell_DP_donohue(_m_dot_per_shell(geom, cond), cond.T, cond.p, AS, geom) \
        * geom.get('n_series', 1)


def _dp_shell_bell_delaware(AS, geom, G, cond):
    return shell_bell_delaware_DP(_m_dot_per_shell(geom, cond), cond.h, cond.p, AS, geom) \
        * geom.get('n_series', 1)


def _dp_tube_and_fins(AS, geom, G, cond):
    return DP_tube_and_fins(AS, geom, cond.p, cond.h, cond.m_dot)


DP_CORRELATIONS = {
    **{n: _make_pipe_1p(n) for n in _PIPE_1P},
    **{n: _make_pipe_2p(n) for n in _PIPE_2P},
    "Choi": _dp_pipe_choi,
    "Shell_Kern_DP": _dp_shell_kern,
    "Shell_Donohue_DP": _dp_shell_donohue,
    "Shell_Bell_Delaware_DP": _dp_shell_bell_delaware,
    "Tube_And_Fins_DP": _dp_tube_and_fins,
}

DP_ALIASES = {
    "Muller_Steinhagen_Heck_DP": "MSH",
    "Choi_DP": "Choi",
    "Cheng_CO2_DP": "Cheng-CO2",
}


# =============================================================================
# DISPATCHERS -- the two functions the heat exchanger model calls
# =============================================================================

def _resolve(name, registry, aliases, kind):
    name = aliases.get(name, name)
    try:
        return registry[name]
    except KeyError:
        raise ValueError(f"Unknown {kind} correlation {name!r}. Available: "
                         f"{sorted(registry) + sorted(aliases)}") from None


def check_correlation_name(name, kind):
    """Raise a clear error early if a correlation name is unknown."""
    if kind == "htc":
        _resolve(name, HTC_CORRELATIONS, HTC_ALIASES, "heat transfer")
    else:
        _resolve(name, DP_CORRELATIONS, DP_ALIASES, "pressure drop")


def heat_transfer_coefficient(AS, geom, G, correlation, cond):
    """Heat transfer coefficient [W/(m^2 K)] of one side of one cell."""
    f = _resolve(correlation, HTC_CORRELATIONS, HTC_ALIASES, "heat transfer")
    return float(f(AS, geom, G, cond))


def pressure_drop(AS, geom, G, correlation, cond):
    """Pressure drop [Pa] over the full flow length of one side."""
    f = _resolve(correlation, DP_CORRELATIONS, DP_ALIASES, "pressure drop")
    return float(f(AS, geom, G, cond))
