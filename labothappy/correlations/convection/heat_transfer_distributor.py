"""
Single entry point for all heat transfer coefficient correlations.

    htc = heat_transfer_coefficient(AS, geom, G, correlation, **extra)

AS          : CoolProp AbstractState, already set at the mean state of the cell (h, p).
              It is left in the same state on return.
geom        : geometry of ONE side, built by the component model
                'D' hydraulic diameter [m], 'L' flow length [m]
              + correlation-specific keys (chevron_angle, plate_length, ...)
G           : mass flux [kg/(m^2 s)]
correlation : name, see available_correlations()
**extra     : operating data that AS does not hold. Only the correlations that
              need them read them (a clear error is raised if one is missing):
                m_dot     : mass flow rate of the side, whole heat exchanger [kg/s]
                T_wall    : wall temperature [K]
                x         : vapour quality [-] (default: quality of AS)
                q         : heat flux [W/m^2]
                Q         : heat rate of the cell [W]
                DT_lm     : (corrected) log mean temperature difference of the cell [K]
                htc_other : heat transfer coefficient of the other side [W/(m^2 K)]
                h_min     : enthalpy at the cell boundary of lower enthalpy [J/kg]
                h_max     : enthalpy at the cell boundary of higher enthalpy [J/kg]

Returns the heat transfer coefficient [W/(m^2 K)].
"""

from contextlib import contextmanager
from functools import lru_cache

import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

from labothappy.correlations.convection.pipe_htc import (
    gnielinski_pipe_htc, Meshram, Liu_sCO2, Cheng_sCO2,
    horizontal_flow_boiling, flow_boiling_gungor_winterton, choi_boiling, thome_condensation)
from labothappy.correlations.convection.plate_htc import (
    water_plate_HTC, martin_holger_plate_HTC, martin_BPHEX_HTC, muley_manglik_BPHEX_HTC,
    thonon_plate_HTC, kumar_plate_HTC, amalfi_plate_HTC, han_boiling_BPHEX_HTC,
    han_cond_BPHEX_HTC, shah_condensation_plate_HTC)
from labothappy.correlations.convection.shell_and_tube_htc import shell_htc_kern, shell_bell_delaware_htc
from labothappy.correlations.convection.fins_htc import htc_tube_and_fins
from labothappy.correlations.convection.tube_bank_htc import ext_tube_film_condens
from labothappy.correlations.convection.printed_circuit_htc import PCHE_Lee, PCHE_conv
from labothappy.correlations.properties.thermal_conductivity import conducticity_R1233zd

X_MIN, X_MAX = 1e-4, 1.0 - 1e-4   # Quality clipping for the two-phase correlations [-]


# =============================================================================
# 0) PROPERTY HELPERS
# =============================================================================

def _fluid(AS):
    """Fluid name usable by PropsSI. The incompressible backend has no fluid_names()
    and PropsSI needs the 'INCOMP::' prefix for those fluids."""
    try:
        return AS.fluid_names()[0]
    except ValueError:
        return "INCOMP::" + AS.name()


def _need(extra, name, correlation):
    """Value of an extra argument, with a clear error if the model did not give it."""
    if extra.get(name) is None:
        raise ValueError(f"Correlation {correlation!r} needs the extra argument {name!r}.")
    return extra[name]


@contextmanager
def _restored(AS):
    """Give the state back to the caller after a correlation changed it."""
    h, p = AS.hmass(), AS.p()
    try:
        yield AS
    finally:
        AS.update(CP.HmassP_INPUTS, h, p)


def _conductivity(AS, T, p):
    """Thermal conductivity, with the R1233zd(E) fit that CoolProp lacks."""
    if _fluid(AS) == 'R1233zd(E)':
        return conducticity_R1233zd(T, p)
    return AS.conductivity()


def _bulk(AS):
    """Transport properties at the bulk state of the cell."""
    T, p = AS.T(), AS.p()
    try:
        mu, cp, rho = AS.viscosity(), AS.cpmass(), AS.rhomass()
        k = _conductivity(AS, T, p)
    except ValueError:  # Tabular backends sometimes fail on transport properties: use (p, T)
        mu, cp, rho = PropsSI(("V", "CPMASS", "D"), "P", p, "T", T, _fluid(AS))
        k = conducticity_R1233zd(T, p) if _fluid(AS) == 'R1233zd(E)' else PropsSI("L", "P", p, "T", T, _fluid(AS))
    return {"mu": mu, "k": k, "Pr": mu * cp / k, "cp": cp, "rho": rho, "T": T, "p": p}


def _mu_wall(AS, T_wall):
    """Viscosity at the wall temperature. Without T_wall, or if the wall state
    cannot be computed, the bulk viscosity is used (viscosity ratio = 1)."""
    mu_bulk = _bulk(AS)["mu"]
    if T_wall is None:
        return mu_bulk
    with _restored(AS):
        for dT in (0.0, -1.0, 1.0):
            try:
                AS.update(CP.PT_INPUTS, AS.p(), T_wall + dT)
                return AS.viscosity()
            except ValueError:
                continue
    return mu_bulk


@lru_cache(maxsize=2048)
def _saturated_cached(fluid, p):
    AS = CP.AbstractState("HEOS", fluid)
    AS.update(CP.PQ_INPUTS, p, 0)
    mu_l, rho_l, h_l, cp_l, T_sat = AS.viscosity(), AS.rhomass(), AS.hmass(), AS.cpmass(), AS.T()
    k_l = conducticity_R1233zd(T_sat, p) if fluid == 'R1233zd(E)' else AS.conductivity()
    AS.update(CP.PQ_INPUTS, p, 1)
    rho_v, h_v = AS.rhomass(), AS.hmass()
    return {"mu_l": mu_l, "rho_l": rho_l, "k_l": k_l, "Pr_l": mu_l * cp_l / k_l,
            "rho_v": rho_v, "i_fg": h_v - h_l, "T_sat": T_sat}


def _saturated(AS):
    """Saturated liquid and vapour properties at the cell pressure (cached per 1 Pa)."""
    return _saturated_cached(_fluid(AS), float(round(AS.p())))


def _x(AS, extra):
    """Vapour quality: the one given by the model, otherwise the one of AS."""
    x = extra.get('x')
    if x is None or np.isnan(x):
        x = AS.Q()
    return min(X_MAX, max(X_MIN, x))


# =============================================================================
# 1) SINGLE PHASE - channels and pipes
# =============================================================================

def _htc_gnielinski(AS, geom, G, **extra):
    b = _bulk(AS)
    return gnielinski_pipe_htc(b["mu"], b["Pr"], _mu_wall(AS, extra.get('T_wall')), b["k"], G, geom['D'], geom['L'])[0]

def _htc_pche_lee(AS, geom, G, **extra):
    b = _bulk(AS)
    return PCHE_Lee(geom['alpha'], geom['D_c'], G, b["k"], geom['L'], b["mu"], b["Pr"], b["rho"])

def _htc_pche_conv(AS, geom, G, **extra):
    b = _bulk(AS)
    return PCHE_conv(geom['alpha'], geom['D_c'], G, b["k"], geom['L'], b["mu"],
                     _mu_wall(AS, extra.get('T_wall')), b["Pr"], b["T"], geom['type_channel'])


# =============================================================================
# 2) SINGLE PHASE - plates
# =============================================================================

def _htc_water_plate(AS, geom, G, **extra):
    b = _bulk(AS)
    return water_plate_HTC(b["mu"], b["Pr"], b["k"], G, geom['D'])

def _htc_martin_holger_plate(AS, geom, G, **extra):
    b = _bulk(AS)
    m_dot = _need(extra, 'm_dot', 'martin_holger_plate_HTC')
    return martin_holger_plate_HTC(b["mu"], b["Pr"], b["k"], m_dot, geom['n_channels'], b["T"], b["p"],
                                   _fluid(AS), geom['D'], geom['plate_length'], geom['plate_width'],
                                   geom['amplitude'], geom['chevron_angle'])

def _htc_martin_bphex(AS, geom, G, **extra):
    b = _bulk(AS)
    return martin_BPHEX_HTC(b["mu"], _mu_wall(AS, extra.get('T_wall')), b["Pr"], b["k"], G,
                            geom['D'], geom['chevron_angle'])

def _htc_muley_manglik(AS, geom, G, **extra):
    b = _bulk(AS)
    return muley_manglik_BPHEX_HTC(b["mu"], _mu_wall(AS, extra.get('T_wall')), b["Pr"], b["k"], G,
                                   geom['D'], geom['chevron_angle'])

def _htc_thonon_plate(AS, geom, G, **extra):
    b = _bulk(AS)
    return thonon_plate_HTC(b["mu"], b["Pr"], b["k"], G, geom['D'], geom['chevron_angle'])

def _htc_kumar_plate(AS, geom, G, **extra):
    b = _bulk(AS)
    return kumar_plate_HTC(b["mu"], b["Pr"], b["k"], G, geom['D'], _mu_wall(AS, extra.get('T_wall')),
                           geom['chevron_angle'])


# =============================================================================
# 3) SINGLE PHASE - shell side and finned side
# =============================================================================

def _m_dot_per_shell(geom, extra, correlation):
    return _need(extra, 'm_dot', correlation) / geom.get('n_parallel', 1)

def _htc_shell_kern(AS, geom, G, **extra):
    T_wall = _need(extra, 'T_wall', 'Shell_Kern_HTC')
    with _restored(AS):
        return shell_htc_kern(_m_dot_per_shell(geom, extra, 'Shell_Kern_HTC'), T_wall, AS.T(), AS.p(), AS, geom)[0]

def _htc_shell_bell_delaware(AS, geom, G, **extra):
    T_wall = _need(extra, 'T_wall', 'Shell_Bell_Delaware_HTC')
    return shell_bell_delaware_htc(_m_dot_per_shell(geom, extra, 'Shell_Bell_Delaware_HTC'), AS.T(), T_wall,
                                   AS.p(), _fluid(AS), geom)

def _htc_tube_and_fins(AS, geom, G, **extra):
    return htc_tube_and_fins(_fluid(AS), geom, AS.p(), AS.hmass(),
                             _m_dot_per_shell(geom, extra, 'Tube_And_Fins'), geom['Fin_type'])[0]


# =============================================================================
# 4) SUPERCRITICAL
# =============================================================================

def _htc_meshram(AS, geom, G, **extra):
    b = _bulk(AS)
    return Meshram(geom['D'], G, b["k"], b["mu"], b["Pr"])

def _htc_liu(AS, geom, G, **extra):
    b = _bulk(AS)
    T_wall = _need(extra, 'T_wall', 'Liu')
    return Liu_sCO2(G, b["p"], T_wall, b["k"], b["rho"], b["mu"], b["cp"], geom['D'], _fluid(AS))

def _htc_cheng_sco2(AS, geom, G, **extra):
    b = _bulk(AS)
    args = [_need(extra, name, 'Cheng_sCO2') for name in ('q', 'T_wall', 'h_min', 'h_max')]
    q, T_wall, h_min, h_max = args
    return Cheng_sCO2(G, q, T_wall, b["p"], h_min, h_max, b["mu"], b["k"], geom['D'], _fluid(AS))


# =============================================================================
# 5) TWO PHASE - boiling
# =============================================================================

def _htc_han_boiling_bphex(AS, geom, G, **extra):
    s = _saturated(AS)
    DT_lm, Q, htc_other = [_need(extra, name, 'Han_Boiling_BPHEX_HTC') for name in ('DT_lm', 'Q', 'htc_other')]
    return float(np.squeeze(han_boiling_BPHEX_HTC(
        _x(AS, extra), s["mu_l"], s["k_l"], s["Pr_l"], s["rho_l"], s["rho_v"], s["i_fg"], G,
        DT_lm, Q, htc_other, geom['D'], geom['chevron_angle'], geom['corrugation_pitch'])))

def _htc_amalfi_plate(AS, geom, G, **extra):
    m_dot = _need(extra, 'm_dot', 'amalfi_plate_HTC')
    with _restored(AS):
        return amalfi_plate_HTC(geom['D'], geom['plate_length'], geom['plate_width'], geom['amplitude'],
                                geom['chevron_angle'], geom['n_channels'], geom['A'], m_dot, AS.p(), AS)

def _htc_flow_boiling(AS, geom, G, **extra):
    q = _need(extra, 'q', 'Flow_boiling')
    with _restored(AS):
        return horizontal_flow_boiling(AS, G, AS.p(), _x(AS, extra), geom['D'], q)

def _htc_gungor_winterton(AS, geom, G, **extra):
    s = _saturated(AS)
    q = _need(extra, 'q', 'Flow_boiling_gungor_winterton')
    return flow_boiling_gungor_winterton(_fluid(AS), G, AS.p(), _x(AS, extra), geom['D'], q,
                                         s["mu_l"], s["Pr_l"], s["k_l"])

def _htc_choi_boiling(AS, geom, G, **extra):
    q = _need(extra, 'q', 'choi_boiling')
    with _restored(AS):
        return choi_boiling(AS, AS.p(), _x(AS, extra), geom['D'], G, q)


# =============================================================================
# 6) TWO PHASE - condensation
# =============================================================================

def _htc_han_cond_bphex(AS, geom, G, **extra):
    s = _saturated(AS)
    return han_cond_BPHEX_HTC(_x(AS, extra), s["mu_l"], s["k_l"], s["Pr_l"], s["rho_l"], s["rho_v"], G,
                              geom['D'], geom['corrugation_pitch'], geom['chevron_angle'],
                              geom['port_distance'], geom['n_channels'], extra.get('m_dot'), geom.get('t_channel'))

def _htc_shah_condensation_plate(AS, geom, G, **extra):
    m_dot = _need(extra, 'm_dot', 'shah_condensation_plate_HTC')
    with _restored(AS):
        return shah_condensation_plate_HTC(geom['D'], geom['port_distance'], geom['port_width'],
                                           geom['amplitude'], geom['enlargement_factor'], m_dot, AS.p(),
                                           geom['n_channels'], AS)

def _htc_thome_condensation(AS, geom, G, **extra):
    s = _saturated(AS)
    T_wall = _need(extra, 'T_wall', 'Thome_Condensation')
    with _restored(AS):
        return thome_condensation(AS, geom['D'], G, AS.p(), s["T_sat"], T_wall, _x(AS, extra))

def _htc_gnielinski_liquid(AS, geom, G, **extra):
    """Gnielinski with saturated-liquid properties (liquid-only reference for two-phase flow)."""
    s = _saturated(AS)
    return gnielinski_pipe_htc(s["mu_l"], s["Pr_l"], s["mu_l"], s["k_l"], G, geom['D'], geom['L'])[0]

def _htc_ext_tube_film_condens(AS, geom, G, **extra):
    T_wall = _need(extra, 'T_wall', 'ext_tube_film_condens')
    return ext_tube_film_condens(geom['Tube_OD'], _fluid(AS), AS.T(), T_wall, G / AS.rhomass())


# =============================================================================
# REGISTRY + DISTRIBUTOR
# =============================================================================

CORRELATIONS = {
    # name                            : (function,                       regime)
    "Gnielinski":                       (_htc_gnielinski,                "single_phase"),
    "Lee":                              (_htc_pche_lee,                  "single_phase"),
    "PCHE_conv":                        (_htc_pche_conv,                 "single_phase"),
    "water_plate_HTC":                  (_htc_water_plate,               "single_phase"),
    "martin_holger_plate_HTC":          (_htc_martin_holger_plate,       "single_phase"),
    "martin_BPHEX_HTC":                 (_htc_martin_bphex,              "single_phase"),
    "muley_manglik_BPHEX_HTC":          (_htc_muley_manglik,             "single_phase"),
    "thonon_plate_HTC":                 (_htc_thonon_plate,              "single_phase"),
    "kumar_plate_HTC":                  (_htc_kumar_plate,               "single_phase"),
    "Shell_Kern_HTC":                   (_htc_shell_kern,                "single_phase"),
    "Shell_Bell_Delaware_HTC":          (_htc_shell_bell_delaware,       "single_phase"),
    "Tube_And_Fins":                    (_htc_tube_and_fins,             "single_phase"),
    "Meshram":                          (_htc_meshram,                   "supercritical"),
    "Liu":                              (_htc_liu,                       "supercritical"),
    "Cheng_sCO2":                       (_htc_cheng_sco2,                "supercritical"),
    "Han_Boiling_BPHEX_HTC":            (_htc_han_boiling_bphex,         "two_phase"),
    "amalfi_plate_HTC":                 (_htc_amalfi_plate,              "two_phase"),
    "Flow_boiling":                     (_htc_flow_boiling,              "two_phase"),
    "Flow_boiling_gungor_winterton":    (_htc_gungor_winterton,          "two_phase"),
    "choi_boiling":                     (_htc_choi_boiling,              "two_phase"),
    "Han_cond_BPHEX":                   (_htc_han_cond_bphex,            "two_phase"),
    "shah_condensation_plate_HTC":      (_htc_shah_condensation_plate,   "two_phase"),
    "Thome_Condensation":               (_htc_thome_condensation,        "two_phase"),
    "gnielinski_pipe_htc":              (_htc_gnielinski_liquid,         "two_phase"),
    "ext_tube_film_condens":            (_htc_ext_tube_film_condens,     "two_phase"),
}

# Other names accepted for the same correlations (used in existing scripts)
ALIASES = {
    "Liu_sCO2": "Liu",
    "PCHE_Lee": "PCHE_conv",   # The original model called PCHE_conv under the name 'PCHE_Lee'
}


def available_correlations(regime=None):
    """Names of the available correlations, optionally filtered by regime
    ('single_phase', 'two_phase' or 'supercritical')."""
    return [name for name, (_, reg) in CORRELATIONS.items() if regime is None or reg == regime]


def heat_transfer_coefficient(AS, geom, G, correlation, **extra):
    """Heat transfer coefficient [W/(m^2 K)] of one side of one cell."""
    name = ALIASES.get(correlation, correlation)
    try:
        func, _ = CORRELATIONS[name]
    except KeyError:
        raise ValueError(f"Unknown heat transfer correlation {correlation!r}. "
                         f"Available: {available_correlations()}") from None
    return float(func(AS, geom, G, **extra))
