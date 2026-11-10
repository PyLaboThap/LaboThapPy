# -*- coding: utf-8 -*-
"""
In-pipe (tube-side) HTC dispatcher, same style as the BPHX / STHE dispatchers:
flow state (SP / EVAP / COND) -> correlation, with the (AS, geom, cond) formalism.

Automatic choice = fluid family priority list, filtered by what each correlation can actually
do with the data available (required cond / geom keys, validity window, fluid data).
The first feasible correlation of the list is taken.
"""

import warnings
import numpy as np

from labothappy.correlations.convection.pipe_HTC_correlation_bank import (
    pipe_gnielinski_SP_HTC, pipe_dittus_boelter_SP_HTC, pipe_meshram_HTC,
    pipe_liu_sCO2_HTC, pipe_cheng_sCO2_HTC,
    pipe_steiner_taborek_boiling_HTC, pipe_vdi_boiling_HTC,
    pipe_gungor_winterton_boiling_HTC, pipe_choi_boiling_HTC,
    pipe_thome_condensing_HTC, pipe_dobson_chato_condensing_HTC,
    _tube_ID, _tube_G, _A0_VDI)

#%% FLUID FAMILIES

_WATER = {'water', 'h2o'}
_AMMONIA = {'ammonia', 'nh3'}
_CO2 = {'co2', 'carbondioxide'}
_HYDROCARBONS = {'methane', 'ethane', 'propane', 'n-butane', 'isobutane', 'n-pentane', 'isopentane',
                 'cyclopentane', 'n-propane', 'isohexane', 'n-hexane', 'cyclohexane', 'n-heptane', 'benzene', 'toluene',
                 'propylene', 'ethylene', 'n-octane', 'n-decane', 'n-dodecane'}


def fluid_family(AS):
    """'water', 'ammonia', 'co2', 'refrigerant' (R134a, R1233zd(E), R245fa, ...), 'hydrocarbon' or 'other'."""
    name = AS.fluid_names()[0].lower()
    if name in _WATER:
        return 'water'
    if name in _AMMONIA:
        return 'ammonia'
    if name in _CO2:
        return 'co2'
    if name in _HYDROCARBONS:
        return 'hydrocarbon'
    if (name.startswith('r') and len(name) > 1 and name[1].isdigit()) or name.startswith(('hfo', 'novec')):
        return 'refrigerant'
    return 'other'


def _D(geom):
    return _tube_ID(geom)

def _has(cond, *keys):
    return all(k in cond for k in keys)

#%% SINGLE PHASE

_PIPE_SP_HTC = {
    'Gnielinski': pipe_gnielinski_SP_HTC,
    'Dittus_Boelter': pipe_dittus_boelter_SP_HTC,
    'Meshram': pipe_meshram_HTC,
    'Liu': pipe_liu_sCO2_HTC,
    'Cheng': pipe_cheng_sCO2_HTC,
}


def _is_supercritical(AS, cond):
    p = cond.get('p')
    return p is not None and p > AS.p_critical()


def select_pipe_HTC(AS, geom, cond):
    """Automatic choice of the single-phase in-tube HTC (heuristic).
      1. supercritical CO2 cooling, inside the validity window of a dedicated correlation:
           Cheng (D ~ 10 mm, G 497-1346, q 97-400 kW/m2, p 7.5-23.5 MPa; needs T_w, q, h_in/h_out), else
           Liu   (D 4-10.7 mm, p 7.5-8.5 MPa, G 74-200; needs T_w, rho, cp)
      2. otherwise Gnielinski (laminar Stephan-Preusser / transition Gnielinski / Sieder-Tate;
         needs Tube_L), or Dittus-Boelter when Tube_L is not known.
    Supercritical CO2 outside both windows falls back on Gnielinski, with a warning (the property
    variation is not captured). Meshram (PCHE) is never automatic: give it explicitly."""
    D = _D(geom)
    G = _tube_G(geom, cond)
    if _is_supercritical(AS, cond) and fluid_family(AS) == 'co2':
        p, q = cond['p'], cond.get('q')
        if q is not None and _has(cond, 'T_w') and 8e-3 <= D <= 12e-3 and 496.7 <= G <= 1346.2 \
                and 97.4e3 <= q <= 400.3e3 and 7.53e6 <= p <= 23.51e6 and ('h_in' in cond or 'h' in cond):
            return 'Cheng'
        if _has(cond, 'T_w', 'rho', 'cp') and 4e-3 <= D <= 10.7e-3 and 7.5e6 <= p <= 8.5e6 and 74 <= G <= 200:
            return 'Liu'
        warnings.warn("Supercritical CO2 outside the Cheng and Liu validity windows (or wall data missing): "
                      "Gnielinski used, property variations not captured.")
    return 'Gnielinski' if 'Tube_L' in geom else 'Dittus_Boelter'


def pipe_SP_HTC(AS, geom, cond, choice=None):
    """Single-phase in-tube heat transfer coefficient [W/m2/K]. choice: None/'auto' or a name of
    _PIPE_SP_HTC.  cond : G (or m_dot), mu, Pr, k (+ mu_w, heating, ...: see the correlations)."""
    if choice in (None, 'auto'):
        choice = select_pipe_HTC(AS, geom, cond)
    elif choice not in _PIPE_SP_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_PIPE_SP_HTC)} or None/'auto'.")
    return _PIPE_SP_HTC[choice](AS, geom, cond)

#%% TWO-PHASE : BOILING

_PIPE_EVAP_HTC = {
    'Steiner_Taborek': pipe_steiner_taborek_boiling_HTC,
    'VDI': pipe_vdi_boiling_HTC,
    'Gungor_Winterton': pipe_gungor_winterton_boiling_HTC,
    'Choi': pipe_choi_boiling_HTC,
}

# Priority by fluid family (first feasible wins). Reasons:
#  - water       : Gungor-Winterton is built on a water + refrigerant database; VDI / Steiner-Taborek need a_0.
#  - ammonia     : Steiner-Taborek / VDI (reference nucleate boiling data exist); GW has no ammonia in its base.
#  - co2         : Choi for mini-tubes (D <= 3.5 mm, R22/R134a/CO2 mini-channel data), else VDI, Steiner-Taborek.
#  - refrigerant : Gungor-Winterton needs no fluid-specific constant (works for R1233zd(E), R245fa, ...).
#  - hydrocarbon : Steiner-Taborek (own h_nb_o: cyclopentane estimate in the bank), then VDI.
_EVAP_PRIORITY = {
    'water':       ['Gungor_Winterton', 'VDI', 'Steiner_Taborek'],
    'ammonia':     ['Steiner_Taborek', 'VDI'],
    'co2':         ['Choi', 'VDI', 'Steiner_Taborek'],
    'refrigerant': ['Gungor_Winterton', 'Steiner_Taborek', 'VDI'],
    'hydrocarbon': ['Steiner_Taborek', 'VDI', 'Gungor_Winterton'],
    'other':       ['Steiner_Taborek', 'Gungor_Winterton', 'VDI'],
}


def _a0_known(AS, cond):
    return 'a_0' in cond or AS.fluid_names()[0] in _A0_VDI


def _evap_feasible(name, AS, geom, cond):
    if not _has(cond, 'p', 'x', 'q'):
        return False
    fluid = AS.fluid_names()[0]
    if name == 'VDI':
        return _a0_known(AS, cond)
    if name == 'Steiner_Taborek':
        return 'Tube_L' in geom and ('h_nb_o' in cond or fluid == 'Cyclopentane' or _a0_known(AS, cond))
    if name == 'Choi':
        return _D(geom) <= 3.5e-3                          # mini-tube data only
    return True                                           # Gungor_Winterton


def select_pipe_evaporation_HTC(AS, geom, cond):
    """Automatic choice of the in-tube evaporation HTC: priority list of the fluid family,
    first correlation whose data (cond / geom keys, a_0, tube size) are available."""
    fam = fluid_family(AS)
    for name in _EVAP_PRIORITY[fam]:
        if _evap_feasible(name, AS, geom, cond):
            return name
    for name in _PIPE_EVAP_HTC:                           # last resort: anything that can run
        if _evap_feasible(name, AS, geom, cond):
            warnings.warn(f"Evaporation HTC: no priority correlation feasible for family {fam!r}: {name} used.")
            return name
    raise ValueError("Evaporation HTC: no feasible correlation (need cond p, x, q and, for "
                     "Steiner-Taborek / VDI, Tube_L and a_0 / h_nb_o for this fluid).")


def pipe_boiling_HTC(AS, geom, cond, choice=None):
    """Evaporation in-tube heat transfer coefficient [W/m2/K]. choice: None/'auto' or a name of _PIPE_EVAP_HTC.
    cond : G (or m_dot), p, x, q [W/m2]; a_0 or h_nb_o for fluids missing from the VDI table.
    geom : Tube_ID (or Tube_OD, Tube_t), Tube_L."""
    if choice in (None, 'auto'):
        choice = select_pipe_evaporation_HTC(AS, geom, cond)
    elif choice not in _PIPE_EVAP_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_PIPE_EVAP_HTC)} or None/'auto'.")
    return _PIPE_EVAP_HTC[choice](AS, geom, cond)

#%% TWO-PHASE : CONDENSING

_PIPE_COND_HTC = {
    'Thome': pipe_thome_condensing_HTC,
    'Dobson_Chato': pipe_dobson_chato_condensing_HTC,
}


def _thome_available():
    try:
        from labothappy.correlations.void_fraction.void_fraction_old import void_fraction  # noqa: F401
        return True
    except ImportError:
        return False


def select_pipe_condensation_HTC(AS, geom, cond):
    """Automatic choice of the in-tube condensation HTC (heuristic); both need cond T_w.
      - small tube (D < 8 mm) and refrigerant / CO2 / water family -> Dobson-Chato
        (data: 3-7 mm tubes, G 25-800 kg/m2/s, R12 / R22 / R134a / R32-R125)
      - otherwise                                                  -> Thome (flow-regime based, wide
        refrigerant database, 3-21 mm tubes; fluid enters through its properties only)
    Thome falls back on Dobson-Chato if the void fraction module cannot be imported."""
    if 'T_w' not in cond:
        raise KeyError("Condensation HTC: cond needs 'T_w'.")
    D = _D(geom)
    fam = fluid_family(AS)
    if D < 8e-3 and fam in ('refrigerant', 'co2', 'water'):
        return 'Dobson_Chato'
    return 'Thome' if _thome_available() else 'Dobson_Chato'


def pipe_condensing_HTC(AS, geom, cond, choice=None):
    """Condensation in-tube heat transfer coefficient [W/m2/K]. choice: None/'auto' or a name of _PIPE_COND_HTC.
    cond : G (or m_dot), p, x, T_w [K].    geom : Tube_ID (or Tube_OD, Tube_t)."""
    if choice in (None, 'auto'):
        choice = select_pipe_condensation_HTC(AS, geom, cond)
    elif choice not in _PIPE_COND_HTC:
        raise ValueError(f"Unknown correlation {choice!r}: {sorted(_PIPE_COND_HTC)} or None/'auto'.")
    return _PIPE_COND_HTC[choice](AS, geom, cond)

#%% GLOBAL DISPATCHER

_STATE_DISPATCH = {
    'SP': pipe_SP_HTC,
    'EVAP': pipe_boiling_HTC,
    'COND': pipe_condensing_HTC,
}

def pipe_HTC(AS, geom, cond, flow_state, choices=None):
    """HTC inside the tubes of a shell & tube heat exchanger.
    flow_state : 'SP', 'EVAP' or 'COND'
    choices    : dict {'SP': name, 'EVAP': name, 'COND': name}.
                 None, empty, 'auto' (or empty / 'auto' for a state) = automatic selection."""
    if flow_state not in _STATE_DISPATCH:
        raise ValueError(f"Unknown flow state {flow_state!r}: expected one of {sorted(_STATE_DISPATCH)}.")

    if not choices or (isinstance(choices, str) and choices.lower() == 'auto'):
        choices = {}
    elif not isinstance(choices, dict):
        raise TypeError(f"choices must be None, 'auto' or a dict like {{'SP': 'Gnielinski'}}, "
                        f"got {type(choices).__name__}.")

    unknown = set(choices) - set(_STATE_DISPATCH)
    if unknown:
        raise ValueError(f"Unknown keys in choices: {sorted(unknown)} (expected {sorted(_STATE_DISPATCH)}).")

    choice = choices.get(flow_state)
    if not choice or (isinstance(choice, str) and choice.lower() == 'auto'):
        choice = None
    return _STATE_DISPATCH[flow_state](AS, geom, cond, choice)

