# -*- coding: utf-8 -*-
"""
PCHE HTC dispatcher: HX type -> flow state (SP / EVAP / COND) -> correlation, (AS, geom, cond) formalism.

geom['type_channel'] : 'Straight' (default) or 'Zigzag'.

SP    : PCHE correlations  'Straight', 'Lee' (zigzag), 'Bennett_Chen' (zigzag, sCO2 only, explicit choice)
        or any pipe correlation name ('Gnielinski', 'Dittus_Boelter', 'Cheng', 'Liu', ...), evaluated
        on the equivalent tube (D = D_h, G of one channel).
        auto : Straight -> 'Straight' ; Zigzag -> 'Lee'.
EVAP / COND : curved pipe correlations (HCHE dispatcher) on the equivalent tube: 'Zhao' / 'Nariai' (EVAP),
        'Zakeralhoseini' (COND); a pipe correlation name falls back to the straight pipe_HTC.
Every non-PCHE path uses the curved pipe functions (HCHE dispatcher): SP auto name -> 'Sigalotti'.
geom['R_curv'] : radius of curvature of the channel bends [m]; absent -> straight tube (R -> infinity).

The equivalent tube is read from geom: Tube_ID (= D_h), Tube_L (= developed channel length), n_tubes,
Tube_pass; cond['G'] is the mass flux of one channel (the pipe bank uses it first).
"""
import warnings
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import (
    PCHE_straight_HTC, PCHE_bennett_chen_HTC, PCHE_lee_zigzag_HTC)
from labothappy.correlations.convection.curved_pipe_HTC import curved_pipe_HTC   

_PCHE_SP_HTC = {'Straight': PCHE_straight_HTC, 'Lee': PCHE_lee_zigzag_HTC, 'Bennett_Chen': PCHE_bennett_chen_HTC}

# Pipe correlations established for tubes of several mm: warn below this hydraulic diameter [m]
_PIPE_MIN_D = {'Cheng': 8e-3, 'Liu': 4e-3, 'Gungor_Winterton': 3e-3, 'Friedel': 0.0}
_CONFINEMENT_D = 3e-3     # two-phase confinement effects below ~3 mm


def select_PCHE_HTC(geom, cond):
    """Straight channel -> Straight; zigzag -> Lee (Bennett_Chen only on explicit request)."""
    return 'Lee' if str(geom.get('type_channel', 'Straight')).lower() == 'zigzag' else 'Straight'


def _equivalent_tube(geom):
    g = dict(geom)
    Dh = g.get('Dh') or g['Tube_ID']
    g.update(Tube_ID=Dh, Dh=Dh, Tube_pass=1)
    g.setdefault('Tube_L', g.get('L_c', g.get('L')))
    g.setdefault('L', g['Tube_L'])
    if 'R_curv' not in g:
        if str(g.get('type_channel', 'Straight')).lower() == 'zigzag':
            warnings.warn("PCHE: zigzag channel without geom['R_curv']: curved pipe correlations run as a straight tube.")
        g['R_curv'] = 1e6
    return g


def _pipe(AS, geom, cond, state, name):
    g = _equivalent_tube(geom)
    if name not in (None, 'auto') and g['Tube_ID'] < _PIPE_MIN_D.get(name, 0.0):
        warnings.warn(f"PCHE: pipe correlation {name!r} used with D_h = {g['Tube_ID']*1e3:.2f} mm, "
                      "below the tube diameters it was established for.")
    if state != 'SP' and g['Tube_ID'] < _CONFINEMENT_D:
        warnings.warn(f"PCHE: two-phase flow in a D_h = {g['Tube_ID']*1e3:.2f} mm channel: confinement "
                      "effects are not in the tube correlations.")
    return curved_pipe_HTC(AS, g, {**cond, 'side': 'tube'}, state, {state: name})


def PCHE_SP_HTC(AS, geom, cond, choice=None):
    if choice in (None, 'auto'):
        choice = select_PCHE_HTC(geom, cond)
    if choice in _PCHE_SP_HTC:
        return _PCHE_SP_HTC[choice](AS, geom, cond)
    return _pipe(AS, geom, cond, 'SP', choice)        # 'Sigalotti' (curved pipe) or a pipe correlation name


def PCHE_HTC(AS, geom, cond, flow_state, choices=None):
    choice = (choices or {}).get(flow_state)
    if flow_state == 'SP':
        return PCHE_SP_HTC(AS, geom, cond, choice)
    if flow_state in ('EVAP', 'COND'):
        return _pipe(AS, geom, cond, flow_state, choice)
    raise ValueError(f"PCHE HTC: flow_state {flow_state!r}: expected 'SP', 'EVAP' or 'COND'.")