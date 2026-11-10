# -*- coding: utf-8 -*-
"""
curved pipe (helical coil) HTC dispatcher, (AS, geom, cond) formalism. cond['side']: 'tube' (inside the coil) or 'shell'.

tube  SP   : 'Sigalotti' (curved pipe, default)  or any pipe correlation name (straight tube)
      EVAP : 'Zhao', 'Nariai' or pipe names; needs cond q.
             auto: on the liquid-only Reynolds number Re_lo = G D / mu_l (saturated liquid)
                   Re_lo <  Re_cr -> 'Zhao'     (laminar liquid)
                   Re_lo >= Re_cr -> 'Nariai'   (turbulent liquid)
             Re_cr = 12730 (D / 2R)^0.2 (Ito) if geom has R_curv (or D_coil), else 2300.
      COND : 'Zakeralhoseini' (default; R1234yf, extrapolated to other fluids with a warning) or pipe names
shell side (outside the coil): NO correlation provided -> NotImplementedError.
"""
import warnings
from labothappy.correlations.convection.curved_pipe_HTC_correlation_bank import (
    curved_pipe_sigalotti_SP_HTC, curved_pipe_zakeralhoseini_condensing_HTC, curved_pipe_zhao_boiling_HTC, curved_pipe_nariai_boiling_HTC,
    _D, _R_curv)
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props
from labothappy.correlations.convection.pipe_HTC import pipe_HTC

_HCHE = {'SP': {'Sigalotti': curved_pipe_sigalotti_SP_HTC},
         'EVAP': {'Zhao': curved_pipe_zhao_boiling_HTC, 'Nariai': curved_pipe_nariai_boiling_HTC},
         'COND': {'Zakeralhoseini': curved_pipe_zakeralhoseini_condensing_HTC}}
_AUTO = {'SP': 'Sigalotti', 'COND': 'Zakeralhoseini'}      # EVAP: see select_evap_HTC


def select_evap_HTC(AS, geom, cond):
    """Laminar liquid-only flow -> 'Zhao'; turbulent -> 'Nariai'."""
    D = _D(geom)
    Re_lo = cond['G'] * D / _sat_props(AS, cond['p'], cond)['mu_l']
    try:
        Re_cr = 12730 * (D / (2 * _R_curv(geom))) ** 0.2
    except KeyError:
        Re_cr = 2300.0
    return 'Zhao' if Re_lo < Re_cr else 'Nariai'


def _equivalent_tube(geom):
    g = dict(geom)
    D = g.get('Dh') or g.get('Tube_ID') or g.get('D')
    g.update(Tube_ID=D, Dh=D, Tube_pass=1)
    g.setdefault('Tube_L', g.get('L'))
    g.setdefault('n_tubes', 1)
    return g


def curved_pipe_HTC(AS, geom, cond, flow_state, choices=None):
    side = cond.get('side', geom.get('side', 'tube'))
    if side != 'tube':
        raise NotImplementedError("HCHE HTC: shell side (outside the coil) has no correlation yet.")
    if flow_state not in _HCHE:
        raise ValueError(f"HCHE HTC: flow_state {flow_state!r}: expected 'SP', 'EVAP' or 'COND'.")
    choice = (choices or {}).get(flow_state)
    if choice in (None, 'auto'):
        choice = select_evap_HTC(AS, geom, cond) if flow_state == 'EVAP' else _AUTO[flow_state]
    if choice in _HCHE[flow_state]:
        if choice == 'Zakeralhoseini':
            warnings.warn("HCHE: Zakeralhoseini condensation correlation established for R1234yf.")
        return _HCHE[flow_state][choice](AS, geom, cond)
    warnings.warn(f"HCHE: pipe correlation {choice!r} on the straight equivalent tube (no curvature effect).")
    return pipe_HTC(AS, _equivalent_tube(geom), cond, flow_state, {flow_state: choice})