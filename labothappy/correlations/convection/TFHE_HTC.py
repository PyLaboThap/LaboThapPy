# -*- coding: utf-8 -*-
"""
Tube & fins (TFHE) HTC dispatcher: HX type -> side ('fin' / 'tube') -> flow state -> correlation,
(AS, geom, cond) formalism.

side 'fin'  : finned side (gas): 'Pismennyi' (= None / 'auto'), routed on geom['Fin_type'] ('Annular' or 'Square'),
              or a forced shape 'Pismennyi_Annular' / 'Pismennyi_Square'. Two-phase / 'vapor-wet' cells (EVAP / COND)
              reuse the same correlation (phase change ignored, warning).
              Returns the REDUCED coefficient (referred to the total outer area geom['A'], fin
              efficiency included).
side 'tube' : SP / EVAP / COND inside the tubes: pipe correlations (pipe_HTC).
"""
import warnings
from labothappy.correlations.convection.TFHE_HTC_correlation_bank import (
    TFHE_pismennyi_HTC, TFHE_pismennyi_annular_HTC, TFHE_pismennyi_square_HTC)
from labothappy.correlations.convection.pipe_HTC_correlation_bank import _sat_props, _clamp_x
from labothappy.correlations.convection.pipe_HTC import pipe_HTC

_TFHE_FIN_HTC = {'Pismennyi': TFHE_pismennyi_HTC,                       # routed on geom['Fin_type']
                 'Pismennyi_Annular': TFHE_pismennyi_annular_HTC,       # forced fin shape
                 'Pismennyi_Square': TFHE_pismennyi_square_HTC}


def select_fin_HTC(geom, cond):
    """Automatic choice: Pis'mennyi, routed on geom['Fin_type'] ('Annular' or 'Square')."""
    return 'Pismennyi'


def _homogeneous_props(AS, cond):
    """Bulk properties of a TWO-PHASE cell for the gas correlation (crude homogeneous model, quality
    cond['x']): rho = 1/(x/rho_v + (1-x)/rho_l), mu = McAdams, k and Pr = mass-weighted.
    cond['mu_l'], ['mu_v'], ['k_l'], ['k_v'] override CoolProp (see _sat_props)."""
    x = _clamp_x(cond.get('x', 0.5))
    s = _sat_props(AS, cond['p'], cond)
    return dict(rho=1.0 / (x / s['rho_v'] + (1 - x) / s['rho_l']),
                mu=1.0 / (x / s['mu_v'] + (1 - x) / s['mu_l']),
                k=x * s['k_v'] + (1 - x) * s['k_l'],
                Pr=x * s['Pr_v'] + (1 - x) * s['Pr_l'])


def _ensure_gas_props(AS, cond, flow_state):
    """cond with rho, mu, k, Pr of the cell. Missing (not provided by the model's _cond_dict): single-phase or
    'vapor-wet' cell -> read from AS (set at the mean state of the cell); two-phase cell -> homogeneous mixture."""
    if all(k in cond for k in ('rho', 'mu', 'k', 'Pr')):
        return cond
    if flow_state == 'SP' or cond.get('phase') == 'vapor-wet':
        return {**cond, 'rho': AS.rhomass(), 'mu': AS.viscosity(), 'k': AS.conductivity(), 'Pr': AS.Prandtl()}
    return {**cond, **_homogeneous_props(AS, cond)}


def TFHE_fin_HTC(AS, geom, cond, flow_state, choices=None):
    """Finned side. SP: gas correlation. 'EVAP' / 'COND' (two-phase or 'vapor-wet' cells): NO phase-change
    correlation exists on the finned side, so the SAME correlation is evaluated with the bulk properties of
    the cell (a 'vapor-wet' cell already holds them; a two-phase cell gets homogeneous-mixture properties)
    and a warning is raised: the phase change is ignored (sensible convection only)."""
    if flow_state not in ('SP', 'EVAP', 'COND'):
        raise ValueError(f"TFHE finned side: flow state {flow_state!r}, expected 'SP', 'EVAP' or 'COND'.")
    if flow_state != 'SP':
        warnings.warn(f"TFHE finned side: no {flow_state} correlation, the gas correlation is used for a "
                      f"{cond.get('phase')!r} cell (phase change ignored).", stacklevel=2)
    cond = _ensure_gas_props(AS, cond, flow_state)
    choice = (choices or {}).get(flow_state) or (choices or {}).get('SP')
    if choice in (None, 'auto'):
        choice = select_fin_HTC(geom, cond)
    elif choice in ('Pismennyi_Annular', 'Pismennyi_Square') and \
            str(geom.get('Fin_type', '')).lower() != choice.split('_')[1].lower():
        warnings.warn(f"TFHE: {choice!r} requested but geom['Fin_type'] = {geom.get('Fin_type')!r}.", stacklevel=2)
    if choice not in _TFHE_FIN_HTC:
        raise ValueError(f"Unknown finned-side correlation {choice!r}: {sorted(_TFHE_FIN_HTC)}, 'Pismennyi' or None/'auto'.")
    return _TFHE_FIN_HTC[choice](AS, geom, cond)


def TFHE_HTC(AS, geom, cond, flow_state, choices=None):
    side = cond.get('side', geom.get('side'))
    if side == 'fin':
        return TFHE_fin_HTC(AS, geom, cond, flow_state, choices)
    if side == 'tube':
        return pipe_HTC(AS, geom, cond, flow_state, choices)
    raise ValueError(f"TFHE: cond['side'] (or geom['side']) = {side!r}: expected 'fin' or 'tube'.")