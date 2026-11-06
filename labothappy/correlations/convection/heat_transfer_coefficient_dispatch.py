# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 16:12:55 2026

@author: basil
"""

# Adapt the import path of BPHE_HTC to where the plate dispatchers live.
from labothappy.correlations.convection.BPHX_HTC import BPHE_HTC

# # The previous registry-based function is kept for the heat exchanger types that are
# # not yet handled by the three-level scheme (HTX type -> flow state -> correlation).
# def _heat_transfer_coefficient_legacy(AS, geom, G, correlation, cond):
#     f = _resolve(correlation, HTC_CORRELATIONS, HTC_ALIASES, "heat transfer")
#     return float(f(AS, geom, G, cond))


# def _plate_cond(AS, G, cond, single_phase):
#     """Dict of cell conditions read by the BPHX correlations, built from the
#     CellConditions of the model. AS must be set at the cell mean state (h, p)."""
#     d = {'G': G, 'p': cond.p, 'x': cond.x, 'q': cond.q, 'm_dot': cond.m_dot}
#     if single_phase:
#         b = bulk_properties(AS, cond)                       # needs AS at (h, p)
#         d.update(mu=b['mu'], k=b['k'], Pr=b['Pr'],
#                  mu_w=wall_viscosity(AS, cond, b['mu']))    # moves AS: keep it last
#     return d


from dataclasses import fields

def heat_transfer_coefficient(AS, geom, cond, correlation):
    """cond is a dictionary built by the model (_cond_dict): it holds G and flow_state,
    and every field of CellConditions by name."""
    if geom.get('HTX_Type') == 'BPHE':
        state = cond.get('flow_state')                       # was: getattr(cond, 'flow_state', '')
        if state not in ('SP', 'EVAP', 'COND'):
            raise ValueError(f"cond['flow_state'] = {state!r}: expected 'SP', 'EVAP' or 'COND'.")
        choices = correlation if isinstance(correlation, dict) else {state: correlation}
        return float(BPHE_HTC(AS, geom, cond, state, choices))   # the dict goes straight through

    # cell = CellConditions(**{f.name: cond[f.name] for f in fields(CellConditions)})
    # return _heat_transfer_coefficient_legacy(AS, geom, cond['G'], correlation, cell)

