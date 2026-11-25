# -*- coding: utf-8 -*-
"""
PCHE HTC dispatcher: HX type -> flow state (SP / EVAP / COND) -> correlation, (AS, geom, cond) formalism.

geom['type_channel'] : 'Straight' (default), 'Zigzag' or 'S_fin' (S-shaped fins).

SP    : PCHE correlations  'Straight', 'Lee' (zigzag), 'Bennett_Chen' (zigzag, sCO2 only, explicit choice)
        or any pipe correlation name ('Gnielinski', 'Dittus_Boelter', 'Cheng', 'Liu', ...), evaluated
        on the equivalent tube (D = D_h, G of one channel).
        auto : Straight -> 'Straight' ; Zigzag -> 'Jiang' (since 2026-11-23; before: 'Lee').
        'Jiang' (zigzag, Jiang et al. 2018 30 / 40 deg fits, turbulent data 3000-20600,
        extrapolated with a warning below): auto for zigzag. 'Lee': explicit choice only.
        'Ngo' (zigzag fins 52 deg, sCO2, 3500 < Re < 22000, Ngo et al. 2007): explicit choice only.
        'Ngo_S_fin' (S-shaped fins 52 deg, sCO2, 3500 < Re < 23000, Ngo et al. 2007):
        auto for geom['type_channel'] = 'S_fin'.
        'Kim_No_2012' (zigzag, laminar Re < 2500, helium, angle 5-45 deg, Kim & No 2012; optional
        geom['P_zz']), 'Kim_No_2011' (zigzag 15 deg, laminar, helium / water),
        'Liu_straight' (straight, sCO2 cooling, 3600 < Re < 36500, needs cond T_w, p, cp):
        explicit choice only.
        'Kwon' (straight, nitrogen gas, 8500 < Re < 17000, Kwon et al. 2018): explicit choice only.
        Zigzag, explicit choice only (2026-11-23): 'Katz' (37 deg, sCO2 / helium, Katz 2021),
        'Kim_SG_2016' (32.5 / 40 deg, sCO2 CFD), 'Baik' (32.5 deg, CO2 15000-85000 / water
        50-200, modified diameter D_c / sqrt 2), 'Saeed' (40 deg, sCO2 Pr 2-13 / water),
        'Chen_2016' (15 deg, helium 1400-3558). Branch chosen from the fluid of AS.
EVAP / COND : PCHE two-phase correlations (explicit choice only, Kwon et al. 2018, straight
        channel Dh 1.83 mm): COND 'Kwon' (CO2 at 0.7 MPa, Akers-type), EVAP 'Kwon_film_boiling'
        (liquid nitrogen, FILM boiling, needs cond T_w).
        Otherwise, WITHOUT geom['R_curv'] (straight channel, zigzag without R_curv): straight pipe
        correlations (pipe_HTC, auto and names as pipe_HTC) on the equivalent tube
        (since 2026-11-23).
        WITH geom['R_curv']: curved pipe correlations (HCHE dispatcher) on the equivalent tube:
        'Zhao' / 'Nariai' (EVAP), 'Zakeralhoseini' (COND); a pipe correlation name falls back to the
        straight pipe_HTC.
Every non-PCHE path uses the curved pipe functions (HCHE dispatcher): SP auto name -> 'Sigalotti'.
geom['R_curv'] : radius of curvature of the channel bends [m]; absent -> straight tube (R -> infinity).

The equivalent tube is read from geom: D_t_i (= D_h), L_t (= developed channel length), n_t,
n_passes; cond['G'] is the mass flux of one channel (the pipe bank uses it first).
"""
import warnings
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import (
    PCHE_straight_HTC, PCHE_bennett_chen_HTC, PCHE_lee_zigzag_HTC, PCHE_jiang_zigzag_HTC)
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import PCHE_ngo_zigzag_HTC
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import PCHE_ngo_S_fin_HTC
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import (
    PCHE_kim_no_2012_zigzag_HTC, PCHE_kim_no_2011_zigzag_HTC, PCHE_liu_straight_HTC)
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import (
    PCHE_kwon_straight_HTC, PCHE_kwon_condensing_HTC, PCHE_kwon_film_boiling_HTC)
from labothappy.correlations.convection.PCHE_HTC_correlation_bank import (
    PCHE_katz_zigzag_HTC, PCHE_kim_sg_zigzag_HTC, PCHE_baik_zigzag_HTC, PCHE_saeed_zigzag_HTC,
    PCHE_chen_zigzag_HTC)
from labothappy.correlations.convection.curved_pipe_HTC import curved_pipe_HTC
from labothappy.correlations.convection.pipe_HTC import pipe_HTC

_PCHE_SP_HTC = {'Straight': PCHE_straight_HTC, 'Lee': PCHE_lee_zigzag_HTC, 'Bennett_Chen': PCHE_bennett_chen_HTC}
_PCHE_SP_HTC['Jiang'] = PCHE_jiang_zigzag_HTC   # zigzag turbulent, 30 / 40 deg fits (explicit choice only)
_PCHE_SP_HTC['Ngo'] = PCHE_ngo_zigzag_HTC       # zigzag fins 52 deg, sCO2, turbulent (explicit choice only)
_PCHE_SP_HTC['Ngo_S_fin'] = PCHE_ngo_S_fin_HTC  # S-shaped fins 52 deg, sCO2, turbulent (auto for 'S_fin')
_PCHE_SP_HTC['Kim_No_2012'] = PCHE_kim_no_2012_zigzag_HTC  # zigzag laminar, helium (explicit only)
_PCHE_SP_HTC['Kim_No_2011'] = PCHE_kim_no_2011_zigzag_HTC  # zigzag 15 deg laminar (explicit only)
_PCHE_SP_HTC['Liu_straight'] = PCHE_liu_straight_HTC       # straight, sCO2 cooling (explicit only)
_PCHE_SP_HTC['Kwon'] = PCHE_kwon_straight_HTC              # straight, nitrogen gas (explicit only)
# Zigzag correlations added 2026-11-23 (explicit only; fluid-dependent fits pick the branch from AS)
_PCHE_SP_HTC['Katz'] = PCHE_katz_zigzag_HTC                # 37 deg, sCO2 / helium
_PCHE_SP_HTC['Kim_SG_2016'] = PCHE_kim_sg_zigzag_HTC       # 32.5 / 40 deg, sCO2 (CFD)
_PCHE_SP_HTC['Baik'] = PCHE_baik_zigzag_HTC                # 32.5 deg, CO2 turbulent / water lam.
_PCHE_SP_HTC['Saeed'] = PCHE_saeed_zigzag_HTC              # 40 deg, sCO2 precooler / water (CFD)
_PCHE_SP_HTC['Chen_2016'] = PCHE_chen_zigzag_HTC           # 15 deg, helium, Re 1400-3558

# Two-phase PCHE correlations (explicit choice only; other names go to the curved pipe dispatcher)
_PCHE_COND_HTC = {'Kwon': PCHE_kwon_condensing_HTC}                 # CO2 at 0.7 MPa, straight
_PCHE_EVAP_HTC = {'Kwon_film_boiling': PCHE_kwon_film_boiling_HTC}  # liquid N2, film boiling

# Pipe correlations established for tubes of several mm: warn below this hydraulic diameter [m]
_PIPE_MIN_D = {'Cheng': 8e-3, 'Liu': 4e-3, 'Gungor_Winterton': 3e-3, 'Friedel': 0.0}
_CONFINEMENT_D = 3e-3     # two-phase confinement effects below ~3 mm

# Two-phase names that exist only in the curved pipe dispatcher (curved_pipe_HTC.py _HCHE)
_CURVED_2P_NAMES = {'EVAP': ('Zhao', 'Nariai'), 'COND': ('Zakeralhoseini',)}


def select_PCHE_HTC(geom, cond):
    """Straight channel -> Straight; zigzag -> Jiang (since 2026-11-23, was 'Lee': on real data,
    Katz 2021 and Alvarez 2021, the Lee laminar branch under-predicted UA by up to 50 %; Jiang is
    extrapolated with a warning below Re 3000); S-shaped fins -> Ngo_S_fin (only S-fin correlation).
    Bennett_Chen and Lee only on explicit request."""
    type_channel = str(geom.get('type_channel', 'Straight')).lower()
    if type_channel == 's_fin':
        return 'Ngo_S_fin'
    return 'Jiang' if type_channel == 'zigzag' else 'Straight'


def _equivalent_tube(geom):
    g = dict(geom)
    Dh = g.get('Dh') or g['D_t_i']
    g.update(D_t_i=Dh, Dh=Dh, n_passes=1)
    g.setdefault('L_t', g.get('L_ch', g.get('L')))
    g.setdefault('L', g['L_t'])
    if 'R_curv' not in g:
        if str(g.get('type_channel', 'Straight')).lower() == 'zigzag':
            warnings.warn("PCHE: zigzag channel without geom['R_curv']: curved pipe correlations run as a straight tube.")
        g['R_curv'] = 1e6
    return g


def _pipe(AS, geom, cond, state, name):
    g = _equivalent_tube(geom)
    if name not in (None, 'auto') and g['D_t_i'] < _PIPE_MIN_D.get(name, 0.0):
        warnings.warn(f"PCHE: pipe correlation {name!r} used with D_h = {g['D_t_i']*1e3:.2f} mm, "
                      "below the tube diameters it was established for.")
    if state != 'SP' and g['D_t_i'] < _CONFINEMENT_D:
        warnings.warn(f"PCHE: two-phase flow in a D_h = {g['D_t_i']*1e3:.2f} mm channel: confinement "
                      "effects are not in the tube correlations.")
    # Two phase without a radius of curvature (straight channel, or zigzag without R_curv):
    # STRAIGHT pipe correlations (pipe_HTC, auto as pipe_HTC). The curved two-phase correlations
    # (e.g. Zakeralhoseini, Mori-Nakayama liquid Nu) are not valid at R_curv -> infinity
    # (2026-11-23, proposed_changes E20).
    # Curved-only names given explicitly stay on the curved path.
    curved_only = _CURVED_2P_NAMES.get(state, ())
    if state != 'SP' and 'R_curv' not in geom and name not in curved_only:
        return pipe_HTC(AS, g, {**cond, 'side': 'tube'}, state, {state: name})
    if state != 'SP' and 'R_curv' not in geom:
        warnings.warn(f"PCHE: curved-pipe correlation {name!r} on a channel without geom['R_curv'] "
                      "(evaluated with R_curv = 1e6 m, outside its validity).")
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
    if flow_state == 'COND' and choice in _PCHE_COND_HTC:
        return _PCHE_COND_HTC[choice](AS, geom, cond)
    if flow_state == 'EVAP' and choice in _PCHE_EVAP_HTC:
        return _PCHE_EVAP_HTC[choice](AS, geom, cond)
    if flow_state in ('EVAP', 'COND'):
        return _pipe(AS, geom, cond, flow_state, choice)
    raise ValueError(f"PCHE HTC: flow_state {flow_state!r}: expected 'SP', 'EVAP' or 'COND'.")