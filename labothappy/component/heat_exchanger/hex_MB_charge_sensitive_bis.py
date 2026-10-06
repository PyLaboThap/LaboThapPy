"""
Moving-boundary heat exchanger model with charge estimation (clean version).

Based on I. Bell et al., "A generalized moving-boundary algorithm to predict
the heat transfer rate of counterflow heat exchangers for any phase
configuration", Applied Thermal Engineering, 2015.

Differences with hex_MB_charge_sensitive.py (same public interface):

* All heat transfer and pressure drop correlations are called through one
  uniform interface, the same as the pipe correlations:
      heat_transfer_coefficient(AS, geom, G, correlation=..., cond=...)
      pressure_drop(AS, geom, G, correlation=..., cond=...)
  (see toolbox/heat_exchangers/hex_MB_charge_sensitive/hex_MB_correlations.py)
* Hot and cold streams share the same code through a small `_Side` object.
* The residual function handed to the root finder depends on Q only. Pressure
  drops and multi-pass cell positions are updated in an outer fixed-point
  loop around the root finder instead of inside it.
* An oversized exchanger (residual still positive at the pinch-limited heat
  rate) returns that heat rate instead of failing.
* Distributed pressure drops are evaluated once per outer iteration instead
  of at every residual evaluation.
"""

import warnings

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
from scipy.optimize import brentq

from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.heat_connector import HeatConnector

from labothappy.correlations.heat_exchanger.f_lmtd2 import f_lmtd2, F_shell_and_tube
from labothappy.correlations.properties.two_phase import compute_two_phase_density
from labothappy.correlations.void_fraction.void_fraction import (
    compute_void_fraction, void_fraction_homogeneous, void_fraction_zivi, void_fraction_fauske,
    void_fraction_armand_treschev, void_fraction_cioncolini_thome)

from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import determine_cell_overlap
from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.compute_LMTD_multipass import determine_LMTD_multipass
from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.hex_MB_correlations import (
    CellConditions, heat_transfer_coefficient, pressure_drop, check_correlation_name)

# Bounds used by the external pinch analysis (kept from the original model)
T_HOT_OUT_MIN = 218.0         # [K]
T_COLD_OUT_MAX = 273.15 + 481  # [K]

# Map cell phase -> key of the correlation dictionaries and of user-defined values
PHASE_TO_KEY = {'liquid': '1P', 'vapor': '1P', 'transcritical': 'SC',
                'two-phase': '2P', 'vapor-wet': '2P'}
DENSITY_ONLY_VOID_FRACTION = {    # void-fraction models needing only (rho_l, rho_v, x)
    'Homogeneous': void_fraction_homogeneous, 'Zivi': void_fraction_zivi,
    'Fauske': void_fraction_fauske, 'Armand-Treschev': void_fraction_armand_treschev,
    'Cioncolini-Thome': void_fraction_cioncolini_thome}
PHASE_TO_UD = {'liquid': 'Liquid', 'vapor': 'Vapor', 'two-phase': 'Two-Phase',
               'vapor-wet': 'Vapor-wet', 'transcritical': 'Transcritical'}


class _Side:
    """Settings and state of one fluid stream ('H' = hot, 'C' = cold).

    Enthalpy-ordered vectors (hvec, pvec, Tvec, ...) go from low to high
    enthalpy on both sides: index 0 is the cold inlet / hot outlet end.
    """

    def __init__(self, name):
        self.name = name
        self.label = 'hot' if name == 'H' else 'cold'
        # settings
        self.htc_type = None      # 'Correlation' or 'User-Defined'
        self.htc_corr = {}        # {'1P': name, '2P': name, 'SC': name}
        self.htc_user = {}        # {'Liquid': value, ...}
        self.dp_corr = {}
        self.dp_user = 0.0
        # CoolProp states (re-created only if fluid or backend changes)
        self._AS_key = None
        self.AS = self.AS_sat = self.AS_heos = None
        # pressure-drop profile: cumulative fraction of DP vs. normalised
        # enthalpy measured from the side inlet (linear by default)
        self.DP = 0.0
        self.dp_s = np.array([0.0, 1.0])
        self.dp_f = np.array([0.0, 1.0])

    @property
    def two_phase_possible(self):
        return not (self.incomp or self.SC)

    def correlation(self, phase, kind):
        """Correlation name for a cell phase ('SC' falls back to '1P')."""
        table = self.htc_corr if kind == 'htc' else self.dp_corr
        key = PHASE_TO_KEY[phase]
        # A 'vapor-wet' cell has a superheated bulk: its pressure drop is single phase
        if kind == 'dp' and phase == 'vapor-wet':
            key = '1P'
        name = table.get(key) or (table.get('1P') if key == 'SC' else None)
        if name is None:
            raise ValueError(f"No '{key}' {'heat transfer' if kind == 'htc' else 'pressure drop'} "
                             f"correlation given for the {self.label} side (cell phase: {phase}).")
        return name


class HexMBChargeSensitive(BaseComponent):
    """
    **Component**: Heat exchanger (moving boundary, charge sensitive)

    **Model**: Moving-boundary algorithm of Bell et al. (2015). The heat
    exchanger is split in cells bounded by the discretisation nodes and by
    the phase-change points of both fluids. For a heat rate Q, each cell k
    needs a fraction w_k = UA_required,k / UA_available,k of the total
    area. The solution is the Q for which sum(w) = 1. Pressure drops and the
    fluid charge (with a void-fraction model in two-phase cells) follow from
    the converged cells.

    **Geometries**: 'Plate', 'Shell&Tube', 'Tube&Fins', 'PCHE'.

    **Assumptions**:
        - Steady state, no heat loss to the ambient.
        - Pressure drop and volume of a cell proportional to its area fraction.
        - Uniform wall temperature per cell (mean of the four cell temperatures).

    **Connectors**: su_H, su_C, ex_H, ex_C (MassConnector), Q (HeatConnector).

    **Parameters**:
        Flow_Type : 'CounterFlow', 'CrossFlow', 'Shell&Tube' or 'ParallelFlow'
        htc_type  : set by set_htc ('Correlation' or 'User-Defined')
        n_disc    : minimum number of cells
        Geometry  : depends on HTX_Type (see get_required_parameters)
        Optional  : n_series, n_parallel (default 1), AS_Type ('HEOS' or
                    tabular 'BICUBIC&HEOS' by default), roughness [m],
                    void_fraction_model (default 'Zivi'),
                    Q_rtol (1e-6), DP_rtol (1e-3), max_outer_iter (30)

    **Inputs**: P_su_H, T_su_H (or h_su_H), m_dot_H, fluid_H and the same for C.

    **Outputs**: Q_dot [W], ex_H / ex_C states, DP_h / DP_c [Pa],
        Mvec_h / Mvec_c (charge per cell) [kg], M_h / M_c (total charge) [kg].
    """

    HTX_TYPES = ('Plate', 'Shell&Tube', 'Tube&Fins', 'PCHE')

    def __init__(self, HTX_Type):
        super().__init__()
        if HTX_Type not in self.HTX_TYPES:
            raise ValueError(f"Heat exchanger types implemented for this model are: {self.HTX_TYPES}.")
        self.HTX_Type = HTX_Type

        self.su_H = MassConnector()
        self.su_C = MassConnector()
        self.ex_H = MassConnector()
        self.ex_C = MassConnector()
        self.Q = HeatConnector()

        self.H = _Side('H')
        self.C = _Side('C')
        self.params['DP_type'] = None
        self.eval = 0

    # =========================================================================
    # INPUTS, PARAMETERS AND CORRELATION SETTINGS
    # =========================================================================

    def get_required_inputs(self):
        return ['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']

    def get_required_parameters(self):
        general = ['Flow_Type', 'htc_type', 'n_disc']
        corr_1p = {self.H.htc_corr.get('1P'), self.C.htc_corr.get('1P')}

        if self.HTX_Type == 'Plate':
            geometry = ['A_c', 'A_h', 'h', 'l', 'l_v',
                        'C_CS', 'C_Dh', 'C_V_tot', 'C_canal_t', 'C_n_canals',
                        'H_CS', 'H_Dh', 'H_V_tot', 'H_canal_t', 'H_n_canals',
                        'casing_t', 'chevron_angle', 'fooling',
                        'n_plates', 'plate_cond', 'plate_pitch_co', 't_plates', 'w']
        elif self.HTX_Type == 'Shell&Tube':
            geometry = ['Baffle_cut', 'Shell_ID', 'Tube_L', 'Tube_OD', 'Tube_pass', 'Tube_t',
                        'central_spacing', 'foul_s', 'foul_t', 'n_series', 'n_parallel', 'n_tubes',
                        'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']
            if "Shell_Bell_Delaware_HTC" in corr_1p:
                geometry += ['D_OTL', 'N_strips', 'Tubesheet_t', 'clear_BS', 'clear_TB',
                             'inlet_spacing', 'outlet_spacing']
        elif self.HTX_Type == 'Tube&Fins':
            geometry = ['A_flow', 'Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type', 'Finned_tube_flag',
                        'Tube_L', 'Tube_OD', 'Tube_cond', 'Tube_t', 'fouling', 'h', 'k_fin',
                        'Tube_pass', 'n_rows', 'n_series', 'n_parallel', 'n_tubes', 'pitch',
                        'pitch_ratio', 'tube_arrang', 'w', 'Fin_Side']
        else:  # PCHE
            geometry = ['alpha', 'D_c', 'H_V_tot', 'C_V_tot', 'k_cond', 'L_c', 'N_c', 'N_p', 'R_p',
                        't_2', 't_3']
        return general + geometry

    def set_htc(self, htc_type="Correlation", Corr_H=None, Corr_C=None, UD_H_HTC=None, UD_C_HTC=None):
        """
        htc_type : 'User-Defined' (constant values per phase, given in UD_H_HTC /
                   UD_C_HTC with keys 'Liquid', 'Vapor', 'Two-Phase', 'Vapor-wet',
                   'Transcritical') or anything else = correlations, given in
                   Corr_H / Corr_C with keys '1P', '2P' and optionally 'SC'.
        """
        self.params['htc_type'] = htc_type
        self.Corr_H, self.Corr_C = Corr_H, Corr_C
        self.UD_H_HTC, self.UD_C_HTC = UD_H_HTC, UD_C_HTC

        for side, corr, ud in ((self.H, Corr_H, UD_H_HTC), (self.C, Corr_C, UD_C_HTC)):
            if htc_type == "User-Defined":
                side.htc_type = "User-Defined"
                side.htc_user = dict(ud)
            else:
                side.htc_type = "Correlation"
                side.htc_corr = {k: v for k, v in dict(corr).items() if v is not None}
                for name in side.htc_corr.values():
                    check_correlation_name(name, "htc")

    def set_DP(self, DP_type=None, Corr_H=None, Corr_C=None, UD_H_DP=None, UD_C_DP=None):
        """
        DP_type : None (no pressure drop), 'User-Defined' (UD_H_DP / UD_C_DP in Pa),
                  'Correlation_Global' (one evaluation at the inlet state) or
                  'Correlation_Disc' (cell by cell, iterated with the heat balance).
        """
        allowed = (None, "User-Defined", "Correlation_Global", "Correlation_Disc")
        if DP_type not in allowed:
            raise ValueError(f"DP_type must be one of {allowed}.")
        self.params['DP_type'] = DP_type
        for side, corr, ud in ((self.H, Corr_H, UD_H_DP), (self.C, Corr_C, UD_C_DP)):
            side.dp_user = float(ud or 0.0)
            side.dp_corr = {}
            if DP_type in ("Correlation_Global", "Correlation_Disc"):
                side.dp_corr = {k: v for k, v in dict(corr).items() if v is not None}
                for name in side.dp_corr.values():
                    check_correlation_name(name, "dp")

    # =========================================================================
    # SET-UP: GEOMETRY AND FLUIDS
    # =========================================================================

    def _setup_geometry(self):
        """Side geometry dictionaries, mass fluxes, areas, volumes and the
        constant terms of the overall conductance UA."""
        p = self.params
        p.setdefault('n_series', 1)
        p.setdefault('n_parallel', 1)
        n_s, n_p = p['n_series'], p['n_parallel']
        K = p.get('roughness', p.get('Tube_roughness', 0.0))
        m_h, m_c = self.su_H.m_dot, self.su_C.m_dot
        H, C = self.H, self.C

        self._UA_fh, self._UA_fc, self._R_extra = 1.0, 1.0, 0.0  # UA = 1/(fh/(ah Ah) + fc/(ac Ac) + R)

        if self.HTX_Type == 'Plate':
            for side, pre, m, A in ((H, 'H', m_h, p['A_h']), (C, 'C', m_c, p['A_c'])):
                side.geom = {**p, 'D': p[pre + '_Dh'], 'L': p['l'], 'K': K, 'theta': 0.0, 'A': A,
                             'n_canals': p[pre + '_n_canals'], 'canal_t': p[pre + '_canal_t']}
                side.G = (m / p[pre + '_n_canals']) / p[pre + '_CS']
                side.V = p[pre + '_V_tot']
            self._UA_fh = 1.0 + p['fooling']

        elif self.HTX_Type in ('Shell&Tube', 'Tube&Fins'):
            D_in = p['Tube_OD'] - 2 * p['Tube_t']
            A_in_one_tube = np.pi * D_in ** 2 / 4
            L_flow = p['Tube_L'] * p['Tube_pass'] * n_s
            G_of = lambda m: (p['Tube_pass'] / n_p) * m / (A_in_one_tube * p['n_tubes'])

            if self.HTX_Type == 'Shell&Tube':
                n_tot = p['n_tubes'] * n_s * n_p
                A_in = np.pi * D_in * p['Tube_L'] * n_tot
                A_out = np.pi * p['Tube_OD'] * p['Tube_L'] * n_tot
                V_tube = A_in_one_tube * p['Tube_L'] * n_tot
                V_shell = (np.pi / 4 * p['Shell_ID'] ** 2 - np.pi / 4 * p['Tube_OD'] ** 2 * p['n_tubes']) \
                    * p['Tube_L'] * n_s * n_p
                if 'central_spacing' in p:
                    p['cross_passes'] = np.round(p['Tube_L'] / p['central_spacing']) - 1
                p['A_eff'], p['T_V_tot'], p['S_V_tot'] = A_out, V_tube, V_shell
                R_fs = p['foul_s'] / A_out if p.get('foul_s') else 0.0
                R_ft = p['foul_t'] / A_in if p.get('foul_t') else 0.0
                self._R_extra = R_fs + R_ft
                outer_is_hot = p['Shell_Side'] == 'H'
                A_outer, A_inner, V_outer, V_inner = A_out, A_in, V_shell, V_tube
            else:
                self._setup_fin_geometry(D_in)
                A_outer, A_inner = p['A_out_tot'], p['A_in_tot']
                V_outer, V_inner = p['B_V_tot'], p['T_V_tot']
                R_cond = np.log(p['Tube_OD'] / D_in) / (2 * np.pi * p['Tube_cond'] * p['Tube_L']
                                                        * p['n_tubes'] * n_s * n_p)
                self._R_extra = R_cond
                outer_is_hot = p['Fin_Side'] == 'H'

            for side, m, is_outer in ((H, m_h, outer_is_hot), (C, m_c, not outer_is_hot)):
                side.geom = {**p, 'D': D_in, 'L': L_flow, 'K': K, 'theta': 0.0,
                             'A': A_outer if is_outer else A_inner}
                side.G = G_of(m)
                side.V = V_outer if is_outer else V_inner

        else:  # PCHE
            D_h = np.pi * p['D_c'] / (2 + np.pi)
            A_channel = np.pi * p['D_c'] ** 2 / 8
            N = p['N_c'] * p['N_p']
            frac = {'H': p['R_p'] / (1 + p['R_p']), 'C': 1 / (1 + p['R_p'])}
            for side, m, pre in ((H, m_h, 'H'), (C, m_c, 'C')):
                A = frac[pre] * N * (np.pi / 2) * p['D_c'] * p['L_c'] * n_s * n_p
                side.geom = {**p, 'D': D_h, 'L': p['L_c'] * n_s, 'K': K, 'theta': 0.0, 'A': A}
                side.G = m / (A_channel * N * frac[pre]) / n_p
                side.V = p[pre + '_V_tot']

        self.A_h, self.A_c = H.geom['A'], C.geom['A']
        self.G_h, self.G_c = H.G, C.G

    def _setup_fin_geometry(self, D_in):
        p = self.params
        p['pitch'] = p['pitch_V'] = p['pitch_H'] = p['pitch_ratio'] * p['Tube_OD']
        if p['Fin_type'] == "Annular":
            p['N_fins'] = p['Tube_L'] * p['Fin_per_m'] - 1
            r_f, r_t = p['Fin_OD'] / 2, p['Tube_OD'] / 2
            A_fin_tip = 2 * np.pi * r_f * p['Fin_t'] * p['N_fins'] * p['n_tubes']
            A_bare_tube = 2 * np.pi * r_t * p['n_tubes'] * (p['Tube_L'] - p['N_fins'] * p['Fin_t'])
            A_fin_faces = 2 * np.pi * (r_f ** 2 - r_t ** 2) * p['N_fins'] * p['n_tubes']
            p['A_out_tot'] = A_fin_tip + A_bare_tube + A_fin_faces
        elif p['Fin_type'] == "Square":
            p['N_fins'] = p['Tube_L'] * p['Fin_per_m']
            p['Fin_spacing'] = (p['Tube_L'] - p['N_fins'] * p['Fin_t']) / p['N_fins']
            A_r = 2 * (p['Fin_OD'] ** 2 - 0.785 * p['Tube_OD'] ** 2 + 2 * p['Fin_OD'] * p['Fin_t']) \
                * (p['Tube_L'] / p['Fin_spacing']) * p['n_tubes']
            L_t = p['Tube_L'] - p['N_fins'] * p['Fin_t']
            A_t = np.pi * p['Tube_OD'] * (p['Tube_L'] * (1 - p['Fin_t'] / p['Fin_spacing']) * p['n_tubes'] + L_t)
            p['A_out_tot'] = A_r + A_t
        else:
            raise ValueError("Fin geometry is not 'Annular' nor 'Square'")
        p['A_in_tot'] = np.pi * D_in * p['Tube_L'] * p['n_tubes']
        p['B_V_tot'] = p['Tube_L'] * p['w'] * p['h']
        p['T_V_tot'] = np.pi * D_in ** 2 / 4 * p['n_tubes'] * p['Tube_pass'] * p['Tube_L']

    def _setup_fluids(self):
        """CoolProp states, inlet states, supercritical flags and saturation data."""
        for side, su in ((self.H, self.su_H), (self.C, self.su_C)):
            side.fluid = su.fluid
            side.incomp = su.AS.backend_name() == 'IncompressibleBackend'
            if side.incomp:
                backend = "INCOMP"
            else:
                backend = "HEOS" if self.params.get('AS_Type') == 'HEOS' else "BICUBIC&HEOS"
            if side._AS_key != (backend, side.fluid):
                side.AS = CP.AbstractState(backend, side.fluid)
                side.AS_sat = CP.AbstractState(backend, side.fluid)
                side.AS_heos = None if side.incomp else CP.AbstractState("HEOS", side.fluid)
                side._AS_key = (backend, side.fluid)

            side.m_dot, side.h_su, side.p_su = su.m_dot, su.h, su.p
            AS = self._state(side, side.h_su, side.p_su)
            side.T_su = AS.T()
            side.SC = (not side.incomp) and side.p_su >= AS.p_critical()
            if side.two_phase_possible:
                side.AS_sat.update(CP.PQ_INPUTS, side.p_su, 0)
                side.h_bub_su = side.AS_sat.hmass()
                side.AS_sat.update(CP.PQ_INPUTS, side.p_su, 1)
                side.h_dew_su = side.AS_sat.hmass()
                side.x_su = (side.h_su - side.h_bub_su) / (side.h_dew_su - side.h_bub_su)
            else:
                side.x_su = np.nan

        # Names kept from the original model
        H, C = self.H, self.C
        self.AS_H, self.AS_C = H.AS, C.AS
        self.mdot_h, self.h_hi, self.p_hi, self.T_hi = H.m_dot, H.h_su, H.p_su, H.T_su
        self.mdot_c, self.h_ci, self.p_ci, self.T_ci = C.m_dot, C.h_su, C.p_su, C.T_su
        self.SC_h, self.SC_c = H.SC, C.SC
        self.h_incomp_flag, self.c_incomp_flag = int(H.incomp), int(C.incomp)

    def _state(self, side, h, p):
        """Return an AbstractState set at (h, p). Falls back on the full EOS
        when the tabular backend fails (typically very close to saturation)."""
        try:
            side.AS.update(CP.HmassP_INPUTS, h, p)
            return side.AS
        except ValueError:
            if side.AS_heos is None:
                raise
            side.AS_heos.update(CP.HmassP_INPUTS, h, p)
            return side.AS_heos

    # =========================================================================
    # CELL BOUNDARIES
    # =========================================================================

    def _pressure(self, side, s):
        """Pressure at normalised enthalpy s measured from the side inlet."""
        if side.DP == 0.0:
            return np.full_like(np.asarray(s, dtype=float), side.p_su)
        return side.p_su - side.DP * np.interp(s, side.dp_s, side.dp_f)

    def _local_saturation_enthalpies(self, side, h_a, h_b):
        """Bubble and dew enthalpies at the local pressure (pressure drop
        included), strictly inside the enthalpy range ]h_a, h_b[ of the side."""
        if not side.two_phase_possible:
            return []
        out = []
        for Q_flag in (0, 1):
            h_sat = side.h_bub_su if Q_flag == 0 else side.h_dew_su
            if side.DP > 0.0:
                for _ in range(4):   # h_sat = h_sat(p(h_sat)): fixed point, converges in 2-3 steps
                    s = (h_sat - h_a) / (h_b - h_a) if side.name == 'C' else (h_b - h_sat) / (h_b - h_a)
                    side.AS_sat.update(CP.PQ_INPUTS, float(self._pressure(side, min(1.0, max(0.0, s)))), Q_flag)
                    h_sat = side.AS_sat.hmass()
            if h_a < h_sat < h_b:
                out.append(h_sat)
        return out

    def calculate_cell_boundaries(self, Q):
        """Enthalpy, pressure, temperature and quality at every cell boundary for
        the heat rate Q. Nodes are: n_disc+1 evenly spaced nodes, the phase-change
        points of both fluids, and their images on the other fluid (energy balance)."""
        H, C = self.H, self.C
        self.h_co = C.h_su + Q / C.m_dot
        self.h_ho = H.h_su - Q / H.m_dot

        n_nodes = max(int(self.params['n_disc']), 1) + 1
        q_regular = np.linspace(0.0, Q, n_nodes)
        q_sat = [C.m_dot * (h - C.h_su) for h in self._local_saturation_enthalpies(C, C.h_su, self.h_co)]
        q_sat += [H.m_dot * (h - self.h_ho) for h in self._local_saturation_enthalpies(H, self.h_ho, H.h_su)]

        # Only exact duplicates are merged: dropping regular nodes that are merely
        # close to a phase-change node changes the number of cells abruptly and
        # makes the residual jump; a tiny cell contributes a tiny w instead.
        if q_sat:
            tol = 1e-9 * Q
            keep = np.ones(n_nodes, dtype=bool)
            keep[1:-1] = np.min(np.abs(q_regular[1:-1, None] - np.array(q_sat)[None, :]), axis=1) > tol
            q_nodes = np.sort(np.concatenate([q_regular[keep], q_sat]))
        else:
            q_nodes = q_regular

        self.hvec_c = C.h_su + q_nodes / C.m_dot
        self.hvec_h = self.h_ho + q_nodes / H.m_dot
        self.hnorm_c = q_nodes / Q
        self.hnorm_h = self.hnorm_c.copy()
        self.Qvec_c = C.m_dot * np.diff(self.hvec_c)
        self.Qvec_h = H.m_dot * np.diff(self.hvec_h)

        # Pressures (cold inlet at s=0 / index 0, hot inlet at s=0 / index -1)
        self.pvec_c = self._pressure(C, self.hnorm_c)
        self.pvec_h = self._pressure(H, 1.0 - self.hnorm_h)

        for side, hvec, pvec in ((C, self.hvec_c, self.pvec_c), (H, self.hvec_h, self.pvec_h)):
            self._boundary_states(side, hvec, pvec)
        self.Tvec_c, self.svec_c, self.x_vec_c = C.Tvec, C.svec, C.x_flag
        self.Tvec_h, self.svec_h, self.x_vec_h = H.Tvec, H.svec, H.x_flag
        self.Tvec_sat_pure_c, self.Tvec_sat_pure_h = C.Tsat, H.Tsat

        self.DT_pinch = np.min(self.Tvec_h - self.Tvec_c)
        self.DT_ho_ci = self.Tvec_h[0] - self.Tvec_c[0]
        self.DT_hi_co = self.Tvec_h[-1] - self.Tvec_c[-1]

    def _boundary_states(self, side, hvec, pvec):
        n = len(hvec)
        side.hvec, side.pvec = hvec, pvec
        side.Tvec, side.svec = np.empty(n), np.empty(n)
        for i in range(n):
            AS = self._state(side, hvec[i], pvec[i])
            side.Tvec[i], side.svec[i] = AS.T(), AS.smass()

        side.Tsat = np.zeros(n)
        if not side.two_phase_possible:
            side.x = np.full(n, np.nan)
            side.x_flag = np.full(n, 3.0)
            return
        h_l, h_v, T_l, T_v = np.empty(n), np.empty(n), np.empty(n), np.empty(n)
        constant_p = np.all(pvec == pvec[0])
        for i in range(n):
            if constant_p and i > 0:
                h_l[i], h_v[i], T_l[i], T_v[i] = h_l[0], h_v[0], T_l[0], T_v[0]
                continue
            side.AS_sat.update(CP.PQ_INPUTS, pvec[i], 0)
            h_l[i], T_l[i] = side.AS_sat.hmass(), side.AS_sat.T()
            side.AS_sat.update(CP.PQ_INPUTS, pvec[i], 1)
            h_v[i], T_v[i] = side.AS_sat.hmass(), side.AS_sat.T()
        side.h_l, side.h_v = h_l, h_v
        side.Tsat = 0.5 * (T_l + T_v)
        side.x = (hvec - h_l) / (h_v - h_l)
        # Original convention: -2 subcooled, 2 superheated, 3 supercritical/incompressible
        side.x_flag = np.where(side.x < 0, -2.0, np.where(side.x > 1, 2.0, side.x))

    # =========================================================================
    # MAXIMUM HEAT RATE (PINCH ANALYSIS)
    # =========================================================================

    def _compute_Qmax(self):
        """External pinch (outlet at the other inlet temperature), then internal
        pinch at the phase-change points, then a final no-crossing check."""
        H, C = self.H, self.C
        p_ho, p_co = H.p_su - H.DP, C.p_su - C.DP

        def h_at(side, p, T):
            for dT in (0.0, 1e-3, -1e-3):
                try:
                    side.AS.update(CP.PT_INPUTS, p, T + dT)
                    return side.AS.hmass()
                except ValueError:
                    continue
            side.AS_heos.update(CP.PT_INPUTS, p, T)
            return side.AS_heos.hmass()

        Qmax_h = H.m_dot * (H.h_su - h_at(H, p_ho, max(C.T_su, T_HOT_OUT_MIN)))
        Qmax_c = C.m_dot * (h_at(C, p_co, min(H.T_su, T_COLD_OUT_MAX)) - C.h_su)
        Qmax = self.Qmax_ext = min(Qmax_h, Qmax_c)
        self.Qmax_int = None

        # Internal pinch at the hot dew point / cold bubble point (Bell et al. eqs. 9-16)
        self.calculate_cell_boundaries(Qmax)
        crossing = self.Tvec_c > self.Tvec_h + 1e-9
        if H.two_phase_possible and np.any(crossing):
            i = np.where(crossing & np.isclose(self.x_vec_h, 1.0, atol=1e-9))[0]
            if i.size:
                i = i[0]
                h_c_pinch = h_at(C, self.pvec_c[i], self.Tvec_h[i])
                Qmax = self.Qmax_int = C.m_dot * (h_c_pinch - C.h_su) + H.m_dot * (H.h_su - self.hvec_h[i])
                self.calculate_cell_boundaries(Qmax)
        crossing = self.Tvec_c > self.Tvec_h + 1e-9
        if C.two_phase_possible and np.any(crossing):
            i = np.where(crossing & np.isclose(self.x_vec_c, 0.0, atol=1e-9))[0]
            if i.size:
                i = i[0]
                h_h_pinch = h_at(H, self.pvec_h[i], self.Tvec_c[i])
                Qmax = self.Qmax_int = C.m_dot * (self.hvec_c[i] - C.h_su) + H.m_dot * (H.h_su - h_h_pinch)
                self.calculate_cell_boundaries(Qmax)

        # Safety net: a temperature crossing can remain (e.g. single-phase pinch
        # caused by a varying cp). The smallest temperature difference is
        # continuous in Q, so its zero is found with Brent's method.
        if self.DT_pinch < -1e-3:
            def dT_min(Q):
                self.calculate_cell_boundaries(Q)
                return self.DT_pinch
            Qmax = self.Qmax_int = brentq(dT_min, 1e-6 * Qmax, Qmax, xtol=1e-7 * Qmax)

        self.Qmax = Qmax
        return Qmax

    # =========================================================================
    # RESIDUAL FUNCTION  (depends on Q only)
    # =========================================================================

    def objective_function(self, Q):
        """Return 1 - sum(w) for the heat rate Q (w_k: area fraction of cell k)."""
        self.calculate_cell_boundaries(Q)
        self._cell_states()
        self._correction_factors()
        self._compute_htc()
        self._compute_UA()

        if self._multipass:
            self.w = self._w_geom_for_grid()      # cell positions used by the multi-pass routine
            self.LMTD_matrix, self.LMTD = determine_LMTD_multipass(self)
        self.UA_req = self.Qvec_h / (self.F * self.LMTD)
        self.w = self.UA_req / self.UA_avail
        self.w_sum = float(np.sum(self.w))
        self.w_cumsum = np.cumsum(self.w)
        self.eval += 1
        self._last_Q = Q
        return 1.0 - self.w_sum

    def _cell_states(self):
        """Phase, mean temperatures, pressures and qualities of every cell."""
        H, C = self.H, self.C
        Thi, Tho = self.Tvec_h[1:], self.Tvec_h[:-1]
        Tci, Tco = self.Tvec_c[:-1], self.Tvec_c[1:]
        self.n_cells = n = len(Thi)
        T_wall = 0.25 * (Thi + Tho + Tci + Tco)

        for side, T_in, T_out in ((H, Thi, Tho), (C, Tci, Tco)):
            side.h_mean = 0.5 * (side.hvec[1:] + side.hvec[:-1])
            side.p_mean = 0.5 * (side.pvec[1:] + side.pvec[:-1])
            side.T_mean = 0.5 * (T_in + T_out)
            side.T_wall = T_wall.copy()
            side.Tsat_mean = 0.5 * (side.Tsat[1:] + side.Tsat[:-1])
            if side.incomp:
                side.phases = np.full(n, 'liquid', dtype=object)
                side.x_mean = np.full(n, np.nan)
            elif side.SC:
                side.phases = np.full(n, 'transcritical', dtype=object)
                side.x_mean = np.full(n, np.nan)
            else:
                x_cell = (side.h_mean - 0.5 * (side.h_l[1:] + side.h_l[:-1])) \
                    / (0.5 * (side.h_v[1:] + side.h_v[:-1]) - 0.5 * (side.h_l[1:] + side.h_l[:-1]))
                side.phases = np.where(x_cell < 0, 'liquid', np.where(x_cell > 1, 'vapor', 'two-phase')).astype(object)
                xc = np.clip(side.x, 0.0, 1.0)
                side.x_mean = 0.5 * (xc[1:] + xc[:-1])

        # Hot vapour cooled by a wall below its saturation temperature: condensation
        # on the wall ('vapor-wet'); the wall temperature seen by the fluid is T_sat.
        if H.two_phase_possible:
            vap = H.phases == 'vapor'
            H.T_wall[vap] = np.maximum(T_wall[vap], H.Tsat_mean[vap])
            wet = vap & (T_wall <= H.Tsat_mean)
            H.phases[wet] = 'vapor-wet'
            H.x_mean[wet] = 1.0

        self.phases_h, self.phases_c = H.phases, C.phases

        # Counter-flow LMTD of every cell (also used by the Han boiling correlation)
        DTA = np.maximum(Thi - Tco, 1e-6)
        DTB = np.maximum(Tho - Tci, 1e-6)
        with np.errstate(divide='ignore', invalid='ignore'):
            lmtd = (DTA - DTB) / np.log(DTA / DTB)
        self.LMTD = np.where(np.abs(DTA - DTB) < 1e-9 * np.maximum(DTA, DTB), DTA, lmtd)

    def _correction_factors(self):
        """LMTD correction factor F per cell for non counter-flow arrangements."""
        n = self.n_cells
        self.F = np.ones(n)
        flow = self.params['Flow_Type']
        if flow == "CounterFlow" or (flow == 'Shell&Tube' and self.params.get('Tube_pass', 1) == 1):
            return
        Thi, Tho = self.Tvec_h[1:], self.Tvec_h[:-1]
        Tci, Tco = self.Tvec_c[:-1], self.Tvec_c[1:]
        for k in range(n):
            dTc, dTh = Tco[k] - Tci[k], Thi[k] - Tho[k]
            if dTc <= 1e-9 or dTh <= 1e-9:
                continue                                  # one side isothermal: F = 1
            R = dTh / dTc
            P = dTc / (Thi[k] - Tci[k])
            if flow == 'Shell&Tube':
                # (the original rule "P > 0.99 -> F = 0.01" was removed: it made the
                # residual jump by orders of magnitude where the exact F is still ~1)
                if R > 10:
                    F = 1.0
                elif self.params['Tube_pass'] % 2 == 0:
                    F = F_shell_and_tube(R, P, self.params['n_series'])
                    if not (np.isfinite(F) and F > 0):      # outside the formula's range
                        F = self._f_lmtd2(R, P)
                else:
                    F = self._f_lmtd2(R, P)
            else:
                F = self._f_lmtd2(R, P)
            self.F[k] = max(F, 0.0)

    def _f_lmtd2(self, R, P):
        """f_lmtd2 with two guards: C_r consistent with R (C_min/C_max of the
        cell), and R snapped to 1 when it is within round-off of 1 (the general
        formula of f_lmtd2 is 0/0 there and returns F = 0)."""
        if abs(R - 1.0) < 1e-6:
            R = 1.0
        F = f_lmtd2(R, P, self.params, min(R, 1.0 / R))
        return F if (np.isfinite(F) and 0.0 < F <= 1.0 + 1e-9) else 1.0

    def _cell_conditions(self, side, k, alpha_other=np.nan):
        n = self.n_cells
        return CellConditions(
            fluid=side.fluid, m_dot=side.m_dot, p=side.p_mean[k], h=side.h_mean[k],
            T=side.T_mean[k], T_wall=side.T_wall[k], x=side.x_mean[k], T_sat=side.Tsat_mean[k],
            h_in=side.hvec[k], h_out=side.hvec[k + 1], p_su=side.p_su,
            Q=self.Qvec_h[k], q=self.Qvec_h[k] * n / side.geom['A'],
            DT_lm=self.F[k] * self.LMTD[k], alpha_other=alpha_other)

    def _compute_htc(self):
        """Heat transfer coefficient of every cell, hot side first (the cold
        side may need the hot-side value, e.g. Han boiling correlation)."""
        n = self.n_cells
        for side in (self.H, self.C):
            side.alpha = np.empty(n)
            for k in range(n):
                phase = side.phases[k]
                if side.htc_type == "User-Defined":
                    side.alpha[k] = side.htc_user[PHASE_TO_UD[phase]]
                    continue
                other = self.H.alpha[k] if side is self.C else np.nan
                cond = self._cell_conditions(side, k, alpha_other=other)
                AS = self._state(side, cond.h, cond.p)
                try:
                    side.alpha[k] = heat_transfer_coefficient(
                        AS, side.geom, side.G, side.correlation(phase, 'htc'), cond)
                except Exception as e:
                    raise RuntimeError(f"{side.label} side, cell {k} ({phase}, T={cond.T:.2f} K, "
                                       f"p={cond.p:.0f} Pa): {type(e).__name__}: {e}") from e
        self.alpha_h, self.alpha_c = self.H.alpha, self.C.alpha

    def _compute_UA(self):
        """Available conductance of every cell if it covered the whole exchanger."""
        a_h, a_c = self.H.alpha, self.C.alpha
        if not self._multipass:
            self.UA_avail = 1.0 / (self._UA_fh / (a_h * self.A_h) + self._UA_fc / (a_c * self.A_c)
                                   + self._R_extra)
            return
        # Multi-pass shell & tube: every tube cell exchanges with the shell cells it
        # overlaps. Matrices are indexed [tube cell, shell cell].
        shell_hot = self.params['Shell_Side'] == 'H'
        a_tube, A_tube = (a_c, self.A_c) if shell_hot else (a_h, self.A_h)
        a_shell, A_shell = (a_h, self.A_h) if shell_hot else (a_c, self.A_c)
        self.overlap_matrix = determine_cell_overlap(self._w_geom_for_grid(), self.params['Tube_pass'],
                                                     self.params['Shell_Side'])
        UA_pair = 1.0 / (1.0 / (a_tube[:, None] * A_tube) + 1.0 / (a_shell[None, :] * A_shell)
                         + self._R_extra)
        self.UA_matrix = UA_pair * self.overlap_matrix
        self.UA_avail = self.UA_matrix.sum(axis=1) / self.overlap_matrix.sum(axis=1)

    # ---- multi-pass: cell positions used for the tube/shell overlap -----------

    @property
    def _multipass(self):
        return self.HTX_Type == 'Shell&Tube' and self.params.get('Tube_pass', 1) > 1

    def _w_geom_for_grid(self):
        """Cell positions assumed for the tube/shell overlap: length proportional
        to the heat of the cell (the original model assumed equal cells, which is
        the same thing on a regular grid).

        Note: making these positions self-consistent with the computed area
        fractions was tested (fixed point with Anderson acceleration) and gave
        LARGER errors against the exact 1-shell/2-pass eps-NTU solution (up to
        +12% instead of -1.6%/+2.9%): one vector of cell positions cannot describe
        both sides of a multi-pass exchanger. The positions are therefore kept
        fixed, which also keeps the residual a function of Q only."""
        return np.maximum(np.diff(self.hnorm_c), 1e-12)

    # =========================================================================
    # PRESSURE DROPS
    # =========================================================================

    def _init_pressure_drops(self):
        DP_type = self.params.get('DP_type')
        for side in (self.H, self.C):
            side.dp_s, side.dp_f = np.array([0.0, 1.0]), np.array([0.0, 1.0])
            if DP_type is None:
                side.DP = 0.0
            elif DP_type == "User-Defined":
                side.DP = side.dp_user
            else:
                side.DP = self._inlet_pressure_drop(side)
            side.DP = self._limit_DP(side, side.DP)

    def _inlet_pressure_drop(self, side):
        """Pressure drop of the whole side evaluated at the inlet state."""
        if side.SC:
            phase = 'transcritical'
        elif side.two_phase_possible and 0.002 < side.x_su < 0.999:
            phase = 'two-phase'
        else:
            phase = 'liquid'
        cond = CellConditions(fluid=side.fluid, m_dot=side.m_dot, p=side.p_su, h=side.h_su,
                              T=side.T_su, T_wall=0.5 * (self.H.T_su + self.C.T_su),
                              x=side.x_su, p_su=side.p_su, h_in=side.h_su, h_out=side.h_su)
        AS = self._state(side, side.h_su, side.p_su)
        return pressure_drop(AS, side.geom, side.G, side.correlation(phase, 'dp'), cond)

    def _distributed_pressure_drops(self):
        """Cell-by-cell pressure drops for the current (converged) cells.
        Returns the new totals and updates the pressure profiles."""
        frac = self.w / np.sum(self.w)          # length fraction of every cell
        new = {}
        for side in (self.H, self.C):
            DPvec = np.empty(self.n_cells)
            for k in range(self.n_cells):
                cond = self._cell_conditions(side, k)
                AS = self._state(side, cond.h, cond.p)
                DPvec[k] = frac[k] * pressure_drop(AS, side.geom, side.G,
                                                   side.correlation(side.phases[k], 'dp'), cond)
            side.DPvec = DPvec
            # cumulative profile from the side inlet (cold: index 0, hot: last index)
            s_cells = self.hnorm_c if side.name == 'C' else (1.0 - self.hnorm_h)[::-1]
            dp_from_inlet = DPvec if side.name == 'C' else DPvec[::-1]
            total = float(np.sum(DPvec))
            cum = np.concatenate([[0.0], np.cumsum(dp_from_inlet)])
            new[side.name] = (total, s_cells, cum / total if total > 0 else s_cells)
        self.DPvec_h, self.DPvec_c = self.H.DPvec, self.C.DPvec
        return new

    def _limit_DP(self, side, DP):
        if DP > 0.8 * side.p_su:
            warnings.warn(f"{side.label} side pressure drop ({DP:.0f} Pa) above 80% of the inlet "
                          f"pressure; limited to that value.", stacklevel=3)
            return 0.8 * side.p_su
        return max(DP, 0.0)

    # =========================================================================
    # SOLVE
    # =========================================================================

    def solve(self, only_external=False, and_solve=True):
        if not self.check_calculable():
            raise ValueError("Component not calculable: check the inputs.")
        self._setup_geometry()          # before the parameter check: it derives some parameters
        self.check_parametrized()
        self._setup_fluids()
        H, C = self.H, self.C
        if not (H.T_su - C.T_su > 1e-2 and H.m_dot > 0 and C.m_dot > 0):
            raise ValueError("Hot and cold temperatures seem to be reversed or a flow rate is not positive.")

        rtol_dp = self.params.get('DP_rtol', 1e-3)
        max_outer = int(self.params.get('max_outer_iter', 30))
        disc = self.params.get('DP_type') == "Correlation_Disc"
        self._init_pressure_drops()
        self.eval = 0
        Q = None
        self.outer_iterations = 0

        # Outer fixed point on the distributed pressure drops (one pass otherwise):
        # solve Q with the current pressure profile, recompute the cell pressure
        # drops from the converged cells, repeat until the totals stop changing.
        for it in range(max_outer):
            self.outer_iterations = it + 1
            self._compute_Qmax()
            if only_external or not and_solve:
                self.Q_dot = self.Qmax
                return self.Qmax
            Q = self._solve_heat_rate(Q)
            if not disc:
                break
            new = self._distributed_pressure_drops()
            updates, converged = [], True
            for side in (H, C):
                DP_new, s, f = new[side.name]
                DP_new = self._limit_DP(side, DP_new)
                converged &= abs(DP_new - side.DP) <= max(1.0, rtol_dp * DP_new)
                updates.append((side, DP_new, s, f))
            if converged:
                break
            for side, DP_new, s, f in updates:
                side.DP, side.dp_s, side.dp_f = DP_new, s, f
        else:
            warnings.warn(f"Distributed pressure drops not converged after {max_outer} outer "
                          f"iterations.", stacklevel=2)

        self._finalize(Q)
        return Q

    def _solve_heat_rate(self, Q_guess=None):
        """Find Q such that sum(w) = 1, within ]0, Qmax[."""
        Qmax = self.Qmax
        rtol = self.params.get('Q_rtol', 1e-6)
        # Residual values are cached within one solve: brentq re-evaluates the
        # bracket ends, which were just computed to check the bracket.
        self._f_cache = {}

        def f(Q):
            if Q not in self._f_cache:
                self._f_cache[Q] = self.objective_function(Q)
            return self._f_cache[Q]
        self._f = f
        Q_hi = Qmax * (1.0 - 1e-5)

        # Warm start: narrow bracket around the previous solution
        if Q_guess is not None and Q_guess < Q_hi:
            a, b = 0.98 * Q_guess, min(1.02 * Q_guess, Q_hi)
            fa, fb = f(a), f(b)
            if fa > 0 > fb:
                return self._brent(a, b, rtol)

        f_hi = f(Q_hi)
        if f_hi >= 0.0:
            # The exchanger is large enough to reach the pinch-limited heat rate
            self.pinched = True
            if self._last_Q != Q_hi:
                self.objective_function(Q_hi)
            return Q_hi
        self.pinched = False
        for frac in (1e-2, 1e-4, 1e-6):
            Q_lo = frac * Qmax
            if f(Q_lo) > 0.0:
                return self._brent(Q_lo, Q_hi, rtol)
        raise RuntimeError("No heat rate found with sum(w) < 1, even at 1e-6*Qmax: "
                           "check the geometry and the heat transfer coefficients.")

    def _brent(self, a, b, rtol):
        Q, self.results = brentq(self._f, a, b, xtol=1e-9 * b, rtol=rtol, full_output=True)
        # leave the stored cell data at the solution (re-evaluate only if needed)
        res = self._f_cache[Q] if self._last_Q == Q else self.objective_function(Q)
        if not abs(res) < 1e-3:
            warnings.warn(f"Residual 1-sum(w) = {res:.3g} at the root found (Q = {Q:.6g} W): "
                          "the residual function is discontinuous there.", stacklevel=3)
        return Q

    def _finalize(self, Q):
        H, C = self.H, self.C
        self.Q_dot = Q
        self.Q.set_Q_dot(Q)
        self.epsilon_th = Q / self.Qmax
        self.residual = 1.0 - float(np.sum(self.w))
        self.DP_h = float(self.pvec_h[-1] - self.pvec_h[0])
        self.DP_c = float(self.pvec_c[0] - self.pvec_c[-1])
        self.p_ho, self.p_co = self.pvec_h[0], self.pvec_c[-1]
        if not hasattr(self, 'DPvec_h') or self.params.get('DP_type') != "Correlation_Disc":
            self.DPvec_h = np.abs(np.diff(self.pvec_h))
            self.DPvec_c = np.abs(np.diff(self.pvec_c))
        self.Avec_h = self.w * self.A_h
        self.Avec_c = self.w * self.A_c
        self._compute_charge()

        self.ex_H.reset()
        self.ex_H.set_fluid(H.fluid)
        self.ex_H.set_p(self.pvec_h[0])
        self.ex_H.set_h(self.hvec_h[0])
        self.ex_H.set_m_dot(H.m_dot)
        self.ex_C.reset()
        self.ex_C.set_fluid(C.fluid)
        self.ex_C.set_p(self.pvec_c[-1])
        self.ex_C.set_h(self.hvec_c[-1])
        self.ex_C.set_m_dot(C.m_dot)
        self.H_su, self.C_su, self.H_ex, self.C_ex = self.su_H, self.su_C, self.ex_H, self.ex_C
        self.solved = True

    def _compute_charge(self):
        """Fluid mass in every cell: cell volume times the mean of the boundary
        densities (void-fraction model in two-phase)."""
        model = self.params.get('void_fraction_model', 'Zivi')
        frac = self.w / np.sum(self.w)
        for side in (self.H, self.C):
            n = len(side.hvec)
            D, eps = np.empty(n), np.full(n, -1.0)
            for i in range(n):
                AS = self._state(side, side.hvec[i], side.pvec[i])
                x = side.x[i] if side.two_phase_possible else np.nan
                if side.two_phase_possible and 0.0 < x < 1.0:
                    side.AS_sat.update(CP.PQ_INPUTS, side.pvec[i], 0)
                    rho_l = side.AS_sat.rhomass()
                    side.AS_sat.update(CP.PQ_INPUTS, side.pvec[i], 1)
                    rho_v = side.AS_sat.rhomass()
                    if model in DENSITY_ONLY_VOID_FRACTION:
                        # only densities needed: no transport property (e.g. R1233zd(E)
                        # has no viscosity model) and no extra AbstractState
                        eps[i] = DENSITY_ONLY_VOID_FRACTION[model](rho_l, rho_v, x)
                    else:
                        A_cross = np.pi * side.geom['D'] ** 2 / 4    # so that m_dot/A_cross = G
                        eps[i] = compute_void_fraction(AS, side.geom, side.G * A_cross,
                                                       void_fraction_model=model)
                    D[i] = compute_two_phase_density(x=x, rho_l=rho_l, rho_g=rho_v, alpha=eps[i])
                else:
                    D[i] = AS.rhomass()
            side.Dvec, side.eps_void = D, eps
            side.Vvec = side.V * frac
            side.Mvec = side.Vvec * 0.5 * (D[1:] + D[:-1])
        H, C = self.H, self.C
        self.Dvec_h, self.Dvec_c = H.Dvec, C.Dvec
        self.eps_void_h, self.eps_void_c = H.eps_void, C.eps_void
        self.Vvec_h, self.Vvec_c = H.Vvec, C.Vvec
        self.Mvec_h, self.Mvec_c = H.Mvec, C.Mvec
        self.M_h, self.M_c = float(np.sum(H.Mvec)), float(np.sum(C.Mvec))

    # =========================================================================
    # PLOTS AND PRINTS
    # =========================================================================

    def plot_objective_function(self, N=100):
        """Plot 1 - sum(w) against Q (then restore the solved state)."""
        Q = np.linspace(1e-3 * self.Qmax, self.Qmax * (1 - 1e-5), N)
        r = np.array([self.objective_function(q) for q in Q])
        fig, ax = plt.subplots()
        ax.plot(Q, r)
        ax.axhline(0, color='k', lw=0.8)
        ax.set_xlabel('Q [W]')
        ax.set_ylabel('1 - sum(w) [-]')
        ax.grid(True)
        if getattr(self, 'solved', False):
            self.objective_function(self.Q_dot)

    def plot_ph_pair(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for side, hvec, pvec in (('H', self.hvec_h, self.pvec_h), ('C', self.hvec_c, self.pvec_c)):
                fluid = self.H.fluid if side == 'H' else self.C.fluid
                diagram = PropertyPlot(fluid, "PH", unit_system="EUR")
                plt.plot(0.001 * np.asarray(hvec), np.asarray(pvec) * 1e-5, 's-')
                diagram.calc_isolines()
                diagram.title(f"{'Hot' if side == 'H' else 'Cold'} side P-h diagram. Fluid: {fluid}")
                diagram.show()

    def plot_Ts_pair(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for side, svec, Tvec in (('H', self.svec_h, self.Tvec_h), ('C', self.svec_c, self.Tvec_c)):
                fluid = self.H.fluid if side == 'H' else self.C.fluid
                diagram = PropertyPlot(fluid, "TS", unit_system="EUR")
                plt.plot(0.001 * svec, Tvec - 273.15, 's-')
                diagram.calc_isolines()
                diagram.title(f"{'Hot' if side == 'H' else 'Cold'} side T-s diagram. Fluid: {fluid}")
                diagram.show()

    def plot_cells(self, fName='', dpi=400):
        plt.figure(figsize=(4, 3))
        plt.plot(self.hnorm_h, self.Tvec_h, 'rs-')
        plt.plot(self.hnorm_c, self.Tvec_c, 'bs-')
        plt.xlim(0, 1)
        plt.ylabel('T [K]')
        plt.xlabel(r'$\hat h$ [-]')
        plt.grid(True)
        plt.tight_layout(pad=0.2)
        if fName != '':
            plt.savefig(fName, dpi=dpi)
        plt.show()

    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        for name in ('su_C', 'su_H', 'ex_C', 'ex_H'):
            c = getattr(self, name)
            print(f"  - {name}: fluid={c.fluid}, T={c.T}, p={c.p}, m_dot={c.m_dot}")
        print("======================")
