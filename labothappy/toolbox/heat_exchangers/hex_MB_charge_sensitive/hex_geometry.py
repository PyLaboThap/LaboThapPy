# -*- coding: utf-8 -*-
"""
Geometry resolution of the heat exchangers of HexMBChargeSensitive, one common formalism.

For every HTX_Type a class derived from HXGeometry turns the user parameters into

    GeometryResult
        .H, .C            SideGeometry of the hot / cold side
                              .geom    dictionary handed to heat_transfer_coefficient / pressure_drop
                              .A       heat transfer area of the side [m2]
                              .A_flow  total free flow area crossed by ONE stream [m2]  ->  G = m_dot / A_flow
                              .V_ch    channel (tube, shell, ...) volume [m3]
                              .V_port  manifold / port volume [m3]
        .UA_f             (f_h, f_c): UA = 1 / (f_h/(a_h A_h) + f_c/(a_c A_c) + R_extra)
        .R_extra          wall conduction + fouling resistance [K/W], referred to the same UA
        .contact          None (ideal counter-flow, cell k of one side faces cell k of the other)
                          or a ContactMap: which cells of the two streams are separated by the same wall.

The ContactMap is the common formalism of the multi-pass shell & tube and of the spiral exchanger:

    M = contact.overlap(model)           raw overlap [row cell, column cell] (row side = contact.ref)
    LMTD_matrix, LMTD = contact.lmtd(model)
    contact.cell_radius(model, 'H'|'C')  mean radius of curvature of the cells (curved channels) or None
    contact.pinch(model, Q_hi)           heat rate limited by the internal pinch of the contact map, or None

The routines behind are those of the toolbox (cell_overlap_MBHX.py, compute_LMTD_multipass.py). They read the
model (model.w = length fraction of the cells, same on both streams, cell order of the model: index 0 =
cold inlet / hot outlet; model.UA_matrix, model.Tvec_h/c, model.params). The model then uses
UA_avail[i] = sum_j M_ij UA_pair_ij / sum_j M_ij.

Spiral exchanger (SPHE) - what is solved
----------------------------------------
(parameter names, radial_limits and theta_end are those of the previous model, read by the toolbox routines)
Two plates A and B of thickness t are wound on a core. In a radial cut, from the core to the outside:
    A | channel H (gap b_H) | B | channel C (gap b_C) | A | channel H | ...
Plate centre-lines are Archimedean spirals r_A(phi) = r0 + c phi, r_B = r_A + b_H + t, pitch per turn
p = b_H + b_C + 2 t = 2 pi c. The channels have the same angular range [0, 2 pi n_turns].
    hot channel H(phi) touches  C(phi)         through plate B   (phi in [0, Phi])
                         and    C(phi - 2 pi)  through plate A   (phi in [2 pi, Phi])
Walls without a counterpart (inner face of the first hot turn, outer face of the last cold turn) are adiabatic
and are NOT counted in A. The unknown geometry is solved from ONE of
    n_turns | L_plate (developed length of one plate) | D_ext | A (heat transfer area).
"""
from dataclasses import dataclass, field
import warnings
import numpy as np
from scipy.optimize import brentq


# =============================================================================
# RESULT CONTAINERS
# =============================================================================

@dataclass
class SideGeometry:
    geom: dict
    A: float
    A_flow: float
    V_ch: float
    V_port: float = 0.0

    @property
    def V(self):
        return self.V_ch + self.V_port


@dataclass
class GeometryResult:
    H: SideGeometry
    C: SideGeometry
    R_extra: float = 0.0
    UA_f: tuple = (1.0, 1.0)
    contact: object = None
    derived: dict = field(default_factory=dict)     # derived quantities, also written in params


# =============================================================================
# CONTACT MAPS
# =============================================================================

class ContactMap:
    """Which cells of the stream 'ref' (rows) are separated from which cells of the other stream (columns)
    by the same wall."""
    ref = 'H'

    def overlap(self, model):
        raise NotImplementedError

    def lmtd(self, model):
        """(LMTD_matrix, LMTD_vector[row cell])"""
        raise NotImplementedError

    def cell_radius(self, model, channel):
        return None

    def pinch(self, model, Q_hi):
        return None


class STHEMultipassContact(ContactMap):
    """Multi-pass shell & tube (toolbox: determine_cell_overlap_STHE, compute_LMTD_STHE). Rows = tube cells."""

    def __init__(self, n_passes, shell_side):
        self.n_passes, self.shell_side = n_passes, shell_side
        self.ref = 'C' if shell_side == 'H' else 'H'            # the tube side

    def overlap(self, model):
        from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import (
            determine_cell_overlap_STHE)
        return determine_cell_overlap_STHE(model.w, self.n_passes, self.shell_side)

    def lmtd(self, model):
        from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.compute_LMTD_multipass import (
            compute_LMTD_STHE)
        return compute_LMTD_STHE(model)


class SpiralContact(ContactMap):
    """Spiral plate exchanger (toolbox: determine_cell_overlap_SPHE, compute_LMTD_SPHE,
    determine_cell_mean_curvature_radius_SPHE, internal_pinching_SPHE). Rows = hot cells.
    Both streams are numbered from the core (flow_from = 'inner'): cell 0 = cold inlet = hot outlet at the core,
    i.e. the cold stream enters at the core and the hot one at the periphery."""
    ref = 'H'

    def overlap(self, model):
        from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import (
            determine_cell_overlap_SPHE)
        return determine_cell_overlap_SPHE(model)

    def lmtd(self, model):
        from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.compute_LMTD_multipass import (
            compute_LMTD_SPHE)
        LMTD_matrix, LMTD_hot, LMTD_cold, UA = compute_LMTD_SPHE(model)
        return LMTD_matrix, np.maximum(LMTD_hot, 1e-6)

    def cell_radius(self, model, channel):
        from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.cell_overlap_MBHX import (
            determine_cell_mean_curvature_radius_SPHE)
        return determine_cell_mean_curvature_radius_SPHE(model, channel='hot' if channel == 'H' else 'cold')[1]

    def pinch(self, model, Q_hi):
        from labothappy.toolbox.heat_exchangers.hex_MB_charge_sensitive.compute_LMTD_multipass import (
            internal_pinching_SPHE)
        return internal_pinching_SPHE(model, stream='hot', Q_hi=Q_hi)


# =============================================================================
# SPIRAL HELPERS
# =============================================================================

def spiral_arc(a, c, phi):
    """Arc length of the Archimedean spiral r = a + c phi from 0 to phi (closed form)."""
    F = lambda r: (r * np.sqrt(r * r + c * c) + c * c * np.arcsinh(r / c)) / (2.0 * c)
    return F(a + c * np.asarray(phi, float)) - F(a)


def spiral_geometry(b_H, b_C, t, H, r0_in, inner_ch, N):
    """Double spiral of the model. Radial offsets from the spiral r0_in + a*theta (a = pitch / 2 pi), at one angle:
        sheet1 (-t/2, t/2) | inner channel | sheet2 | outer channel | (sheet1 of the next turn at +pitch)
    The INNER channel runs over theta in [0, 2 pi N], the OUTER one over [0, 2 pi N - 2 pi].
        sheet2 separates inner(theta) from outer(theta)            theta in [0, Theta - 2 pi]
        sheet1 separates inner(theta) from outer(theta - 2 pi)     theta in [2 pi, Theta]
    Returns a dict (lengths of channels and plates, exchange area of ONE unit, outer diameter)."""
    p = b_H + b_C + 2.0 * t
    a = p / (2.0 * np.pi)
    b_in, b_out = (b_H, b_C) if inner_ch == 'hot' else (b_C, b_H)
    off = {'sheet1': 0.0, 'inner': 0.5 * t + 0.5 * b_in, 'sheet2': t + b_in, 'outer': 1.5 * t + b_in + 0.5 * b_out}
    Theta = 2.0 * np.pi * N
    end = {'sheet1': Theta, 'inner': Theta, 'sheet2': Theta, 'outer': Theta - 2.0 * np.pi}
    arc = lambda el, th0=0.0, th1=None: (
        0.0 if (end[el] if th1 is None else th1) <= th0 else
        float(spiral_arc(r0_in + off[el] + a * th0, a, (end[el] if th1 is None else th1) - th0)))
    L_in, L_out = arc('inner'), arc('outer')
    L_c1 = arc('sheet1', 2.0 * np.pi, Theta)                  # plate sheet1 in contact
    L_c2 = arc('sheet2', 0.0, Theta - 2.0 * np.pi)            # plate sheet2 in contact
    return dict(pitch=p, a=a, b_in=b_in, b_out=b_out, off=off, Theta=Theta, end=end, r0_in=r0_in,
                L_inner=L_in, L_outer=L_out, L_sheet1=arc('sheet1'), L_sheet2=arc('sheet2'),
                L_contact_1=L_c1, L_contact_2=L_c2, A_exchange=H * (L_c1 + L_c2), H=H,
                D_ext=2.0 * (r0_in + a * Theta + 1.5 * t + b_in), inner_ch=inner_ch, t=t)


def solve_spiral_turns(p, target, value):
    """Number of turns N giving 'value' for target 'D_ext' (closed form) or 'A_sp' (one unit, root)."""
    t, bH, bC = p['t'], p['b_H'], p['b_C']
    pitch = bH + bC + 2.0 * t
    b_in = bH if p['inner_ch'] == 'hot' else bC
    if target == 'D_ext':
        return (value / 2.0 - p['r0_in'] - 1.5 * t - b_in) / pitch
    f = lambda N: spiral_geometry(bH, bC, t, p['H'], p['r0_in'], p['inner_ch'], N)['A_exchange'] - value
    if f(1.0) > 0:
        raise ValueError(f"Spiral: A = {value:.4g} m2 is below the area of a one-turn spiral: check r0_in, H.")
    if f(5000.0) < 0:
        raise ValueError(f"Spiral: A = {value:.4g} m2 would need more than 5000 turns.")
    return brentq(f, 1.0, 5000.0, xtol=1e-10, rtol=1e-12)


def spiral_radial_limits(sg):
    """radial_limits and theta_end as the toolbox routines read them (see spiral_geometry)."""
    t, b_in, p = sg['t'], sg['b_in'], sg['pitch']
    sheet1 = (-t / 2, t / 2)
    inner = (t / 2, t / 2 + b_in)
    sheet2 = (t / 2 + b_in, t / 2 + b_in + t)
    outer = (t / 2 + b_in + t, p - t / 2)
    hot, cold = (inner, outer) if sg['inner_ch'] == 'hot' else (outer, inner)
    Theta = sg['Theta']
    theta_end = {'sheet1': Theta, 'hot': Theta, 'sheet2': Theta, 'cold': Theta}
    theta_end['cold' if sg['inner_ch'] == 'hot' else 'hot'] = Theta - 2.0 * np.pi
    return {'sheet1': sheet1, 'hot': hot, 'sheet2': sheet2, 'cold': cold}, theta_end


def spiral_mean_curvature_radius(r_center, a, theta_end, n=2000):
    """Arc-length weighted mean radius of curvature of the centre line r = r_center + a theta."""
    th = np.linspace(0.0, theta_end, n + 1)
    mid = 0.5 * (th[1:] + th[:-1])
    r = r_center + a * mid
    rho = (r ** 2 + a ** 2) ** 1.5 / (r ** 2 + 2 * a ** 2)
    ds = np.diff(spiral_arc(r_center, a, th))
    return float(np.sum(rho * ds) / np.sum(ds))


# =============================================================================
# GEOMETRIES
# =============================================================================

def wall_resistance(htx, p, A=None, L=None, n_tot=1):
    """Wall conduction resistance [K/W], referred to the same area as the film resistances (UA = 1/(1/(aA)+...+R)).
    Applied ONLY if p['k_wall'] (wall conductivity [W/m/K]) is set; 0 otherwise.
        BPHE : plate,       R = t_p / (k A)      A heat transfer area of the exchanger
        SPHE : sheet,       R = t / (k A)             A = A_exchange x n_parallel
        PCHE : plate wall,  R = t_wall / (k A)        t_wall = p['t_wall'], default p['t_2'] (check that t_2 is the
                                                      wall between a hot and a cold channel)
        SCHE : wall,        R = t_wall / (k A_wall)   t_wall = p['t_wall'] (default: smallest distance between channels)
        STHE, TFHE : tube,  R = ln(OD / ID) / (2 pi k L n_tot)   L = L_t, n_tot = n_t x n_series x n_parallel"""
    k = p.get('k_wall')
    if not k:
        return 0.0
    if htx in ('STHE', 'TFHE'):
        return np.log(p['D_t_o'] / (p['D_t_o'] - 2 * p['t_t'])) / (2 * np.pi * k * L * n_tot)
    t = {'BPHE': lambda: p['t_p'], 'SPHE': lambda: p['t'], 'SCHE': lambda: p['t_wall'],
         'PCHE': lambda: p.get('t_wall', p.get('t_3'))}[htx]()
    if t is None:
        raise KeyError(f"{htx}: wall thickness missing for the wall resistance (k_wall is set).")
    return t / (k * A)


def _sizing_parameter(sizing, p):
    """Spirals (SPHE, SCHE): ONE of the sizing parameters (e.g. D_ext, n_turns) is needed.
    Returns [the one given by the user] (values written by a previous solve are recognised
    and ignored, as in solve), or [sizing[0]] if none is given (reported as missing)."""
    p = p or {}
    written = p.get('_spiral_written', {})
    given = [k for k in sizing if p.get(k) is not None and p[k] != written.get(k)]
    if given:
        return given[:1]
    return [sizing[0]]


class HXGeometry:
    HTX_TYPE = None

    def required_parameters(self, p=None):
        """Geometry parameters needed by solve, for the geometry options set in p (parameter
        dictionary of the model; None = options not set yet). The parameters needed by the
        correlations are added by the model (HexMB.get_required_parameters)."""
        raise NotImplementedError

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        """p : parameter dictionary (derived quantities are written back in it)."""
        raise NotImplementedError

    @staticmethod
    def _side_geom(p, base, htx):
        g = {**p, **base}
        g['HTX_Type'] = htx
        return g


class BPHEGeometry(HXGeometry):
    HTX_TYPE = 'BPHE'

    def required_parameters(self, p=None):
        return ['A_C', 'A_H', 'h', 'l', 'l_v', 'w', 'n_p', 't_p', 'n_ch_C', 'n_ch_H',
                'beta', 'R_f']

    @staticmethod
    def _wavelength_from_phi(phi, a):
        """Corrugation wavelength Lambda [m] of a sinusoidal profile of amplitude a [m] whose area enlargement
        factor is phi (Martin: X = 2 pi a / Lambda, phi = 1/6 (1 + sqrt(1 + X^2) + 4 sqrt(1 + X^2/2)))."""
        if not 1.0 < phi < 30.0:
            raise ValueError(f"Enlargement factor phi = {phi:.4g} outside ]1, 30[: check the heat transfer "
                             "area, the number of plates and the projected plate dimensions.")
        f = lambda X: (1 + np.sqrt(1 + X ** 2) + 4 * np.sqrt(1 + X ** 2 / 2)) / 6 - phi
        X = brentq(f, 1e-9, 50.0)
        return 2 * np.pi * a / X

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        p['n_ch'] = p['n_ch_H'] + p['n_ch_C']
        p['b'] = (p['h'] - p['n_p'] * p['t_p']) / p['n_ch']       # channel gap
        p['a'] = p['b'] / 2                                                  # corrugation amplitude
        n_active = p['n_p'] - 2                                               # end plates do not exchange
        A_projected = n_active * p['l_v'] * p['w']
        p['phi'] = p['A_H'] / A_projected
        p['Lambda'] = self._wavelength_from_phi(p['phi'], p['a'])
        p['X'] = 2 * np.pi * p['a'] / p['Lambda']
        p['Dh_C'] = p['Dh_H'] = 2 * p['b']
        p['A_flow_C'] = p['b'] * p['w'] * p['n_ch_C']
        p['A_flow_H'] = p['b'] * p['w'] * p['n_ch_H']

        D_port = p.get('D_port', 0)
        if D_port == 0 and p['l'] > p['l_v']:
            D_port = 0.7 * (p['l'] - p['l_v'])
            if p.get('l_h'):
                D_port = min(D_port, 0.7 * (p['w'] - p['l_h']))
            warnings.warn(f"D_port not given: estimated at {D_port*1e3:.0f} mm from l and l_v "
                          "(used for the port volume only).", stacklevel=3)
        A_port = np.pi * D_port ** 2 / 4 if D_port else 0.0
        p['D_port'] = D_port

        A_plate = p['l_v'] * p['w']
        for pre in ('H', 'C'):
            if p.get('V_' + pre) is not None:
                p['V_ch_' + pre], p['V_port_' + pre] = p['V_' + pre], 0.0
            else:
                p['V_ch_' + pre] = p['n_ch_' + pre] * p['b'] * (A_plate - 4 * A_port)
                p['V_port_' + pre] = 2 * A_port * p['h'] + p.get('V_p_' + pre, 0.0)

        sides = {}
        for pre, A in (('H', p['A_H']), ('C', p['A_C'])):
            g = self._side_geom(p, {'Dh': p['Dh_' + pre], 'L': p['l'], 'K': K, 'gamma': 0.0, 'A': A,
                                    'n_ch_side': p['n_ch_' + pre], 'beta': p['beta']}, 'BPHE')
            sides[pre] = SideGeometry(g, A, p['A_flow_' + pre], p['V_ch_' + pre], p['V_port_' + pre])
        return GeometryResult(sides['H'], sides['C'], R_extra=wall_resistance('BPHE', p, A=sides['H'].A))


class _TubeGeometry(HXGeometry):
    """Common part of STHE and TFHE: tubes of n_t x n_s x n_p, n_passes passes."""

    def _tube_basics(self, p, n_s, n_p):
        D_t_i = p['D_t_o'] - 2 * p['t_t']
        A_one = np.pi * D_t_i ** 2 / 4
        L_flow = p['L_t'] * p['n_passes'] * n_s
        A_flow = A_one * p['n_t'] * n_p / p['n_passes']        # G = (n_passes/n_p) m / (A_one n_t)
        return D_t_i, A_one, L_flow, A_flow


class STHEGeometry(_TubeGeometry):
    HTX_TYPE = 'STHE'

    def required_parameters(self, p=None):
        # The Bell-Delaware leakage and bypass parameters (D_OTL, N_strips, ...) are added by the
        # model when a Bell_Delaware correlation is chosen (tables of STHE_HTC / STHE_DP).
        return ['BC', 'D_s', 'L_t', 'D_t_o', 'n_passes', 't_t', 'B',
                'R_f_S', 'R_f_T', 'n_series', 'n_parallel', 'n_t', 'PR', 'k_wall',
                'theta', 'shell_side']

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        D_t_i, A_one, L_flow, A_flow = self._tube_basics(p, n_s, n_p)
        n_tot = p['n_t'] * n_s * n_p
        A_in = np.pi * D_t_i * p['L_t'] * n_tot
        A_out = np.pi * p['D_t_o'] * p['L_t'] * n_tot
        V_tube = A_one * p['L_t'] * n_tot
        V_shell = (np.pi / 4 * p['D_s'] ** 2 - np.pi / 4 * p['D_t_o'] ** 2 * p['n_t']) \
            * p['L_t'] * n_s * n_p
        if 'B' in p:
            p['n_B'] = np.round(p['L_t'] / p['B']) - 1
        p['V_T'], p['V_S'] = V_tube, V_shell
        R_fs = p['R_f_S'] / A_out if p.get('R_f_S') else 0.0
        R_ft = p['R_f_T'] / A_in if p.get('R_f_T') else 0.0
        outer_is_hot = p['shell_side'] == 'H'
        sides = {}
        for pre, is_outer in (('H', outer_is_hot), ('C', not outer_is_hot)):
            A = A_out if is_outer else A_in
            g = self._side_geom(p, {'Dh': D_t_i, 'L': L_flow, 'K': K, 'gamma': 0.0, 'A': A,
                                    'side': 'shell' if is_outer else 'tube'}, 'STHE')
            sides[pre] = SideGeometry(g, A, A_flow, V_shell if is_outer else V_tube)
        contact = STHEMultipassContact(p['n_passes'], p['shell_side']) if p.get('n_passes', 1) > 1 else None
        R_wall = wall_resistance('STHE', p, L=p['L_t'], n_tot=n_tot)
        return GeometryResult(sides['H'], sides['C'], R_extra=R_fs + R_ft + R_wall, contact=contact)


class TFHEGeometry(_TubeGeometry):
    HTX_TYPE = 'TFHE'

    def required_parameters(self, p=None):
        return ['A_flow', 'D_f_o', 'f_density', 't_f', 'fin_type', 'finned', 'L_t', 'D_t_o',
                'k_wall', 't_t', 'R_f', 'h', 'k_f', 'n_passes', 'n_rows', 'n_series', 'n_parallel',
                'n_t', 'P', 'PR', 'theta', 'w', 'fin_side']

    def _fins(self, p, D_t_i):
        p['P'] = p['P_V'] = p['P_H'] = p['PR'] * p['D_t_o']
        if p['fin_type'] == "Annular":
            p['n_f'] = p['L_t'] * p['f_density'] - 1
            r_f, r_t = p['D_f_o'] / 2, p['D_t_o'] / 2
            A_fin_tip = 2 * np.pi * r_f * p['t_f'] * p['n_f'] * p['n_t']
            A_bare_tube = 2 * np.pi * r_t * p['n_t'] * (p['L_t'] - p['n_f'] * p['t_f'])
            A_fin_faces = 2 * np.pi * (r_f ** 2 - r_t ** 2) * p['n_f'] * p['n_t']
            p['A_f'] = A_fin_tip + A_bare_tube + A_fin_faces
        elif p['fin_type'] == "Square":
            p['n_f'] = p['L_t'] * p['f_density']
            p['L_f'] = (p['L_t'] - p['n_f'] * p['t_f']) / p['n_f']
            A_r = 2 * (p['D_f_o'] ** 2 - 0.785 * p['D_t_o'] ** 2 + 2 * p['D_f_o'] * p['t_f']) \
                * (p['L_t'] / p['L_f']) * p['n_t']
            L_bare = p['L_t'] - p['n_f'] * p['t_f']                  # bare tube length [m]
            A_bare = np.pi * p['D_t_o'] * (p['L_t'] * (1 - p['t_f'] / p['L_f']) * p['n_t'] + L_bare)
            p['A_f'] = A_r + A_bare
        else:
            raise ValueError("Fin geometry is not 'Annular' nor 'Square'")
        p['A_t'] = np.pi * D_t_i * p['L_t'] * p['n_t']
        p['V_F'] = p['L_t'] * p['w'] * p['h']
        p['V_T'] = np.pi * D_t_i ** 2 / 4 * p['n_t'] * p['n_passes'] * p['L_t']

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        D_t_i, A_one, L_flow, A_flow = self._tube_basics(p, n_s, n_p)
        self._fins(p, D_t_i)
        # tube conduction (k_wall is a required parameter of the finned tubes)
        R_cond = wall_resistance('TFHE', p, L=p['L_t'], n_tot=p['n_t'] * n_s * n_p)
        outer_is_hot = p['fin_side'] == 'H'
        sides = {}
        for pre, is_outer in (('H', outer_is_hot), ('C', not outer_is_hot)):
            A = p['A_f'] if is_outer else p['A_t']
            g = self._side_geom(p, {'Dh': D_t_i, 'L': L_flow, 'K': K, 'gamma': 0.0, 'A': A,
                                    'side': 'fin' if is_outer else 'tube'}, 'TFHE')
            sides[pre] = SideGeometry(g, A, A_flow, p['V_F'] if is_outer else p['V_T'])
        return GeometryResult(sides['H'], sides['C'], R_extra=R_cond)


class PCHEGeometry(HXGeometry):
    """Printed circuit heat exchanger: n_p plates of n_ch_p etched channels, hot / cold plate ratio R_p.

    Channel cross-section, p['ch_section']:
        'semicircular' (default, previous behaviour): etched half-circle of diameter D_ch [m]
            A_ch = pi D_ch^2 / 8,   P_wet = (pi/2 + 1) D_ch,   D_h = pi D_ch / (2 + pi)
        'rectangular': passage between two etched fins (Ngo et al. 2007, Table 1 and Fig. 3):
            g_f [m] = fin gap (width of the passage, measured normal to the fins)
            d_f [m] = fin depth (etching depth)
            A_ch = g_f d_f,        P_wet = 2 (g_f + d_f),  D_h = 4 A_ch / P_wet
            (Ngo: g_f 1.31 mm, d_f 0.94 mm -> D_h 1.095 mm, paper 1.09 mm; 96 channels -> A_c
            11.82e-5 m2 = Table 1.)
    Channel path, p['type_channel'] (used by the correlations only): 'Straight', 'Zigzag', 'S_fin'.
    L_ch [m] = developed length of one channel (along the zigzag / S path).
    Heat transfer area of a side = N_side * P_wet * L_ch (all walls of the passage, cover plate
    included, as in the semicircular case)."""
    HTX_TYPE = 'PCHE'

    def required_parameters(self, p=None):
        p = p or {}
        required = ['k_wall', 'L_ch', 'n_ch_p', 'n_p', 'R_p', 't_2', 't_3']
        # channel section: semicircle of diameter D_ch, or rectangle g_f x d_f
        if str(p.get('ch_section', 'semicircular')).lower() == 'rectangular':
            required += ['g_f', 'd_f']
        else:
            required += ['D_ch']
        # channel path: the bending angle alpha is needed for zigzag and S-fin channels
        # (a correlation that needs alpha for straight channels adds it itself)
        if str(p.get('type_channel', 'Straight')).lower() in ('zigzag', 's_fin'):
            required += ['alpha']
        return required

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        section = str(p.get('ch_section', 'semicircular')).lower()
        if section == 'semicircular':
            if p.get('D_ch') is None:
                raise KeyError("PCHE: ch_section 'semicircular' needs D_ch [m].")
            A_ch = np.pi * p['D_ch'] ** 2 / 8
            P_wet = (np.pi / 2 + 1) * p['D_ch']
        elif section == 'rectangular':
            if p.get('g_f') is None or p.get('d_f') is None:
                raise KeyError("PCHE: ch_section 'rectangular' needs g_f (fin gap) and "
                               "d_f (fin depth) [m].")
            A_ch = p['g_f'] * p['d_f']
            P_wet = 2 * (p['g_f'] + p['d_f'])
        else:
            raise ValueError(f"PCHE: ch_section {p.get('ch_section')!r}: expected "
                             "'semicircular' or 'rectangular'.")
        D_h = 4 * A_ch / P_wet        # = pi D_ch / (2 + pi) for the semicircle
        N = p['n_ch_p'] * p['n_p']
        frac = {'H': p['R_p'] / (1 + p['R_p']), 'C': 1 / (1 + p['R_p'])}
        sides = {}
        for pre in ('H', 'C'):
            N_side = N * frac[pre]
            A = N_side * P_wet * p['L_ch'] * n_s * n_p
            g = self._side_geom(p, {'L': p['L_ch'] * n_s, 'K': K, 'gamma': 0.0, 'A': A, 'Dh': D_h,
                                    'D_t_i': D_h, 'n_passes': 1, 'L_t': p['L_ch'] * n_s,
                                    'n_t': N_side * n_p}, 'PCHE')
            sides[pre] = SideGeometry(g, A, A_ch * N_side * n_p, N_side * A_ch * p['L_ch'] * n_s * n_p)
        return GeometryResult(sides['H'], sides['C'], R_extra=wall_resistance('PCHE', p, A=sides['H'].A))
    
    
class SPHEGeometry(HXGeometry):
    """Spiral plate heat exchanger (double Archimedean spiral, see spiral_geometry). Parameter names are those
    of the previous model.
        H               plate height (axial width of the channels) [m]
        b_H, b_C        channel gaps [m]            t   plate thickness [m]
        r0_in           radius of the centre line of sheet1 at theta = 0 [m]
        inner_ch        'hot' or 'cold': the stream of the inner channel (the outer one is one turn shorter)
        ONE of  D_ext (outer diameter [m]) | n_turns (number of turns) | A_sp (exchange area of one unit [m2])
        k_wall          sheet conductivity (wall resistance, optional)    R_f_H, R_f_C   [m2 K/W] (optional)
        n_parallel      identical spirals in parallel (area and volume x n, flow split)
    Written in the parameters: the two of n_turns, D_ext, A_sp that were not given (recognised and ignored at
    the next solve), L_H, L_C, L_p_1, L_p_2, L_contact_1, L_contact_2, radial_limits, theta_end, Dh_H, Dh_C,
    H_V_tot, C_V_tot.
    Cells are numbered from the core on both streams: the cold stream enters at the core, the hot one at the
    periphery (counter-current). n_series > 1 is not available (the toolbox overlap handles one spiral)."""
    HTX_TYPE = 'SPHE'
    SIZING = ('D_ext', 'n_turns', 'A_sp')

    def required_parameters(self, p=None):
        return ['H', 'b_H', 'b_C', 't', 'r0_in', 'inner_ch'] + _sizing_parameter(self.SIZING, p)

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        if n_s != 1:
            raise NotImplementedError("SPHE: n_series > 1 is not implemented (spirals in parallel: n_parallel).")
        written = p.get('_spiral_written', {})            # derived values written by a previous solve
        given = [k for k in self.SIZING if p.get(k) is not None and p[k] != written.get(k)]
        if len(given) != 1:
            raise ValueError(f"SPHE: give exactly one of {self.SIZING} (got {given or 'none'}).")
        for k in ('H', 'b_H', 'b_C', 't', 'r0_in', 'inner_ch'):
            if p.get(k) is None:
                raise KeyError(f"SPHE: parameter '{k}' is missing.")
        if p['inner_ch'] not in ('hot', 'cold'):
            raise ValueError("SPHE: inner_ch must be 'hot' or 'cold'.")
        key = given[0]
        N = p['n_turns'] if key == 'n_turns' else solve_spiral_turns(p, key, p[key])
        if N < 1.0:
            raise ValueError(f"SPHE: {N:.2f} turns < 1: the outer diameter is too small for one full turn.")
        if N < 2.0:
            warnings.warn(f"SPHE: only {N:.2f} turns: the wall-contact map is not representative of a spiral.",
                          stacklevel=3)
        sg = spiral_geometry(p['b_H'], p['b_C'], p['t'], p['H'], p['r0_in'], p['inner_ch'], N)
        hot_in = p['inner_ch'] == 'hot'
        L_H, L_C = (sg['L_inner'], sg['L_outer']) if hot_in else (sg['L_outer'], sg['L_inner'])
        p['radial_limits'], p['theta_end'] = spiral_radial_limits(sg)     # read by the toolbox routines
        # the two sizing values that were not given are written (the given one is kept as it is) and
        # recorded in '_spiral_written', so that the next solve recognises and ignores them
        derived = {'n_turns': N, 'D_ext': sg['D_ext'], 'A_sp': sg['A_exchange']}
        del derived[key]
        p.update(derived)
        p['_spiral_written'] = derived
        p.update(L_H=L_H, L_C=L_C, L_p_1=sg['L_sheet1'], L_p_2=sg['L_sheet2'],
                 L_contact_1=sg['L_contact_1'], L_contact_2=sg['L_contact_2'])
        A = sg['A_exchange'] * n_p
        R_wall = wall_resistance('SPHE', p, A=A)
        R_foul = sum(p[k] / A for k in ('R_f_H', 'R_f_C') if p.get(k))

        sides = {}
        for pre, b, L, el in (('H', p['b_H'], L_H, 'hot'), ('C', p['b_C'], L_C, 'cold')):
            Dh = 2.0 * b * p['H'] / (b + p['H'])
            A_ch = b * p['H']
            p['Dh_' + pre], p[pre + '_V_tot'] = Dh, L * A_ch * n_p
            lo, hi = p['radial_limits'][el]
            R_mean = spiral_mean_curvature_radius(p['r0_in'] + 0.5 * (lo + hi), sg['a'], p['theta_end'][el])
            g = self._side_geom(p, {'Dh': Dh, 'L': L, 'K': K, 'gamma': 0.0, 'A': A,
                                    'R_curv': R_mean, 'b': b, 'H': p['H'], 'D_t_i': Dh,
                                    'L_t': L, 'n_passes': 1, 'n_t': n_p, 'A_ch': A_ch}, 'SPHE')
            sides[pre] = SideGeometry(g, A, A_ch * n_p, A_ch * L * n_p)
        return GeometryResult(sides['H'], sides['C'], R_extra=R_wall + R_foul, contact=SpiralContact(),
                              derived=dict(n_turns=N, D_ext=sg['D_ext']))


class SCHEGeometry(HXGeometry):
    """Spiral channel heat exchanger: elliptical channels in a double Archimedean spiral, hot and cold channels
    staggered in stacked layers (ported from the previous model).
        H                       height of the stack [m]            r0_in   radius at theta = 0 [m]
        w_ch_C / h_ch_C   width / height of the cold channels [m];  w_ch_H / h_ch_H the same for hot
        P_H / P_V       radial / vertical centre-to-centre distance of two channels of the SAME fluid [m]
        inner_ch           'hot' or 'cold' (the outer channel is one turn shorter)
        ONE of D_ext | n_turns
        exchange_frac (1)       fraction of the channel perimeter exchanging heat
        k_wall                  wall conductivity, optional: R_wall = t_wall / (k_wall A_wall);  t_wall defaults to
                                the minimum wall between neighbouring channels
    Cells are numbered from the core on both streams (cold enters at the core, hot at the periphery).
    Aliases written for the toolbox routines: t, b_C, b_H, radial_limits, theta_end.
    n_p = number of channel layers (int(H / P_V))."""
    HTX_TYPE = 'SCHE'
    SIZING = ('D_ext', 'n_turns')

    def required_parameters(self, p=None):
        return ['H', 'w_ch_C', 'h_ch_C', 'w_ch_H', 'h_ch_H', 'P_H', 'P_V', 'r0_in',
                'inner_ch'] + _sizing_parameter(self.SIZING, p)

    @staticmethod
    def _ellipse(w, h):
        a_, b_ = w / 2.0, h / 2.0
        A = np.pi * a_ * b_
        Per = np.pi * (3 * (a_ + b_) - np.sqrt((3 * a_ + b_) * (a_ + 3 * b_)))          # Ramanujan
        return A, Per, 4.0 * A / Per

    def solve(self, p, m_h, m_c, n_s=1, n_p=1, K=0.0):
        if n_s != 1:
            raise NotImplementedError("SCHE: n_series > 1 is not implemented (spirals in parallel: n_parallel).")
        for k in ('H', 'w_ch_C', 'h_ch_C', 'w_ch_H', 'h_ch_H', 'P_H', 'P_V', 'r0_in',
                  'inner_ch'):
            if p.get(k) is None:
                raise KeyError(f"SCHE: parameter '{k}' is missing.")
        written = p.get('_spiral_written', {})
        given = [k for k in self.SIZING if p.get(k) is not None and p[k] != written.get(k)]
        if len(given) != 1:
            raise ValueError(f"SCHE: give exactly one of {self.SIZING} (got {given or 'none'}).")
        Cw, Ch, Hw, Hh = p['w_ch_C'], p['h_ch_C'], p['w_ch_H'], p['h_ch_H']
        pr, pv, H, r0 = p['P_H'], p['P_V'], p['H'], p['r0_in']
        a = pr / (2.0 * np.pi)

        C_A, C_Per, p['Dh_C'] = self._ellipse(Cw, Ch)
        H_A, H_Per, p['Dh_H'] = self._ellipse(Hw, Hh)
        if min(pr - Cw, pr - Hw) <= 0:
            raise ValueError("SCHE: P_H <= channel width: two channels of the same fluid touch.")
        if min(pv - Ch, pv - Hh) <= 0:
            raise ValueError("SCHE: P_V <= channel height: two channels of the same fluid touch.")
        frac_ex = p.setdefault('exchange_frac', 1.0)
        P_ex_c, P_ex_h = frac_ex * C_Per, frac_ex * H_Per
        P_ex = np.sqrt(P_ex_h * P_ex_c)

        # minimum wall between a cold channel and the diagonal hot neighbour (sampling of the two ellipses)
        phi = np.linspace(0.0, 2 * np.pi, 361)[:-1]
        xc, zc = 0.5 * Cw * np.cos(phi), 0.5 * Ch * np.sin(phi)
        xh, zh = 0.5 * pr + 0.5 * Hw * np.cos(phi), 0.5 * pv + 0.5 * Hh * np.sin(phi)
        inside = (((xc - 0.5 * pr) / (0.5 * Hw)) ** 2 + ((zc - 0.5 * pv) / (0.5 * Hh)) ** 2 < 1).any() \
            or ((xh / (0.5 * Cw)) ** 2 + (zh / (0.5 * Ch)) ** 2 < 1).any()
        if inside:
            raise ValueError("SCHE: neighbouring hot and cold channels overlap (P_H, P_V too small).")
        t_wall_min = float(np.sqrt((xc[:, None] - xh[None, :]) ** 2 + (zc[:, None] - zh[None, :]) ** 2).min())

        # radial positions at theta = 0 (plan projection), centres shifted by pitch/2
        b_inner, b_outer = (Hw, Cw) if p['inner_ch'] == 'hot' else (Cw, Hw)
        c_in = 0.5 * b_inner
        c_out = c_in + 0.5 * pr
        inner = (c_in - 0.5 * b_inner, c_in + 0.5 * b_inner)
        outer = (c_out - 0.5 * b_outer, c_out + 0.5 * b_outer)
        g = max(0.5 * pr - 0.5 * (Cw + Hw), 0.0)
        sheet2 = (c_in + 0.25 * pr - 0.5 * g, c_in + 0.25 * pr + 0.5 * g)
        sheet1 = (c_in - 0.25 * pr - 0.5 * g, c_in - 0.25 * pr + 0.5 * g)
        hot, cold = (inner, outer) if p['inner_ch'] == 'hot' else (outer, inner)
        p['radial_limits'] = RL = {'sheet1': sheet1, 'hot': hot, 'sheet2': sheet2, 'cold': cold}

        key = given[0]
        outer_offset = max(inner[1], sheet2[1])
        N = p['n_turns'] if key == 'n_turns' else (p['D_ext'] / 2.0 - r0 - outer_offset) / pr
        if N < 1.0:
            raise ValueError(f"SCHE: {N:.2f} turns < 1: the outer diameter is too small for one full turn.")
        Theta = 2.0 * np.pi * N
        outer_el = 'cold' if p['inner_ch'] == 'hot' else 'hot'
        p['theta_end'] = TE = {'sheet1': Theta, 'hot': Theta, 'sheet2': Theta, 'cold': Theta}
        TE[outer_el] = Theta - 2.0 * np.pi

        def length(el, th0=0.0, th1=None):
            th1 = TE[el] if th1 is None else th1
            if th1 <= th0:
                return 0.0
            return float(spiral_arc(r0 + a * th0 + np.mean(RL[el]), a, th1 - th0))

        L_H, L_C = length('hot'), length('cold')
        p['L_H'], p['L_C'] = L_H, L_C
        p['L_p_1'], p['L_p_2'] = length('sheet1'), length('sheet2')

        n_ch = int(H // pv)
        if n_ch < 1:
            raise ValueError("SCHE: H < P_V: no channel layer.")
        p['n_p'], p['t_plate'] = n_ch, 0.5 * pv           # n_p = number of channel layers

        L_contact = length('sheet1', 2.0 * np.pi) + length('sheet2', 0.0, TE['sheet2'] - 2.0 * np.pi)
        A_h = L_contact * n_ch * 0.5 * P_ex_h * n_p
        A_c = L_contact * n_ch * 0.5 * P_ex_c * n_p
        A_wall = L_contact * n_ch * 0.5 * P_ex * n_p
        p['H_contact'] = n_ch * 0.5 * P_ex
        p['t_wall'] = p.get('t_wall', t_wall_min)
        R_wall = wall_resistance('SCHE', p, A=A_wall)

        # aliases read by the toolbox routines (t is a numerical alias here, not a thickness)
        p['t'], p['b_C'], p['b_H'] = 0.5 * (pr - Cw - Hw), Cw, Hw
        # the sizing value that was not given is written (the given one is kept as it is) and recorded
        # in '_spiral_written', so that the next solve recognises and ignores it
        derived = {'n_turns': N, 'D_ext': 2.0 * (r0 + a * Theta + outer_offset)}
        del derived[key]
        p.update(derived)
        p['_spiral_written'] = derived
        p['A_flow_H'], p['A_flow_C'] = H_A * n_ch, C_A * n_ch
        p['H_V_tot'], p['C_V_tot'] = L_H * H_A * n_ch * n_p, L_C * C_A * n_ch * n_p

        sides = {}
        for pre, el, L, A_ch, Dh, A in (('H', 'hot', L_H, H_A, p['Dh_H'], A_h), ('C', 'cold', L_C, C_A, p['Dh_C'], A_c)):
            lo, hi = RL[el]
            R_mean = spiral_mean_curvature_radius(r0 + 0.5 * (lo + hi), a, TE[el])
            gm = self._side_geom(p, {'Dh': Dh, 'L': L, 'K': K, 'gamma': 0.0, 'A': A, 'R_curv': R_mean,
                                     'D_t_i': Dh, 'L_t': L, 'n_passes': 1, 'n_t': n_ch * n_p,
                                     'A_ch': A_ch}, 'SCHE')
            sides[pre] = SideGeometry(gm, A, A_ch * n_ch * n_p, A_ch * n_ch * L * n_p)
        return GeometryResult(sides['H'], sides['C'], R_extra=R_wall, contact=SpiralContact(),
                              derived=dict(n_turns=N, D_ext=p['D_ext']))


GEOMETRIES = {cls.HTX_TYPE: cls for cls in (BPHEGeometry, STHEGeometry, TFHEGeometry, PCHEGeometry, SPHEGeometry,
                                             SCHEGeometry)}


def make_geometry(HTX_Type):
    if HTX_Type not in GEOMETRIES:
        raise ValueError(f"No geometry for {HTX_Type!r}: choose from {tuple(GEOMETRIES)}.")
    return GEOMETRIES[HTX_Type]()