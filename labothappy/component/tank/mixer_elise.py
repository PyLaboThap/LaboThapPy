"""
Created on Tue Nov 17 2026

@author: Elise

Mixer (junction) counterpart of TankSpliter (Basile, modified by Elise).
Same optional geometry, read in the flow direction:

    n inlet tubes  ->  transition (converging/diverging)  ->  one outlet body

If the geometry is not provided, the component is an ideal adiabatic mixer
(enthalpy mass-weighted, outlet pressure from `P_mix_model`, no charge computed).
"""

import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

from labothappy.correlations.pressure_drop.dp_corr_distributor import compute_dp
from labothappy.correlations.properties.two_phase import compute_two_phase_density
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction, VOID_FRACTION_MODELS

from labothappy.connector.mass_connector import MassConnector
from labothappy.component.base_component import BaseComponent


class TankMixer(BaseComponent):
    """
    **Component**: Mixer

    **Model**: Adiabatic mixing of n streams, with optional pressure drop and charge

    **Description**:

        The `Mixer` component merges n incoming streams into one outlet stream.
        Mass and energy are conserved (adiabatic, no work):

            m_dot_ex = sum(m_dot_su_i)
            h_ex     = sum(m_dot_su_i * h_su_i) / m_dot_ex

        Geometry (all optional) describes the component as three parts in series:

            1. Inlet tubes      : n straight tubes, inner diameter D_su[i], length L_su[i]
            2. Transition       : converging (or diverging) section of length L_trans,
                                  from the n inlet tubes to the outlet body
            3. Outlet body      : straight tube, inner diameter D_ex, length L_ex

        NOTE: with respect to TankSpliter, the roles are swapped. To reuse a splitter
        geometry for the mirrored component: D_su(mixer) = D_ex(splitter),
        L_su(mixer) = L_ex(splitter), D_ex(mixer) = D_su(splitter),
        L_ex(mixer) = L_su(splitter). L_trans, V_trans, V_int are unchanged.

        - If D_su, L_su, D_ex and L_ex are given, pressure drops are computed.
        - If, additionally, L_trans or V_trans is given (or the total internal
          volume V_int alone), the refrigerant charge is computed.
        - Otherwise, the outlet pressure is obtained from the inlet pressures
          with `P_mix_model` and no charge is computed.

    **Junction pressure**:

        At the junction, all branches share one pressure. With imposed inlet pressures
        and imposed mass flows, the problem is over-determined: each branch "arrives"
        at its own pressure P_su_i - dp_branch_i. These branch pressures are combined
        into one junction pressure with `P_mix_model`:
            'min'           : lowest branch pressure (default, conservative)
            'mass_weighted' : mass-flow weighted average
        The spread between branch pressures is reported in DP_detail['branch_imbalance']:
        in a consistent system model it should tend to zero (the upstream components
        adapt their pressures/flows).

    **Assumptions**:

        - Steady-state, adiabatic, no shaft work.
        - Inlet tubes: properties at each inlet state. Transition and outlet body:
          properties of the mixed stream (h_ex) at the junction pressure
          (pressure drops are assumed small compared to the pressure level).
        - Mixing is assumed to happen in the transition.
        - Two-phase flow: homogeneous model for the pressure drop (McAdams viscosity),
          void-fraction model (Zivi by default) for the charge.

    **Connectors**:

        su_1, ..., su_n (MassConnector): Inlet streams (n = n_inlets).
        ex (MassConnector): Outlet stream.

    **Parameters** (all optional, set via `set_parameters(...)`):

        D_su (float or list): Inner diameter of the inlet tubes [m] (one value or one per inlet).
        L_su (float or list): Length of the inlet tubes [m] (one value or one per inlet).
        D_ex (float): Inner diameter of the outlet body [m].
        L_ex (float): Length of the outlet body [m].
        L_trans (float): Length of the transition section [m]. If missing, the transition
                         is treated as a sudden contraction/expansion for the pressure drop.
        V_trans (float): Internal volume of the transition [m^3]. If missing, estimated as a
                         cone frustum between the equivalent inlet diameter and D_ex.
        V_int (float): Total internal volume [m^3]. Used for the charge only when the detailed
                       geometry is not available.
        roughness (float): Absolute wall roughness [m]. Default 1.5e-6 (drawn tube).
        K_extra (float): Additional local loss coefficient for the flow combining, based on the
                         outlet body (combined flow) velocity. Default 0.
        include_acceleration (bool): If True, the reversible static-pressure change due to the
                         velocity change is added. Default False.
        void_fraction_model (str): 'Zivi' (default), 'Homogeneous' or any model of
                                   VOID_FRACTION_MODELS that does not need the mass flux.
        dp_corr (dict): Pressure-drop correlation names per regime.
                        Default {'single-phase': 'Churchill', 'two-phase': 'MSH'}.
        P_mix_model (str): 'min' (default) or 'mass_weighted'. See "Junction pressure".

    **Inputs** (on each inlet connector su_i):

        fluid (str): Fluid (must be identical for all inlets).
        P (float): Inlet pressure [Pa].
        T or h (float): Inlet temperature [K] or specific enthalpy [J/kg].
        m_dot (float): Inlet mass flow rate [kg/s].

    **Outputs**:

        ex: fluid, P (junction pressure minus outlet body pressure drop), h (mixed), m_dot (sum).

        DP (np.ndarray or None): Pressure drop from each inlet to the outlet, P_su_i - P_ex [Pa].
        DP_detail (dict or None): Breakdown (inlet tubes, transition, acceleration, outlet body,
                                  junction pressures, branch imbalance).
        charge (float or None): Refrigerant charge in the component [kg].
        charge_detail (dict or None): Charge breakdown per section [kg].
    """

    DP_KEYS = ('D_su', 'L_su', 'D_ex', 'L_ex')
    P_MIX_MODELS = ('min', 'mass_weighted')

    def __init__(self, n_inlets):

        super().__init__()

        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.solved = False

        "Number of inlets"
        self.n_inlets = int(n_inlets)
        if self.n_inlets < 1:
            raise ValueError("'Mixer' needs at least one inlet.")

        "Inputs"
        for i in range(self.n_inlets):
            setattr(self, f"su_{i + 1}", MassConnector())

        "Output"
        self.ex = MassConnector()

        "Optional geometry: default to None (= not available)"
        for name in self.get_required_parameters():
            self.params.setdefault(name, None)

        "Results"
        self._reset_results()

    @property
    def inlets(self):
        return [getattr(self, f"su_{i + 1}") for i in range(self.n_inlets)]

    # ------------------------------------------------------------------ #
    # Framework methods
    # ------------------------------------------------------------------ #

    def get_required_inputs(self):
        # One set of inputs per inlet connector
        inputs = []
        for i in range(1, self.n_inlets + 1):
            inputs += [f'P_su_{i}', f'T_su_{i}', f'm_dot_su_{i}']
        return inputs + ['fluid']

    def get_required_parameters(self):
        # Geometry listed as for the other models. None means "not available".
        return ['D_su', 'L_su', 'D_ex', 'L_ex', 'L_trans']

    def check_calculable(self):
        # Checked directly on the connectors (n inlets).
        self.calculable = all(
            c.fluid is not None and c.p is not None and c.h is not None and c.m_dot is not None
            for c in self.inlets)
        return self.calculable

    def check_parametrized(self):
        # The geometry is optional: the mixer is always parametrized.
        self.parametrized = True
        if self.print_flag:
            missing = [k for k in self.get_required_parameters() if self.params.get(k) is None]
            if missing:
                print(f"Component {self.__class__.__name__}: missing parameters {', '.join(missing)}. "
                      f"Pressure drop computed: {self._can_compute_dp()}. The charge may not be computed.")
        return self.parametrized

    def solve(self):

        if not self.check_parametrized():
            raise ValueError("'Mixer' is not fully parametrized.")
        if not self.check_calculable():
            raise ValueError("'Mixer' is not fully calculable.")

        self._reset_results()
        inlets = self.inlets

        fluid = inlets[0].fluid
        if any(c.fluid != fluid for c in inlets):
            raise ValueError("'Mixer' all inlets must carry the same fluid.")

        p_su = np.array([c.p for c in inlets], dtype=float)
        h_su = np.array([c.h for c in inlets], dtype=float)
        m_su = np.array([c.m_dot for c in inlets], dtype=float)
        if np.any(m_su < 0):
            raise ValueError("'Mixer' negative inlet mass flow rate (reverse flow is not modelled).")
        m_dot = m_su.sum()
        if m_dot <= 0:
            raise ValueError("'Mixer' total inlet mass flow rate must be positive.")

        "1) Mass and energy balances"
        h_ex = float(np.dot(m_su, h_su) / m_dot)

        "2) Pressures (optional pressure drop)"
        if self._can_compute_dp():
            p_ex, p_junction = self._compute_pressure_drops(fluid, p_su, h_su, m_su, h_ex)
            if p_ex <= 0:
                raise ValueError("'Mixer' computed pressure drop exceeds inlet pressure.")
            self.DP = p_su - p_ex
        else:
            p_ex = p_junction = self._mix_pressure(p_su, m_su)

        "3) Outlet state"
        self.ex.reset()
        self.ex.set_fluid(fluid)
        self.ex.set_p(p_ex)
        self.ex.set_h(h_ex)
        self.ex.set_m_dot(m_dot)

        "4) Refrigerant charge (optional)"
        if self._can_compute_charge():
            self.charge = self._compute_charge(fluid, p_su, h_su, p_junction, p_ex, h_ex)

        self.solved = True
        return

    # ------------------------------------------------------------------ #
    # Parameter helpers
    # ------------------------------------------------------------------ #

    def _reset_results(self):
        self.DP = None
        self.DP_detail = None
        self.charge = None
        self.charge_detail = None

    def _per_inlet(self, value, name):
        """Broadcast a scalar to all inlets, or check a list has one value per inlet."""
        n = self.n_inlets
        arr = np.atleast_1d(np.asarray(value, dtype=float))
        if arr.size == 1:
            return np.full(n, arr[0])
        if arr.size != n:
            raise ValueError(f"'Mixer' parameter '{name}' must be a scalar or have {n} values.")
        return arr

    def _can_compute_dp(self):
        return all(self.params.get(k) is not None for k in self.DP_KEYS)

    def _can_compute_charge(self):
        detailed = self._can_compute_dp() and (
            self.params.get('L_trans') is not None or self.params.get('V_trans') is not None)
        return detailed or self.params.get('V_int') is not None

    def _geometry(self):
        D_su = self._per_inlet(self.params.get('D_su'), 'D_su')
        D_ex = float(self.params.get('D_ex'))
        A_su = np.pi * D_su**2 / 4
        A_ex = np.pi * D_ex**2 / 4
        D_eq_su = np.sqrt(4 * A_su.sum() / np.pi)  # diameter of one tube with the total inlet area
        return D_su, D_ex, A_su, A_ex, D_eq_su

    def _mix_pressure(self, p, m):
        """Combine branch pressures into one junction pressure."""
        model = self.params.get('P_mix_model') or 'min'
        if model == 'min':
            return float(np.min(p))
        if model == 'mass_weighted':
            return float(np.dot(m, p) / m.sum())
        raise ValueError(f"'Mixer' unknown P_mix_model '{model}'. Available: {self.P_MIX_MODELS}.")

    # ------------------------------------------------------------------ #
    # Fluid properties
    # ------------------------------------------------------------------ #

    @staticmethod
    def _fluid_state(fluid, p, h, void_fraction_model):
        """
        Properties needed for pressure drop and charge.
        Two-phase: void-fraction weighted density (void_fraction_model), McAdams viscosity.
        """
        x = None
        try:
            q = PropsSI('Q', 'P', p, 'H', h, fluid)
            if 0.0 <= q <= 1.0:
                x = q
        except ValueError:
            pass  # e.g. incompressible fluids: treated as single phase

        if x is None:
            return {'two_phase': False,
                    'rho': PropsSI('D', 'P', p, 'H', h, fluid),
                    'mu': PropsSI('V', 'P', p, 'H', h, fluid)}

        rho_l = PropsSI('D', 'P', p, 'Q', 0, fluid)
        rho_v = PropsSI('D', 'P', p, 'Q', 1, fluid)
        AS = CP.AbstractState('HEOS', fluid)
        AS.update(CP.HmassP_INPUTS, h, p)
        alpha = compute_void_fraction(AS, {}, void_fraction_model=void_fraction_model)
        rho = compute_two_phase_density(x, rho_l, rho_v, alpha)
        mu_l = PropsSI('V', 'P', p, 'Q', 0, fluid)
        mu_v = PropsSI('V', 'P', p, 'Q', 1, fluid)
        return {'two_phase': True, 'x': x,
                'rho_l': rho_l, 'rho_v': rho_v,
                'rho': rho,
                'mu': 1 / (x / mu_v + (1 - x) / mu_l)}

    # ------------------------------------------------------------------ #
    # Pressure drop
    # ------------------------------------------------------------------ #

    @staticmethod
    def _transition_K(sigma, theta):
        """
        Crane TP-410 coefficients for gradual/sudden contraction or expansion.
        sigma = A_downstream / A_upstream (here A_outlet_body / A_inlets_total),
        theta = total included angle [rad].
        Returns (K, basis): basis = 'downstream' or 'upstream', i.e. the velocity
        K refers to (always the smaller cross-section).
        Ref: Crane, Technical Paper No. 410M, 1979.
        """
        if sigma < 1:  # contraction -> downstream (small) velocity
            if theta <= np.pi / 4:
                K = 0.8 * np.sin(theta / 2) * (1 - sigma)
            else:
                K = 0.5 * np.sqrt(np.sin(theta / 2)) * (1 - sigma)
            return K, 'downstream'

        beta2 = 1 / sigma  # expansion -> upstream (small) velocity
        if theta <= np.pi / 4:
            K = 2.6 * np.sin(theta / 2) * (1 - beta2)**2
        else:
            K = (1 - beta2)**2
        return K, 'upstream'

    def _compute_pressure_drops(self, fluid, p_su, h_su, m_su, h_ex):
        n = self.n_inlets
        m_dot = m_su.sum()

        D_su, D_ex, A_su, A_ex, D_eq_su = self._geometry()
        L_su = self._per_inlet(self.params.get('L_su'), 'L_su')
        L_ex = float(self.params.get('L_ex'))
        L_trans = self.params.get('L_trans')
        eps = float(self.params.get('roughness', 1.5e-6))
        K_extra = float(self.params.get('K_extra', 0.0))
        include_acc = bool(self.params.get('include_acceleration', False))
        dp_corr = self.params.get('dp_corr') or {'single-phase': 'Churchill', 'two-phase': 'MSH'}
        void_fraction_model = self.params.get('void_fraction_model', 'Zivi')

        "Fluid states: each inlet at its own state, mixed stream at the reference junction pressure"
        st_su = [self._fluid_state(fluid, p_su[i], h_su[i], void_fraction_model) for i in range(n)]
        rho_su = np.array([s['rho'] for s in st_su])
        mu_su = np.array([s['mu'] for s in st_su])

        p_ref = self._mix_pressure(p_su, m_su)  # small dp assumption for the mixed-stream properties
        st_ex = self._fluid_state(fluid, p_ref, h_ex, void_fraction_model)
        rho_ex, mu_ex = st_ex['rho'], st_ex['mu']

        "Velocities"
        v_su = m_su / (rho_su * A_su)      # Vector with the different inlet velocities
        v_ex = m_dot / (rho_ex * A_ex)     # Velocity in the outlet body

        "1) Inlet tubes: friction (one tube per inlet, properties of that inlet)"
        Re_su = rho_su * v_su * D_su / mu_su   # Only reported in DP_detail
        dp_su = np.zeros(n)
        for i, c in enumerate(self.inlets):
            if m_su[i] == 0:
                continue  # Branch closed: no flow, no friction
            corr = dp_corr['two-phase' if st_su[i]['two_phase'] else 'single-phase']
            geom_su = {'A_flow': A_su[i],
                       'D_h': D_su[i],
                       'roughness': eps,
                       'chevron_angle': None}
            dp_su[i] = compute_dp(c.AS, m_su[i], geom_su, L_su[i], corr, void_fraction_model)

        "2) Transition: local loss (contraction/expansion) + optional combining-flow coefficient"
        sigma = A_ex / A_su.sum()   # Area ratio (outlet body area / total inlet area)
        if L_trans is None or L_trans <= 0:
            theta = np.pi  # sudden change of section
        else:
            theta = 2 * np.arctan(abs(D_eq_su - D_ex) / (2 * L_trans))
        K, basis = self._transition_K(sigma, theta)
        if basis == 'downstream':
            dp_trans = np.full(n, K * rho_ex * v_ex**2 / 2)   # Same for all branches
        else:
            dp_trans = K * rho_su * v_su**2 / 2               # Each branch with its own velocity
        dp_trans = dp_trans + K_extra * rho_ex * v_ex**2 / 2

        "3) Optional reversible static-pressure change due to acceleration (per branch)"
        dp_acc = (rho_ex * v_ex**2 - rho_su * v_su**2) / 2 if include_acc else np.zeros(n)

        "4) Junction pressure"
        p_branch = p_su - dp_su - dp_trans - dp_acc     # Pressure at which each branch reaches the junction
        p_junction = self._mix_pressure(p_branch, m_su)

        "5) Outlet body: friction of the mixed stream"
        mixed = MassConnector()          # Temporary connector at the mixed state (gives the AS)
        mixed.set_fluid(fluid)
        mixed.set_p(p_junction)
        mixed.set_h(h_ex)
        corr_ex = dp_corr['two-phase' if st_ex['two_phase'] else 'single-phase']
        geom_ex = {'A_flow': A_ex,
                   'D_h': D_ex,
                   'roughness': eps,
                   'chevron_angle': None}
        dp_ex = compute_dp(mixed.AS, m_dot, geom_ex, L_ex, corr_ex, void_fraction_model)
        Re_ex = rho_ex * v_ex * D_ex / mu_ex   # Only reported in DP_detail

        p_ex = p_junction - dp_ex

        self.DP_detail = {
            'inlet_tubes': dp_su,
            'transition': dp_trans,
            'acceleration': dp_acc,
            'outlet_body': dp_ex,
            'branch_pressures': p_branch,
            'junction_pressure': p_junction,
            'branch_imbalance': float(p_branch.max() - p_branch.min()),
            'K_transition': K,
            'theta_deg': np.degrees(theta),
            'area_ratio': sigma,
            'Re_su': Re_su,
            'Re_ex': Re_ex,
            'two_phase_su': [s['two_phase'] for s in st_su],
            'two_phase_ex': st_ex['two_phase'],
        }
        return p_ex, p_junction

    # ------------------------------------------------------------------ #
    # Refrigerant charge
    # ------------------------------------------------------------------ #

    def _compute_charge(self, fluid, p_su, h_su, p_junction, p_ex, h_ex):
        void_fraction_model = self.params.get('void_fraction_model', 'Zivi')

        "Only a total volume is known: mixed state at the outlet"
        if not self._can_compute_dp() or (self.params.get('L_trans') is None
                                          and self.params.get('V_trans') is None):
            V_int = float(self.params.get('V_int'))
            rho = self._fluid_state(fluid, p_ex, h_ex, void_fraction_model)['rho']
            self.charge_detail = {'total_volume': rho * V_int}
            return rho * V_int

        "Detailed geometry"
        D_su, D_ex, A_su, A_ex, D_eq_su = self._geometry()
        L_su = self._per_inlet(self.params.get('L_su'), 'L_su')
        L_ex = float(self.params.get('L_ex'))

        V_su = A_su * L_su
        V_trans = self.params.get('V_trans')
        if V_trans is None:
            L_trans = float(self.params.get('L_trans'))
            V_trans = np.pi * L_trans / 12 * (D_eq_su**2 + D_eq_su * D_ex + D_ex**2)  # cone frustum
        V_ex = A_ex * L_ex

        rho_su = np.array([self._fluid_state(fluid, p_su[i], h_su[i], void_fraction_model)['rho']
                           for i in range(self.n_inlets)])
        rho_trans = self._fluid_state(fluid, p_junction, h_ex, void_fraction_model)['rho']
        rho_ex = self._fluid_state(fluid, (p_junction + p_ex) / 2, h_ex, void_fraction_model)['rho']

        charge_su = rho_su * V_su
        charge_trans = rho_trans * float(V_trans)
        charge_ex = rho_ex * V_ex

        self.charge_detail = {
            'inlet_tubes': charge_su,
            'transition': charge_trans,
            'outlet_body': charge_ex,
            'volumes': {'inlet_tubes': V_su, 'transition': float(V_trans), 'outlet_body': V_ex},
        }
        return float(charge_su.sum()) + charge_trans + charge_ex