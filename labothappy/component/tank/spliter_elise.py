# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile

Extended with optional pressure-drop and refrigerant-charge computation.
If the geometry is not provided, the component behaves exactly as before
(isobaric, isenthalpic splitting, no charge computed).
"""

import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

from labothappy.correlations.pressure_drop.dp_corr_distributor import compute_dp
from labothappy.correlations.properties.two_phase import compute_two_phase_density
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction, VOID_FRACTION_MODELS

from labothappy.connector.mass_connector import MassConnector
from labothappy.component.base_component import BaseComponent


class TankSpliter(BaseComponent):
    """
    **Component**: Spliter

    **Model**: Mass flow splitting model, with optional pressure drop and charge

    **Description**:

        The `Spliter` component divides an incoming mass flow into multiple outlet
        streams based on specified mass flow fractions. Enthalpy is conserved
        (adiabatic, no work).

        Geometry (all optional) describes the component as three parts in series:

            1. Inlet body       : straight tube, inner diameter D_su, length L_su
            2. Transition       : converging (or diverging) section of length L_trans,
                                  from D_su to the n outlet tubes
            3. Outlet tubes     : n straight tubes, inner diameter D_ex[i], length L_ex[i]

        - If D_su, L_su, D_ex and L_ex are given, pressure drops are computed
          (one value per outlet) and each outlet gets its own pressure.
        - If, additionally, L_trans or V_trans is given (or the total internal
          volume V_int alone), the refrigerant charge is computed.
        - Otherwise, the outlets keep the inlet pressure and no charge is computed
          (original behaviour).

    **Assumptions**:

        - Steady-state, adiabatic, no shaft work.
        - Fluid properties evaluated at inlet conditions for the pressure drop
          (pressure drops are assumed small compared to the pressure level).
        - Two-phase flow: homogeneous model for the pressure drop (McAdams viscosity),
          void-fraction model (Zivi by default) for the charge.
        - Flow maldistribution between outlets is NOT predicted: the split is imposed
          by `outlet_repartition`.

    **Connectors**:

        su (MassConnector): Inlet supply stream.
        ex_1, ..., ex_n (MassConnector): Outlet streams (n = len(outlet_repartition)).

    **Parameters** (all optional, set via `set_parameters(...)` or as attributes of `geom`):

        D_su (float): Inner diameter of the inlet body [m].
        L_su (float): Length of the inlet body [m].
        D_ex (float or list): Inner diameter of the outlet tubes [m] (one value or one per outlet).
        L_ex (float or list): Length of the outlet tubes [m] (one value or one per outlet).
        L_trans (float): Length of the transition section [m]. If missing, the transition
                         is treated as a sudden contraction/expansion for the pressure drop.
        V_trans (float): Internal volume of the transition [m^3] (e.g. from CAD). If missing,
                         estimated as a cone frustum between D_su and the equivalent outlet diameter.
        V_int (float): Total internal volume [m^3]. Used for the charge only when the detailed
                       geometry is not available.
        roughness (float): Absolute wall roughness [m]. Default 1.5e-6 (drawn tube).
        K_extra (float): Additional local loss coefficient (based on outlet velocity) for
                         effects not covered by the correlations (e.g. flow dividing, inner
                         walls). Default 0.
        include_acceleration (bool): If True, the reversible static-pressure change due to the
                         velocity change (rho/2 (v_ex^2 - v_su^2)) is added. Default False.
        void_fraction_model (str): 'Zivi' (default), 'Homogeneous' or any model of VOID_FRACTION_MODELS
                                   that does not need the mass flux (not in NEEDS_FLOW).
        dp_corr (dict): Pressure-drop correlation names (keys of DP_CORRELATIONS in
                        dp_corr_distributor) for the inlet body, per regime.
                        Default {'single-phase': 'Churchill', 'two-phase': 'MSH'}.

    **Inputs**:

        fluid (str): Inlet fluid.
        T_su (float): Inlet temperature [K].
        h_su (float): Inlet specific enthalpy [J/kg].
        P_su (float): Inlet pressure [Pa].
        m_dot (float): Inlet mass flow rate [kg/s].

    **Outputs**:

        For each outlet connector (ex_1, ..., ex_n):
            - fluid: Same as inlet fluid.
            - P: Inlet pressure minus pressure drop of that branch [Pa]
                 (equal to inlet pressure if no geometry).
            - h: Same as inlet enthalpy [J/kg].
            - m_dot: su_m_dot * outlet_repartition[i] [kg/s].

        DP (np.ndarray or None): Pressure drop per outlet [Pa].
        DP_detail (dict or None): Breakdown (inlet body, transition, outlet tubes, acceleration).
        charge (float or None): Refrigerant charge in the component [kg].
        charge_detail (dict or None): Charge breakdown per section [kg].
    """

    DP_KEYS = ('D_su', 'L_su', 'D_ex', 'L_ex')

    def __init__(self, outlet_repartition):
 
        super().__init__()
 
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.solved = False
 
        "Flow split (fixes the number of outlets)"
        self.outlet_repartition = outlet_repartition
 
        "Input"
        self.su = MassConnector()
 
        "Outputs"
        for i in range(len(self.outlet_repartition)):
            outlet_num = i + 1
            setattr(self, f"ex_{outlet_num}", MassConnector())
 
        "Optional geometry: default to None (= not available)"
        for name in self.get_required_parameters():
            self.params.setdefault(name, None)
 
        "Results"
        self._reset_results()
    # ------------------------------------------------------------------ #
    # Framework methods
    # ------------------------------------------------------------------ #
 
    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su', 'T_su', 'm_dot', 'fluid']
 
    def get_required_parameters(self):
        # Geometry listed as for the other models. These parameters are
        # initialised to None in __init__: a parameter left to None means
        # "not available" and the corresponding computation is skipped.
        return ['D_su', 'L_su', 'D_ex', 'L_ex', 'L_trans']

    def check_parametrized(self):
        # The geometry is optional: the splitter is always parametrized.
        # Missing geometry only means that the pressure drops and/or the charge are not computed.
        self.parametrized = True
        if self.print_flag:
            missing = [k for k in self.get_required_parameters() if self.params.get(k) is None]
            if missing:
                print(f"Component {self.__class__.__name__}: missing parameters {', '.join(missing)}. "
                      f"Pressure drop computed: {self._can_compute_dp()}. The charge may not be computed.")
        return self.parametrized


    def solve(self):

        if not self.check_parametrized():
            raise ValueError("'Spliter' is not fully parametrized.")
        if not self.check_calculable():
            raise ValueError("'Spliter' is not fully calculable.")

        tolerance = 1e-3
        if abs(sum(self.outlet_repartition) - 1) >= tolerance:
            raise ValueError(f"'Spliter' total outlet repartition is too far from 1 (tolerance : {tolerance}).")

        self._reset_results()   # Used to clean if there were previous values associated to the results
        n = len(self.outlet_repartition)

        "1) Pressure drop (optional)"
        p_ex = np.full(n, float(self.su.p))
        if self._can_compute_dp():
            self.DP = self._compute_pressure_drops()
            p_ex = self.su.p - self.DP
            if np.any(p_ex <= 0):
                raise ValueError("'Spliter' computed pressure drop exceeds inlet pressure.")

        "2) Outlet states"
        for i in range(n):
            connector = getattr(self, f"ex_{i + 1}")
            connector.reset()
            connector.set_fluid(self.su.fluid)
            connector.set_p(p_ex[i])
            connector.set_h(self.su.h)
            connector.set_m_dot(self.su.m_dot * self.outlet_repartition[i])

        "3) Refrigerant charge (optional)"
        if self._can_compute_charge():
            self.charge = self._compute_charge(p_ex)

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

    def _per_outlet(self, value, name):
        """Broadcast a scalar to all outlets, or check a list has one value per outlet."""
        n = len(self.outlet_repartition)
        arr = np.atleast_1d(np.asarray(value, dtype=float))
        if arr.size == 1:
            return np.full(n, arr[0])
        if arr.size != n:
            raise ValueError(f"'Spliter' parameter '{name}' must be a scalar or have {n} values.")
        return arr

    def _can_compute_dp(self):
        return all(self.params.get(k) is not None for k in self.DP_KEYS)

    def _can_compute_charge(self):
        detailed = self._can_compute_dp() and (
            self.params.get('L_trans') is not None or self.params.get('V_trans') is not None)
        return detailed or self.params.get('V_int') is not None

    def _geometry(self):
        D_su = float(self.params.get('D_su'))
        D_ex = self._per_outlet(self.params.get('D_ex'), 'D_ex')
        A_su = np.pi * D_su**2 / 4
        A_ex = np.pi * D_ex**2 / 4
        D_eq_ex = np.sqrt(4 * A_ex.sum() / np.pi)  # diameter of one tube with the total outlet area
        return D_su, D_ex, A_su, A_ex, D_eq_ex

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
        AS = CP.AbstractState('HEOS', fluid)   # Local state at (p, h): the connector's AS is not modified
        AS.update(CP.HmassP_INPUTS, h, p)
        alpha = compute_void_fraction(AS, {}, void_fraction_model=void_fraction_model)   # No geometry: models of NEEDS_FLOW not available
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
        sigma = A_outlets_total / A_inlet, theta = total included angle [rad].
        Returns (K, basis) where basis tells which velocity K refers to
        (always the smaller cross-section).
        Ref: Flow of Fluids Through Valves, Fittings, and Pipe, Crane Valves North America, Technical Paper No. 410M. 1979.
        Note : Used in EES and the fluid python library
        """
        if sigma < 1:  # contraction -> outlet (small) velocity
            if theta <= np.pi / 4:
                K = 0.8 * np.sin(theta / 2) * (1 - sigma)
            else:
                K = 0.5 * np.sqrt(np.sin(theta / 2)) * (1 - sigma)
            return K, 'outlet'

        beta2 = 1 / sigma  # expansion -> inlet (small) velocity
        if theta <= np.pi / 4:
            K = 2.6 * np.sin(theta / 2) * (1 - beta2)**2
        else:
            K = (1 - beta2)**2
        return K, 'inlet'

    def _compute_pressure_drops(self):
        fluid, p, h, m_dot = self.su.fluid, self.su.p, self.su.h, self.su.m_dot
        rep = np.asarray(self.outlet_repartition, dtype=float)

        D_su, D_ex, A_su, A_ex, D_eq_ex = self._geometry()
        L_su = float(self.params.get('L_su'))
        L_ex = self._per_outlet(self.params.get('L_ex'), 'L_ex')
        L_trans = self.params.get('L_trans')
        eps = float(self.params.get('roughness', 1.5e-6))
        K_extra = float(self.params.get('K_extra', 0.0))
        include_acc = bool(self.params.get('include_acceleration', False))
        dp_corr = self.params.get('dp_corr') or {'single-phase': 'Churchill', 'two-phase': 'MSH'}
        void_fraction_model = self.params.get('void_fraction_model', 'Zivi')

        st = self._fluid_state(fluid, p, h, void_fraction_model)
        rho, mu = st['rho'], st['mu']

        "Velocities"
        m_dot_ex = m_dot * rep             # Vector with the different outlet mass flow rates
        v_su = m_dot / (rho * A_su)        # Velocity in the inlet body
        v_ex = m_dot_ex / (rho * A_ex)     # Vector with the different outlet velocities

        "Friction correlation: same for inlet body and outlet tubes (properties at inlet conditions)"
        corr = dp_corr['two-phase' if st['two_phase'] else 'single-phase']

        "1) Inlet body: friction (correlation from dp_corr_distributor)"
        Re_su = rho * v_su * D_su / mu      # Only reported in DP_detail
        geom_su = {'A_flow': A_su,          # Flow section of the inlet body [m^2]
                   'D_h': D_su,             # Hydraulic diameter [m]
                   'roughness': eps,        # Absolute roughness [m]
                   'chevron_angle': None}   # Not a plate: plate-only correlations raise an error
        dp_su = compute_dp(self.su.AS, m_dot, geom_su, L_su, corr, void_fraction_model)

        "2) Transition: local loss (contraction/expansion) + optional extra coefficient"
        # Pressure drop due to the transition is computed based on the area ratio and the included angle.
        sigma = A_ex.sum() / A_su   # Area ratio (total outlet area / inlet area)
        if L_trans is None or L_trans <= 0:
            theta = np.pi  # sudden change of section
        else:
            theta = 2 * np.arctan(abs(D_su - D_eq_ex) / (2 * L_trans))
        K, basis = self._transition_K(sigma, theta,)
        v_ref = v_ex if basis == 'outlet' else np.full_like(v_ex, v_su)     # np.full_like to broadcast to the same shape as v_ex
        dp_trans = K * rho * v_ref**2 / 2 + K_extra * rho * v_ex**2 / 2

        "3) Outlet tubes: friction (correlation from dp_corr_distributor, one tube per outlet)"
        Re_ex = rho * v_ex * D_ex / mu      # Only reported in DP_detail
        dp_ex = np.zeros(len(D_ex))
        for i in range(len(D_ex)):
            geom_ex = {'A_flow': A_ex[i],         # Flow section of outlet tube i [m^2]
                       'D_h': D_ex[i],            # Hydraulic diameter [m]
                       'roughness': eps,          # Absolute roughness [m]
                       'chevron_angle': None}     # Not a plate: plate-only correlations raise an error
            dp_ex[i] = compute_dp(self.su.AS, m_dot_ex[i], geom_ex, L_ex[i], corr, void_fraction_model)

        "4) Optional reversible static-pressure change due to acceleration"
        dp_acc = rho / 2 * (v_ex**2 - v_su**2) if include_acc else np.zeros_like(v_ex)

        self.DP_detail = {
            'inlet_body': dp_su,
            'transition': dp_trans,
            'outlet_tubes': dp_ex,
            'acceleration': dp_acc,
            'K_transition': K,
            'theta_deg': np.degrees(theta),
            'area_ratio': sigma,
            'Re_su': Re_su,
            'Re_ex': Re_ex,
            'two_phase': st['two_phase'],
        }
        return dp_su + dp_trans + dp_ex + dp_acc

    # ------------------------------------------------------------------ #
    # Refrigerant charge
    # ------------------------------------------------------------------ #

    def _compute_charge(self, p_ex):
        void_fraction_model = self.params.get('void_fraction_model', 'Zivi')
        st_su = self._fluid_state(self.su.fluid, self.su.p, self.su.h, void_fraction_model)
        rho_su = st_su['rho']

        "Only a total volume is known"
        if not self._can_compute_dp() or (self.params.get('L_trans') is None
                                          and self.params.get('V_trans') is None):
            V_int = float(self.params.get('V_int'))
            self.charge_detail = {'total_volume': rho_su * V_int}
            return rho_su * V_int

        "Detailed geometry"
        D_su, D_ex, A_su, A_ex, D_eq_ex = self._geometry()
        L_su = float(self.params.get('L_su'))
        L_ex = self._per_outlet(self.params.get('L_ex'), 'L_ex')
        L_trans = self.params.get('L_trans')

        V_su = A_su * L_su
        V_trans = self.params.get('V_trans')
        if V_trans is None:
            L_trans = float(self.params.get('L_trans'))
            V_trans = np.pi * L_trans / 12 * (D_su**2 + D_su * D_eq_ex + D_eq_ex**2)  # cone frustum
        V_ex = A_ex * L_ex

        st_ex = [self._fluid_state(self.su.fluid, pi, self.su.h, void_fraction_model) for pi in p_ex]
        rho_ex = np.array([st['rho'] for st in st_ex])

        charge_su = rho_su * V_su
        charge_trans = rho_su * float(V_trans)
        charge_ex = rho_ex * V_ex

        self.charge_detail = {
            'inlet_body': charge_su,
            'transition': charge_trans,
            'outlet_tubes': charge_ex,
            'volumes': {'inlet_body': V_su, 'transition': float(V_trans), 'outlet_tubes': V_ex},
        }
        return charge_su + charge_trans + float(charge_ex.sum())