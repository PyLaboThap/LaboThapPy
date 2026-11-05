# -*- coding: utf-8 -*-
"""
Tee junction (diverging "té") with pressure drop and refrigerant charge.

Pressure-drop coefficients identical to the EES Minor Losses library:
    K_tee_straight(RelRough) = 20 * f_T     (flow through the run)
    K_tee_branched(RelRough) = 60 * f_T     (flow through the branch)
where f_T is the fully turbulent Darcy friction factor.

    - Crane Co., "Flow of Fluids Through Valves, Fittings, and Pipe",
      Technical Paper No. 410M, 1979, p. A-29 (reference cited by EES).
    - F-Chart Software, EES Minor Losses library, help pages
      https://fchartsoftware.com/ees/minor_losses/hs6010.htm (K_tee_straight)
      https://fchartsoftware.com/ees/minor_losses/hs6020.htm (K_tee_branched)

f_T: Churchill (1977) correlation evaluated at Re = 1e8. This reproduces the
EES example values to the 4 digits printed in the EES help (K_tee_straight=0.1206,
K_tee_branched=0.3617, K_45deg_elbow=0.09644, K_90deg_elbow=0.1808,
K_ReturnBend=0.3014 for RelRough=0). EES does not document how it evaluates
f_T: this choice was found by matching those published values.
    - Churchill, S.W., "Friction-factor equation spans all fluid-flow regimes",
      Chemical Engineering 84(24), 91-92, 1977.

Reducing tee (D_branch != D): Crane approach, i.e. tee coefficient of the run size
plus a sudden contraction/enlargement (same Crane TP-410 equations as in
TankSpliter._transition_K, theta = 180 deg).
"""

import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

from labothappy.correlations.properties.two_phase import compute_two_phase_density
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction

from labothappy.connector.mass_connector import MassConnector
from labothappy.component.base_component import BaseComponent


# ============================================================================
# EES / CRANE TP-410 CORRELATIONS
# ============================================================================

RE_FULLY_TURBULENT = 1e8   # Reynolds number used to evaluate f_T (matches EES outputs)
K_RUN_FACTOR = 20.0        # Crane TP-410 (1979) p. A-29: standard tee, flow through run
K_BRANCH_FACTOR = 60.0     # Crane TP-410 (1979) p. A-29: standard tee, flow through branch


def friction_factor_churchill(Re, rel_rough):
    """
    Darcy friction factor, Churchill (1977), valid for all regimes.
        f = 8 * [ (8/Re)^12 + (A + B)^(-3/2) ]^(1/12)
        A = [ 2.457 * ln( 1 / ((7/Re)^0.9 + 0.27 * eps/D) ) ]^16
        B = (37530 / Re)^16
    """
    A = (2.457 * np.log(1.0 / ((7.0 / Re)**0.9 + 0.27 * rel_rough)))**16
    B = (37530.0 / Re)**16
    return 8.0 * ((8.0 / Re)**12 + (A + B)**(-1.5))**(1.0 / 12.0)


def friction_factor_fully_turbulent(rel_rough):
    """f_T: Darcy friction factor in the fully turbulent zone (Churchill at Re = 1e8)."""
    return friction_factor_churchill(RE_FULLY_TURBULENT, rel_rough)


def K_tee_straight(rel_rough):
    """EES K_tee_straight(RelRough): standard tee, flow straight through the run."""
    return K_RUN_FACTOR * friction_factor_fully_turbulent(rel_rough)


def K_tee_branched(rel_rough):
    """EES K_tee_branched(RelRough): standard tee, flow through the branch."""
    return K_BRANCH_FACTOR * friction_factor_fully_turbulent(rel_rough)


def K_sudden_area_change(beta):
    """
    Crane TP-410 sudden contraction / enlargement, referred to the velocity in the
    SMALLER pipe (beta = D_small / D_large).
        contraction : K = 0.5 * (1 - beta^2)
        enlargement : K = (1 - beta^2)^2
    Returns (K_contraction, K_enlargement).
    """
    return 0.5 * (1.0 - beta**2), (1.0 - beta**2)**2


# ============================================================================
# LABOTHAPPY COMPONENT
# ============================================================================

class TeeJunction(BaseComponent):
    """
    **Component**: Tee junction (diverging "té") – single and two-phase

    **Model**:
        One inlet (su) splits into a run outlet (ex_1 = ex_run, straight through)
        and a branch outlet (ex_2 = ex_branch, 90 deg). Local pressure losses with
        the EES / Crane TP-410 coefficients:

            dP_run    = K_run    * rho * v_ref_run^2    / 2,  K_run    = 20 f_T
            dP_branch = K_branch * rho * v_ref_branch^2 / 2,  K_branch = 60 f_T

        Two-phase: homogeneous model for the pressure drop (rho = rho_h),
        void-fraction model (Zivi by default) for the charge.

    **Assumptions**:
        - Steady-state, adiabatic, no work, no elevation change (h_ex = h_su).
        - Turbulent flow (the Crane K values are fully turbulent values).
        - Properties evaluated at inlet conditions for the pressure drop.
        - Flow split imposed by `outlet_repartition` (not predicted).

    **Connectors**:
        su (MassConnector)                 : inlet.
        ex_1 = ex_run (MassConnector)      : run outlet (straight through).
        ex_2 = ex_branch (MassConnector)   : branch outlet (90 deg).

    **Parameters**:
        D (float)                : Inner diameter of the tee (run) [m]. Required.
        D_branch (float)         : Inner diameter of the branch [m]. Default D.
                                   If different from D, a sudden contraction/enlargement
                                   is added on the branch path (Crane).
        roughness (float)        : Absolute roughness [m]. Default 0 (smooth, as in the EES examples).
        velocity_basis (str)     : 'combined' (default): both K use the inlet (combined-flow)
                                   velocity, as in the EES examples with one tee velocity.
                                   'own': each K uses the velocity of its own outlet stream.
        K_extra_run (float)      : Additional loss coefficient on the run (own velocity). Default 0.
        K_extra_branch (float)   : Additional loss coefficient on the branch (own velocity). Default 0.
        include_acceleration (bool): Add the reversible static-pressure change
                                   rho/2 (v_ex^2 - v_su^2). Default False.
        L_run (float)            : Internal length of the run (inlet -> run outlet) [m]. For the charge.
        L_branch (float)         : Internal length of the branch [m]. For the charge.
        V_int (float)            : Total internal volume [m^3]. Used for the charge if L_run/L_branch missing.
        void_fraction_model (str): 'Zivi' (default), 'Homogeneous', ... (not in NEEDS_FLOW).

    **Inputs**:
        P_su (float)  : Inlet pressure [Pa].
        h_su (float)  : Inlet specific enthalpy [J/kg].
        m_dot (float) : Inlet mass flow rate [kg/s].
        fluid (str)   : Fluid name.

    **Outputs**:
        ex_run / ex_branch: P, h (= h_su), m_dot.
        DP (np.ndarray)      : [dP_run, dP_branch] [Pa].
        DP_detail (dict)     : Breakdown (tee, reducer, extra, acceleration, K, f_T, Re...).
        charge (float/None)  : Refrigerant charge [kg].
        charge_detail (dict) : Charge breakdown [kg].
    """

    OUTLETS = ('run', 'branch')

    def __init__(self, outlet_repartition=(0.5, 0.5)):

        super().__init__()

        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.solved = False

        "Flow split: [fraction to the run, fraction to the branch]"
        if len(outlet_repartition) != 2:
            raise ValueError("'TeeJunction' outlet_repartition must have 2 values: [run, branch].")
        self.outlet_repartition = list(outlet_repartition)

        "Input"
        self.su = MassConnector()

        "Outputs (ex_1/ex_2 as in the splitter, with explicit aliases)"
        self.ex_1 = MassConnector()
        self.ex_2 = MassConnector()
        self.ex_run = self.ex_1
        self.ex_branch = self.ex_2

        "Results"
        self._reset_results()

    # ------------------------------------------------------------------ #
    # Framework methods
    # ------------------------------------------------------------------ #

    def get_required_inputs(self):
        return ['P_su', 'h_su', 'm_dot', 'fluid']

    def get_required_parameters(self):
        return ['D']

    def solve(self):

        if not self.check_parametrized():
            raise ValueError("'TeeJunction' is not fully parametrized.")
        if not self.check_calculable():
            raise ValueError("'TeeJunction' is not fully calculable.")

        tolerance = 1e-3
        if abs(sum(self.outlet_repartition) - 1) >= tolerance:
            raise ValueError(f"'TeeJunction' total outlet repartition is too far from 1 (tolerance : {tolerance}).")
        if min(self.outlet_repartition) < 0:
            raise ValueError("'TeeJunction' outlet repartition must be positive.")

        self._reset_results()

        "1) Pressure drop"
        self.DP = self._compute_pressure_drops()
        p_ex = self.su.p - self.DP
        if np.any(p_ex <= 0):
            raise ValueError("'TeeJunction' computed pressure drop exceeds inlet pressure.")

        "2) Outlet states"
        for i, connector in enumerate((self.ex_1, self.ex_2)):
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
    # Helpers
    # ------------------------------------------------------------------ #

    def _reset_results(self):
        self.DP = None
        self.DP_detail = None
        self.charge = None
        self.charge_detail = None
        self.Re = None
        self.low_Re_warning = False

    def _can_compute_charge(self):
        detailed = self.params.get('L_run') is not None and self.params.get('L_branch') is not None
        return detailed or self.params.get('V_int') is not None

    def _geometry(self):
        D = float(self.params['D'])
        D_b = self.params.get('D_branch')
        D_b = D if D_b is None else float(D_b)
        A = np.pi * D**2 / 4
        A_b = np.pi * D_b**2 / 4
        return D, D_b, A, A_b

    # ------------------------------------------------------------------ #
    # Fluid properties (same approach as TankSpliter)
    # ------------------------------------------------------------------ #

    @staticmethod
    def _fluid_state(fluid, p, h, void_fraction_model):
        """
        Single phase: rho, mu at (p, h).
        Two-phase   : rho_dp  = homogeneous density (pressure drop),
                      rho     = void-fraction weighted density (charge),
                      mu      = McAdams viscosity.
        """
        x = None
        try:
            q = PropsSI('Q', 'P', p, 'H', h, fluid)
            if 0.0 <= q <= 1.0:
                x = q
        except ValueError:
            pass  # e.g. incompressible fluids: treated as single phase

        if x is None:
            rho = PropsSI('D', 'P', p, 'H', h, fluid)
            return {'two_phase': False, 'rho': rho, 'rho_dp': rho,
                    'mu': PropsSI('V', 'P', p, 'H', h, fluid)}

        rho_l = PropsSI('D', 'P', p, 'Q', 0, fluid)
        rho_v = PropsSI('D', 'P', p, 'Q', 1, fluid)
        mu_l = PropsSI('V', 'P', p, 'Q', 0, fluid)
        mu_v = PropsSI('V', 'P', p, 'Q', 1, fluid)

        AS = CP.AbstractState('HEOS', fluid)   # Local state: the connector's AS is not modified
        AS.update(CP.HmassP_INPUTS, h, p)
        alpha = compute_void_fraction(AS, {}, void_fraction_model=void_fraction_model)

        return {'two_phase': True, 'x': x,
                'rho_l': rho_l, 'rho_v': rho_v, 'alpha': alpha,
                'rho': compute_two_phase_density(x, rho_l, rho_v, alpha),
                'rho_dp': 1.0 / (x / rho_v + (1.0 - x) / rho_l),
                'mu': 1.0 / (x / mu_v + (1.0 - x) / mu_l)}

    # ------------------------------------------------------------------ #
    # Pressure drop
    # ------------------------------------------------------------------ #

    def _compute_pressure_drops(self):
        fluid, p, h, m_dot = self.su.fluid, self.su.p, self.su.h, self.su.m_dot
        rep = np.asarray(self.outlet_repartition, dtype=float)

        D, D_b, A, A_b = self._geometry()
        eps = float(self.params.get('roughness') or 0.0)
        basis = self.params.get('velocity_basis') or 'combined'
        K_extra = np.array([float(self.params.get('K_extra_run') or 0.0),
                            float(self.params.get('K_extra_branch') or 0.0)])
        include_acc = bool(self.params.get('include_acceleration', False))
        void_fraction_model = self.params.get('void_fraction_model') or 'Zivi'

        if basis not in ('combined', 'own'):
            raise ValueError("'TeeJunction' velocity_basis must be 'combined' or 'own'.")

        st = self._fluid_state(fluid, p, h, void_fraction_model)
        rho, mu = st['rho_dp'], st['mu']

        "Velocities"
        m_dot_ex = m_dot * rep                          # [m_run, m_branch]
        v_su = m_dot / (rho * A)                        # combined-flow velocity in the run
        v_own = m_dot_ex / (rho * np.array([A, A_b]))   # [v_run, v_branch]
        v_ref = np.full(2, v_su) if basis == 'combined' else v_own

        "1) Tee coefficients (EES K_tee_straight / K_tee_branched, Crane TP-410 p. A-29)"
        rel_rough = eps / D
        f_T = friction_factor_fully_turbulent(rel_rough)
        K_tee = np.array([K_RUN_FACTOR * f_T, K_BRANCH_FACTOR * f_T])
        dp_tee = K_tee * rho * v_ref**2 / 2

        "2) Reducing tee: sudden contraction/enlargement on the branch path (Crane)"
        dp_red = np.zeros(2)
        K_red = 0.0
        if not np.isclose(D_b, D):
            beta = min(D, D_b) / max(D, D_b)
            K_contr, K_enl = K_sudden_area_change(beta)
            if D_b < D:    # run -> smaller branch: contraction, branch velocity
                K_red = K_contr
                v_small = v_own[1]
            else:          # run -> larger branch: enlargement, velocity of the branch flow at diameter D
                K_red = K_enl
                v_small = m_dot_ex[1] / (rho * A)
            dp_red[1] = K_red * rho * v_small**2 / 2
            if self.print_flag:
                print("TeeJunction: D_branch != D, the EES/Crane tee coefficients are for line-size tees. "
                      "A sudden contraction/enlargement is added on the branch path.")

        "3) Additional user coefficients (own velocities)"
        dp_extra = K_extra * rho * v_own**2 / 2

        "4) Optional reversible static-pressure change due to acceleration"
        dp_acc = rho / 2 * (v_own**2 - v_su**2) if include_acc else np.zeros(2)

        "Reynolds number (validity check: Crane K values are for turbulent flow)"
        self.Re = rho * v_su * D / mu
        Re_own = rho * v_own * np.array([D, D_b]) / mu
        if self.Re < 4000:
            self.low_Re_warning = True
            if self.print_flag:
                print(f"TeeJunction: Re = {self.Re:.0f} < 4000. The Crane/EES tee coefficients are "
                      f"fully turbulent values and underestimate the losses in laminar/transitional flow.")

        self.DP_detail = {
            'tee': dp_tee,
            'reducer': dp_red,
            'extra': dp_extra,
            'acceleration': dp_acc,
            'K_run': K_tee[0],
            'K_branch': K_tee[1],
            'K_reducer': K_red,
            'f_T': f_T,
            'rel_rough': rel_rough,
            'velocity_basis': basis,
            'v_su': v_su,
            'v_ex': v_own,
            'rho_dp': rho,
            'Re_su': self.Re,
            'Re_ex': Re_own,
            'two_phase': st['two_phase'],
        }
        return dp_tee + dp_red + dp_extra + dp_acc

    # ------------------------------------------------------------------ #
    # Refrigerant charge
    # ------------------------------------------------------------------ #

    def _compute_charge(self, p_ex):
        fluid, h = self.su.fluid, self.su.h
        void_fraction_model = self.params.get('void_fraction_model') or 'Zivi'
        rho_su = self._fluid_state(fluid, self.su.p, h, void_fraction_model)['rho']

        "Only a total volume is known"
        if self.params.get('L_run') is None or self.params.get('L_branch') is None:
            V_int = float(self.params.get('V_int'))
            self.charge_detail = {'total_volume': rho_su * V_int, 'volumes': {'total': V_int}}
            return rho_su * V_int

        "Detailed geometry: run (inlet -> run outlet) + branch"
        D, D_b, A, A_b = self._geometry()
        V_run = A * float(self.params['L_run'])
        V_branch = A_b * float(self.params['L_branch'])

        rho_run = self._fluid_state(fluid, p_ex[0], h, void_fraction_model)['rho']
        rho_branch = self._fluid_state(fluid, p_ex[1], h, void_fraction_model)['rho']

        charge_run = V_run * (rho_su + rho_run) / 2     # mean density along the run (as in CurvedElbow)
        charge_branch = V_branch * (rho_su + rho_branch) / 2

        self.charge_detail = {
            'run': charge_run,
            'branch': charge_branch,
            'volumes': {'run': V_run, 'branch': V_branch},
        }
        return charge_run + charge_branch

    # ------------------------------------------------------------------ #
    # Results
    # ------------------------------------------------------------------ #

    def print_results(self):
        d = self.DP_detail
        print("=" * 70)
        print("TEE JUNCTION RESULTS")
        print("=" * 70)
        print(f"Diameter D (run)   : {self.params['D']*1e3:.2f} mm")
        print(f"Split [run, branch]: {self.outlet_repartition}")
        print(f"f_T (Re = 1e8)     : {d['f_T']:.5f}   (eps/D = {d['rel_rough']:.2e})")
        print(f"K_run  = 20 f_T    : {d['K_run']:.4f}")
        print(f"K_branch = 60 f_T  : {d['K_branch']:.4f}")
        print(f"Velocity basis     : {d['velocity_basis']}")
        print(f"Inlet velocity     : {d['v_su']:.2f} m/s")
        print(f"Reynolds (inlet)   : {d['Re_su']:.0f}" + ("   [WARNING: low Re]" if self.low_Re_warning else ""))
        print(f"Two-phase          : {d['two_phase']}  (rho for dP = {d['rho_dp']:.2f} kg/m3)")
        for i, name in enumerate(self.OUTLETS):
            print(f"\n--- {name.upper()} ---")
            print(f"m_dot              : {self.su.m_dot * self.outlet_repartition[i]:.4f} kg/s")
            print(f"dP tee             : {d['tee'][i]:.2f} Pa")
            print(f"dP reducer         : {d['reducer'][i]:.2f} Pa")
            print(f"dP extra           : {d['extra'][i]:.2f} Pa")
            print(f"dP acceleration    : {d['acceleration'][i]:.2f} Pa")
            print(f"dP total           : {self.DP[i]:.2f} Pa")
        if self.charge is not None:
            print(f"\nRefrigerant charge : {self.charge:.6f} kg")
        print("=" * 70)