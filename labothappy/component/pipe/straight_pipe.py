import math
import CoolProp.CoolProp as CP
import numpy as np

from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector

from labothappy.correlations.pressure_drop.straight_pipe_DP import pressure_drop_pipe_single_phase
from labothappy.correlations.pressure_drop.straight_pipe_DP import pressure_drop_pipe_two_phase
from labothappy.correlations.properties.two_phase import compute_two_phase_density
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction

PI = math.pi
EPS = 1e-12

class StraightPipe(BaseComponent):
    """
    **Component**: Straight Pipe

    **Model**:
        Friction pressure loss in a straight circular pipe.
        Supports both single-phase and two-phase flows.

    **Description**:
        Single-phase: ΔP = f · (L/D) · (ρ·v²/2)
        Two-phase: Uses Friedel or MSH correlation with quality-dependent properties.

    **Assumptions**:
        - Adiabatic (no heat transfer)
        - Fully developed flow
        - Circular cross-section

    **Parameters**
    ----------
    D : float
        Pipe inner diameter [m]
    L : float
        Pipe length [m]
    K : float, optional
        Absolute surface roughness [m]. Default: 0 (smooth pipe)
    theta : float, optional
        Pipe inclination angle [deg]. Default: 0 (horizontal)
    two_phase_correlation : str, optional
        'friedel' or 'msh' (Müller-Steinhagen & Heck). Default: 'friedel'

    **Inputs**
    ----------
    P_su : float
        Inlet pressure [Pa]
    h_su : float
        Inlet specific enthalpy [J/kg]
    m_dot : float
        Mass flow rate [kg/s]
    fluid : str
        Fluid name (e.g., 'Water', 'Air', 'R134a')

    **Outputs**
    ----------
    P_ex : float
        Outlet pressure [Pa]
    h_ex : float
        Outlet specific enthalpy [J/kg] (same as inlet, adiabatic)
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()  # Inlet
        self.ex = MassConnector()  # Outlet

    def get_required_inputs(self):
        return ['P_su', 'h_su', 'm_dot', 'fluid']

    def get_required_parameters(self):
        return ['D', 'L', 'K', 'theta']

    def solve(self):
        """
        Solve for outlet pressure and enthalpy.
        
        Detects single-phase vs two-phase and applies appropriate correlation.
        """
        self.check_calculable()
        self.check_parametrized()

        # ====================================================================
        # Get fluid properties at inlet
        # ====================================================================
        self.AS = CP.AbstractState('HEOS', self.su.fluid)
        self.AS.update(CP.HmassP_INPUTS, self.su.h, self.su.p)

        # Check if two-phase
        x = self.AS.Q()  # Returns quality if two-phase, else 0 or 1
        
        if x > EPS and x < 1.0 - EPS:
            # ================================================================
            # TWO-PHASE FLOW
            # ================================================================
            self._solve_two_phase(x)
        else:
            # ================================================================
            # SINGLE-PHASE FLOW
            # ================================================================
            self._solve_single_phase()

        self.solved = True

    def _solve_single_phase(self):
        """Solve single-phase pressure drop."""
        rho_su = self.AS.rhomass()
        A_cross = PI * self.params['D']**2 / 4  # Cross-sectional area [m²]

        self.dP = pressure_drop_pipe_single_phase(self.AS, self.params, self.su.m_dot)

        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_h(self.su.h)
        self.ex.set_p(self.su.p - self.dP)

        self.rho_ex = self.ex.D

        # Charge inventory using mean density
        rho_mean = (rho_su + self.rho_ex) / 2.0
        self.m_charge = A_cross * self.params['L'] * rho_mean

        self.velocity = self.su.m_dot / (rho_su * A_cross)
        self.quality = None


    def _solve_two_phase(self, x):
        """Solve two-phase pressure drop using Friedel."""

        # Get saturation properties at inlet pressure (for charge inventory)
        AS_sat_l = CP.AbstractState('HEOS', self.su.fluid)
        AS_sat_l.update(CP.PQ_INPUTS, self.su.p, 0.0)  # Quality = 0 (saturated liquid)

        AS_sat_g = CP.AbstractState('HEOS', self.su.fluid)
        AS_sat_g.update(CP.PQ_INPUTS, self.su.p, 1.0)  # Quality = 1 (saturated vapor)

        rho_l = AS_sat_l.rhomass()
        rho_g = AS_sat_g.rhomass()

        # Compute void fraction and mass inventory at the inlet state, before
        # pressure_drop_pipe_two_phase advances self.AS to the outlet state.
        alpha = compute_void_fraction(self.AS, self.params, self.su.m_dot, void_fraction_model='Zivi')
        rho_tp = compute_two_phase_density(x, rho_l, rho_g, alpha)
        self.m_charge = rho_tp * self.params['L'] * (PI * self.params['D']**2 / 4)

        self.dP = pressure_drop_pipe_two_phase(self.AS, self.params, self.su.m_dot)

        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_h(self.su.h)
        self.ex.set_p(self.su.p - self.dP)

        self.quality = x
        self.void_fraction = alpha



    def print_results(self):
        """Print summary of straight pipe analysis."""
        print("=" * 60)
        print("STRAIGHT PIPE RESULTS")
        print("=" * 60)
        print(f"Diameter D        : {self.params['D']*1e3:.2f} mm")
        print(f"Length L          : {self.params['L']:.2f} m")
        
        if self.quality is None:
            # Single-phase
            print(f"\nSINGLE-PHASE FLOW")
            print(f"Velocity          : {self.velocity:.2f} m/s")
            print(f"Mass inventory    : {self.m_charge:.4f} kg")
        else:
            # Two-phase
            print(f"\nTWO-PHASE FLOW")
            print(f"Quality x         : {self.quality:.4f}")
            print(f"Void fraction α   : {self.void_fraction:.4f}")
            print(f"Mass inventory    : {self.m_charge:.4f} kg")
        
        print(f"Pressure drop ΔP  : {self.dP:.2f} Pa")
        print("=" * 60)
