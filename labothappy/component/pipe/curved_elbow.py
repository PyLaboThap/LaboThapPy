import math
import CoolProp.CoolProp as CP
import numpy as np

from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector

from correlations. pressure_drop.local_losses.curved_elbow import pressure_drop_curved_elbow
from correlations.properties.dimensionless import compute_reynolds
from correlations.properties.void_fraction import compute_void_fraction

# ============================================================================
# MODELICA DISSIPATION LIBRARY PORT: Curved Bends
# ============================================================================
# Sources:
#   - Idelchik, I.E.: HANDBOOK OF HYDRAULIC RESISTANCE, 3rd edition, 2006
#   - Miller, D.S.: INTERNAL FLOW SYSTEMS, 2nd edition, 1984
#   - VDI-Waermeatlas, 9th edition, Springer-Verlag, 2002, Section Lac 6
# ============================================================================

PI = math.pi
EPS = 1e-12

# ============================================================================
# LABOTHAPPY COMPONENT
# ============================================================================

class CurvedElbow(BaseComponent):
    """
    **Component**: Elbow (Curved Bend) – Single and Two-Phase

    **Model**:
        Local pressure loss through a curved pipe bend using Idelchik correlation.
        For two-phase flows, uses homogeneous density approach.

    **Assumptions**:
        - Adiabatic
        - No elevation change
        - Constant diameter
        - For two-phase: homogeneous flow model (continuous phase distribution)

    **Parameters**
    ----------
    D : float
        Pipe inner diameter [m]
    delta : float
        Bend angle [deg]
    R0 : float
        Curvature radius [m]
    K : float, optional
        Absolute surface roughness [m]. Default: 0 (smooth pipe)

    **Inputs**
    ----------
    P_su : float
        Inlet pressure [Pa]
    h_su : float
        Inlet specific enthalpy [J/kg]
    m_dot : float
        Mass flow rate [kg/s]
    fluid : str
        Fluid name (e.g., 'R410A', 'Water')

    **Outputs**
    ----------
    P_ex : float
        Outlet pressure [Pa]
    h_ex : float
        Outlet specific enthalpy [J/kg]
    m_charge : float
        Refrigerant mass inventory [kg]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()
        self.ex = MassConnector()

    def get_required_inputs(self):
        return ['P_su', 'h_su', 'm_dot', 'fluid']

    def get_required_parameters(self):
        return ['D', 'delta', 'R0']

    def solve(self):
        """Main solve method: detect phase and dispatch to appropriate solver."""
        self.check_calculable()
        self.check_parametrized()

        # Get inlet state
        self.AS = CP.AbstractState('HEOS', self.su.fluid)
        self.AS.update(CP.HmassP_INPUTS, self.su.h, self.su.p)
        x = self.AS.Q()

        # Dispatch based on phase
        if x > EPS and x < 1.0 - EPS:
            self._solve_two_phase(x)
        else:
            self._solve_single_phase()

    def _solve_single_phase(self):
        """Solve single-phase pressure drop and charge."""
        rho_su = self.AS.rhomass()
        mu = self.AS.viscosity()

        A_cross = PI * self.params['D']**2 / 4 # Cross-sectional area [m²]
        v = self.su.m_dot / (rho_su * A_cross) # Mean velocity [m/s]

        # Bend geometry
        delta_rad = self.params['delta'] * np.pi / 180
        R0 = self.params['R0']
        K = self.params.get('K', 0.0)

        # Pressure drop
        self.dP = pressure_drop_curved_elbow(self.params['D'], R0, delta_rad, K, rho_su, mu, self.su.m_dot)

        # Set outlet
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_h(self.su.h)
        self.ex.set_p(self.su.p - self.dP)

        # Get outlet density
        rho_ex = self.ex.D

        # Charge using mean density
        rho_mean = (rho_su + rho_ex) / 2.0
        L = abs(delta_rad) * R0  # Arc length
        self.m_charge = A_cross * L * rho_mean

        # Store diagnostics
        self.is_two_phase = False       
        self.rho_l = None
        self.rho_g = None
        self.rho_h = None
        self.velocity = v
        self.Re = compute_reynolds(rho_su, v, self.params['D'], mu)

        self.solved = True

    def _solve_two_phase(self, x):
        """Solve two-phase pressure drop and charge using homogeneous model."""
        # Get saturation properties
        AS_sat_l = CP.AbstractState('HEOS', self.su.fluid)
        AS_sat_l.update(CP.PQ_INPUTS, self.su.p, 0.0)
        
        AS_sat_g = CP.AbstractState('HEOS', self.su.fluid)
        AS_sat_g.update(CP.PQ_INPUTS, self.su.p, 1.0)
        
        rho_l = AS_sat_l.rhomass()
        rho_g = AS_sat_g.rhomass()
        mu = AS_sat_l.viscosity()

        # Homogeneous density
        alpha = compute_void_fraction(x, rho_l, rho_g)
        rho_h = 1.0 / (alpha / rho_g + (1.0 - alpha) / rho_l)

        A_cross = PI * self.params['D']**2 / 4 # Cross-sectional area [m²]
        v = self.su.m_dot / (rho_h * A_cross) # Mean velocity [m/s]

        # Bend geometry
        delta_rad = self.params['delta'] * np.pi / 180
        R0 = self.params['R0']
        K = self.params.get('K', 0.0)

        # Pressure drop
        self.dP = pressure_drop_curved_elbow(self.params['D'], R0, delta_rad, K, rho_h, mu, self.su.m_dot)

        # Set outlet
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_h(self.su.h)
        self.ex.set_p(self.su.p - self.dP)

        # Charge using homogeneous density (constant along short bend)
        L = abs(delta_rad) * R0  # Arc length
        self.m_charge = A_cross * L * rho_h

        # Store diagnostics
        self.is_two_phase = True
        self.rho_l = rho_l
        self.rho_g = rho_g
        self.rho_h = rho_h
        self.velocity = v
        self.Re = compute_reynolds(rho_h, v, self.params['D'], mu)

        self.solved = True

    def print_results(self):
        """Print summary of bend analysis."""
        print("=" * 70)
        print("ELBOW (CURVED BEND) RESULTS")
        print("=" * 70)
        print(f"Diameter D        : {self.params['D']*1e3:.2f} mm")
        print(f"Bend angle        : {self.params['angle']:.1f}°")
        print(f"Curvature ratio   : {self.params['R0_D']:.2f}")
        print(f"Velocity          : {self.velocity:.2f} m/s")
        print(f"Reynolds number   : {self.Re:.0f}")
        
        if self.is_two_phase:
            print(f"\n--- TWO-PHASE (HOMOGENEOUS MODEL) ---")
            print(f"ρ_liquid          : {self.rho_l:.2f} kg/m³")
            print(f"ρ_vapor           : {self.rho_g:.4f} kg/m³")
            print(f"ρ_homogeneous     : {self.rho_h:.2f} kg/m³")
        
        print(f"\nPressure drop     : {self.deltaP:.2f} Pa")
        print(f"Refrigerant charge: {self.m_charge:.6f} kg")
        print("=" * 70)