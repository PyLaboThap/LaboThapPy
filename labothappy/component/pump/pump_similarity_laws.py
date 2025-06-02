"""
Created on Wed Nov 29 10:38:08 2023

@author: Elise Neven
@email: elise.neven@uliege.be
"""

# Internal imports
import __init__

from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector

# External imports
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import root_scalar

GRAVITY = 9.81  # m/s², Gravitational acceleration constant

class PumpSimilarityLaws(BaseComponent):
    
    """
    Component: Pump with a characteristic curve associated
    
    Model: Model using characteristic curves for head and power and similarity laws
    to calculate the pump performance at different speeds.

    Description:
        Simulates a pump using head and power performance curves (vs flow and speed).
        Uses the similarity laws to calculate the pump performance at different speeds.

    Assumptions:
        - Steady-state
        - Incompressible fluid approximation for head calculation

    Connectors:
        - su: Suction side (MassConnector)
        - ex: Exhaust side (MassConnector)
        - W_mec: Work connector (WorkConnector)

    Required Parameters:
        - curves_head: Dict = List[(flow, head)] dicionary of head curves
        - curves_power: Dict = List[(flow, power)] dictionary of power curves
        - curves_fluid: Fluid used for the characteristic curves
        - speed_ref: Reference speed for the curves (RPM)

    Required Inputs:
        - P_su: Suction pressure [Pa]
        - T_su: Suction temperature [K]
        - P_ex: Exhaust pressure [Pa]
        - N_rot: Rotational speed [RPM]
        - fluid: Actual working fluid

    Outputs:
        - V_dot_pp: Volumetric flow rate [m³/h] 
        - m_dot_pp: Mass flow rate [kg/s]
        - W_dot_pp: Shaft power required by the pump [W]
    """

    def __init__(self):
        """
        Initializes the class and sets up connectors for mass flow and work.
        """
        super().__init__()
        self.su = MassConnector()  # Suction connector
        self.ex = MassConnector()  # Exhaust connector
        self.W_mec = WorkConnector()  # Mechanical work connector

    def get_required_inputs(self):
        """
        Returns a list of required input parameters for the class.
        """
        return ['P_su', 'T_su', 'P_ex', 'N_rot', 'fluid']

    def get_required_parameters(self):
        """
        Returns a list of required parameters for the class.
        """
        return ['curves_head', 'curves_power', 'speed_ref', 'curves_fluid']
    
    def _prepare_interpolators(self):
        """
        Prepares polynomial interpolation functions for head and power curves
        based on the point from the characteristic curve provided for a single-speed pump.
        This method returns the fitted polynomial functions for head and power curves at teh reference speed.
        """
        # Extract the head and power curves from the input parameters
        head_curve = self.params['curves_head']
        power_curve = self.params['curves_power']

        # Extract flow and corresponding head values for the head curve
        flows_head, heads = zip(*sorted(head_curve))

        # Extract flow and corresponding power values for the power curve
        flows_power, powers = zip(*sorted(power_curve))

        self.flows_np = np.array(flows_head)  # Store flow values

        # Fit polynomial curves to both the head and power data (degree 3 or lower)
        deg = min(3, len(flows_head) - 1)  # Avoid overfitting if there are too few data points
        self.head_poly = np.poly1d(np.polyfit(flows_head, heads, deg))  # Polynomial curve for head at one speed of reference
        self.power_poly = np.poly1d(np.polyfit(flows_power, powers, deg))  # Polynomial curve for power at one speed of reference

    def _similarity_fluid(self):
        """
        Calculates head and power for a different fluid by adjusting based on density differences
        between the curve fluid and the actual operating fluid.
        """
        fluid_actual = self.inputs['fluid']  # Get the actual fluid type
        fluid_curve = self.params['curves_fluid']  # Reference fluid used for the characteristic curves

        # Calculate the density ratio between the actual fluid and the curve fluid
        rho_actual = CP.PropsSI('D', 'T', self.inputs['T_su'], 'P', self.inputs['P_su'], fluid_actual)
        rho_curve = CP.PropsSI('D', 'T', 20 + 273, 'P', 2e5, fluid_curve)
        density_ratio = rho_curve / rho_actual  # The ratio of the reference fluid's density to the actual fluid

        # Calculate pressure difference (Delta P) for both fluids
        DeltaP_actual = (self.inputs['P_ex'] - self.inputs['P_su'])  # Pressure difference across the pump
        DeltaP_curve = DeltaP_actual * density_ratio  # Pressure difference for the curve fluid (adjusted by density)

        # Calculate head values for both the actual and curve fluids
        head_actual = DeltaP_actual / (1000 * GRAVITY)  # Head for the actual fluid (converted to meters)
        head_curve = DeltaP_curve / (1000 * GRAVITY)  # Head for the curve fluid (converted to meters)

        return head_curve, head_actual
    
    def _flow_from_head(self, head_target):
        """
        Given a desired head at a given speed, solve for the corresponding flow using 
        the inverse of the head curve (using the affinity laws).
        """
        # Define a residual function to minimize: difference between calculated and target head
        def residual(flow):
            head_estimated = self.head_poly(flow)
            return head_estimated - head_target

        # Set flow bounds: between minimum and maximum flow values from the fitted curve
        flow_min, flow_max = min(self.flows_np), max(self.flows_np)

        # Use root-finding to solve for the flow corresponding to the target head
        res = root_scalar(residual, bracket=(flow_min, flow_max), method='bisect')

        if not res.converged:
            raise ValueError("Could not determine flow from head value.")  # If the root-finding did not converge

        return res.root  # Return the calculated flow value
    
    def _similarity_flow(self, flow_curve):
        """
        Calculate flow using the affinity laws.
        This method scales the flow curve to the new speed.
        """
        flow = flow_curve * (self.inputs['N_rot'] / self.params['speed_ref'])  # Apply the affinity law for flow scaling
        return flow

    def _similarity_power(self, flow_curve):
        """
        Calculate power using the affinity laws.
        This method scales the power curve to the new speed, and adjusts for fluid differences.
        """
        fluid = self.inputs['fluid']  # Get the fluid type
        fluid_curve = self.params['curves_fluid']  # Reference fluid for the curves

        # Apply the affinity law for power scaling
        power = self.power_poly(flow_curve) * (self.inputs['N_rot'] / self.params['speed_ref']) ** 3

        # Density correction: scale power based on the actual fluid density
        rho_actual = CP.PropsSI('D', 'T', self.inputs['T_su'], 'P', self.inputs['P_su'], fluid)
        rho_curve = CP.PropsSI('D', 'T', 20 + 273, 'P', 2e5, fluid_curve)  # Assume reference conditions for the curve fluid
        density_ratio = rho_actual / rho_curve  # Ratio of densities for scaling

        # Scale the power for the actual fluid density
        power = power * density_ratio  # Power in kW (after density correction)

        return power

    def solve(self):
        """
        Solves for the power, flow, and head using the similarity laws, based on the provided inputs.
        """
        # Check if the component is calculable and parametrized
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized): # If the component is not calculable and/or not parametrized
            self.solved = False
            print("PumpSimilarityLaws could not be solved. It is not calculable and/or not parametrized")
            self.print_setup()
            return
        self._prepare_interpolators()  # Prepare the interpolators for head and power curves

        # Compute the head values for both the actual fluid and the curve fluid
        head_curve, head_actual = self._similarity_fluid()

        # Calculate the flow corresponding to the computed head
        flow_curve = self._flow_from_head(head_curve)
        self.flow_actual = self._similarity_flow(flow_curve)  # Scale flow to the new speed

        # Calculate the power based on the scaled flow and speed
        self.power_actual = self._similarity_power(flow_curve)
        self.update_connectors()  # Update the connectors with the calculated values


    def update_connectors(self):
        self.su.set_m_dot(self.flow_actual*self.su.D)  # Update the mass flow rate for the suction connector
        self.ex.set_m_dot(self.su.m_dot)  # Update the mass flow rate for the exhaust connector
        self.W_mec.set_W_dot(self.power_actual)  # Update the work connector with the calculated power
        #Calculate the exhaust enthalpy
        self.h_ex = self.su.h + (self.power_actual / self.su.m_dot)  # Calculate the exhaust enthalpy based on the power and mass flow rate
        self.ex.set_h(self.h_ex)  # Set the exhaust enthalpy for the exhaust connector
        self.ex.set_fluid(self.su.fluid)  # Set the fluid type for the exhaust connector

    def plot_characteristic_curves(self, speeds_to_plot=None, flow_range=None, n_points=100):
        """
        Plots the characteristic curves for head and power versus flow for different speeds,
        using the similarity laws, including density correction.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        # Default to plotting the curve at the current operating speed if none is specified
        if speeds_to_plot is None:
            speeds_to_plot = [self.inputs['N_rot']]

        # Flow range to plot
        if flow_range is None:
            flow_min, flow_max = min(self.flows_np), max(self.flows_np)
        else:
            flow_min, flow_max = flow_range

        base_flows = np.linspace(flow_min, flow_max, n_points)

        # Retrieve fluid densities
        rho_actual = CP.PropsSI('D', 'T', self.inputs['T_su'], 'P', self.inputs['P_su'], self.inputs['fluid'])
        rho_curve = CP.PropsSI('D', 'T', 20 + 273, 'P', 2e5, self.params['curves_fluid'])
        density_ratio = rho_actual / rho_curve

        # Start plotting
        plt.figure(figsize=(12, 5))

        # ---- Head vs Flow Plot ----
        plt.subplot(1, 2, 1)
        
        for speed in speeds_to_plot:
            factor = speed / self.params['speed_ref']
            scaled_flows = base_flows * factor
            heads = [self.head_poly(flow) * factor ** 2 for flow in base_flows]
            plt.plot(scaled_flows, heads, label=f"{speed} RPM")
        plt.title("Head vs Flow Rate")
        plt.xlabel("Flow Rate [m³/h]")
        plt.ylabel("Head [m]")
        plt.grid(True)
        plt.legend()

        # ---- Power vs Flow Plot ----
        plt.subplot(1, 2, 2)
        for speed in speeds_to_plot:
            factor = speed / self.params['speed_ref']
            scaled_flows = base_flows * factor
            powers = [
                self.power_poly(flow) * factor ** 3 * density_ratio
                for flow in base_flows
            ]
            plt.plot(scaled_flows, powers, label=f"{speed} RPM")
            
        plt.title("Power vs Flow Rate (density-corrected)")
        plt.xlabel("Flow Rate [m³/h]")
        plt.ylabel("Power [W]")
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()