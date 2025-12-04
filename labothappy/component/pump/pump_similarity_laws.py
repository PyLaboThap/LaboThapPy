# Internal imports
import CoolProp.CoolProp as CP

# External imports
import numpy as np
from scipy.optimize import root_scalar

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

GRAVITY = 9.81  # m/s², Gravitational acceleration constant


class PumpSimilarityLaws(BaseComponent):
    """
    **Component**: Pump with a characteristic curve associated
    
    **Model**: Model using characteristic curves for head and power and similarity laws
    to calculate the pump performance at different speeds.

    **Description**:

        Simulates a pump using head and power performance curves (vs flow and speed).
        Uses the similarity laws to calculate the pump performance at different speeds.

    **Assumptions**:

        - Steady-state
        - Incompressible fluid approximation for head calculation

    **Connectors**:
        su (MassConnector): Supply (inlet) side of the turbine.

        ex (MassConnector): Exhaust (outlet) side of the turbine.

        W (WorkConnector): Shaft power output from the turbine.

    **Parameters**:

        curves_head: Dict = List[(flow, head)] dicionary of head curves ([m^3/h], [m])

        curves_power: Dict = List[(flow, power)] dictionary of power curves ([m^3/h], [W])

        curves_fluid: Fluid used for the characteristic curves 

        curves_rho : Density of the fluid used for the characteristic curves (kg/m^3)

        speed_ref: Reference speed for the curves (RPM)

        mode: Operation mode of the pump ('P_N', 'M_N', 'P_M')

    **Inputs**:

        *mode = P_N:*

        P_su: Suction pressure [Pa]

        T_su: Suction temperature [K]

        P_ex: Exhaust pressure [Pa]

        N_rot: Rotational speed [RPM]

        fluid: Actual working fluid

        *mode = P_M:*

        P_su: Suction pressure [Pa]

        T_su: Suction temperature [K]

        P_ex: Exhaust pressure [Pa]

        m_dot: Mass flow rate [kg/s]

        fluid: Actual working fluid

        *mode = M_N:*

        P_su: Suction pressure [Pa]

        T_su: Suction temperature [K]

        N_rot: Rotational speed [RPM]

        m_dot: Mass flow rate [kg/s]

        fluid: Actual working fluid

    **Outputs**:

        V_dot: Volumetric flow rate [m³/h]

        W_dot: Shaft power required by the pump [W]

        m_dot: Mass flow rate [kg/s] or P_ex: Exhaust pressure [Pa] or N_rot: Rotational speed [RPM], depending on the mode


    **Notes**:
        - Three modes of operation are available:
            - 'P_N': Given suction and exhaust pressures and rotational speed, calculates flow and power
            - 'M_N': Given suction pressure, temperature, rotational speed, and mass flow rate, calculates exhaust pressure and power
            - 'P_M': Given suction and exhaust pressures and mass flow rate, calculates rotational speed and power
    """

    def __init__(self):
        # Initialize the base component class 
        super().__init__()
        self.su = MassConnector()  # Suction connector
        self.ex = MassConnector()  # Exhaust connector
        self.W = WorkConnector()  # Mechanical work connector

    def get_required_inputs(self):
        # Returns a list of required input variable names. Used to check if the model has enough data to run.

        if self.params["mode"] == "P_N": # If the mode is 'P_N', the rotational speed and exhaust pressure are required while the mass flow rate of the expander is calculated
            return ["P_su", "T_su", "P_ex", "N_rot", "fluid"]

        elif self.params["mode"] == "P_M": # If the mode is 'P_M', the mass flow rate and exhaust pressure are required while the rotational speed of the expander is calculated
            return ["P_su", "T_su", "P_ex", "m_dot", "fluid"]

        elif self.params["mode"] == "M_N": # If the mode is 'M_N', the rotational speed and mass flow rate are required while the exhaust pressure of the expander is calculated
            return ["P_su", "T_su", "N_rot", "m_dot", "fluid"]

        else:
            raise ValueError("Specify a mode for the pump : 'P_N', 'P_M', 'M_N'")

    def get_required_parameters(self):
        # Returns a list of required parameters needed for model execution.

        return ["curves_head", "curves_power", "speed_ref", "curves_rho", "curves_fluid", "mode"]

    def _prepare_interpolators(self):
        """
        Prepares polynomial interpolation functions for head and power curves
        based on the point from the characteristic curve provided for a single-speed pump.
        This method returns the fitted polynomial functions for head and power curves at teh reference speed.
        """
        # Extract the head and power curves from the input parameters
        head_curve = self.params["curves_head"]
        power_curve = self.params["curves_power"]

        # Extract flow and corresponding head values for the head curve
        flows_head, heads = zip(*sorted(head_curve))

        # Extract flow and corresponding power values for the power curve
        flows_power, powers = zip(*sorted(power_curve))

        self.flows_np = np.array(flows_head)  # Store flow values

        # Fit polynomial curves to both the head and power data (degree 3 or lower)
        deg = min(
            3, len(flows_head) - 1
        )  # Avoid overfitting if there are too few data points
        self.head_poly = np.poly1d(
            np.polyfit(flows_head, heads, deg)
        )  # Polynomial curve for head at one speed of reference
        self.power_poly = np.poly1d(
            np.polyfit(flows_power, powers, deg)
        )  # Polynomial curve for power at one speed of reference

    def _similarity_fluid(self):
        """
        Calculates head for a different fluid by adjusting based on density differences
        between the curve fluid and the actual operating fluid.
        """

        # Calculate the density ratio between the actual fluid and the curve fluid

        # rho_actual = CP.PropsSI("D", "T", self.inputs["T_su"], "P", self.inputs["P_su"], self.su.fluid)
        rho_actual = self.su.D  # Density of the actual fluid at suction conditions
        rho_curve = self.params["curves_rho"]
        density_ratio = rho_curve / rho_actual # The ratio of the reference fluid's density to the actual fluid
        # Calculate pressure difference (Delta P) for both fluids
        DeltaP_actual = self.inputs["P_ex"] - self.inputs["P_su"] # Pressure difference across the pump
        DeltaP_curve = DeltaP_actual * density_ratio # Pressure difference for the curve fluid (adjusted by density)

        # Calculate head values for both the actual and curve fluids
        head_actual = DeltaP_actual / (rho_actual * GRAVITY) # Head for the actual fluid (converted to meters)
        head_curve = DeltaP_curve / (rho_curve * GRAVITY) # Head for the curve fluid (converted to meters)

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
        res = root_scalar(
            residual, x0=flow_min, method="secant"
        )  # bracket=(flow_min, flow_max), method="bisect") #to be checked

        if not res.converged:
            raise ValueError(
                "Could not determine flow from head value."
            )  # If the root-finding did not converge

        return res.root  # Return the calculated flow value

    def _flow_in_curve_from_m_dot(self):
        """
        Calculate flow using the affinity laws.
        This method scales the flow curve to the new speed.
        """
        V_dot = (self.su.m_dot * 3600) / self.su.D  # m^3/h

        flow = V_dot * (
            self.inputs["N_rot"] / self.params["speed_ref"]
        )  # Apply the affinity law for flow scaling
        return flow
    
    def _similarity_flow(self, flow_curve):
        """
        Calculate flow using the affinity laws.
        This method scales the flow curve to the new speed.
        """
        flow = flow_curve * (self.inputs["N_rot"] / self.params["speed_ref"])  # Apply the affinity law for flow scaling
        return flow
    
    def _similarity_speed(self, flow_curve):
        """
        Calculate speed using the affinity laws.
        This method scales the speed based on the flow curve.
        """

        # Apply the affinity law for speed scaling
        flow_actual = self.inputs["m_dot"] * 3600 / self.su.D  # m^3/h
        speed = self.params["speed_ref"] * (flow_actual / flow_curve)

        return speed

    def _similarity_power(self, flow_curve, speed_actual):
        """
        Calculate power using the affinity laws.
        This method scales the power curve to the new speed, and adjusts for fluid differences.
        """

        # Apply the affinity law for power scaling
        power = (
            self.power_poly(flow_curve)
            * (speed_actual/ self.params["speed_ref"]) ** 3
        )

        # Density correction: scale power based on the actual fluid density
        self.AS_actual.update(CP.HmassP_INPUTS, self.su.h, self.inputs["P_su"])
        rho_actual = self.AS_actual.rhomass()
        self.AS_curve.update(CP.PT_INPUTS, 2e5, 20 + 273.15)
        rho_curve = self.AS_curve.rhomass()
        density_ratio = rho_actual / rho_curve  # Ratio of densities for scaling

        # Scale the power for the actual fluid density
        power = power * density_ratio  # Power in W (after density correction)

        return power

    def solve(self):
        """
        Solves for the power, flow, and head using the similarity laws, based on the provided inputs.
        """
        # Check if the component is calculable and parametrized
        self.check_calculable()
        self.check_parametrized()

        self.AS_curve = CP.AbstractState("HEOS", self.params["curves_fluid"])
        self.AS_actual = CP.AbstractState("HEOS", self.su.fluid)

        if (not self.calculable):  # If the component is not calculable and/or not parametrized
            self.solved = False
            print("PumpSimilarityLaws could not be solved. It is not calculable.")
            self.print_setup()
            return

        elif (not self.parametrized):  # If the component is not calculable and/or not parametrized
            self.solved = False
            print("PumpSimilarityLaws could not be solved. It is not parametrized.")
            self.print_setup()
            return

        self._prepare_interpolators()  # Prepare the interpolators for head and power curves

        if self.params["mode"] == "P_N":
            # Compute the head values for both the actual fluid and the curve fluid
            head_curve, head_actual = self._similarity_fluid()

            # Calculate the flow corresponding to the computed head
            flow_curve = self._flow_from_head(head_curve) # Flow in the curve fluid at reference speed corresponding to the head
            self.flow_actual = self._similarity_flow(flow_curve)  # Scale flow to the new speed

            # Calculate the power based on the scaled flow and speed
            self.power_actual = self._similarity_power(flow_curve, self.inputs['N_rot'])  # Power in watts

        elif self.params["mode"] == "P_M":
            # Compute the head values for both the actual fluid and the curve fluid
            head_curve, head_actual = self._similarity_fluid()

            # Calculate the flow corresponding to the computed head
            flow_curve = self._flow_from_head(head_curve) # Flow in the curve fluid at reference speed corresponding to the head
            self.speed_actual = self._similarity_speed(flow_curve)  # Scale speed to the new flow

            # Calculate the power based on the scaled flow and speed
            self.power_actual = self._similarity_power(flow_curve, self.speed_actual)  # Power in watts

        elif self.params["mode"] == "M_N":
            rho_curve = self.params["curves_rho"]
            flow_curve = self._flow_in_curve_from_m_dot()
            head_curve = self.head_poly(flow_curve)

            DP_curve = (
                head_curve
                * GRAVITY
                * rho_curve
                * (self.inputs["N_rot"] / self.params["speed_ref"]) ** 2
            )
            density_ratio = (rho_curve / self.su.D)  # The ratio of the reference fluid's density to the actual fluid

            self.DP_actual = DP_curve / density_ratio
            self.power_actual = self._similarity_power(flow_curve, self.inputs['N_rot'])  # Power in watts

        self.update_connectors()  # Update the connectors with the calculated values

    def update_connectors(self):
        if self.params["mode"] == "P_N":
            self.su.set_m_dot(self.flow_actual * self.su.D / 3600)  # Update the mass flow rate for the suction connector
            self.su.set_V_dot(self.flow_actual)  # Update the volumetric flow rate for the suction connector
            self.ex.set_m_dot(self.su.m_dot)  # Update the mass flow rate for the exhaust connector
            self.W.set_N_rot(self.inputs["N_rot"])  # Update the rotational speed for the work connector
            self.ex.set_p(self.inputs["P_ex"])  # Update the exhaust pressure for the exhaust connector

        elif self.params["mode"] == "P_M":
            self.ex.set_m_dot(self.su.m_dot)  # Update the mass flow rate for the exhaust connector
            self.su.set_V_dot(self.su.m_dot * 3600 / self.su.D)  # Update the mass flow rate for the suction connector
            self.ex.set_m_dot(self.su.m_dot)  # Update the mass flow rate for the exhaust connector
            self.W.set_N_rot(self.speed_actual)  # Update the rotational speed for the work connector
            self.ex.set_p(self.inputs["P_ex"])  # Update the exhaust pressure for the exhaust connector

        elif self.params["mode"] == "M_N":
            self.ex.set_m_dot(self.su.m_dot)  # Update the mass flow rate for the exhaust connector
            self.su.set_V_dot(self.su.m_dot * 3600 / self.su.D)  # Update the mass flow rate for the suction connector
            self.ex.set_V_dot(self.su.m_dot * 3600 / self.su.D)  # Update the mass flow rate for the exhaust connector
            self.ex.set_p(self.su.p + self.DP_actual)  # Update the work connector with the calculated power

        self.W.set_W_dot_el(self.power_actual)  # Update the work connector with the calculated power

        # Calculate the exhaust enthalpy
        self.h_ex = self.su.h + (self.power_actual / self.su.m_dot)  # Calculate the exhaust enthalpy based on the power and mass flow rate
        self.ex.set_h(self.h_ex)  # Set the exhaust enthalpy for the exhaust connector
        self.ex.set_fluid(self.su.fluid)  # Set the fluid type for the exhaust connector

        self.solved = True

    def print_results(self):
        """
        Prints the results of the pump simulation, including flow rate, power, and other relevant parameters.
        """
        if not self.solved:
            print("PumpSimilarityLaws has not been solved yet.")
            return

        print("---- Pump Similarity Laws Results ----")
        print(f"Mode: {self.params['mode']}")
        print(f"Suction Pressure (P_su): {self.su.p:.2f} Pa")
        print(f"Suction Temperature (T_su): {self.su.T:.2f} K")
        print(f"Exhaust Pressure (P_ex): {self.ex.p:.2f} Pa")
        print(f"Mass Flow Rate (m_dot): {self.su.m_dot:.4f} kg/s")
        print(f"Volumetric Flow Rate (V_dot): {self.su.V_dot:.4f} m³/h")
        print(f"Rotational Speed (N_rot): {self.W.N_rot} RPM")
        print(f"Electrical consumption (W_dot_el): {self.W.W_dot_el:.2f} W")
        print("--------------------------------------")

    def plot_characteristic_curves(
        self, speeds_to_plot=None, flow_range=None, n_points=100
    ):
        """
        Plots the characteristic curves for head and power versus flow for different speeds,
        using the similarity laws, including density correction.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        # Default to plotting the curve at the current operating speed if none is specified
        if speeds_to_plot is None:
            speeds_to_plot = [self.inputs["N_rot"]]

        # Flow range to plot
        if flow_range is None:
            flow_min, flow_max = min(self.flows_np), max(self.flows_np)
        else:
            flow_min, flow_max = flow_range

        base_flows = np.linspace(flow_min, flow_max, n_points)

        # Retrieve fluid densities
        rho_actual = self.su.D
        rho_curve = self.params["curves_rho"]
        density_ratio = rho_actual / rho_curve

        # Start plotting
        plt.figure(figsize=(12, 5))

        # ---- Head vs Flow Plot ----
        plt.subplot(1, 2, 1)

        for speed in speeds_to_plot:
            factor = speed / self.params["speed_ref"]
            scaled_flows = base_flows * factor
            heads = [self.head_poly(flow) * factor**2 for flow in base_flows]
            plt.plot(scaled_flows, heads, label=f"{speed} RPM")

        plt.title("Head vs Flow Rate")
        plt.xlabel("Flow Rate [m³/h]")
        plt.ylabel("Head [m]")
        plt.grid(True)
        plt.legend()

        # ---- Power vs Flow Plot ----
        plt.subplot(1, 2, 2)
        for speed in speeds_to_plot:
            factor = speed / self.params["speed_ref"]
            scaled_flows = base_flows * factor
            powers = [
                self.power_poly(flow) * factor**3 * density_ratio for flow in base_flows
            ]
            plt.plot(scaled_flows, powers, label=f"{speed} RPM")

        plt.title("Power vs Flow Rate (density-corrected)")
        plt.xlabel("Flow Rate [m³/h]")
        plt.ylabel("Power [W]")
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()
