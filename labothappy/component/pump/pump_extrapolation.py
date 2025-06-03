"External modules"

from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d

"Internal modules"
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector


class PumpExtrapolationModel(BaseComponent):
    """
    Component: Pump with extrapolated performance curves

    Model: Empirical extrapolation-based pump model

    **Description**:

        This model simulates the behavior of a dynamic pump (e.g., centrifugal) based on extrapolated performance curves.
        The curves (flowrate, head, isentropic efficiency, and NPSH required) are scaled with rotational speed to represent
        off-design operation. The model computes the outlet state and power consumption accordingly.

    **Assumptions**:

        - Steady-state operation.
        - Speed-based performance extrapolation is valid.
        - Incompressible working fluid (liquid).
        - One-dimensional flow behavior.

    **Connectors**:

        su (MassConnector): Suction (inlet) side of the pump.

        ex (MassConnector): Exhaust (outlet) side of the pump.

        W_pp (WorkConnector): Mechanical power input to the pump shaft.

    **Parameters**:

        Omega_rated: Rated pump speed [rpm]

        min_flowrate: Minimum valid flowrate [m³/h]

        rated_flowrate: Rated flowrate [m³/h]

        max_flowrate: Maximum valid flowrate [m³/h]

        PI_rated: Rated power input [W]

        D_p: Hydraulic diameter [m]

        V_dot_curve: Volumetric flowrate curve [m³/h]

        Delta_H_curve: Corresponding head rise curve [m]

        eta_is_curve: Isentropic efficiency curve [-]

        NPSH_r_curve: Net Positive Suction Head required curve [m]

        eta_m: Mechanical efficiency of the pump [-]

        eta_max_motor: Maximum electrical efficiency of the motor [%]

        W_dot_el_rated: Rated electric power consumption [W]

    **Inputs**:

        su_p: Inlet pressure [Pa]

        su_T: Inlet temperature [K]

        ex_p: Outlet pressure [Pa]

        su_fluid: Working fluid [-]

        Omega_pp: Actual rotational speed [rpm]

    **Outputs**:

        h_ex: Outlet specific enthalpy [J/kg]

        eta_is: Isentropic efficiency [-]

        NPSH_r: Required NPSH [m]

        W_dot_wf: Work done on fluid [W]

        W_dot_el: Electric power consumption [W]

        V_dot: Volumetric flow rate [m³/h]

    **Notes**:

        The model extrapolates the performance curves to handle off-design speed operation. Flow outside the valid range
        triggers a warning. Electrical efficiency is estimated using empirical polynomial fits as a function of power ratio.
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()  # Inlet side connector
        self.ex = MassConnector()  # Outlet side connector
        self.W_pp = WorkConnector()  # Mechanical shaft power connector

    def get_required_inputs(self):
        """
        Return a list of required input variable names.
        Used to check if the model has enough data to run.
        """
        self.sync_inputs()  # Ensure connector values are reflected in the inputs dict
        return ["su_p", "su_T", "ex_p", "su_fluid", "Omega_pp"]

    def sync_inputs(self):
        """
        Synchronize the internal `inputs` dictionary with values from the connectors.
        """
        if self.su.fluid is not None:
            self.inputs["su_fluid"] = self.su.fluid
        if self.su.T is not None:
            self.inputs["su_T"] = self.su.T
        elif self.su.h is not None:
            self.inputs["su_h"] = self.su.h
        if self.su.p is not None:
            self.inputs["su_p"] = self.su.p
        if self.ex.p is not None:
            self.inputs["ex_p"] = self.ex.p
        if self.W_pp.N is not None:
            self.inputs["Omega_pp"] = self.W_pp.N * 60

    def set_inputs(self, **kwargs):
        """
        Set inputs via keyword arguments and update connector states accordingly.
        """
        self.inputs.update(kwargs)

        # Push inputs into connector objects
        if "su_fluid" in self.inputs:
            self.su.set_fluid(self.inputs["su_fluid"])
        if "su_T" in self.inputs:
            self.su.set_T(self.inputs["su_T"])
        elif "su_h" in self.inputs:
            self.su.set_h(self.inputs["su_h"])
        if "su_p" in self.inputs:
            self.su.set_p(self.inputs["su_p"])
        if "ex_p" in self.inputs:
            self.ex.set_p(self.inputs["ex_p"])
        if "Omega_pp" in self.inputs:
            self.W_pp.set_N(self.inputs["Omega_pp"] / 60)
            self.Omega_pp = self.inputs["Omega_pp"]

    def get_required_parameters(self):
        """
        Return a list of required parameters needed for model execution.
        """
        return [
            "Omega_rated",
            "min_flowrate",
            "rated_flowrate",
            "max_flowrate",
            "PI_rated",
            "D_p",
            "V_dot_curve",
            "Delta_H_curve",
            "eta_is_curve",
            "NPSH_r_curve",
            "eta_m",
            "eta_max_motor",
            "W_dot_el_rated",
        ]

    def print_setup(self):
        """
        Print current configuration of inputs, connectors, and parameters.
        Useful for debugging or inspection.
        """
        print("=== Pump Setup ===")
        print("Connectors:")
        print(
            f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}"
        )
        print(
            f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}"
        )
        print("\nInputs:")
        for input_name in self.get_required_inputs():
            if input_name in self.input_values:
                print(f"  - {input_name}: {self.input_values[input_name]}")
            else:
                print(f"  - {input_name}: Not set")
        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")
        print("======================")

    def eta_el(self, P_ratio):
        """
        Estimate electrical efficiency based on power ratio using piecewise polynomial fits.

        Parameters:
            P_ratio (float): Ratio of actual to rated electric power [%]

        Returns:
            eta_el (float): Estimated electrical efficiency (fraction)
            
        References:
            Kostas
        """
        coefs_20p = [
            78.74503721229180,
            1.54402269709448,
            -0.04662069008665,
            0.00069559243591,
            -0.00000499382422,
            0.00000001349770,
        ]
        coefs_m20 = [
            0.82025554862776,
            4.78234707015054,
            0.73842411551209,
            -0.06398392686793,
            0.00134594665523,
        ]

        eta_ratio = 0
        if P_ratio >= 20:
            for i, c in enumerate(coefs_20p):
                eta_ratio += c * P_ratio**i
        else:
            for i, c in enumerate(coefs_m20):
                eta_ratio += c * P_ratio**i

        return eta_ratio * self.params["eta_max_motor"]

    def solve(self):
        """
        Main solving routine that extrapolates pump curves, computes thermodynamic and performance outputs.
        """
        self.check_calculable()
        self.check_parametrized()

        g = 9.81  # Gravity [m/s²]

        if not self.calculable or not self.parametrized:
            print("Component is not calculable or not parametrized")
            return

        if self.su.p > self.ex.p or self.Omega_pp < 0:
            print(
                "Supply pressure is higher than exhaust pressure or rotation speed is negative"
            )
            return

        self.V_dot_flag = 0  # Flag if flow is outside curve range

        # Scale curves with respect to current speed
        self.V_dot_curve = self.params["V_dot_curve"] * (
            self.Omega_pp / self.params["Omega_rated"]
        )
        self.DH_curve = (
            self.params["Delta_H_curve"]
            * (self.Omega_pp / self.params["Omega_rated"]) ** 2
        )
        self.eta_is_curve = self.params["eta_is_curve"]
        self.NPSH_r_curve = (
            self.params["NPSH_r_curve"]
            * (self.Omega_pp / self.params["Omega_rated"]) ** 2
        )

        # Create interpolation functions for extrapolated performance
        DH_V = interp1d(
            self.DH_curve, self.V_dot_curve, kind="linear", fill_value="extrapolate"
        )
        V_eta = interp1d(
            self.V_dot_curve, self.eta_is_curve, kind="linear", fill_value="extrapolate"
        )
        V_NPSH_r = interp1d(
            self.V_dot_curve, self.NPSH_r_curve, kind="linear", fill_value="extrapolate"
        )

        # Compute head difference [m]
        DP = self.ex.p - self.su.p
        self.DH = DP / (g * self.su.D)

        # Compute volumetric flowrate [m³/h]
        self.V_dot = DH_V(self.DH)

        # Flag if extrapolation is outside defined curves
        if self.V_dot > self.V_dot_curve[-1] or self.V_dot < self.V_dot_curve[0]:
            self.V_dot_flag = 1

        # Mass flow rate [kg/s]
        self.m_dot = self.su.D * self.V_dot / 3600
        self.su.set_m_dot(self.m_dot)
        self.ex.set_m_dot(self.m_dot)

        # Isentropic efficiency
        self.eta_is = V_eta(self.V_dot)

        # NPSH required
        self.NPSH_r = V_NPSH_r(self.V_dot)

        # Isentropic outlet enthalpy
        h_ex_s = PropsSI("H", "P", self.ex.p, "S", self.su.s, self.su.fluid)

        # Actual outlet enthalpy
        h_ex = self.su.h + (h_ex_s - self.su.h) / self.eta_is
        self.ex.set_h(h_ex)

        # Hydraulic power [W]
        self.W_dot_hyd = self.su.m_dot * (h_ex_s - self.su.h)

        # Actual shaft work on fluid [W]
        Delta_h = self.ex.h - self.su.h
        self.W_dot_wf = self.su.m_dot * Delta_h

        # Electrical efficiency and electric power consumption
        P_ratio = (
            100
            * (self.W_dot_wf / (self.params["eta_m"] * self.eta_is))
            / self.params["W_dot_el_rated"]
        )
        self.eta_el = self.eta_el(P_ratio) / 100

        self.W_dot_el = self.W_dot_wf / (self.params["eta_m"] * self.eta_el)

        # Set output to connector
        self.W_pp.set_W_dot(self.W_dot_wf)

        if self.V_dot_flag:
            print("Flowrate outside possible range. Impossible operation.")

        self.defined = True
