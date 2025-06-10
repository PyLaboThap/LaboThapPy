"External modules"

import numpy as np
from CoolProp.CoolProp import PropsSI

"Internal modules"
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector


class Turb_polyn_eff(BaseComponent):
    """
    Component: Turbine with choked flow and polynomial efficiency maps

    Model: Steady-state turbine with polynomial-based isentropic and electrical efficiencies

    **Description**:

        This model simulates a small-scale radial or axial flow turbine operating under steady-state,
        choked-flow conditions. It estimates isentropic efficiency from a 3D polynomial map based on
        condenser pressure and mass flow rate, and electrical efficiency from a polynomial fit versus
        power output ratio. The outlet enthalpy, shaft power, and electric output are computed from
        energy balances and property calculations using CoolProp.

    **Assumptions**:

        - Steady-state, one-dimensional flow.
        - Flow is always choked at turbine throat (critical condition).
        - CoolProp is used for accurate fluid property evaluation.
        - Electrical efficiency is a function of output-to-rated power ratio.
        - No heat losses to surroundings.

    **Connectors**:

        su (MassConnector): Supply (inlet) side of the turbine.

        ex (MassConnector): Exhaust (outlet) side of the turbine.

        W_turb (WorkConnector): Shaft power output from the turbine.

    **Parameters**:

        D_inlet: Turbine throat diameter [m]

        N_turb_rated: Rated rotational speed [rpm]

        turb_voltage: Generator voltage rating [V]

        turb_phases: Number of electrical phases [-]

        eta_max_motor: Maximum generator electrical efficiency [%]

        W_dot_el_rated: Rated electric output [W]

        eta_m: Mechanical efficiency of the turbine shaft [-]

        eta_is_coefs: Coefficients for 3D isentropic efficiency polynomial [-]

        eta_is_coefs_red: Reduced polynomial coefficients for validation or simplified use [-]

        A_th: Throat flow area [m²]

    **Inputs**:

        su_p: Inlet pressure [Pa]

        su_T or su_h: Inlet temperature [K] or enthalpy [J/kg]

        ex_p: Outlet pressure [Pa]

        su_fluid: Working fluid [-]

        N_rot: Actual shaft rotational speed [rpm]

    **Outputs**:

        h_ex: Outlet enthalpy [J/kg]

        eta_is: Isentropic efficiency [-]

        W_dot: Shaft work output [W]

        eta_el: Electrical efficiency [-]

        W_dot_el: Electrical output [W]

        m_dot: Mass flow rate through turbine [kg/s]

    **Notes**:

        - Flow rate is computed from sonic velocity, area, and inlet conditions (choked assumption).
        - Isentropic enthalpy drop is calculated using inlet entropy and outlet pressure.
        - Electrical efficiency uses separate polynomials above and below 20% load.
        - No dynamic behavior is included; suitable for steady-state energy system simulations.
    """

    def __init__(self):
        super().__init__()
        # Define connectors: suction (su), exhaust (ex), and work (W_turb)
        self.su = MassConnector()
        self.ex = MassConnector()
        self.W_turb = WorkConnector()

    def get_required_inputs(self):
        """
        Determines and sets the required input values from connected components.
        Returns the list of all required input labels.
        """

        if self.inputs == {}:
            # Populate inputs dictionary with known inlet conditions
            if self.su.T is not None:
                self.inputs["su_T"] = self.su.T
            elif self.su.h is not None:
                self.inputs["su_h"] = self.su.h
            if self.su.p is not None:
                self.inputs["su_p"] = self.su.p
            if self.ex.p is not None:
                self.inputs["ex_p"] = self.ex.p
            if self.W_turb.N is not None:
                self.inputs["N_rot"] = self.W_turb.N
            if self.su.fluid is not None:
                self.inputs["su_fluid"] = self.su.fluid

        if self.inputs != {}:
            # Set connector values from the inputs dictionary
            self.su.set_fluid(self.inputs["su_fluid"])
            if "su_T" in self.inputs:
                self.su.set_T(self.inputs["su_T"])
            elif "su_h" in self.inputs:
                self.su.set_h(self.inputs["su_h"])
            if "su_p" in self.inputs:
                self.su.set_p(self.inputs["su_p"])
            if "ex_p" in self.inputs:
                self.ex.set_p(self.inputs["ex_p"])
            if "N_rot" in self.inputs:
                self.W_turb.set_N(self.inputs["N_rot"])
            if "su_fluid" in self.inputs:
                self.su.set_fluid(self.inputs["su_fluid"])

        return ["su_p", "su_T", "ex_p", "N_rot", "su_fluid"]

    def get_required_parameters(self):
        """
        Returns a list of parameter names that must be defined to perform calculations.
        """
        return [
            "D_inlet",  # Inlet hydraulic diameter
            "N_turb_rated",  # Rated turbine speed
            "turb_voltage",  # Turbine electrical voltage
            "turb_phases",  # Number of turbine phases
            "eta_max_motor",  # Max motor (generator) efficiency
            "W_dot_el_rated",  # Rated electrical output
            "eta_m",  # Mechanical efficiency
            "eta_is_coefs",  # Polynomial coefficients for isentropic efficiency
            "eta_is_coefs_red",  # (Unused here but listed)
            "A_th",  # Turbine throat area
        ]

    def print_setup(self):
        """
        Prints the current configuration of the turbine component, including connectors, inputs, and parameters.
        """
        print("=== Turbine Setup ===")
        print("Connectors:")
        print(
            f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}"
        )
        print(
            f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}"
        )
        print(f"  - W_dot: speed={self.W_turb.N}")

        print("\nInputs:")
        for input in self.get_required_inputs():
            if input in self.inputs:
                print(f"  - {input}: {self.inputs[input]}")
            else:
                print(f"  - {input}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")

        print("======================")

    def eta_el_fun(self, P_ratio):
        """
        Computes the electrical efficiency as a function of normalized power output (P_ratio).
        Uses two separate polynomial fits depending on the magnitude of P_ratio.
        """

        # Polynomial coefficients for high and low power ranges
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

        # Use appropriate polynomial based on power ratio
        if P_ratio >= 20:
            for i in range(len(coefs_20p)):
                eta_ratio += coefs_20p[i] * P_ratio**i
        else:
            for i in range(len(coefs_m20)):
                eta_ratio += coefs_m20[i] * P_ratio**i

        return eta_ratio * self.params["eta_max_motor"]

    def eta_is_turb(self):
        """
        Computes isentropic efficiency using a 3rd-degree polynomial function of
        exhaust pressure and mass flow rate.
        """

        a = self.params["eta_is_coefs"]
        p_cond = self.ex.p * 1e-5  # Convert Pa to bar
        m_dot = self.su.m_dot

        eta = (
            a[0]
            + a[1] * p_cond
            + a[2] * m_dot
            + a[3] * p_cond**2
            + a[4] * p_cond * m_dot
            + a[5] * m_dot**2
            + a[6] * p_cond**3
            + a[7] * p_cond**2 * m_dot
            + a[8] * p_cond * m_dot**2
            + a[9] * m_dot**3
        )

        return eta

    def solve(self):
        """
        Solves the turbine model: calculates outlet state, flow rate, mechanical and electrical work output.
        """

        self.check_calculable()
        self.check_parametrized()

        if self.calculable and self.parametrized:

            # Calculate speed of sound at inlet using ideal gas approximation
            R = 8.3144
            R_M = R / PropsSI("M", self.su.fluid)
            cp_in, cv_in = PropsSI(
                ("C", "CVMASS"), "T", self.su.T, "P", self.su.p, self.su.fluid
            )
            gamma = cp_in / cv_in
            self.a = np.sqrt(gamma * R_M * self.su.T)

            # Mass flow rate using choked flow assumption
            self.m_dot = self.a * self.su.D * self.params["A_th"]
            self.su.set_m_dot(self.m_dot)

            # Calculate isentropic efficiency from the empirical map
            self.eta_is = self.eta_is_turb()

            # Compute isentropic enthalpy at exhaust pressure and entropy
            h_ex_is = PropsSI("H", "P", self.ex.p, "S", self.su.s, self.su.fluid)

            # Actual enthalpy based on isentropic efficiency
            h_ex = self.su.h - (self.eta_is / 100) * (self.su.h - h_ex_is)

            # Set exhaust state
            self.ex.set_h(h_ex)
            self.ex.set_m_dot(self.su.m_dot)

            # Compute shaft work
            self.W_dot = (self.su.h - h_ex) * self.su.m_dot

            # Electrical efficiency and electrical output
            self.eta_el = (
                self.eta_el_fun(100 * self.W_dot / self.params["W_dot_el_rated"]) / 100
            )
            self.W_el = self.W_dot * self.eta_el

            # Set mechanical output in the work connector
            self.W_turb.set_W_dot(self.W_dot)

            self.defined = True

        else:
            if not self.calculable:
                print("Input of the component not completely known")
            if not self.parametrized:
                print("Parameters of the component not completely known")
