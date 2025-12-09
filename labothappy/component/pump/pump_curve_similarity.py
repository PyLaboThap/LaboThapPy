# External imports
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import CoolProp.CoolProp as CP

# Internal imports
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

class PumpCurveSimilarity(BaseComponent):
    """
    **Component**: Pump with a characteristic curve associated
    
    **Model**: Model using characteristic curves for head and power and similarity laws
    to calculate the pump performance at different speeds.

    **Description**:

        Simulates a pump using head, efficiency and NPSHr performance curves (vs flow and speed).
        Uses the similarity laws to calculate the pump performance at different speeds.

    **Assumptions**:

        - Steady-state
        - Incompressible fluid approximation for head calculation

    **Connectors**:

        su (MassConnector): Supply (inlet) side of the turbine.

        ex (MassConnector): Exhaust (outlet) side of the turbine.

        W (WorkConnector): Shaft power output from the turbine.

    **Parameters**:

        V_dot_curve: Array of volumetric flowrate curve [m³/h]

        Delta_H_curve: Corresponding array from head rise curve [m]

        eta_is_curve: Corresponding array from isentropic efficiency curve [-]

        NPSH_r_curve: Corresponding array from NPSH required curve [m]

        N_rot_rated: Rated pump speed [rpm]

        mode : Operation mode ("P_N", "P_M", or "M_N")
            - 'P_N': Given suction and exhaust pressures and rotational speed, calculates flow and power
            - 'M_N': Given suction pressure, temperature, rotational speed, and mass flow rate, calculates exhaust pressure and power
            - 'P_M': Given suction and exhaust pressures and mass flow rate, calculates rotational speed

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

        return ["V_dot_curve", "Delta_H_curve", "eta_is_curve", "NPSH_r_curve", "N_rot_rated", "mode"]

    def solve(self):
        """
        Main solving routine that extrapolates pump curves, computes thermodynamic and performance outputs.
        """
        self.check_calculable()
        self.check_parametrized()

        if not self.calculable or not self.parametrized:
            print("Component is not calculable or not parametrized")
            return

        if self.params["mode"] == "P_N": # Given P_su, T_su, P_ex, N_rot
            if self.su.p > self.ex.p:
                print("Supply pressure is higher than exhaust pressure")
                return
            self.solve_PN()

        elif self.params["mode"] == "P_M": # Given P_su, T_su, P_ex, m_dot
            if self.su.p > self.ex.p:
                print("Supply pressure is higher than exhaust pressure")
                return
            self.solve_PM()

        elif self.params["mode"] == "M_N": # Given P_su, T_su, N_rot, m_dot
            self.solve_MN()

    def solve_PN(self):
        self.AS = CP.AbstractState("HEOS", self.su.fluid)

        g = 9.81  # Gravity [m/s²]

        self.V_dot_flag = 0  # Flag if flow is outside curve range

        # Scale curves with respect to current speed
        self.V_dot_curve = self.params["V_dot_curve"] * (self.W.N_rot / self.params["N_rot_rated"])

        self.Delta_H_curve = (self.params["Delta_H_curve"] * (self.W.N_rot / self.params["N_rot_rated"]) ** 2)

        self.eta_is_curve = self.params["eta_is_curve"]

        self.NPSH_r_curve = (self.params["NPSH_r_curve"] * (self.W.N_rot / self.params["N_rot_rated"]) ** 2)

        # Create interpolation functions for extrapolated performance
        Delta_H_V = interp1d(self.Delta_H_curve, self.V_dot_curve, kind="linear", fill_value="extrapolate")

        V_eta_is = interp1d(self.V_dot_curve, self.eta_is_curve, kind="linear", fill_value="extrapolate")

        V_NPSH_r = interp1d(self.V_dot_curve, self.NPSH_r_curve, kind="linear", fill_value="extrapolate")

        # Compute head difference [m]
        DP = self.ex.p - self.su.p
        self.Delta_H = DP / (g * self.su.D)

        # Compute volumetric flowrate [m³/h]
        self.V_dot = Delta_H_V(self.Delta_H)

        # Flag if extrapolation is outside defined curves
        if self.V_dot > self.V_dot_curve[-1] or self.V_dot < self.V_dot_curve[0]:
            self.V_dot_flag = 1

        # Mass flow rate [kg/s]
        self.m_dot = self.su.D * self.V_dot / 3600
        self.su.set_m_dot(self.m_dot)
        self.ex.set_m_dot(self.m_dot)

        # Isentropic efficiency
        self.eta_is = V_eta_is(self.V_dot)

        # NPSH required
        self.NPSH_r = V_NPSH_r(self.V_dot)

        # Isentropic outlet enthalpy
        self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        h_ex_s = self.AS.hmass()

        # Actual outlet enthalpy
        h_ex = self.su.h + (h_ex_s - self.su.h) / self.eta_is
        self.ex.set_h(h_ex)

        # Ideal hydraulic power
        self.W_dot_hyd = g * self.Delta_H * self.m_dot # W

        # Shaft power
        self.W_dot_shaft = self.su.m_dot * (self.ex.h - self.su.h)  # W
        self.W.set_W_dot(self.W_dot_shaft)

        # Check for impossible operation
        if self.V_dot_flag:
            print("Flowrate outside possible range. Impossible operation.")
            return

    def solve_PM(self):
        self.AS = CP.AbstractState("HEOS", self.su.fluid)
        g = 9.81 # Gravity [m/s²]

        # Known inputs
        self.m_dot = self.su.m_dot

        # Convert mass flow → volumetric flow in m³/h
        self.V_dot = self.m_dot / self.su.D * 3600

        # Required head from pressure rise
        DP = self.ex.p - self.su.p
        H_required = DP / (g * self.su.D)

        # Rated curves
        V_dot_rated = self.params["V_dot_curve"]
        H_rated = self.params["Delta_H_curve"]
        N_rated = self.params["N_rot_rated"]
        NPSH_r_rated = self.params["NPSH_r_curve"]

        # Interpolator for the head curve at rated speed
        H_V_interp = interp1d(V_dot_rated, H_rated,
                            kind="linear", fill_value="extrapolate")

        # Define the equation for N that must be zero:
        # H_curve(Q,N) - H_required = 0
        def head_error(N):
            if N <= 0:
                return 1e9   # avoid non-physical speeds

            # Scale flow relative to rated speed using similarity laws
            V_dot_equiv = self.V_dot * (N_rated / N)

            # Get rated head at this equivalent flow
            H_at_rated = H_V_interp(V_dot_equiv)

            # Apply scaling law
            H_at_speed = H_at_rated * (N / N_rated)**2

            return H_at_speed - H_required

        # Solve for N
        sol = root_scalar(head_error, bracket=[1, 5*N_rated], method="bisect")
        N_rot = sol.root
        self.W.set_N_rot(N_rot)

        # Now compute pump performance at this speed

        # Rebuild scaled curves
        self.V_dot_curve = V_dot_rated * (N_rot / N_rated)
        self.Delta_H_curve = H_rated * (N_rot / N_rated)**2
        self.NPSH_r_curve = NPSH_r_rated * (N_rot / N_rated)**2
        self.eta_is_curve = self.params["eta_is_curve"]

        # Interpolators
        eta_interp = interp1d(self.V_dot_curve, self.eta_is_curve,
                            kind="linear", fill_value="extrapolate")
        NPSH_r_interp = interp1d(self.V_dot_curve, self.NPSH_r_curve,
                            kind="linear", fill_value="extrapolate")
        
        # NPSH
        self.NPSH_r = NPSH_r_interp(self.V_dot)

        # Efficiency
        self.eta_is = eta_interp(self.V_dot)

        # Update mass flow
        self.V_dot = self.V_dot

        # Thermodynamics
        self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        h_ex_s = self.AS.hmass()
        h_ex = self.su.h + (h_ex_s - self.su.h) / self.eta_is

        self.ex.set_h(h_ex)

        # Power
        # Ideal hydraulic power
        self.W_dot_hyd = g * H_required * self.m_dot  # W

        # Actual shaft power
        self.W_dot_shaft = self.m_dot * (h_ex - self.su.h)
        self.W.set_W_dot(self.W_dot_shaft)

    def solve_MN(self):
        self.AS = CP.AbstractState("HEOS", self.su.fluid)
        g = 9.81 # Gravity [m/s²]
        self.V_dot_flag = 0  # Flag if flow is outside curve range

        # Known inputs
        self.m_dot = self.su.m_dot
        rho = self.su.D
        N = self.W.N_rot
        N_rated = self.params["N_rot_rated"]

        # Convert mass flow to volumetric flow (m³/h)
        V_dot = self.m_dot / rho * 3600
        self.V_dot = V_dot

        # Rated curves
        V_dot_rated = self.params["V_dot_curve"]
        H_rated = self.params["Delta_H_curve"]
        eta_rated = self.params["eta_is_curve"]
        NPSH_r_rated = self.params["NPSH_r_curve"]

        # Interpolators for rated curves
        H_interp = interp1d(V_dot_rated, H_rated, kind="linear", fill_value="extrapolate")
        eta_interp = interp1d(V_dot_rated, eta_rated, kind="linear", fill_value="extrapolate")
        NPSH_r_interp = interp1d(V_dot_rated, NPSH_r_rated, kind="linear", fill_value="extrapolate")

        # 1. Map actual Q to equivalent Q at rated speed
        V_dot_eq = V_dot * (N_rated / N)

        # 2. Rated head at Q_eq
        H_at_rated = H_interp(V_dot_eq)

        # 3. Head at actual speed
        H = H_at_rated * (N / N_rated)**2
        self.Delta_H = H

        # 4. Compute exhaust pressure
        P_ex = self.su.p + rho * g * H
        self.ex.set_p(P_ex)

        # Flag out-of-range flows
        if V_dot < V_dot_rated[0] or V_dot > V_dot_rated[-1]:
            self.V_dot_flag = 1

        # Efficiency at this speed (same scaling method)
        eta = eta_interp(V_dot)
        self.eta_is = eta

        # NPSH required
        self.NPSH_r = NPSH_r_interp(V_dot)

        # 5. Thermodynamic outlet enthalpy
        # Compute isentropic enthalpy
        self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        h_ex_s = self.AS.hmass()
        
        # Actual enthalpy rise
        h_ex = self.su.h + (h_ex_s - self.su.h)/eta
        self.ex.set_h(h_ex)

        # 6. Power
        self.W_dot_hyd = g * H * self.m_dot  # Ideal hydraulic power [W]

        self.W_dot_shaft = self.m_dot * (h_ex - self.su.h)
        self.W.set_W_dot(self.W_dot_shaft)

    def print_results(self):
        """
        Prints the results of the pump simulation.
        """
        print("Pump Results:")
        print(f"  Volumetric Flow Rate (V_dot): {self.V_dot:.3f} m³/h")
        print(f"  Mass Flow Rate (m_dot): {self.m_dot:.3f} kg/s")
        print(f"  Outlet Pressure (P_ex): {self.ex.p:.3f} Pa")
        print(f"  Rotational Speed (N_rot): {self.W.N_rot:.3f} RPM")
        print(f"  Shaft Power (W_dot): {self.W_dot_shaft:.3f} W")
        print(f"  Ideal Hydraulic Power (W_dot_hyd): {self.W_dot_hyd:.3f} W")
        print(f"  Isentropic Efficiency (eta_is): {self.eta_is:.3f}")
        print(f"  NPSH Required: {self.NPSH_r:.3f} m")

    def plot_characteristic_curves(self, speeds_to_plot=None, rho_curve=None, n_points=200):
        """
        Plot characteristic curves (H, Power, Efficiency, NPSH_r) vs Flow (Q)
        for different rotational speeds using pump similarity laws and
        density correction according to the current working fluid.

        Uses the same nomenclature and curve definitions as the solver:
        - V_dot_curve
        - Delta_H_curve
        - eta_is_curve
        - NPSH_curve
        - N_rot_rated
        - density_curve
        """

        import matplotlib.pyplot as plt
        import numpy as np
        from scipy.interpolate import interp1d

        # ---------------------------------------------------------
        # 1. Default speeds: use current operating speed
        # ---------------------------------------------------------
        if speeds_to_plot is None:
            speeds_to_plot = [self.inputs["N_rot"]]

        # ---------------------------------------------------------
        # 2. Load manufacturer curves (reference curves)
        # ---------------------------------------------------------
        Q_base = np.array(self.params["V_dot_curve"])       # m³/h
        H_base = np.array(self.params["Delta_H_curve"])     # m
        eta_base = np.array(self.params["eta_is_curve"])    # -
        NPSH_base = np.array(self.params["NPSH_r_curve"])     # m

        # Interpolators for plotting on uniform flow grid
        H_of_Q = interp1d(Q_base, H_base, fill_value="extrapolate")
        eta_of_Q = interp1d(Q_base, eta_base, fill_value="extrapolate")
        NPSH_of_Q = interp1d(Q_base, NPSH_base, fill_value="extrapolate")

        # ---------------------------------------------------------
        # 3. Density correction (power only)
        # ---------------------------------------------------------
        rho_actual = self.su.D                         # actual fluid density
        rho_curve = rho_curve       # reference fluid density
        density_ratio = rho_actual / rho_curve

        # ---------------------------------------------------------
        # 4. Flow range for uniform plotting
        # ---------------------------------------------------------
        Q_min = np.min(Q_base)
        Q_max = np.max(Q_base)
        Q_grid = np.linspace(Q_min, Q_max, n_points)

        # ---------------------------------------------------------
        # 5. Create 4 plots
        # ---------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        axH, axP, axE, axN = axes.flatten()

        # ---------------------------------------------------------
        # 6. Loop over speeds and scale curves
        # ---------------------------------------------------------
        for N in speeds_to_plot:

            N_ref = self.params["N_rot_rated"]
            s = N / N_ref   # speed ratio

            # Similarity laws:
            # Q  ~ N
            # H  ~ N²
            # P  ~ ρ N³

            Q_scaled = Q_grid * s
            H_scaled = H_of_Q(Q_grid) * s**2
            eta_scaled = eta_of_Q(Q_grid)
            NPSH_scaled = NPSH_of_Q(Q_grid) * s**2

            # Hydraulic power at reference density
            P_h = rho_curve * 9.81 * (Q_grid / 3600) * H_of_Q(Q_grid)

            # Shaft power with similarity + density correction
            P_scaled = (P_h / eta_scaled) * s**3 * density_ratio

            # --- Plot Head ---
            axH.plot(Q_scaled, H_scaled, label=f"{N:.0f} RPM")

            # --- Plot Power ---
            axP.plot(Q_scaled, P_scaled, label=f"{N:.0f} RPM")

            # --- Plot Efficiency ---
            axE.plot(Q_scaled, eta_scaled, label=f"{N:.0f} RPM")

            # --- Plot NPSH_r ---
            axN.plot(Q_scaled, NPSH_scaled, label=f"{N:.0f} RPM")

        # ---------------------------------------------------------
        # 7. Format axes
        # ---------------------------------------------------------
        axH.set_title("Head vs Flow Rate")
        axH.set_xlabel("Flow Rate [m³/h]")
        axH.set_ylabel("Head [m]")
        axH.grid(True)
        axH.legend()

        axP.set_title("Power vs Flow Rate")
        axP.set_xlabel("Flow Rate [m³/h]")
        axP.set_ylabel("Shaft Power [W]")
        axP.grid(True)
        axP.legend()

        axE.set_title("Isentropic Efficiency vs Flow")
        axE.set_xlabel("Flow Rate [m³/h]")
        axE.set_ylabel("η_is [-]")
        axE.grid(True)
        axE.legend()

        axN.set_title("NPSH Required vs Flow")
        axN.set_xlabel("Flow Rate [m³/h]")
        axN.set_ylabel("NPSH_r [m]")
        axN.grid(True)
        axN.legend()

        plt.tight_layout()
        plt.show()
