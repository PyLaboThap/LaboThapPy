# External imports
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import CoolProp.CoolProp as CP

# Internal imports
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector

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

        V_dot_curve: Volumetric flowrate curve [m³/h]

        Delta_H_curve: Corresponding head rise curve [m]

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

        # return ["P_su", "T_su", "P_ex", "N_rot", "fluid"]
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

        return ["V_dot_curve", "Delta_H_curve", "eta_is_curve", "NPSH_curve", "N_rot_rated", "mode"]


    def solve(self):
        """
        Main solving routine that extrapolates pump curves, computes thermodynamic and performance outputs.
        """
        self.check_calculable()
        self.check_parametrized()

        if not self.calculable or not self.parametrized:
            print("Component is not calculable or not parametrized")
            return

        if self.params["mode"] == "P_N":
            if self.su.p > self.ex.p:
                print("Supply pressure is higher than exhaust pressure")
                return
            self.solve_PN()
        elif self.params["mode"] == "P_M":
            if self.su.p > self.ex.p:
                print("Supply pressure is higher than exhaust pressure")
                return
            self.solve_PM()
        elif self.params["mode"] == "M_N":
            self.solve_MN()

    def solve_PN(self):
        self.AS = CP.AbstractState("HEOS", self.su.fluid)

        g = 9.81  # Gravity [m/s²]

        self.V_dot_flag = 0  # Flag if flow is outside curve range

        # Scale curves with respect to current speed
        self.V_dot_curve = self.params["V_dot_curve"] * (self.W.N_rot / self.params["N_rot_rated"])

        self.Delta_H_curve = (self.params["Delta_H_curve"] * (self.W.N_rot / self.params["N_rot_rated"]) ** 2)

        self.eta_is_curve = self.params["eta_is_curve"]

        self.NPSH_r_curve = (self.params["NPSH_curve"] * (self.W.N_rot / self.params["N_rot_rated"]) ** 2)

        # Create interpolation functions for extrapolated performance
        Delta_H_V = interp1d(self.Delta_H_curve, self.V_dot_curve, kind="linear", fill_value="extrapolate")

        V_eta_is = interp1d(self.V_dot_curve, self.eta_is_curve, kind="linear", fill_value="extrapolate")

        V_NPSH = interp1d(self.V_dot_curve, self.NPSH_r_curve, kind="linear", fill_value="extrapolate")

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
        self.NPSH_r = V_NPSH(self.V_dot)

        # Isentropic outlet enthalpy
        self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        h_ex_s = self.AS.hmass()

        # Actual outlet enthalpy
        h_ex = self.su.h + (h_ex_s - self.su.h) / self.eta_is
        self.ex.set_h(h_ex)

        # Ideal hydraulic power
        W_dot_hyd = g * self.Delta_H * self.m_dot # W

        # Shaft power
        self.W_dot_wf = self.su.m_dot * (self.ex.h - self.su.h)  # W
        self.W.set_W_dot(self.W_dot_wf)
        eta_is_bis = W_dot_hyd / self.W_dot_wf

        # Check for impossible operation
        if self.V_dot_flag:
            print("Flowrate outside possible range. Impossible operation.")
            return


    def solve_PM(self):
        self.AS = CP.AbstractState("HEOS", self.su.fluid)
        g = 9.81

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
        NPSH_rated = self.params["NPSH_curve"]

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
        self.NPSH_r_curve = NPSH_rated * (N_rot / N_rated)**2
        self.eta_is_curve = self.params["eta_is_curve"]

        # Interpolators
        eta_interp = interp1d(self.V_dot_curve, self.eta_is_curve,
                            kind="linear", fill_value="extrapolate")
        NPSH_interp = interp1d(self.V_dot_curve, self.NPSH_r_curve,
                            kind="linear", fill_value="extrapolate")
        
        # NPSH
        NPSH = NPSH_interp(self.V_dot)
        self.NPSH_r = NPSH

        # Efficiency
        eta_is = eta_interp(self.V_dot)
        self.eta_is = eta_is

        # Update mass flow
        self.V_dot = self.V_dot

        # Thermodynamics
        self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        h_ex_s = self.AS.hmass()
        h_ex = self.su.h + (h_ex_s - self.su.h) / eta_is

        self.ex.set_h(h_ex)

        # Power
        self.W_dot_wf = self.m_dot * (h_ex - self.su.h)
        self.W.set_W_dot(self.W_dot_wf)

    def solve_MN(self):
        self.AS = CP.AbstractState("HEOS", self.su.fluid)
        g = 9.81

        # Known inputs
        m_dot = self.su.m_dot
        rho = self.su.D
        N = self.W.N_rot
        N_rated = self.params["N_rot_rated"]

        # Convert mass flow to volumetric flow (m³/h)
        V_dot = m_dot / rho * 3600
        self.V_dot = V_dot

        # Rated curves
        V_dot_rated = self.params["V_dot_curve"]
        H_rated = self.params["Delta_H_curve"]
        eta_rated = self.params["eta_is_curve"]
        NPSH_rated = self.params["NPSH_curve"]

        # Interpolators for rated curves
        H_interp = interp1d(V_dot_rated, H_rated, kind="linear", fill_value="extrapolate")
        eta_interp = interp1d(V_dot_rated, eta_rated, kind="linear", fill_value="extrapolate")
        NPSH_interp = interp1d(V_dot_rated, NPSH_rated, kind="linear", fill_value="extrapolate")

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
        self.NPSH_r = NPSH_interp(V_dot)

        # 5. Thermodynamic outlet enthalpy
        # Compute isentropic enthalpy
        self.AS.update(CP.PSmass_INPUTS, self.ex.p, self.su.s)
        h_ex_s = self.AS.hmass()
        
        # Actual enthalpy rise
        h_ex = self.su.h + (h_ex_s - self.su.h)/eta
        self.ex.set_h(h_ex)

        # 6. Shaft power
        self.W_dot_wf = m_dot * (h_ex - self.su.h)
        self.W.set_W_dot(self.W_dot_wf)





    def print_results(self):
        """
        Prints the results of the pump simulation.
        """
        print("Pump Results:")
        print(f"  Volumetric Flow Rate (V_dot): {self.V_dot:.3f} m³/h")
        print(f"  Mass Flow Rate (m_dot): {self.m_dot:.3f} kg/s")
        print(f"  Outlet Pressure (P_ex): {self.ex.p:.3f} Pa")
        print(f"  Rotational Speed (N_rot): {self.W.N_rot:.3f} RPM")
        print(f"  Shaft Power (W_dot): {self.W_dot_wf:.3f} W")
        print(f"  Isentropic Efficiency (eta_is): {self.eta_is:.3f}")
        print(f"  NPSH Required: {self.NPSH_r:.3f} m")

