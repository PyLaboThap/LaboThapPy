from labothappy.connector.mass_connector import MassConnector
from labothappy.correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model
from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
from labothappy.toolbox.turbomachinery.mean_line_axial_turbine_mapping import map_plot, map_plot_clean, plot_power_eta_vs_mdot, filter_sparse_by_proximity

from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution
import pyswarms as ps

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

class AxialTurbineMeanLine(BaseComponent):
    """
    **Component**: Mean-Line 1D Axial Turbine Model
 
    **Model**: Steady-state multi-stage mean-line turbine model (1D) with Aungier loss model.
 
    **Description**:
 
        This model simulates a large-scale axial flow turbine operating under steady-state conditions.
        It computes the isentropic efficiency, the shaft power and the outlet pressure based on inlet
        conditions (flow rate included) and an imposed rotational speed. Performance maps can be
        approximated from this model with a built-in method (CFD is however advised for more precision).
        Parameters for this model can be generated from the AxialTurbineMeanLineDesign sizing model.
 
    **Assumptions**:
 
        - Steady-state, one-dimensional flow.
        - Losses expressed using the Aungier model.
        - CoolProp is used for accurate fluid property evaluation.
        - No heat losses to surroundings.
 
    **Connectors**:
 
        su (MassConnector): Supply (inlet) side of the turbine.
 
        ex (MassConnector): Exhaust (outlet) side of the turbine.
 
        W (WorkConnector): Shaft power output from the turbine.
 
    **Parameters**:
 
        r_m: Turbine mean radius [m]
 
        nStages: Number of stages [-]
 
        damping: Damping factor for stage iterations [-]
 
        delta_tip: Blade tip clearance [m]
 
        N_lw: Number of lashing wires [-]
 
        D_lw: Lashing wire diameter [m]
 
        e_blade: Blade roughness [m]
       
        mdot_rated: Rated mass flow rate [kg/s] (For map generation)
 
        DP_rated: Rated pressure ratio [-] (For map generation)
 
        N_rot_rated: Rated rotational speed [rpm] (For map generation)
 
    **Stage Parameters (one value per stage, same parameters for rotor blades with R suffix)**:
 
        h_blade_S: Stator Blade height [m]
 
        chord_S: Stator chord length [m]
 
        xhi_S1: Stator inlet blade angle [rad]
 
        xhi_S2: Stator outlet blade angle [rad]
 
        pitch_S: Stator blade pitch [m]
           
        o_S: Stator throat opening [m]
       
        A_th: Throat flow area [m²]
       
        t_TE_S: Stator blade trailing edge thickness [m]
 
        t_blade_S: Stator blade thickness [m]
 
        n_blade_S: Stator blade number [-]
       
        R_c_S = Stator blade suction side radius of curvature [m]        
 
    **Inputs**:
 
        m_dot: Mass flow rate [kg/s]
       
        P_su: Inlet pressure [Pa]
 
        T_su or h_su: Inlet temperature [K] or enthalpy [J/kg]
 
        fluid: Working fluid [-]
 
        N_rot: Actual shaft rotational speed [rpm]
 
        P_ex: Outlet pressure [Pa] (For map generation)
 
    **Outputs**:
 
        h_ex: Outlet enthalpy [J/kg]
 
        eta_is: Isentropic efficiency [-]
 
        W_dot: Shaft work output [W]
 
        P_ex: Exhaust pressure [Pa] (Except for map generation)
 
    **Notes**:
 
        - Outlet State is the total state.
        - No dynamic behavior is included; suitable for steady-state energy system simulations.
       
    """
    
    def __init__(self, fluid):
        super().__init__()
        
        # Inputs
        self.inputs = {}
        
        # Params
        self.params = {}  

        # Abstract State 
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Blade Dictionnary
        self.stages = []

        # Velocity Triangle Data
        self.Vel_Tri_Last_Stage = {}

        self.su = MassConnector()
        self.ex = MassConnector()
        self.W = WorkConnector()
        
        self.Dh0_stage_guess = 0

    def get_required_inputs(self):
        """
        Returns a list of required input variable names.
        Used to check if the model has enough data to run.
        """
        return ["P_su", "T_su", "m_dot", "N_rot", "fluid"]
    
    def get_map_required_inputs(self):
        """
        Returns a list of required input variable names.
        Used to check if the model has enough data to run.
        """
        return ["P_su", "P_ex", "T_su", "m_dot", "N_rot", "fluid"]
    
    def get_required_parameters(self):
        """
        Returns a list of required parameters needed for model execution.
        """
        return ["r_m", "nStages", "damping", "delta_tip", "N_lw", "D_lw", "e_blade"]
    
    def get_map_required_parameters(self):
        """
        Returns a list of required parameters needed for model execution.
        """
        return ["r_m", "nStages", "mdot_rated", "DP_rated",
            "N_rot_rated", "damping", "delta_tip", "N_lw",
            "D_lw", "e_blade"]

    # ---------------- Stage Sub Class ----------------------------------------------------------------------
    
    class stage(object):
        
        def __init__(self, fluid):
            self.total_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3])
            self.static_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3])
            self.AS = CP.AbstractState('HEOS', fluid)
            
            self.eta_is_R = None
            self.eta_is_S = None
            
            self.A_flow_S = None
            self.A_flow_R = None
            
            self.h_blade_S = None
            self.h_blade_R = None
            
            self.chord_S = None
            self.chord_R = None
            
            self.stage = None
            self.AR = None
            
            self.xhi_S1 = None
            self.xhi_S2 = None
            
            self.xhi_R1 = None
            self.xhi_R2 = None
            
            self.Vel_Tri_R = {}
            self.Vel_Tri_S = {}
            
        def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.total_states.loc[position, 'H'] = self.AS.hmass()
            self.total_states.loc[position, 'S'] = self.AS.smass()
            self.total_states.loc[position, 'P'] = self.AS.p()
            self.total_states.loc[position, 'D'] = self.AS.rhomass()
            try:
                self.total_states.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                self.total_states.loc[position, 'A'] = -1
            self.total_states.loc[position, 'V'] = self.AS.viscosity()        
            
            return
        
        def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            
            self.static_states.loc[position, 'H'] = self.AS.hmass()
            self.static_states.loc[position, 'S'] = self.AS.smass()
            self.static_states.loc[position, 'P'] = self.AS.p()
            self.static_states.loc[position, 'D'] = self.AS.rhomass()
            try:
                self.static_states.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                self.static_states.loc[position, 'A'] = -1
            self.static_states.loc[position, 'V'] = self.AS.viscosity()         

            return

    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_stage_parameters(self, **parameters):
        """
        Assign stage parameters from arrays.
        If a stage doesn't exist, instantiate it.
        """
        # number of stages required = longest array among all parameters
        n_stages = max(len(arr) for arr in parameters.values())
    
        # ensure self.stages exists and has enough elements
        if not hasattr(self, "stages"):
            self.stages = []
        while len(self.stages) < n_stages:
            # ⬇️ replace Stage() with your actual Stage class constructor
            self.stages.append(self.stage(self.fluid))
    
        # assign parameters
        for key, array in parameters.items():
            for i in range(len(array)):
                setattr(self.stages[i], key, array[i])
    
    # # --- 1) Faire un snapshot picklable de la machine ---
    def _snapshot_from_machine(self):
        base_inputs = dict(
            P_su = float(self.inputs['P_su']),
            T_su = float(self.inputs['T_su']),
            P_ex = float(self.inputs['P_ex']),
            m_dot = float(self.inputs['m_dot']),
            N_rot = float(self.inputs['N_rot']),
            fluid = self.fluid
        )
        base_params = dict(self.params)  # copie légère

        stage_params = {}
        import numbers, numpy as np
        def _ok(v): return (v is None) or isinstance(v, (numbers.Number, np.floating, np.integer))
        if self.stages:
            keys = set().union(*(vars(st).keys() for st in self.stages))
            blacklist = {'AS','total_states','static_states','Vel_Tri_S','Vel_Tri_R',
                          'eta_is_R','eta_is_S','M1_S','M2_S','M2_R','M3_R',
                          'Y_vec_S','Y_vec_R','delta_S','delta_R','beta_g_S','beta_g_R'}
            for k in sorted(keys - blacklist):
                vals = [getattr(st, k, None) for st in self.stages]
                if any((_ok(v) for v in vals)) and all((_ok(v) for v in vals)):
                    stage_params[k] = vals
        return base_inputs, base_params, stage_params
    
    def _wegstein_solve(self, f, x0, tol=1e-8, max_iter=1000, q_min=-5.0, q_max=0.0):
        """
        Wegstein iteration to solve x = f(x).
    
        Uses a secant-based mixing weight computed from the last two iterates.
        Falls back to a damped FPI step on the first iteration (no history yet).
    
        Parameters
        ----------
        f         : callable(x) -> x_out  (both in scaled units, e.g. *1e-5)
        x0        : initial guess (numpy array)
        tol       : convergence tolerance on sum of relative residuals
        max_iter  : maximum iterations before raising RuntimeError
        q_min/max : clamping bounds for the secant slope (controls stability)
        """
        x_in      = np.array(x0, dtype=float)
        x_out_prev = None
        x_in_prev  = None
    
        for c in range(max_iter):
            x_out = np.array(f(x_in), dtype=float)
    
            res = np.sum(np.abs((x_in - x_out) / np.where(np.abs(x_out) > 1e-30, x_out, 1.0)))
            
            # print(f"res : {res}")
            
            if res < tol:
                # Run one final evaluation at the converged point
                f(x_out)
                return x_out
    
            if x_in_prev is None:
                # First step: plain damped FPI (no history to build secant)
                x_new = (1.0 - self.params['damping']) * x_in + self.params['damping'] * x_out
            else:
                dx   = x_in  - x_in_prev
                df   = x_out - x_out_prev
                # Component-wise Wegstein weight; fall back to damped FPI where dx≈0
                mask = np.abs(dx) > 1e-30
                q    = np.where(mask, df / dx, 0.0)
                q    = np.clip(q, q_min, q_max)
                w    = np.where(mask, q / (q - 1.0), self.params['damping'])
                x_new = w * x_in + (1.0 - w) * x_out
    
            x_in_prev  = x_in.copy()
            x_out_prev = x_out.copy()
            x_in       = x_new
    
        raise RuntimeError(
            f"Wegstein failed to converge after {max_iter} iterations "
            f"(final residual: {res:.3e})"
        )
    
    # ---------------- Loss Models ------------------------------------------------------------------------

    def stator_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]

        # 1) Guess outlet state
        h_static_out = x[0]*1e5
        p_static_out = x[1]*1e5
                
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        stage.Vel_Tri_S['u'] = u = self.u
        
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][2]*A_flow)
        
        if self.curr_stage_index == 0:
            stage.Vel_Tri_S['alpha1'] = alpha1 = stage.xhi_S1
            stage.Vel_Tri_S['vu1'] = vu1 = vm*np.tan(alpha1)
        else:
            stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
            stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u 
            stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
        
        
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][1]

        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        h02 = h0in
        
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][2]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)

        # 4) Compute cord, aspect ratio, blade pitch and blade number
        
        stage.Re_s = stage.chord_S*(stage.static_states['D'][2]*vm)/(stage.static_states['V'][2])
        stage.AR_S = stage.h_blade_S/stage.chord_S
                
        # 5) Loss model
        
        stage.beta_g_S = np.arcsin(stage.o_S/stage.pitch_S)
        
        stage.Y_vec_S = aungier_loss_model(stage.Vel_Tri_S['alpha1'], stage.Vel_Tri_S['alpha2'], stage.beta_g_S*180/np.pi, stage.xhi_S1, stage.chord_S, 
                               0, self.params['D_lw'], self.params['e_blade'], stage.h_blade_S, stage.static_states['V'][2], 
                               stage.M1_S, stage.M2_S, self.params['N_lw'], stage.R_c_S, stage.static_states['D'][2], stage.pitch_S, stage.t_blade_S, stage.t_TE_S,
                               vm, v2,1)
        
        self.compute_deviation_stator(stage)
        alpha2_calc = stage.xhi_S2 + stage.delta_S
        
        v2_new = vm/np.cos(alpha2_calc)
        
        Y = stage.Y_vec_S['Y_tot']
        p0_out = (stage.total_states['P'][1] + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = stage.total_states['H'][2]-v2_new**2/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
                
        pout_calc = stage.static_states['P'][2]

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        return np.array([hout, pout_calc])*1e-5 # (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def rotor_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]
        
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x*1e5
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)
        
        stage.Vel_Tri_R['u'] = u = self.u
        
        A_flow = stage.h_blade_R*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_R['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][3]*A_flow)
        stage.Vel_Tri_R['vu2'] = vu2 = vm*np.tan(stage.Vel_Tri_R['alpha2'])    
        
        stage.Vel_Tri_R['wu2'] = wu2 = vu2 - u
        stage.Vel_Tri_R['w2'] = w2 = np.sqrt(wu2**2 + vm**2)
        stage.M2_R = w2/stage.static_states['A'][2]
        stage.Vel_Tri_R['beta2'] = np.arctan(wu2/vm)
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][2]
        h0in = hin + (w2**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        h03 = stage.total_states['H'][2]
        stage.Vel_Tri_R['w3'] = w3 = np.sqrt(2*(h03 - h_static_out))
        stage.M3_R = w3/stage.static_states['A'][3]
        stage.Vel_Tri_R['beta3'] = -np.arccos(vm/w3)
                
        # 4) Compute cord, aspect ratio, pitch and blade number
        stage.Re_r = stage.chord_R*(stage.static_states['D'][3]*vm)/(stage.static_states['V'][3])
        stage.AR_R = stage.h_blade_R/stage.chord_R
                
        # 5) Loss model
        
        stage.beta_g_R = np.arcsin(stage.o_R/stage.pitch_R)  # mid-passage metal angle        
                
        stage.Y_vec_R = aungier_loss_model(-stage.Vel_Tri_R['beta2'], -stage.Vel_Tri_R['beta3'], stage.beta_g_R*180/np.pi, -stage.xhi_R1, stage.chord_R, 
                               self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.static_states['V'][3], 
                               stage.M2_R, stage.M3_R, self.params['N_lw'], stage.R_c_R, stage.static_states['D'][3], stage.pitch_R, stage.t_blade_R, stage.t_TE_R,
                               vm, w3,1)

        self.compute_deviation_rotor(stage)
        beta3_calc = stage.xhi_R2 + stage.delta_R

        w3_new = vm/np.cos(beta3_calc)

        Y = stage.Y_vec_R['Y_tot']
                
        p0_out = (stage.total_states['P'][2] + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.total_states['S'][3]
        
        hout = h0in-(w3_new**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 3)
        
        pout_calc = stage.static_states['P'][3]

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][2])
        hout_s = self.AS.hmass()

        stage.eta_is_R = (stage.static_states['H'][2]-stage.static_states['H'][3])/(stage.static_states['H'][2]-hout_s)
        
        return np.array([hout, pout_calc])*1e-5 # (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def last_blade_row_system(self, x):
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x*1e5
        
        stage = self.stages[-1]
        
        stage.Vel_Tri_S['u'] = u = self.u
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][2]*A_flow)
        
        stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
        stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u 
        stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
        
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][1]

        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)    
                
        h02 = h0in
        
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][2]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)
        
        # 5) Estimate pressure losses 
        AR_S = stage.h_blade_S/stage.chord_S
        solidity = (stage.chord_S/stage.pitch_S)
        
        a = 0.0117 # NACA blade - 0.007 : C.4 circular-arc blade
        
        D_e = (np.cos(alpha2)/np.cos(alpha1))*(1.12+a*(alpha1 - stage.xhi_S1)+0.61*np.cos(alpha1)**2 / solidity * (np.tan(alpha1)-np.tan(alpha2)))
        
        P_cst = np.cos(alpha2)/2 * solidity * (v1/v2)**2 # Profile Constant
        
        Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
        # 5.2) Cohen : Endwall losses
        EW_Cst = np.cos((alpha1+alpha2)/2)**3 / np.cos(alpha1)**2  # Endwall Constant

        Yew = 0.02*(solidity/AR_S)/EW_Cst

        # Pressure loss 
        DP_loss = (Yp+Yew)*(v1**2)*stage.static_states['D'][1]/2
        p0_out = stage.total_states['P'][1]-DP_loss
                
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(v2**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
        
        pout_calc = stage.static_states['P'][2]

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][1])
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.static_states['H'][1]-stage.static_states['H'][2])/(stage.static_states['H'][1]-hout_s)

        # print(f"h0in: {h0in}")
        # print(f"h1: {stage.static_states['H'][1]}")
        # print(f"kinetic1: {(self.Vel_Tri_Last_Stage['vu1']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2}")
        # print(f"h2: {stage.static_states['H'][2]}")
        # print(f"kinetic2: {(self.Vel_Tri_Last_Stage['vu2']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2}")

        return np.array([hout, pout_calc])*1e-5 # return (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def compute_deviation_stator(self, stage):
        
        delta_0S = np.arcsin((stage.o_S/stage.pitch_S)*(1+(1-stage.o_S/stage.pitch_S)*(2*stage.beta_g_S/np.pi)**2)) - stage.beta_g_S
        
        if stage.M2_S <= 0.5:
            stage.delta_S = delta_0S
        else:
            X = 2*stage.M2_S-1
            stage.delta_S = delta_0S*(1-10*X**3 + 15*X**4 - 6*X**5)
        
        return 

    def compute_deviation_rotor(self, stage):
                
        delta_0R = np.arcsin((stage.o_R/stage.pitch_R)*(1+(1-stage.o_R/stage.pitch_R)*(2*stage.beta_g_R/np.pi)**2)) - abs(stage.beta_g_R)
        
        if stage.M3_R <= 0.5:
            stage.delta_R = delta_0R
        else:
            X = 2*stage.M3_R-1
            stage.delta_R = delta_0R*(1-10*X**3 + 15*X**4 - 6*X**5)

        return 

    # ---------------- Flow Computations ------------------------------------------------------------------

    def computeVelTriangle(self):

        # Velocities over u
        self.Vel_Tri['vu2OverU'] = (2*(1-self.R) + self.psi)/2
        self.Vel_Tri['vu3OverU'] = (2*(1-self.R) - self.psi)/2
        self.Vel_Tri['vmOverU']  = self.phi
        
        self.Vel_Tri['wu2OverU']  = self.Vel_Tri['vu2OverU'] - 1
        self.Vel_Tri['wu3OverU']  = self.Vel_Tri['vu3OverU'] - 1

        self.Vel_Tri['v2OverU']  = np.sqrt(self.Vel_Tri['vu2OverU']*self.Vel_Tri['vu2OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['w2OverU']  = np.sqrt(self.Vel_Tri['wu2OverU']*self.Vel_Tri['wu2OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['v3OverU']  = np.sqrt(self.Vel_Tri['vu3OverU']*self.Vel_Tri['vu3OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])
        self.Vel_Tri['w3OverU']  = np.sqrt(self.Vel_Tri['wu3OverU']*self.Vel_Tri['wu3OverU']+self.Vel_Tri['vmOverU']*self.Vel_Tri['vmOverU'])

        # Angles in radians
        self.Vel_Tri['alpha1'] = self.Vel_Tri['alpha3'] = np.arctan(self.Vel_Tri['vu3OverU']/self.Vel_Tri['vmOverU'])
        self.Vel_Tri['alpha2'] = np.arctan(self.Vel_Tri['vu2OverU']/self.Vel_Tri['vmOverU'])

        self.Vel_Tri['beta1'] = self.Vel_Tri['beta3'] = np.arctan(self.Vel_Tri['wu3OverU']/self.Vel_Tri['vmOverU'])
        self.Vel_Tri['beta2'] = np.arctan(self.Vel_Tri['wu2OverU']/self.Vel_Tri['vmOverU'])
        
        return 
    
    def computeVelTriangleLastStage(self):

        self.Vel_Tri_Last_Stage['u'] = self.Vel_Tri['u']
        self.Vel_Tri_Last_Stage['vu2'] = 0
        self.Vel_Tri_Last_Stage['vu1'] = self.Vel_Tri['vu3']
        self.Vel_Tri_Last_Stage['vm']  = self.Vel_Tri['vm']
        
        self.Vel_Tri_Last_Stage['wu2'] = self.Vel_Tri_Last_Stage['vu2'] - self.Vel_Tri_Last_Stage['u']
        self.Vel_Tri['v2'] = np.sqrt(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w2'] = np.sqrt(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w3'] = np.sqrt(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)

        # Angles in radians
        self.Vel_Tri_Last_Stage['alpha1'] = self.Vel_Tri['alpha3'] 
        self.Vel_Tri_Last_Stage['alpha2'] = 0

        self.Vel_Tri_Last_Stage['beta1'] = self.Vel_Tri['beta3']
        self.Vel_Tri_Last_Stage['beta2'] = np.arctan(self.Vel_Tri['u']/self.Vel_Tri['vm'])
        
        return 
    
    def computeBladeRow(self,stage_index,row_type):
        stage = self.stages[stage_index]
        
        self.curr_stage_index = stage_index
               
        if row_type == 'S': # Stator
            
            # print("Stator")
            
            if 'P_ex' not in self.inputs:
                RP_1_row = 5**(1/(2*self.nStages)) 
            else:
                RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
            
            if self.Dh0_stage_guess !=0:
                h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2    
            else:
                h_out_guess = stage.static_states['H'][1]*0.99
                
            pout_guess = stage.static_states['P'][1]/RP_1_row
            # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])         
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
            
            x_out = self._wegstein_solve(self.stator_blade_row_system, x0_disc)
            self.stator_blade_row_system(x_out)
                        
            # print(f'Y_S : {stage.Y_vec_S}')

        else: # Rotor

            # print("Rotor")

            if 'P_ex' not in self.inputs:
                RP_1_row = 5**(1/(2*self.nStages)) 
            else:
                RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
                
            if self.Dh0_stage_guess !=0:
                h_out_guess = stage.static_states['H'][2] - self.Dh0_stage_guess/2    
            else:
                h_out_guess = stage.static_states['H'][2]*0.99
            
            pout_guess = stage.static_states['P'][2]/RP_1_row
            # sol = minimize(self.rotor_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])    
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
                        
            x_out = self._wegstein_solve(self.rotor_blade_row_system, x0_disc)
            self.rotor_blade_row_system(x_out)
            
            # print(f'Y_R : {stage.Y_vec_R}')

        return
            
    def computeRepeatingStages(self):
                
        self.nStages = self.params['nStages']
        
        for i in range(int(self.nStages)):
                    
            if i == 0:
                self.computeBladeRow(i, 'S')
                
                self.compute_deviation_stator(self.stages[i])
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']
                
                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

                self.Dh0_stage_guess = self.stages[i].total_states['H'][1] - self.stages[i].total_states['H'][3]
                
            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
                self.computeBladeRow(i, 'S')
                self.stages[i].Vel_Tri_R['alpha2'] = self.stages[i].Vel_Tri_S['alpha2']

                self.computeBladeRow(i, 'R')
                self.stages[i+1].Vel_Tri_S['beta1'] = self.stages[i].Vel_Tri_R['beta3']

        return
    
    def computeLastStage(self):
        stage = self.stages[-1]
        
        stage.static_states.loc[1] = self.stages[-2].static_states.loc[3]
        
        if 'P_ex' not in self.inputs:
            RP_1_row = 5**(1/(2*self.nStages)) 
        else:
            RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
                
        h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2  
        pout_guess = stage.static_states['P'][1]/RP_1_row
        # sol = minimize(self.last_blade_row_system, x0=(h_out_guess,pout_guess), bounds=[(self.stages[-1].static_states['H'][1], h_out_guess), (self.stages[-1].static_states['P'][1], pout_guess)])         
        
        # Initial guess vector
        x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

        res = 1
        x_in = x0_disc
        
        c = 0
            
        while res > 1e-8:

            if c > 1000:
                raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")                

            x_out = self.last_blade_row_system(x_in) 

            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            
            x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                        
            c += 1
        
        self.last_blade_row_system(x_out)
        
        return
    
    def generate_map_m_dot_N_rot(
        self,
        m_dot_grid=None,
        N_rot_grid=None,
        *,
        m_dot_range=None,  # (min, max, n)
        N_rot_range=None,  # (min, max, n)
        fixed_P_su=None,   # override supply pressure for the map (Pa)
        fixed_T_su=None,   # override supply temperature for the map (K)
        fixed_P_ex=None,   # override exhaust/static outlet pressure (Pa)
        per_point_hook=None,  # callable(self) to tweak params before solve
        max_retries=2,
        mach_limit=1.2,
        pressure_tol=0.02,
        verbose=False
        ):
        """
        Build an operation map by sweeping mass flow and speed.
    
        Parameters
        ----------
        m_dot_grid : iterable of float
            Mass-flow values [kg/s]. If None, use m_dot_range.
        N_rot_grid : iterable of float
            Rotational speeds [rpm]. If None, use N_rot_range.
        m_dot_range : (mmin, mmax, n)
            Range spec if m_dot_grid is None.
        N_rot_range : (nmin, nmax, n)
            Range spec if N_rot_grid is None.
        fixed_P_su, fixed_T_su, fixed_P_ex : float or None
            If given, overrides the current inputs for the whole map (Pa, K, Pa).
        per_point_hook : callable(self) -> None
            Called right before `solve()` for each (ṁ, N) point to adjust geometry,
            clearances, or params based on current operating point.
        max_retries : int
            Retries for a point if convergence fails (e.g., tweak damping).
        mach_limit : float
            Warn if any stage exit Mach exceeds this value.
        pressure_tol : float
            Relative tolerance on outlet pressure: warn if |P_ex_calc - P_ex_target| / P_ex_target > pressure_tol.
        verbose : bool
            Print progress.
    
        Returns
        -------
        pandas.DataFrame
            One row per (ṁ, N). Columns include:
            ['m_dot','N_rot','P_su','T_su','P_ex_target','P_ex_calc',
             'PR','W_dot','eta_is','converged','mach_warn','pressure_warn','notes']
        """
        import numpy as _np
        import pandas as _pd
        
        # Build grids if ranges given
        if m_dot_grid is None:
            if m_dot_range is None:
                # default: ±30% around current
                m0 = self.inputs.get('m_dot', 1.0)
                m_dot_grid = _np.linspace(0.7*m0, 1.3*m0, 9)
            else:
                m_dot_grid = _np.linspace(*m_dot_range)
        else:
            m_dot_grid = _np.array(list(m_dot_grid), dtype=float)
    
        if N_rot_grid is None:
            if N_rot_range is None:
                N0 = self.inputs.get('N_rot', 1000.0)
                N_rot_grid = _np.linspace(0.6*N0, 1.2*N0, 9)
            else:
                N_rot_grid = _np.linspace(*N_rot_range)
        else:
            N_rot_grid = _np.array(list(N_rot_grid), dtype=float)
    
        # Cache original boundary conditions to restore later
        _P_su0 = self.inputs.get('P_su', None)
        _T_su0 = self.inputs.get('T_su', None)
        _P_ex0 = self.inputs.get('P_ex', None)
    
        rows = []
        total_pts = len(m_dot_grid) * len(N_rot_grid)
        idx = 0
    
        for N in N_rot_grid:
            for m in m_dot_grid:
                idx += 1
                if verbose:
                    print(f"[{idx}/{total_pts}] N={N:.2f} rpm, ṁ={m:.3f} kg/s")
    
                # Set operating point
                self.set_inputs(
                    m_dot = float(m),
                    P_su  = float(fixed_P_su if fixed_P_su is not None else _P_su0),
                    T_su  = float(fixed_T_su if fixed_T_su is not None else _T_su0),
                    N_rot = float(N),
                    fluid = self.fluid,
                    P_ex  = float(fixed_P_ex if fixed_P_ex is not None else _P_ex0),
                )
    
                # Optionally let user tweak geometry/params per point
                if per_point_hook is not None:
                    try:
                        per_point_hook(self)
                    except Exception as e:
                        # Non-fatal: record the note and proceed
                        hook_note = f"per_point_hook failed: {e}"
                    else:
                        hook_note = ""
    
                # Try to solve with limited retries (e.g., adjust damping on the fly)
                converged = False
                notes = []
                for attempt in range(max_retries + 1):
                    try:
                        # Example adaptive damping: increase damping on retries
                        if attempt > 0 and 'damping' in self.params:
                            self.params['damping'] = min(0.8, self.params['damping'] * 1.5)
    
                        # Re-init the *first* stage inlet each attempt to reduce drift
                        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
    
                        self.solve()
                        converged = True
                        break
                    except Exception as e:
                        notes.append(f"attempt {attempt}: {e}")
    
                # Collect metrics (even if failed; fill NaNs)
                if converged:
                    try:
                        P_su = self.inputs['P_su']
                        T_su = self.inputs['T_su']
                        P_ex_target = self.inputs['P_ex']
                        P_ex_calc = float(self.stages[-1].static_states['P'][2])
                
                        W_dot = getattr(self, 'W_dot', _np.nan)
                        eta   = getattr(self, 'eta_is', _np.nan)
                
                        # >>> INSERT THESE THREE LINES HERE <<<
                        RP_target = P_su / P_ex_target if P_ex_target else _np.nan
                        RP_calc   = P_su / P_ex_calc   if _np.isfinite(P_ex_calc) and P_ex_calc else _np.nan
                        PR = RP_target  # keep old name for compatibility
                        # >>> END INSERT <<<
                
                        # ---- health checks (unchanged) ----
                        machs = []
                        for st in self.stages:
                            for label in ('M2_S','M3_R'):
                                if hasattr(st, label):
                                    val = getattr(st, label)
                                    if val is not None:
                                        machs.append(val)
                        mach_warn = bool(len([x for x in machs if _np.isfinite(x) and x > mach_limit]) > 0)
                
                        pressure_warn = False
                        if P_ex_target and _np.isfinite(P_ex_calc):
                            rel_err = abs(P_ex_calc - P_ex_target)/P_ex_target
                            pressure_warn = rel_err > pressure_tol
                            if pressure_warn:
                                notes.append(f"P_ex mismatch {rel_err:.2%}")
    
                    except Exception as e:
                        # If something unexpected happens while reading results
                        P_su = self.inputs.get('P_su', _np.nan)
                        T_su = self.inputs.get('T_su', _np.nan)
                        P_ex_target = self.inputs.get('P_ex', _np.nan)
                        P_ex_calc = _np.nan
                        PR = _np.nan
                        W_dot = _np.nan
                        eta = _np.nan
                        mach_warn = False
                        pressure_warn = True
                        notes.append(f"post-process error: {e}")
    
                else:
                    # Not converged
                    P_su = self.inputs.get('P_su', _np.nan)
                    T_su = self.inputs.get('T_su', _np.nan)
                    P_ex_target = self.inputs.get('P_ex', _np.nan)
                    P_ex_calc = _np.nan
                    PR = _np.nan
                    W_dot = _np.nan
                    eta = _np.nan
                    mach_warn = False
                    pressure_warn = True
    
                rows.append(dict(
                    m_dot=float(m),
                    N_rot=float(N),
                    P_su=float(P_su) if P_su is not None else _np.nan,
                    T_su=float(T_su) if T_su is not None else _np.nan,
                    P_ex_target=float(P_ex_target) if P_ex_target is not None else _np.nan,
                    P_ex_calc=float(P_ex_calc) if _np.isfinite(P_ex_calc) else _np.nan,
                    PR=float(PR) if _np.isfinite(PR) else _np.nan,
                    W_dot=float(W_dot) if _np.isfinite(W_dot) else _np.nan,
                    eta_is=float(eta) if _np.isfinite(eta) else _np.nan,
                    converged=bool(converged),
                    mach_warn=bool(mach_warn),
                    pressure_warn=bool(pressure_warn),
                    notes="; ".join([hook_note] + notes) if 'hook_note' in locals() else "; ".join(notes)
                ))
    
        # Restore original boundary conditions
        self.set_inputs(
            m_dot=self.inputs['m_dot'],
            P_su=_P_su0,
            T_su=_T_su0,
            N_rot=self.inputs['N_rot'],
            fluid=self.fluid,
            P_ex=_P_ex0
        )
    
        df = _pd.DataFrame(rows)
    
        # Useful sorted ordering
        df.sort_values(by=['N_rot','m_dot'], inplace=True, ignore_index=True)
        return df

    
    def solve(self):
        
        # Check if the component is calculable and parametrized
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized): # If the component is not calculable and/or not parametrized
            self.solved = False
            print("AxialTurbineMeanLine could not be solved. It is not calculable and/or not parametrized")
            self.print_setup()
            return

        if self.su.m_dot is not None:
            self.inputs['m_dot'] = self.su.m_dot
        
        self.omega_rads = 2*np.pi*self.inputs['N_rot']/60
        self.u = self.omega_rads*self.params['r_m']
        
        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
        # self.stages.append(self.stage(self.fluid))
        
        self.computeRepeatingStages()
        
        self.computeLastStage()
        
        hin = self.stages[0].total_states['H'][1]
        hout = self.stages[-1].static_states['H'][2]
        
        self.h_ex = hout
        self.p_ex = self.stages[-1].static_states['P'][2]
        
        self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])

        hout_s = self.AS.hmass()
        
        self.W_dot = self.inputs['m_dot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - hout_s)
        self.update_connectors()

        self.solved = True
                
        return 

    def update_connectors(self):
        """Update the connectors with the calculated values."""
        
        self.ex.reset()
        
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_p(self.p_ex)
        self.ex.set_h(self.h_ex)

        self.W.set_W_dot(self.W_dot)
        self.W.set_N_rot(self.inputs['N_rot'])
        
        return
    
#%%

def _eval_point_from_snapshot(m, N, base_inputs, base_params, stage_params):
    try:
        turb = AxialTurbineMeanLine(base_inputs['fluid'])
        turb.set_inputs(
            m_dot=float(m),
            P_su=base_inputs['P_su'],
            T_su=base_inputs['T_su'],
            N_rot=float(N),
            fluid=base_inputs['fluid'],
            P_ex=base_inputs['P_ex'],
        )
        turb.set_parameters(**base_params)
        
        if stage_params:
            turb.set_stage_parameters(**stage_params)

        # print(f"Computing for N = {float(N)} RPM and m_dot = {float(m)} kg/s \n ")

        turb.solve()

        P_ex_calc = float(turb.stages[-1].static_states['P'][2])
        RP_calc = turb.inputs['P_su'] / P_ex_calc if P_ex_calc else np.nan
        RP_target = turb.inputs['P_su'] / turb.inputs['P_ex'] if turb.inputs.get('P_ex') else np.nan

        # ✅ Ajout du print de debug
        # print(f"[DONE] m={m:.2f}, N={N:.1f}, W={turb.W_dot:.2e}, eta={turb.eta_is:.3f}, RP={RP_calc:.2f}")
        
        if turb.eta_is < 0.3:
            return dict(
                m_dot=float(m), N_rot=float(N),
                P_su=float(base_inputs.get('P_su', np.nan)),
                T_su=float(base_inputs.get('T_su', np.nan)),
                P_ex_target=float(base_inputs.get('P_ex', np.nan)),
                P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
                W_dot=np.nan, eta_is=np.nan, converged=False
            )
            
        else:
            return dict(
                m_dot=float(m),
                N_rot=float(N),
                P_su=float(turb.inputs['P_su']),
                T_su=float(turb.inputs['T_su']),
                P_ex_target=float(turb.inputs['P_ex']),
                P_ex_calc=P_ex_calc,
                RP_target=RP_target,
                RP_calc=RP_calc,
                W_dot=float(getattr(turb, 'W_dot', np.nan)),
                eta_is=float(getattr(turb, 'eta_is', np.nan)),
                converged=True,
                note=""
            )

    except Exception as e:
        return dict(
            m_dot=float(m), N_rot=float(N),
            P_su=float(base_inputs.get('P_su', np.nan)),
            T_su=float(base_inputs.get('T_su', np.nan)),
            P_ex_target=float(base_inputs.get('P_ex', np.nan)),
            P_ex_calc=np.nan, RP_target=np.nan, RP_calc=np.nan,
            W_dot=np.nan, eta_is=np.nan, converged=False, note=str(e)
        )

import os, sys
import pandas as pdproce
from tqdm import tqdm
from joblib import Parallel, delayed
import joblib
from contextlib import contextmanager
from joblib.parallel import BatchCompletionCallBack

# ---- tqdm <-> joblib bridge: update bar when each batch completes ----
@contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar."""
    class TqdmBatchCompletionCallback(BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_cb = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_cb
        tqdm_object.close()

# ----- YOUR multiprocessing map using joblib+loky -----
def generate_map_processes(machine, m_grid, N_grid, max_workers=-1, desc="Operation map"):
    """
    Multiprocessing version using joblib (loky). Compatible with Spyder/Windows.
    - max_workers: -1 => use all cores, or pass an int.
    """
    # Prevent thread oversubscription inside each worker
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")

    base_inputs, base_params, stage_params = machine._snapshot_from_machine()

    tasks = [(m, N) for N in N_grid for m in m_grid]
    total = len(tasks)

    with tqdm(total=total, desc=desc, unit="pt",
              dynamic_ncols=True, miniters=1, mininterval=0,
              ascii=True, file=sys.stdout) as bar, tqdm_joblib(bar):
        results = Parallel(
            n_jobs=max_workers,
            backend="loky",       # processes
            prefer="processes",   # explicit
        )(
            delayed(_eval_point_from_snapshot)(m, N, base_inputs, base_params, stage_params)
            for (m, N) in tasks
        )

    return pd.DataFrame(results).sort_values(['N_rot', 'm_dot'], ignore_index=True)

#%%

if __name__ == "__main__":
    
    case_study = "TCO2_ORC"
    
    if case_study == "Salah_Case":
        Turb_OD = AxialTurbineMeanLine('CO2')
        
        Turb_OD.set_inputs(
              m_dot = 655.18,
              P_su = 25000000.0,
              T_su = 923,
              N_rot = 1506.9946780513785, # 996.4330963327212,
              fluid = 'CO2', 
              P_ex = 100*1e5
            )
        
        Turb_OD.set_parameters(
            r_m = 0.261423771889,
            nStages = 12,
            mdot_rated = 655.18,
            DP_rated = 2.5,
            damping = 0.5,
            delta_tip = 0.0004,
            N_lw = 0,
            D_lw = 0,
            e_blade = 2e-06
            )
        
        Turb_OD.set_stage_parameters(
            h_blade_S = [0.05893535333, 0.06254061127, 0.06644686421, 0.07068270152, 0.07527979324, 0.08027326565, 0.08570212832, 0.09160976073, 0.09804446776, 0.105060115, 0.1127168568, 0.1210819728, 0.1254195395],
            chord_S = [0.008645525688, 0.009066896935, 0.009520522543, 0.01000898242, 0.01053510451, 0.01110199118, 0.0117130487, 0.0123720203, 0.0130830232, 0.01385059029, 0.01467971707, 0.01557591465, 0.01603830059],
            xhi_S1 = [-0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816, -0.6455466816],
            xhi_S2 = [1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265, 1.146678265],
            pitch_S = [0.006966193979, 0.007305716867, 0.007671228935, 0.008064808966, 0.008488735591, 0.008945508564, 0.009437872521, 0.009968843587, 0.01054173924, 0.0111602119, 0.01182828672, 0.01255040431, 0.01292297508],
            o_S = [0.001744270233, 0.001829283609, 0.001920804434, 0.002019353216, 0.002125500504, 0.002239872211, 0.00236315557, 0.002496105791, 0.002639553539, 0.002794413345, 0.00296169307, 0.003142504605, 0.003235792862],
            t_TE_S = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
            t_blade_S = [0.003025933991, 0.003173413927, 0.00333218289, 0.003503143847, 0.003687286579, 0.003885696912, 0.004099567045, 0.004330207105, 0.00457905812, 0.004847706601, 0.005137900974, 0.005451570129, 0.005613405206],
            n_blade_S = [236, 225, 214, 204, 194, 184, 174, 165, 156, 147, 139, 131, 127],
            R_c_S = [0.01352981312, 0.01418923794, 0.01489913922, 0.01566355437, 0.01648690899, 0.01737405813, 0.01833033244, 0.01936158987, 0.02047427367, 0.02167547758, 0.02297301932, 0.02437552349, 0.02509913423],
        
            h_blade_R = [0.06085456338, 0.06461845758, 0.06869830659, 0.07312430505, 0.07792992342, 0.08315231023, 0.08883274986, 0.09501718485, 0.1017568126, 0.1091087686, 0.1171369103, 0.1259127187, None],
            chord_R = [0.008870626618, 0.009309055087, 0.009781089003, 0.01028943789, 0.01083707309, 0.01142725573, 0.01206356808, 0.01274994866, 0.01349073167, 0.01429069144, 0.01515509243, 0.01608974592, None],
            xhi_R1 = [0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, 0.5693595771, None],
            xhi_R2 = [-1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, -1.17766032, None],
            pitch_R = [0.007763880051, 0.008147607851, 0.008560748305, 0.009005672883, 0.009484982197, 0.01000153051, 0.01055845315, 0.01115919723, 0.01180755622, 0.01250770875, 0.01326426248, 0.01408230362, None],
            o_R = [0.00229560867, 0.002409068546, 0.002531225097, 0.002662779514, 0.002804500743, 0.002957232725, 0.003121902507, 0.003299529328, 0.003491234828, 0.00369825454, 0.003921950849, 0.004163827632, None],
            t_TE_R = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, None],
            t_blade_R = [0.003104719316, 0.00325816928, 0.003423381151, 0.003601303262, 0.00379297558, 0.003999539504, 0.004222248827, 0.00446248203, 0.004721756085, 0.005001742002, 0.005304282352, 0.00563141107, None],
            n_blade_R = [212, 202, 192, 182, 173, 164, 156, 147, 139, 131, 124, 117, None],
            R_c_R = [0.01388208476, 0.01456820328, 0.01530691263, 0.01610245309, 0.01695947464, 0.0178830808, 0.01887887765, 0.01995302875, 0.02111231693, 0.02236421375, 0.02371695786, 0.02517964358, None],
            )
        
    elif case_study == "TCO2_ORC":
        Turb_OD = AxialTurbineMeanLine('CO2')
    
        turb_params = {'type': 'Axial Turbine',
          'mdot_rated': 318.437021666738,
          'Wdot_rated': 15149047.61803014,
          'N_rot_rated': 2201.9292383818693,
          'total_to_static_efficiency': 0.8918913361126354,
          'DP_rated': 2.93,
          'n_stages': 8,
          'p0_su': 15309670.5,
          'T0_su': 406.4,
          'p_ex': 5220928,
          'r_m': 0.2424432562137853,
          'delta_tip': 0.0004,
          'N_lw': 0,
          'D_lw': 0,
          'e_blade': 2e-06,
          'stator': {'h_blade_S': [0.027129712191845578,
            0.029658270219229916,
            0.03257968025854356,
            0.035956572944882276,
            0.039861618532234466,
            0.04437960089736192,
            0.049608759632974515,
            0.05566774162922002,
            0.058599071296795195],
          'chord_S': [0.013511359662391169,
            0.014028083649020674,
            0.014637713665581014,
            0.015349906293429992,
            0.016175091304622864,
            0.01712470590035643,
            0.018211169517019312,
            0.019449229232574965,
            0.020047912130238918],
          'xhi_S1': [-0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645,
            -0.6356400014580645],
          'xhi_S2': [1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192,
            1.2113211735665192],
          'pitch_S': [0.013134937068248984,
            0.013637265273227886,
            0.014229911172864795,
            0.014922262318938251,
            0.015724457913053796,
            0.016647616519272383,
            0.01770381157205685,
            0.01890737930002005,
            0.019489383064344484],
          'o_S': [0.003775183728696827,
            0.003919560611208578,
            0.004089896193751755,
            0.0042888886036597795,
            0.004519452004032509,
            0.004784782041861723,
            0.0050883488086367655,
            0.005434272757825279,
            0.005601549626355477],
          't_TE_S': [0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005],
          't_blade_S': [0.0018320275186014946,
            0.0019218083226827988,
            0.0020268089311736032,
            0.0021329017027740434,
            0.002240138294595593,
            0.0023488859174960303,
            0.0024852733652763366,
            0.0026024429397798927,
            0.0026780412271797437],
          'n_blade_S': [116, 112, 107, 102, 97, 92, 86, 81, 78],
          'R_c_S': [0.02161610169772873,
            0.02244278076806165,
            0.023418095226800897,
            0.02455749412199926,
            0.025877663488153915,
            0.02739690109176542,
            0.02913507612488498,
            0.031115781648830443,
            0.03207358239753057]},
          'rotor': {'h_blade_R': [0.028299841921147508,
            0.031006065260253488,
            0.03413278574452733,
            0.03774694200977252,
            0.04192639330728815,
            0.046762238719206164,
            0.05236044258786883,
            0.05884941443434066,
            None],
          'chord_R': [0.013752412553972845,
            0.014312062993939407,
            0.014968897751956336,
            0.01573287734406901,
            0.016614837213121982,
            0.01762675393333128,
            0.018781759644470535,
            0.020095591613924796,
            None],
          'xhi_R1': [0.7347958438124123,
            0.7347958438124123,
            0.7347958438124123,
            0.7347958438124123,
            0.7347958438124123,
            0.7347958438124123,
            0.7347958438124123,
            0.7347958438124123,
            None],
          'xhi_R2': [-1.2018078698570764,
            -1.2018078698570764,
            -1.2018078698570764,
            -1.2018078698570764,
            -1.2018078698570764,
            -1.2018078698570764,
            -1.2018078698570764,
            -1.2018078698570764,
            None],
          'pitch_R': [0.011977859609954216,
            0.012465295132576906,
            0.013037374721346386,
            0.013702773629594827,
            0.014470929140613027,
            0.015352272422314595,
            0.016358241098895984,
            0.017502541767555512,
            None],
          'o_R': [0.0027485430110031138,
            0.0028603941716149692,
            0.002991668489953377,
            0.003144356664497973,
            0.0033206242556974887,
            0.003522864889341437,
            0.0037537031413616907,
            0.004016284245811816,
            None],
          't_TE_R': [0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            0.0005,
            None],
          't_blade_R': [0.0013752412553972846,
            0.0014312062993939408,
            0.0014968897751956338,
            0.001573287734406901,
            0.0016614837213121984,
            0.001762675393333128,
            0.0018781759644470536,
            0.00200955916139248,
            None],
          'n_blade_R': [127, 122, 117, 111, 105, 99, 93, 87, None],
          'R_c_R': [0.022001749326772676,
            0.02289710413397014,
            0.023947938933921296,
            0.02517018902353028,
            0.02658118944829365,
            0.028200100888758198,
            0.0300479327530034,
            0.032149862253444544,
            None]},
          'CAPEX': {'Turbine': 832340.4268232549,
          'Alternator': 475326.0713036128,
          'Installation': 457683.27434440365,
          'Total': 1765349.7724712715}}
        
        
        # turb_params = {'type': 'Axial Turbine',
        #  'mdot_rated': 318.437021666738,
        #  'Wdot_rated': 16525640.389554193,
        #  'N_rot_rated': 2139.2697278709124,
        #  'total_to_static_efficiency': 0.9238504768389768,
        #  'DP_rated': 2.93,
        #  'n_stages': 17,
        #  'p0_su': 15309670.5,
        #  'T0_su': 406.4,
        #  'p_ex': 5220928,
        #  'r_m': 0.2034514959316595,
        #  'delta_tip': 0.0004,
        #  'N_lw': 0,
        #  'D_lw': 0,
        #  'e_blade': 2e-06,
        #  'stator': {'h_blade_S': [0.043675818843015114,
        #    0.04555075447331882,
        #    0.04756421118011338,
        #    0.04972691902234203,
        #    0.05205017510501014,
        #    0.054546596942447775,
        #    0.05722928920039919,
        #    0.060112617286739035,
        #    0.06321209857721537,
        #    0.06654490533621829,
        #    0.07012889817037861,
        #    0.07398381132317491,
        #    0.07813111368635203,
        #    0.08259363826885417,
        #    0.08739739527936775,
        #    0.09257033527958411,
        #    0.09814301040318812,
        #    0.10106509645143522],
        #   'chord_S': [0.014565147357473614,
        #    0.01480943096399099,
        #    0.015075739833271963,
        #    0.01536524874907234,
        #    0.01567912889606207,
        #    0.016018669329296093,
        #    0.0163851244995266,
        #    0.016779838143617463,
        #    0.01720421402838227,
        #    0.017659790605438676,
        #    0.018148061183290848,
        #    0.01867066392643876,
        #    0.019229341748854812,
        #    0.01982586410819994,
        #    0.020462313888188075,
        #    0.021140853881350716,
        #    0.021863805194343753,
        #    0.022242295668683392],
        #   'xhi_S1': [-0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995,
        #    -0.45212668370301995],
        #   'xhi_S2': [1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009,
        #    1.1815777847522009],
        #   'pitch_S': [0.014340200905962212,
        #    0.014580711757621719,
        #    0.01484290771038489,
        #    0.015127945404460429,
        #    0.015436977936556033,
        #    0.015771274453355428,
        #    0.016132070031673943,
        #    0.01652068765548968,
        #    0.016938509411618372,
        #    0.01738704999158629,
        #    0.017867779641003646,
        #    0.01838231122429742,
        #    0.018932360737599515,
        #    0.019519670310785425,
        #    0.020146290659181845,
        #    0.02081435117280855,
        #    0.021526137110781208,
        #    0.02189878211805869],
        #   'o_S': [0.005208284333331462,
        #    0.005295636589336021,
        #    0.005390864758173709,
        #    0.0054943889253871074,
        #    0.005606627889538167,
        #    0.005728042597919949,
        #    0.005859081623831385,
        #    0.0060002254679831415,
        #    0.006151976096920658,
        #    0.006314883638511587,
        #    0.006489482078104392,
        #    0.0066763571994396175,
        #    0.006876132243141885,
        #    0.007089439941471828,
        #    0.00731702510327684,
        #    0.007559661111632391,
        #    0.0078181779604366,
        #    0.007953520635621198],
        #   't_TE_S': [0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005],
        #   't_blade_S': [0.002873334793915352,
        #    0.002929701684898296,
        #    0.0030135521719042804,
        #    0.003069167448835821,
        #    0.003154388975412103,
        #    0.003241018376036457,
        #    0.003329098410875826,
        #    0.003418696869414256,
        #    0.0035099092960128643,
        #    0.0035631878913450514,
        #    0.0036558599177864298,
        #    0.0037504780951033443,
        #    0.0038472621962256357,
        #    0.00399581622006547,
        #    0.004100639839190145,
        #    0.004208772389728993,
        #    0.004320647254195877,
        #    0.004377395538210385],
        #   'n_blade_S': [89,
        #    88,
        #    86,
        #    85,
        #    83,
        #    81,
        #    79,
        #    77,
        #    75,
        #    74,
        #    72,
        #    70,
        #    68,
        #    65,
        #    63,
        #    61,
        #    59,
        #    58],
        #   'R_c_S': [0.021374275566569574,
        #    0.021732760447906438,
        #    0.022123567284125116,
        #    0.02254842006408143,
        #    0.02300903749499203,
        #    0.02350731125823086,
        #    0.024045082253543307,
        #    0.024624322407564403,
        #    0.02524709174055817,
        #    0.02591564792201565,
        #    0.026632182374151173,
        #    0.027399099094575358,
        #    0.028218955800187365,
        #    0.029094349160606794,
        #    0.030028335796503665,
        #    0.03102408959430885,
        #    0.032085016765580995,
        #    0.03264044950507237]},
        #  'rotor': {'h_blade_R': [0.04463901386550254,
        #    0.04658356426213087,
        #    0.0486717482309841,
        #    0.050914854071433466,
        #    0.053324844149276726,
        #    0.055914350111055416,
        #    0.05869720192740975,
        #    0.06168834199415402,
        #    0.06490389260233821,
        #    0.06836171319331613,
        #    0.07208042173763626,
        #    0.07608058437140806,
        #    0.08038461889099627,
        #    0.08501639257628242,
        #    0.09000310695231376,
        #    0.09537399662687153,
        #    0.10116101737079411,
        #    None],
        #   'chord_R': [0.014691692713363517,
        #    0.014947179258483243,
        #    0.015225236215973392,
        #    0.015527083786834149,
        #    0.015853953554919285,
        #    0.01620707828550691,
        #    0.016587774192778914,
        #    0.016997416201988986,
        #    0.01743743839432239,
        #    0.017909417206940486,
        #    0.01841488909821388,
        #    0.018955540396719053,
        #    0.019533173914154756,
        #    0.02014962496120879,
        #    0.020807057041219323,
        #    0.021507713925394085,
        #    0.02225399673609491,
        #    None],
        #   'xhi_R1': [0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    0.44794837042429314,
        #    None],
        #   'xhi_R2': [-1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    -1.2103795569237705,
        #    None],
        #   'pitch_R': [0.014515686298242616,
        #    0.014768112115657546,
        #    0.015042837952000385,
        #    0.015341069390268617,
        #    0.015664023259946775,
        #    0.016012917557790274,
        #    0.0163890527297408,
        #    0.016793787229453454,
        #    0.017228537957827104,
        #    0.017694862466312748,
        #    0.018194278806516355,
        #    0.018728453104805354,
        #    0.01929916657520216,
        #    0.019908232541380466,
        #    0.02055778858791782,
        #    0.02125005160565259,
        #    0.021987393951511107,
        #    None],
        #   'o_R': [0.005304932398056339,
        #    0.005397184453480851,
        #    0.00549758631942464,
        #    0.005606578590715294,
        #    0.005724605972344555,
        #    0.005852113595897372,
        #    0.005989576725006908,
        #    0.0061374918229238915,
        #    0.006296376713148282,
        #    0.006466800621618437,
        #    0.006649318338578578,
        #    0.00684453876997805,
        #    0.007053112881935934,
        #    0.007275703375430132,
        #    0.007513091456491554,
        #    0.007766087314579369,
        #    0.008035557956107168,
        #    None],
        #   't_TE_R': [0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    0.0005,
        #    None],
        #   't_blade_R': [0.0029339168824975563,
        #    0.0029906358249125425,
        #    0.003075651966960518,
        #    0.0031621491077244,
        #    0.003217899079117751,
        #    0.0033057153990275227,
        #    0.0033950293310365115,
        #    0.0034859234193912745,
        #    0.003578509603973971,
        #    0.003672929300044586,
        #    0.0037693633336161607,
        #    0.0038680292077971096,
        #    0.003969185344541105,
        #    0.004073139878085224,
        #    0.004180239496230803,
        #    0.004290889279250552,
        #    0.004405559116686325,
        #    None],
        #   'n_blade_R': [88,
        #    87,
        #    85,
        #    83,
        #    82,
        #    80,
        #    78,
        #    76,
        #    74,
        #    72,
        #    70,
        #    68,
        #    66,
        #    64,
        #    62,
        #    60,
        #    58,
        #    None],
        #   'R_c_R': [0.021559980197087613,
        #    0.021934905330690042,
        #    0.02234295242329613,
        #    0.022785912113324047,
        #    0.023265591743468497,
        #    0.023783800390158066,
        #    0.024342469590639706,
        #    0.024943617040343095,
        #    0.02558934724570997,
        #    0.026281973619812878,
        #    0.02702375090706352,
        #    0.02781715378559627,
        #    0.028664827871902193,
        #    0.029569466474564345,
        #    0.03053424452311857,
        #    0.03156245473976845,
        #    0.03265762075869161,
        #    None]},
        #  'CAPEX': {'Turbine': 871078.6076124497,
        #   'Alternator': 497049.7207146683,
        #   'Installation': 478844.91491449124,
        #   'Total': 1846973.243241609}}
        
        Turb_OD.set_inputs(
            m_dot = turb_params['mdot_rated'], # kg/s
            P_su = turb_params['p0_su'], # Pa 
            T_su = turb_params['T0_su'], # K
            fluid = 'CO2',
            N_rot = turb_params['N_rot_rated'], # 1985.8640516476623, # RPM
            P_ex = turb_params['p_ex'], # 5742510, # Pa
            )

        Turb_OD.set_parameters(
            r_m = turb_params['r_m'],
            nStages = turb_params['n_stages'],
            mdot_rated = turb_params['mdot_rated'],
            DP_rated = turb_params['DP_rated'],
            N_rot_rated = turb_params['N_rot_rated'], # RPM 
            damping = 0.3,
            delta_tip = 0.0004,
            N_lw = 0,
            D_lw = 0,
            e_blade = 2e-06
            )
        
        Turb_OD.set_stage_parameters(
            # --- Stator ---
            h_blade_S  = turb_params['stator']['h_blade_S'],
            chord_S    = turb_params['stator']['chord_S'],
            xhi_S1     = turb_params['stator']['xhi_S1'],
            xhi_S2     = turb_params['stator']['xhi_S2'],
            pitch_S    = turb_params['stator']['pitch_S'],
            o_S        = turb_params['stator']['o_S'],
            t_TE_S     = turb_params['stator']['t_TE_S'],
            t_blade_S  = turb_params['stator']['t_blade_S'],
            n_blade_S  = turb_params['stator']['n_blade_S'],
            R_c_S      = turb_params['stator']['R_c_S'],
            # --- Rotor ---
            h_blade_R  = turb_params['rotor']['h_blade_R'],
            chord_R    = turb_params['rotor']['chord_R'],
            xhi_R1     = turb_params['rotor']['xhi_R1'],
            xhi_R2     = turb_params['rotor']['xhi_R2'],
            pitch_R    = turb_params['rotor']['pitch_R'],
            o_R        = turb_params['rotor']['o_R'],
            t_TE_R     = turb_params['rotor']['t_TE_R'],
            t_blade_R  = turb_params['rotor']['t_blade_R'],
            n_blade_R  = turb_params['rotor']['n_blade_R'],
            R_c_R      = turb_params['rotor']['R_c_R'],
        )
        
        Turb_OD.solve()
        
    df_map = generate_map_processes(
        Turb_OD,
        m_grid=np.linspace(0.6*Turb_OD.params['mdot_rated'], 1.8*Turb_OD.params['mdot_rated'], 50),
        N_grid=np.linspace(0.3*Turb_OD.params['N_rot_rated'], 1.8*Turb_OD.params['N_rot_rated'], 50),
        max_workers=-2
    )
    
    df_clean = filter_sparse_by_proximity(df_map, rp_col=None, group_by='N_rot',
                                          rp_tol_rel=0.6, m_tol_rel=0.6, min_neighbors=2)
    
    fig, ax = map_plot(
        df_clean, rp_col='RP_calc',
        use_grid=True, nx=600, ny=600,
        # triangulation cleaning
        min_circle_ratio=0.01,   # less aggressive
        max_area_factor=50.0,     # allow larger cells before masking
        long_edge_q=1,         # drop triangles that bridge big gaps
        # hole filling + smoothing
        fill_holes=True, hole_method='nearest', hole_smooth_sigma=0.8,
        smooth_sigma=0.6,         # gentle overall blur
        # cosmetics
        show_points=True,        # hide all raw dots (see patch below for "used" only)
        levels=24, focus_high=True, max_iso_speeds=4,
        figsize=(9,6), dpi=220
    )
    plt.show()

    # fig, ax = map_plot_clean(
    #     df_clean,
    #     rp_col='RP_calc',          # or 'RP_calc' if that’s your chosen column
    #     levels=24,
    #     nx=500, ny=500,
    #     min_circle_ratio=0.01,
    #     long_edge_q=0.92,     # tighten to 0.90 if you still see bridging
    #     smooth_sigma=0.6      # small smoothing of the gridded field
    # )
    # plt.show()
    
    _ = plot_power_eta_vs_mdot(df_map, speeds=None, max_lines=5)  # auto-picks up to 5 speeds
    plt.show()
    