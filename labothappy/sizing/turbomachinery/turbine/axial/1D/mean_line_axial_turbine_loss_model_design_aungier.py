#!/usr/bin/python3

# --- loading libraries 

from connector.mass_connector import MassConnector
from correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model

from component.expander.turbine_mean_line_Aungier import AxialTurbineMeanLine
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq, root
import pyswarms as ps

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

#%%

# ---- joblib worker (process-based) ----
import os, numpy as np
from joblib import Parallel, delayed

from contextlib import contextmanager
from tqdm import tqdm
import joblib
from joblib.parallel import BatchCompletionCallBack

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

_SOLVER = None  # cached per-process

def _eval_particle(x, cls, fluid, params, stage_params, inputs, match_p_ex=False):
    """Evaluate one particle using a per-process cached solver."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")

    global _SOLVER
    if _SOLVER is None:
        s = cls(fluid)
        s.set_parameters(**params)
        if stage_params:
            s.set_stage_parameters(**stage_params)
        s.set_inputs(**inputs)
        _SOLVER = s

    x = np.asarray(x, dtype=float)
    if match_p_ex:
        # x here is WITHOUT r_m: [psi, phi, R, Re_min, M_1_st]
        return float(_SOLVER.design_system_match_p(x))
    else:
        # x includes r_m (and M_1_st if you keep it)
        return float(_SOLVER.design_system(x))

#%%

class AxialTurbineMeanLineDesign(object):

    def __init__(self, fluid):
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
        self.Vel_Tri = {}
        self.Vel_Tri_Last_Stage = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None
        self.allowable_positions = []

    def reset(self):

        self.AS = CP.AbstractState('HEOS', self.fluid)
        
        # Blade Dictionnary
        self.stages = []
        
        # Velocity Triangle Data
        self.Vel_Tri = {}
        self.Vel_Tri_Last_Stage = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None
        
        return
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
            
        def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            ts = self.total_states  # shorthand
            ts.loc[position, 'H'] = self.AS.hmass()
            ts.loc[position, 'S'] = self.AS.smass()
            ts.loc[position, 'P'] = self.AS.p()
            ts.loc[position, 'D'] = self.AS.rhomass()
            try:
                ts.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                ts.loc[position, 'A'] = -1
            ts.loc[position, 'V'] = self.AS.viscosity()
            return
        
        def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
            self.AS.update(CP_INPUTS, input_1, input_2)
            ss = self.static_states
            ss.loc[position, 'H'] = self.AS.hmass()
            ss.loc[position, 'S'] = self.AS.smass()
            ss.loc[position, 'P'] = self.AS.p()
            ss.loc[position, 'D'] = self.AS.rhomass()
            try:
                ss.loc[position, 'A'] = self.AS.speed_sound()
            except Exception:
                ss.loc[position, 'A'] = -1
            ss.loc[position, 'V'] = self.AS.viscosity()
            return

    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value
    
    # ---------------- Result Plot Methods ----------------------------------------------------------------

    def plot_geometry(self, fontsize = 16, ticksize = 12):
        
        r_m_line = np.ones(len(self.r_tip))*self.r_m
        
        x = np.linspace(0,len(self.r_tip)-1, len(self.r_tip))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1
        
        plt.figure()
        plt.plot(self.r_tip)
        plt.plot(self.r_hub)
        plt.plot(r_m_line)
                
        plt.axis([-0.5, len(self.r_tip)-0.5, 0, max(self.r_tip)*1.2])
        plt.legend(["$r_{tip}$", "$r_{hub}$", "$r_{m}$"])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.grid()
        plt.ylabel("Length or radius [m]", fontsize= fontsize)
        plt.show()

    def plot_n_blade(self, fontsize = 16, ticksize = 12):
        n_blade_plot = np.array(self.n_blade)

        x = np.linspace(0,len(n_blade_plot)-1, len(n_blade_plot))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1     

        plt.figure()
        plt.plot(x[::2], n_blade_plot[::2], 'o', label="Stator Blades")  # even indices
        plt.plot(x[1::2], n_blade_plot[1::2], 'o', label="Rotor Blades")  # odd indices
        plt.axis([-0.5, len(self.r_tip)-0.5, 0, max(n_blade_plot.flatten())*1.2])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.legend()
        plt.grid()
        plt.ylabel("Blade number [-]", fontsize= fontsize)
        plt.show()

    def plot_radius_verif(self, fontsize = 16, ticksize = 12):
        
        x = np.linspace(0,len(self.r_ratio2)-1, len(self.r_ratio2))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1
            
        plt.figure()
        plt.plot(self.r_ratio2)
        plt.axis([-0.5, len(self.r_ratio2)-0.5, 0, max(self.r_ratio2)*1.2])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.grid()
        plt.ylabel("$\\left[ r_{ext}/r_{hub} \\right]^2$ [-]", fontsize= fontsize)
        plt.show()

        plt.figure()
        plt.plot(self.r_hub_tip)
        plt.axis([-0.5, len(self.r_hub_tip)-0.5, 0, 1])
        plt.xticks(ticks=x, labels=labels, size=ticksize)
        plt.grid()
        plt.ylabel("$\\left[ r_{hub}/r_{tip} \\right]$ [-]", fontsize= fontsize)
        plt.show()

    def plot_Mollier(self, fontsize = 16, ticksize = 12):
        # Thermo Prop
        x = np.linspace(0,len(self.r_tip)-1, len(self.r_tip))
        
        labels = []
        i = 1
        
        while len(labels) < len(x):
            labels.append("S" + str(i))
            labels.append("R" + str(i))
            i += 1
        
        x2 = np.linspace(0,len(self.r_tip), len(self.r_tip)+1)
        labels2 = ['0'] + labels
        
        p = [self.stages[0].static_states['P'][1]]
        s = [self.stages[0].static_states['S'][1]]
        h = [self.stages[0].static_states['H'][1]]
        
        for i in range(self.nStages):
            p.append(self.stages[i].static_states['P'][2])
            p.append(self.stages[i].static_states['P'][3])
        
            s.append(self.stages[i].static_states['S'][2])
            s.append(self.stages[i].static_states['S'][3])
            
            h.append(self.stages[i].static_states['H'][2])
            h.append(self.stages[i].static_states['H'][3])
        
        plt.figure()
        plt.plot(np.array(p)*1e-3)
        plt.axis([-0.5, len(self.r_tip)+0.5, 0, max(np.array(p)*1e-3)*1.2])
        plt.xticks(ticks=x2, labels=labels2, size=ticksize)
        plt.grid()
        plt.ylabel("Oulet Pressure [kPa]", fontsize= fontsize)
        plt.show()
        
        plt.figure()
        plt.plot(s, h)
        plt.plot([s[0], s[0]], [h[0], h[-1]])
        
        # Define entropy range (in J/(kg·K))
        entropy_range = np.linspace(s[0], s[-1], 100)  # Adjust range for your fluid
        
        for P in p:
            enthalpy = [PropsSI('H', 'S', s, 'P', P, self.fluid) for s in entropy_range]  # Enthalpy in kJ/kg
            entropy = entropy_range  # Entropy in kJ/(kg·K)
            plt.plot(entropy, enthalpy, color = 'grey', alpha=0.3, label=f'P = {P/1e5} bar')
        
        plt.ylabel("$Enthalpy$ [J/kg]", fontsize= fontsize)
        plt.xlabel("$Entropy$ [J/(kg x K)]", fontsize= fontsize)
        plt.legend(["real", "isentropic"])
        plt.show()

    # ---------------- Export Geom ------------------------------------------------------------------------
    
    def export_params(self, turb_var_name: str = "Turb_OD"):
        """
        Print a ready-to-copy call to set_parameters() with selected parameters,
        mixing values from self.params, self, and self.inputs as specified.
        """
        # --- keys to pull from self.params (fallback to attributes if missing) ---
        param_keys = ["damping", "delta_tip", "N_lw", "D_lw", "e_blade", "Omega_rated"]
    
        selected = {}
    
        # 1) Derived/off-design style entries from attributes & inputs
        r_m = getattr(self, "r_m", None)
        nStages = getattr(self, "nStages", None)
    
        mdot_rated = None
        if hasattr(self, "inputs"):
            mdot_rated = self.inputs.get("mdot", self.inputs.get("m_dot"))
    
        DP_rated = None
        if hasattr(self, "inputs"):
            # Robust name fallbacks
            p0_su = (self.inputs.get("p0_su", None) or
                     self.inputs.get("P0_su", None) or
                     self.inputs.get("P_su", None))
            p_ex  = (self.inputs.get("p_ex", None) or
                     self.inputs.get("P_ex", None))
            if p0_su is not None and p_ex not in (None, 0):
                DP_rated = p0_su / p_ex
    
        # Insert in a sensible order
        if r_m is not None:        selected["r_m"] = r_m
        if nStages is not None:    selected["nStages"] = nStages
        if mdot_rated is not None: selected["mdot_rated"] = mdot_rated
        if DP_rated is not None:   selected["DP_rated"] = DP_rated
    
        # 2) Design params from self.params (or attributes)
        for k in param_keys:
            if hasattr(self, "params") and k in self.params:
                selected[k] = self.params[k]
            elif hasattr(self, k):
                selected[k] = getattr(self, k)
    
        # --- formatting helpers ---
        def fmt(v):
            if isinstance(v, float):
                return f"{v:.12g}"
            if isinstance(v, (list, tuple)):
                inner = ", ".join(fmt(x) for x in v)
                return f"[{inner}]"
            return repr(v)
    
        args_str = ",\n    ".join(f"{k} = {fmt(v)}" for k, v in selected.items())
    
        print(f"""{turb_var_name}.set_parameters(
        {args_str}
        )""")


    def export_stage_vectors(self):
        """
        Export stage geometry as Python vectors (lists), one per variable.
        Stator and rotor variables are grouped separately.
        """
        # Define stator and rotor attribute lists
        stator_vars = [
            "h_blade_S", "chord_S", "xhi_S1", "xhi_S2",
            "pitch_S", "o_S", "t_TE_S", "t_blade_S", "n_blade_S", "R_c_S"
        ]
        rotor_vars = [
            "h_blade_R", "chord_R", "xhi_R1", "xhi_R2",
            "pitch_R", "o_R", "t_TE_R", "t_blade_R", "n_blade_R", "R_c_R"
        ]
    
        # Helper for formatting values nicely
        def fmt(val):
            if isinstance(val, float):
                return f"{val:.10g}"
            return repr(val)
    
        print("# --- Stage geometry vectors ---")
    
        # Collect and print stator arrays
        for var in stator_vars:
            values = []
            for stage in self.stages:
                values.append(getattr(stage, var, None))
            if any(v is not None for v in values):  # only print if something exists
                arr = [fmt(v) for v in values]
                print(f"{var} = [{', '.join(arr)}],")
    
        print()
    
        # Collect and print rotor arrays
        for var in rotor_vars:
            values = []
            for stage in self.stages:
                values.append(getattr(stage, var, None))
            if any(v is not None for v in values):
                arr = [fmt(v) for v in values]
                print(f"{var} = [{', '.join(arr)}],")
                
    # ---------------- Loss Models ------------------------------------------------------------------------
    def compute_stator_t_max(self, stage):
        # --- constants from stage / design ---
        c   = stage.chord_S
        h   = stage.h_blade_S
        vm  = self.Vel_Tri['vm']
        a1  = self.Vel_Tri['alpha1']     # radians
        a2  = self.Vel_Tri['alpha2']     # radians
        mdot = self.inputs['mdot']
        Nbl  = stage.n_blade_S
    
        # allowable stress (what you had: 2 * 130 MPa)
        sigma_allow = 2 * 130e6
    
        # camber (deg) for B,n interpolation (same as your code)
        xhi = abs(180.0 * (self.Vel_Tri["alpha1"] - self.Vel_Tri["alpha2"]) / np.pi)
    
        # ----- interpolate B and n ONCE -----
        B_values  = np.array([951.083, 904.889, 829.035, 753.251, 665.223, 568.022, 461.654, 359.362, 273.379])
        n_values  = np.array([1.694,   1.605,   1.520,   1.436,   1.347,   1.268,   1.179,   1.094,   1.000])
        xhi_grid  = np.array([40, 50, 60, 70, 80, 90, 100, 110, 120])
    
        # clamp to table range to avoid extrapolation surprises
        xhi_clamped = np.clip(xhi, xhi_grid.min(), xhi_grid.max())
        B = np.interp(xhi_clamped, xhi_grid, B_values)
        n = np.interp(xhi_clamped, xhi_grid, n_values)
    
        # ----- closed-form thickness -----
        K = (h/2.0) * mdot * vm * (np.tan(a1) + np.tan(a2)) / Nbl
    
        # t = ((K*B*c^(n-3)) / (10^n * sigma_allow))^(1/n)
        t = ((K * B * (c**(n - 3.0))) / ( (10.0**n) * sigma_allow ))**(1.0/n)
    
        # optional: sanity bounds
        t_min = 0.1*c   # 0.3 mm
        t_max = 0.35*c   # e.g. 5% chord
        t = float(np.clip(t, t_min, t_max))
    
        # populate stage + return
        stage.t_blade_S    = t
        # compute resulting sigma for reporting
        z = (1.0 / B) * (10.0 * t / c)**n
        stage.sigma_bmax_S = (1.0 / (z * c**3)) * (h/2.0) * mdot * vm * (np.tan(a1) + np.tan(a2)) / Nbl
        return

    def compute_rotor_t_max(self, stage):
        """
        Fast, closed-form rotor thickness from Saravanamuttoo z-correlation.
        Uses rotor metal angles for camber and rotor (relative) angles in load term.
        """
    
        # --- inputs/geometry ---
        c   = stage.chord_R
        h   = stage.h_blade_R
        Nbl = stage.n_blade_R
        mdot = self.inputs['mdot']
        vm   = self.Vel_Tri['vm']
    
        # ***** IMPORTANT *****
        # Use rotor METAL angles for camber (replace with your stored metal angles if named differently)
        # If you only have them in radians:
        beta2_metal = self.Vel_Tri['beta2']   # radians
        beta3_metal = self.Vel_Tri['beta3']   # radians
        xhi_deg = abs(np.degrees(beta2_metal - beta3_metal))  # rotor blade camber in degrees
    
        # use RELATIVE flow angles in the load term. If your stored ones are the same, reuse:
        # guard against tan() singularities at ±90°
        eps = 1e-6  # rad
        def safe_tan(theta):
            # wrap to (-pi/2, pi/2) and clip
            th = ((theta + np.pi/2) % np.pi) - np.pi/2
            th = np.clip(th, -np.pi/2 + eps, np.pi/2 - eps)
            return np.tan(th)
    
        tan_sum = safe_tan(beta2_metal) + safe_tan(beta3_metal)
    
        # Allowable stress (your previous 130 MPa with SF=2)
        sigma_allow = 130e6 * 2
    
        # --- interpolate B, n once (clamped to table) ---
        B_values = np.array([951.083, 904.889, 829.035, 753.251, 665.223, 568.022, 461.654, 359.362, 273.379])
        n_values = np.array([  1.694,   1.605,   1.520,   1.436,   1.347,   1.268,   1.179,   1.094,   1.000])
        xhi_grid = np.array([40, 50, 60, 70, 80, 90, 100, 110, 120], dtype=float)
    
        xhi_clamped = float(np.clip(xhi_deg, xhi_grid[0], xhi_grid[-1]))
        B = np.interp(xhi_clamped, xhi_grid, B_values)
        n = np.interp(xhi_clamped, xhi_grid, n_values)
    
        # --- closed-form thickness from sigma(t) = K*B*c^(n-3) / (10^n * t^n) ---
        # K = (h/2) * mdot * vm * (tan(beta2) + tan(beta3)) / Nbl
        K = abs((h * 0.5) * mdot * vm * tan_sum / Nbl)
    
        # t = ((K * B * c^(n-3)) / (10^n * sigma_allow))^(1/n)
        t = ((K * B * (c**(n - 3.0))) / ((10.0**n) * sigma_allow))**(1.0 / n)
    
        # optional: sanity bounds
        t_min = 0.1*c   # 0.3 mm
        t_max = 0.35*c   # e.g. 5% chord
        t = float(np.clip(t, t_min, t_max))
    
        # write back + compute resulting sigma for reporting
        stage.t_blade_R = t
                
        z = (10.0 * t / c)**n / B
        stage.sigma_bmax_R = (1.0 / (z * c**3)) * (h * 0.5) * mdot * vm * tan_sum / Nbl
        return

    def stator_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]

        # 1) Guess outlet state
        h_static_out = x[0]*1e5
        p_static_out = x[1]*1e5
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu2']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu2']**2)
        stage.M_S = max(v,w)/stage.static_states['A'][2]
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(2*np.pi*self.r_m)
        
        if stage.h_blade_S > self.r_m*0.8:
            raise ValueError('')
        
        # 4) Compute cord, aspect ratio, blade pitch and blade number
                
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        stage.pitch_S = stage.chord_S/self.solidityStator
        stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
        
        self.compute_stator_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha2'] - self.Vel_Tri['alpha1']
        
        stage.R_c_S = 2*np.sin(abs(camber)/2)*stage.chord_S # Geometrical estimation of blade curvature radius
                
        stage.M1_S = self.Vel_Tri['v1']/stage.static_states['A'][1]
        stage.M2_S = self.Vel_Tri['v2']/stage.static_states['A'][2]
        
        stage.beta_g_S = 0.5*(self.Vel_Tri['alpha1'] + self.Vel_Tri['alpha2'])
        stage.o_S = np.sin(stage.beta_g_S)*stage.pitch_S

        stage.t_TE_S = max(self.params['t_TE_o']*stage.pitch_S * np.cos(self.Vel_Tri['alpha2'])/(1+self.params['t_TE_o']), self.params['t_TE_min'])
        
        stage.Y_vec_S = aungier_loss_model(self.Vel_Tri['alpha1'], self.Vel_Tri['alpha2'], stage.beta_g_S*180/np.pi, self.Vel_Tri['alpha1'], stage.chord_S, 
                               0, self.params['D_lw'], self.params['e_blade'], stage.h_blade_S, stage.static_states['V'][2], 
                               stage.M1_S, stage.M2_S, self.params['N_lw'], stage.R_c_S, stage.static_states['D'][2], stage.pitch_S, stage.t_blade_S, stage.t_TE_S,
                               self.Vel_Tri['vm'], self.Vel_Tri['v2'],1)
                
        # stage.Y_vec_S = aungier_loss_model(self.Vel_Tri['beta1'], self.Vel_Tri['beta2'], self.Vel_Tri['beta1'], stage.chord_S, stage.delta_tip, stage.e_blade)    
        Y = stage.Y_vec_S['Y_tot']
                
        # p0_out = stage.total_states['P'][1] - (Y*(stage.total_states['P'][1] - p_static_out))
        p0_out = (stage.total_states['P'][1] + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)/2
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
        
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu3']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)
        stage.M_R = max(v,w)/stage.static_states['A'][3]
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][2]
        h0in = hin + (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][2], 2)            
        
        # 3) Compute A_flow and h_blade based on r_m 
        stage.A_flow_R = self.inputs['mdot']/(stage.static_states['D'][3]*self.Vel_Tri['vm'])
        stage.h_blade_R = stage.A_flow_R/(2*np.pi*self.r_m)
        
        if stage.h_blade_R > self.r_m*0.8:
            raise ValueError('')
        
        # 4) Compute cord, aspect ratio, pitch and blade number
        stage.chord_R = (self.params['Re_min']*stage.static_states['V'][3])/(stage.static_states['D'][3]*self.Vel_Tri['vm'])            
        stage.AR_R = stage.h_blade_R/stage.chord_R
        
        stage.pitch_R = stage.chord_R/self.solidityRotor
        stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
        
        self.compute_rotor_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha3'] - self.Vel_Tri['alpha2']
        
        stage.R_c_R = 2*np.sin(abs(camber)/2)*stage.chord_R # Geometrical estimation of blade curvature radius
                
        stage.M2_R = self.Vel_Tri['w2']/stage.static_states['A'][2]
        stage.M3_R = self.Vel_Tri['w3']/stage.static_states['A'][3]
        
        stage.beta_g_R = 0.5*(self.Vel_Tri['beta2'] + self.Vel_Tri['beta3'])  # mid-passage metal angle
        stage.o_R = abs(np.sin(stage.beta_g_R)*stage.pitch_R)
        
        
        stage.t_TE_R = max(self.params['t_TE_o']*stage.pitch_R * np.cos(self.Vel_Tri['alpha2'])/(1+self.params['t_TE_o']), self.params['t_TE_min'])
        
        stage.Y_vec_R = aungier_loss_model(-self.Vel_Tri['beta2'], -self.Vel_Tri['beta3'], stage.beta_g_R*180/np.pi, -self.Vel_Tri['beta2'], stage.chord_R, 
                               self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.static_states['V'][3], 
                               stage.M2_R, stage.M3_R, self.params['N_lw'], stage.R_c_R, stage.static_states['D'][3], stage.pitch_R, stage.t_blade_R, stage.t_TE_R,
                               self.Vel_Tri['vm'], self.Vel_Tri['w3'],1)

        Y = stage.Y_vec_R['Y_tot']
                
        p0_out = (stage.total_states['P'][2] + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.total_states['S'][3]
        
        hout = h0in-(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)/2
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
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        v = np.sqrt(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu2']**2)
        w = np.sqrt(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['wu2']**2)
        stage.M_S = max(v,w)/stage.static_states['A'][2]
        
        # 2) Compute total inlet state
        hin = stage.static_states['H'][1]
        h0in = hin + (self.Vel_Tri_Last_Stage['vu1']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][1], 1)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_S = self.inputs['mdot']/(stage.static_states['D'][2]*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(2*np.pi*self.r_m)

        if stage.h_blade_S > self.r_m*0.8:
            raise ValueError('')

        # 4) Compute cord, aspect ratio, blade pitch and blade number
                
        stage.chord_S = (self.params['Re_min']*stage.static_states['V'][2])/(stage.static_states['D'][2]*self.Vel_Tri['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        stage.pitch_S = stage.chord_S/self.solidityStator
        stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
        
        self.compute_stator_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha2'] - self.Vel_Tri['alpha1']
        
        stage.R_c_S = 2*np.sin(abs(camber)/2)*stage.chord_S # Geometrical estimation of blade curvature radius
                
        stage.M1_S = self.Vel_Tri['v1']/stage.static_states['A'][1]
        stage.M2_S = self.Vel_Tri['v2']/stage.static_states['A'][2]
        
        stage.beta_g_S = 0.5*(self.Vel_Tri['alpha1'] + self.Vel_Tri['alpha2'])
        stage.o_S = np.sin(stage.beta_g_S)*stage.pitch_S

        stage.t_TE_S = max(self.params['t_TE_o']*stage.pitch_S * np.cos(self.Vel_Tri['alpha2'])/(1+self.params['t_TE_o']), self.params['t_TE_min'])
        
        # 5) Estimate pressure losses 
        # 5.1) Aungier : Profile pressure losses                     
        v_1 = np.sqrt(self.Vel_Tri_Last_Stage["vm"]**2 + self.Vel_Tri_Last_Stage["vu1"]**2)
        v_2 = np.sqrt(self.Vel_Tri_Last_Stage["vm"]**2 + self.Vel_Tri_Last_Stage["vu2"]**2)
        
        a = 0.0117 # NACA blade - 0.007 : C.4 circular-arc blade
        
        alpha = self.Vel_Tri_Last_Stage['alpha1']
        alpha_star = self.Vel_Tri_Last_Stage['alpha1']
        
        D_e = (np.cos(self.Vel_Tri_Last_Stage['alpha2'])/np.cos(self.Vel_Tri_Last_Stage['alpha1']))*(1.12+a*(alpha - alpha_star)+0.61*np.cos(self.Vel_Tri_Last_Stage['alpha1'])**2 / self.solidityStator * (np.tan(self.Vel_Tri_Last_Stage['alpha1'])-np.tan(self.Vel_Tri_Last_Stage['alpha2'])))
        
        P_cst = np.cos(self.Vel_Tri_Last_Stage["alpha2"])/2 * self.solidityStator * (v_1/v_2)**2 # Profile Constant
        
        Yp = 0.004*(1+3.1*(D_e - 1)**2 + 0.4*(D_e-1)**8)/P_cst
    
        # 5.2) Cohen : Endwall losses
        EW_Cst = np.cos((self.Vel_Tri_Last_Stage["alpha1"]+self.Vel_Tri_Last_Stage["alpha2"])/2)**3 / np.cos(self.Vel_Tri_Last_Stage["alpha1"])**2  # Endwall Constant

        Yew = 0.02*(self.solidityStator/stage.AR_S)/EW_Cst

        # Pressure loss 
        DP_loss = (Yp+Yew)*(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu1']**2)*stage.static_states['D'][1]/2
        p0_out = stage.total_states['P'][1]-DP_loss
                
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][2]
        
        hout = h0in-(self.Vel_Tri_Last_Stage['vu2']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2
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
        
            RP_1_row = (self.inputs['p0_su']/self.inputs['p_ex'])**(1/(2*self.nStages))
            h_out_guess = stage.static_states['H'][1] - self.Dh0Stage/2
            pout_guess = stage.static_states['P'][1]/RP_1_row
            # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])         
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
            
            res = 1
            x_in = x0_disc
            
            c = 0
            
            while res > 1e-8:
                
                if c > 100:
                    exit()
                
                # print(f"x_in : {x_in}")
                
                x_out = self.stator_blade_row_system(x_in)

                # print(f"x_out : {x_out}")
                
                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                              
                # print(f"new x_in : {x_in}")

                c += 1
                
            self.stator_blade_row_system(x_out)
            
            stage.xhi_S1 = self.Vel_Tri['alpha1']
            self.compute_deviation_stator(stage)
            stage.xhi_S2 = self.Vel_Tri['alpha2'] - stage.delta_S
            
            # print(f'Y_S : {stage.Y_vec_S}')

        else: # Rotor

            # print("Rotor")

            RP_1_row = (self.inputs['p0_su']/self.inputs['p_ex'])**(1/(2*self.nStages))
            h_out_guess = stage.static_states['H'][2] - self.Dh0Stage/2
            pout_guess = stage.static_states['P'][2]/RP_1_row
            # sol = minimize(self.rotor_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])    
            
            # Initial guess vector
            x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

            res = 1
            x_in = x0_disc
            
            c = 0
            
            while res > 1e-8:

                if c > 100:
                    exit()                

                # print(f"x_in : {x_in}")

                x_out = self.rotor_blade_row_system(x_in) 

                # print(f"x_out : {x_out}")

                res_vec = abs((x_in - x_out)/x_out)
                res = sum(res_vec)
                
                x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                       
                # print(f"new x_in : {x_in}")

                c += 1
            
            self.rotor_blade_row_system(x_out)
                        
            stage.xhi_R1 = self.Vel_Tri['beta2']
            self.compute_deviation_rotor(stage)
            stage.xhi_R2 = self.Vel_Tri['beta3'] - stage.delta_R
                        
            # print(f'Y_R : {stage.Y_vec_R}')

        return
            
    def computeRepeatingStages(self):
        
        for i in range(int(self.nStages)):
            if i == 0:
                self.computeBladeRow(i, 'S')
                self.computeBladeRow(i, 'R')

            else:
                self.stages[i].static_states.loc[1] = self.stages[i-1].static_states.loc[3]
                
                self.computeBladeRow(i, 'S')
                self.computeBladeRow(i, 'R')
            
            self.r_tip.append(self.r_m + self.stages[i].h_blade_S/2)
            self.r_hub.append(self.r_m - self.stages[i].h_blade_S/2)
            self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
            self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
        
            self.r_tip.append(self.r_m + self.stages[i].h_blade_R/2)
            self.r_hub.append(self.r_m - self.stages[i].h_blade_R/2)
            self.r_hub_tip.append(self.r_hub[-1]/self.r_tip[-1])
            self.r_ratio2.append((self.r_tip[-1]/self.r_hub[-1])**2)
            
        return
    
    def computeLastStage(self):
        stage = self.stages[-1]
        stage.static_states.loc[1] = self.stages[-2].static_states.loc[3]
        
        h_out_guess = self.stages[-2].total_states['H'][3]
        pout_guess = self.stages[-2].total_states['P'][3]
        # sol = minimize(self.last_blade_row_system, x0=(h_out_guess,pout_guess), bounds=[(self.stages[-1].static_states['H'][1], h_out_guess), (self.stages[-1].static_states['P'][1], pout_guess)])         
        
        # Initial guess vector
        x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5

        res = 1
        x_in = x0_disc
        
        c = 0
            
        while res > 1e-8:

            if c > 100:
                exit()                

            x_out = self.last_blade_row_system(x_in) 

            res_vec = abs((x_in - x_out)/x_out)
            res = sum(res_vec)
            
            x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out 
                        
            c += 1
        
        self.last_blade_row_system(x_out)
        
        stage.xhi_S1 = self.Vel_Tri['alpha1']
        self.compute_deviation_stator(stage)
        stage.xhi_S2 = self.Vel_Tri['alpha2'] - stage.delta_S
        
        return
    
    def cost_estimation(self):
        
        if self.fluid == 'CO2' or self.fluid == 'CarbonDioxide' or self.fluid == 'R744':
            """
            SCO2 POWER CYCLE COMPONENT COST CORRELATIONS FROM DOE DATA
            SPANNING MULTIPLE SCALES AND APPLICATIONS (2019)
            
            Nathan T. Weiland,  Blake W. Lance, Sandeep R. Pidaparti
            
            Especially good for 10 - 750 
            Based on 2017 CEPCI (chemical plant cost index) for dollars
            """
            
            T_500 = 273.15+500
            
            if self.inputs['T0_su'] > T_500: # K  
                f = 1 + 1.106*1e-4*(self.inputs['T0_su'] - T_500)**2
            else:
                f = 1
        
            W_dot_MW = self.W_dot/1e6
            self.CAPEX_turb = 182600*W_dot_MW**0.5561 * f
        
        else:
            """
            Supplementary Information:
            Techno-economic analysis of recuperated Joule-Brayton
            pumped thermal energy storage (2022)
                        
            Joshua D. McTigue, Pau Farres-Antunez, Kavin Sundarnath Jawaharlal Ayyanathanc, Christos N. Markides, Alexander J. White           
            Based on 1995 CEPCI (chemical plant cost index) for dollars
            """
            rho0_air = PropsSI('D', 'T', 15+273.15, 'P', 101325, 'Air')
            rho_out = self.stages[-1].static_states['D'][2]
            
            W_dot_MW = self.W_dot/1e6
            
            n = 1
            
            if self.eta_is >= 0.92:
                 fact_1 = 1.051 * 266 * self.inputs['mdot'] * np.log(self.inputs['p0_su']/self.inputs['p_ex']) / (0.94 - self.eta_is)
            else:
                 fact_1 = 1.051 * 266 * self.inputs['mdot'] * np.log(self.inputs['p0_su']/self.inputs['p_ex']) / (0.92 - self.eta_is)                
                
            fact_2 = (1+np.exp(0.036*(self.inputs['T0_su']) - 1.207 * 54.4))
            
            fact_3 = (rho_out/rho0_air)**(-n)
            
            self.CAPEX_turb = fact_1*fact_2*fact_3
            
        # Generator Costs
        
        """
        SCO2 POWER CYCLE COMPONENT COST CORRELATIONS FROM DOE DATA
        SPANNING MULTIPLE SCALES AND APPLICATIONS (2019)
        
        Nathan T. Weiland,  Blake W. Lance, Sandeep R. Pidaparti
        
        Especially good for 10 - 750 
        Based on 2017 CEPCI (chemical plant cost index) for dollars
        """
        
        self.eta_alt = 0.97
        self.W_dot_el = self.W_dot*self.eta_alt
        
        W_dot_el_MW = self.W_dot_el/1e6
        
        self.CAPEX_alt = 108900 * W_dot_el_MW**0.5463
        
        self.f_install = 0.35 
        self.CAPEX_install = self.f_install*(self.CAPEX_alt + self.CAPEX_turb)
        
        self.CAPEX_tot = self.CAPEX_turb + self.CAPEX_alt + self.CAPEX_install
            
        return
    
#%%
    def design_system(self, x):
        self.penalty = -1

        self.reset()
        
        if "Omega_choices" in self.params:
            self.psi = x[0]
            self.phi = x[1]
            self.R = x[2]
            self.params['Re_min'] = x[3]
            self.params['Omega'] = self.params['Omega_choices'][np.abs(self.params['Omega_choices'] - x[4]).argmin()]
            self.params['M_1_st'] = x[5]
        else: # Iterate on r_m
            self.psi = x[0]
            self.phi = x[1]
            self.R = x[2]
            self.params['Re_min'] = x[3]
            self.r_m = x[4]
            self.params['M_1_st'] = x[5]

        self.r_tip = []
        self.r_hub = []
        self.r_hub_tip = []
        self.r_ratio2 = []
        
        # First Stator Instanciation
        self.stages.append(self.stage(self.fluid))
        self.stages[0].update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        
        "------------- 1) Isentropic Expansion Calculation -----------------------------------------------" 
        s_in = self.stages[0].total_states['S'][1]
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s_in)
        
        h_is_ex = self.AS.hmass()
        Dh0s = self.stages[0].total_states['H'][1] - h_is_ex
        
        Dh0 = self.inputs['W_dot']/self.inputs['mdot']
        
        self.eta_is = Dh0/Dh0s
        
        "------------- 2) Velocity Triangle Computation (+ Solodity) -------------------------------------" 
        self.computeVelTriangle()
        
        self.solidityStator = 2*np.cos(self.Vel_Tri['alpha2'])/np.cos(self.Vel_Tri['alpha1'])*np.sin(abs(self.Vel_Tri['alpha2']-self.Vel_Tri['alpha1']))/self.params['Zweifel']
        self.solidityRotor  = 2*np.cos(self.Vel_Tri['beta3'])/np.cos(self.Vel_Tri['beta2'])*np.sin(abs(self.Vel_Tri['beta3']-self.Vel_Tri['beta2']))/self.params['Zweifel']
        
        "------------- 3) Guess u from vMax (subsonic flow)  ---------------------------------------------" 
        
        v1st = self.AS.speed_sound() * self.params['M_1_st']
        
        # Assume u based on the maximum speed
        self.Vel_Tri['u'] = v1st / max([self.Vel_Tri['v2OverU'],self.Vel_Tri['w3OverU']])
        
        "------------- 4) Compute number of stage + recompute u  -----------------------------------------" 
        
        # Compute required number of stages based on assumed u
        Dh0Stage = self.psi * self.Vel_Tri['u']**2
        self.nStages = int(round(Dh0/Dh0Stage))
        
        for i in range(self.nStages-1):
            self.stages.append(self.stage(self.fluid))
        
        # Recompute u based on the number of stages to satisfy the work. As r_m is constant, u is contant accross stages
        self.Dh0Stage = Dh0/self.nStages
        self.Vel_Tri['u'] = np.sqrt(self.Dh0Stage/self.psi)

        if "Omega_choices" in self.params:
            self.omega_rads = self.params['Omega']*(2*np.pi)/60
            self.r_m = self.Vel_Tri['u']/(2*self.omega_rads)
        else:
            self.omega_rads = self.Vel_Tri['u']/(2*self.r_m)
            self.params['Omega'] = self.omega_rads*60/(2*np.pi)
        
        "------------- 5) Compute complete velocity triangles and exit losses ----------------------------" 

        # Compute velocity triangle with the value of u
        self.Vel_Tri['vm'] = self.Vel_Tri['vmOverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu2'] = self.Vel_Tri['vu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu3'] = self.Vel_Tri['vu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu2'] = self.Vel_Tri['wu2OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['wu3'] = self.Vel_Tri['wu1'] = self.Vel_Tri['wu3OverU'] * self.Vel_Tri['u']
        self.Vel_Tri['vu1'] = self.Vel_Tri['vu3']

        self.Vel_Tri['v1'] = np.sqrt(self.Vel_Tri['vu1']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['v2'] = np.sqrt(self.Vel_Tri['vu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['v3'] = np.sqrt(self.Vel_Tri['vu3']**2 + self.Vel_Tri['vm']**2)

        self.Vel_Tri['w1'] = np.sqrt(self.Vel_Tri['wu1']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w2'] = np.sqrt(self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)
        self.Vel_Tri['w3'] = np.sqrt(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)

        self.exit_loss = (self.Vel_Tri['vm']**2+self.Vel_Tri['vu3']**2)/2
        self.exit_loss_W = self.exit_loss*self.inputs['mdot']

        "------------- 6) Iterate on r_m and compute repeating stages ------------------------" 

        h_in = self.stages[0].total_states['H'][1] - (self.Vel_Tri['vm']**2)/2
        self.stages[0].update_static_AS(CP.HmassSmass_INPUTS, h_in, s_in, 1)
        
        try:
            self.computeRepeatingStages()        
        except:
            # print("Error in stages")
            self.obj = 10000
            return self.obj
        
        "------------- 7) Compute number of blades per stage ---------------------------" 

        self.n_blade = []

        for stage in self.stages:
              stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
              self.n_blade.append(stage.n_blade_S)

              stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
              self.n_blade.append(stage.n_blade_R)

        "------------- 8) Add last redirecting stator stage -------------------------------------------------------------" 

        self.stages.append(self.stage(self.fluid))
        
        self.computeVelTriangleLastStage()
        
        try:
            self.computeLastStage()
        except:
            # print("Error in last stage")
            self.obj = 10000
            return self.obj

        "------------- 9) Compute main outputs -------------------------------------------------------------" 
        
        hin = self.stages[0].total_states['H'][1]
        hout = self.stages[-1].static_states['H'][2]
        
        self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])

        hout_s = self.AS.hmass()
        
        self.W_dot = self.inputs['mdot']*(hin-hout)
                
        self.eta_is = (hin - hout)/(hin - hout_s)

        penalty_1 = max(self.r_hub_tip[0] - self.params['r_hub_tip_max'],0)*10
        penalty_2 = max(self.params['r_hub_tip_min'] - self.r_hub_tip[-1],0)*10
        
        if abs((self.inputs["p_ex"] - self.stages[-1].static_states['P'][2])/self.inputs["p_ex"]) >= 1e-1:
            penalty_3 = 3*abs((self.inputs["p_ex"] - self.stages[-1].static_states['P'][2])/self.inputs["p_ex"])**2
            self.Pressure_Deviation = self.inputs["p_ex"] - self.stages[-1].static_states['P'][2]
        else:
            penalty_3 = 0
            self.Pressure_Deviation = self.inputs["p_ex"] - self.stages[-1].static_states['P'][2]

        self.penalty = penalty_1 + penalty_2 + penalty_3

        if self.eta_is > 0 and self.eta_is <= 1:
            if self.penalty == 0:
                if self.W_dot > self.inputs['W_dot']:
                    self.inputs['W_dot'] = self.W_dot
                    
            self.obj = -self.eta_is + self.penalty
            # print(f"opt 'success' : {obj}")
        else:
            # print("Bad eta_is")
            self.obj = 10000

        if self.obj < 10000:
            self.allowable_positions.append([self.obj, self.psi, self.phi, self.R, self.params['Re_min'], self.r_m, self.eta_is])
        
        return self.obj

    from scipy.optimize import brentq

    def _pout_minus_target(self, r_m, x_no_rm):
        """
        For a given r_m, evaluate P_out - p_ex at the end of design_system.
        x_no_rm is [psi, phi, R, Re_min, M_1_st] (no r_m here).
        Returns a scalar.
        """
        # Build a full x with r_m inserted in position 4
        x = np.array([x_no_rm[0], x_no_rm[1], x_no_rm[2], x_no_rm[3], r_m, x_no_rm[4]], dtype=float)
        obj = self.design_system(x)
        # If design failed badly, return a large signed residual that pushes the solver away
        if obj >= 1e4 or not hasattr(self, "stages") or not self.stages:
            # keep sign based on the last known deviation if available
            dev = getattr(self, "Pressure_Deviation", np.nan)
            if np.isnan(dev):
                # arbitrary sign; push upward
                return +1e6
            return -float(dev)  # dev = p_ex - p_out => Pout - Pex = -dev
        # Here we have a valid evaluation
        # Pressure_Deviation = p_ex - p_out (you set this in design_system)
        return -float(self.Pressure_Deviation)  # => P_out - p_ex
            
    
    def design_system_match_p(self, x_no_rm):
        """
        Solve for r_m in [r_m_min, r_m_max] so that outlet static pressure ~ p_ex.
        x_no_rm: [psi, phi, R, Re_min, M_1_st]
        Returns objective (same convention as design_system).
        """
        rmin, rmax = self.params["r_m_bounds"]
        # 1) Try to bracket a root quickly
        fmin = self._pout_minus_target(rmin, x_no_rm)
        fmax = self._pout_minus_target(rmax, x_no_rm)
    
        r_star = None
        if np.isfinite(fmin) and np.isfinite(fmax) and fmin * fmax < 0.0:
            # 2) Root exists in the bracket
            try:
                r_star = brentq(lambda r: self._pout_minus_target(r, x_no_rm), rmin, rmax, maxiter=30, xtol=1e-6)
            except Exception:
                r_star = None
    
        if r_star is None:
            # 3) No sign change or brentq failed: pick the best r_m on a coarse grid
            grid = np.linspace(rmin, rmax, 15)
            residuals = []
            for r in grid:
                residuals.append(abs(self._pout_minus_target(r, x_no_rm)))
            r_star = grid[int(np.argmin(residuals))]
    
        # Final consistent call at r_star (populates all fields)
        obj = self._pout_minus_target(r_star, x_no_rm)  # runs design_system inside
        # _pout_minus_target already made a final design_system() call; objective is stored there.
        # Retrieve the actual objective from the last run:
        # (design_system returns obj; to avoid recompute, just take it from self.allowable_positions if needed)
        return getattr(self, "obj", 0.0) if hasattr(self, "obj") else 1e4


    def design(self):
        bounds = (np.array([
            self.params['psi_bounds'][0],
            self.params['phi_bounds'][0],
            self.params['R_bounds'][0],
            self.params['Re_bounds'][0],
            self.params['r_m_bounds'][0],
        ]),
        np.array([
            self.params['psi_bounds'][1],
            self.params['phi_bounds'][1],
            self.params['R_bounds'][1],
            self.params['Re_bounds'][1],
            self.params['r_m_bounds'][1],
        ]))
    
        def objective_wrapper(x):
            # Round inputs if needed for discrete steps
            rounded_x = np.copy(x)
            # (Optional rounding logic goes here)
            return np.array([self.design_system(xi) for xi in rounded_x])
    
        # Initialize the optimizer
        optimizer = ps.single.GlobalBestPSO(
            n_particles=20,
            dimensions=5,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
                
        # Custom stopping logic
        patience = 5
        tol = 1e-3
        max_iter = 1
        no_improve_counter = 0
        best_cost = np.inf
    
        for i in range(max_iter):
            # One iteration step
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
    
            current_best = optimizer.swarm.best_cost
    
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
    
            # print(f"[{i+1}] Best cost: {best_cost:.6f}")
    
            if no_improve_counter >= patience:
                print("Stopping early due to stagnation.")
                break
             
        best_pos = optimizer.swarm.best_pos
             
        self.design_system(best_pos)
        
        self.cost_estimation()
        
        "------------- Print Main Results -------------------------------------------------------------" 
        
        print(f"Parameters : {self.psi, self.phi, self.R, self.params['Re_min'], self.r_m}")
        
        print(f"Turbine mean radius: {self.r_m} [m]")
        print(f"Turbine rotation speed: {self.params['Omega']} [RPM]")
        print(f"Turbine number of stage : {self.nStages} [-]")
        print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")        
        print(f"Turbine Generation : {self.W_dot} [W]")        
        
        return best_pos

#%% 

    # def design_parallel(self, n_jobs=-1, backend="loky", chunksize="auto"):
    #     import numpy as np
    #     import pyswarms as ps
    #     from joblib import Parallel, delayed
        
    #     if "Omega_choices" in self.params:
    #         bounds = (np.array([
    #             self.params['psi_bounds'][0],
    #             self.params['phi_bounds'][0],
    #             self.params['R_bounds'][0],
    #             self.params['Re_bounds'][0],
    #             self.params['Omega_choices'][0],
    #             self.params['M_1st_bounds'][0],
    #         ]), np.array([
    #             self.params['psi_bounds'][1],
    #             self.params['phi_bounds'][1],
    #             self.params['R_bounds'][1],
    #             self.params['Re_bounds'][1],
    #             self.params['Omega_choices'][-1],
    #             self.params['M_1st_bounds'][1],
    #         ]))
    #     else: # Iterate on r_m
    #         bounds = (np.array([
    #             self.params['psi_bounds'][0],
    #             self.params['phi_bounds'][0],
    #             self.params['R_bounds'][0],
    #             self.params['Re_bounds'][0],
    #             self.params['r_m_bounds'][0],
    #             self.params['M_1st_bounds'][0],
    #         ]), np.array([
    #             self.params['psi_bounds'][1],
    #             self.params['phi_bounds'][1],
    #             self.params['R_bounds'][1],
    #             self.params['Re_bounds'][1],
    #             self.params['r_m_bounds'][1],
    #             self.params['M_1st_bounds'][1],
    #         ]))
        
    #     # ---- normalize & snapshot inputs (PRIMITIVES ONLY) ----
    #     inp = dict(self.inputs)  # shallow copy
    #     def pick(*names, default=None):
    #         for n in names:
    #             if n in inp and inp[n] is not None:
    #                 return inp[n]
    #         return default
    
    #     inputs_snapshot = {
    #         # common aliases: lowercase/uppercase and 0 vs no-0
    #         "p0_su": pick("p0_su", "P0_su", "P_su", "p_su"),
    #         "T0_su": pick("T0_su", "t0_su", "T_su", "t_su"),
    #         "p_ex" : pick("p_ex", "P_ex"),
    #         "mdot" : pick("mdot", "m_dot"),
    #         "W_dot": pick("W_dot", "W"),
    #     }
    
    #     # quick sanity check to fail early & clearly
    #     missing = [k for k,v in inputs_snapshot.items() if v is None]
    #     if missing:
    #         raise ValueError(f"Missing required inputs for design: {missing}. "
    #                          f"Found keys: {sorted(self.inputs.keys())}")
    
    #     # ---- snapshot params/stage_params too ----
    #     snapshot = {
    #         "cls": type(self),                        # class object (cloudpickle handles it)
    #         "fluid": self.fluid,
    #         "params": dict(self.params),              # ensure only picklable values
    #         "stage_params": getattr(self, "stage_params", None),
    #         "inputs": inputs_snapshot,
    #     }
    
    #     def objective_wrapper(X):
    #         X = np.asarray(X, dtype=float)
    #         with tqdm_joblib(tqdm(total=len(X), desc="Particles", unit="pt")):
    #             results = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
    #                 delayed(_eval_particle)(xi, snapshot["cls"], snapshot["fluid"],
    #                                         snapshot["params"], snapshot["stage_params"],
    #                                         snapshot["inputs"])
    #                 for xi in X
    #             )
    #         return results

    
    #     optimizer = ps.single.GlobalBestPSO(
    #         n_particles=20, dimensions=6,
    #         options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
    #         bounds=bounds
    #     )
    
    #     patience, tol, max_iter = 5, 1e-3, 40
    #     no_improve, best_cost = 0, float("inf")
    
    #     for i in range(max_iter):
    #         optimizer.optimize(objective_wrapper, iters=1, verbose=False)
    #         cur = optimizer.swarm.best_cost
    #         if cur < best_cost - tol:
    #             best_cost, no_improve = cur, 0
    #         else:
    #             no_improve += 1
    #         print(f"[{i+1}] Best cost: {best_cost:.6f}")
    #         if no_improve >= patience:
    #             print("Stopping early due to stagnation.")
    #             break
    
    #     best_pos = optimizer.swarm.best_pos
    #     self.design_system(best_pos)
    #     self.cost_estimation()
    
    #     print(f"Parameters : {self.psi, self.phi, self.R, self.params['Re_min'], self.r_m}")
    #     print(f"P_in : {self.stages[0].static_states['P'][1]} [Pa]")
    #     print(f"P_out : {self.stages[-1].static_states['P'][2]} [Pa]")
    #     print(f"r_h_tip_in : {self.r_hub_tip[0]}")
    #     print(f"r_h_tip_out : {self.r_hub_tip[-1]}")
    #     print(f"M_1_st : {self.params['M_1_st']}")
    #     print(f"Turbine mean radius: {self.r_m} [m]")
    #     print(f"Turbine rotation speed: {self.params['Omega']} [RPM]")
    #     print(f"Turbine number of stage : {self.nStages} [-]")
    #     print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")
    #     print(f"Turbine Generation : {self.W_dot} [W]")
    #     return best_pos

    def design_parallel(self, n_jobs=-1, backend="loky", chunksize="auto", match_p_ex=True):
        import numpy as np
        import pyswarms as ps
        from joblib import Parallel, delayed
    
        if match_p_ex:
            # No r_m in the swarm; r_m will be solved by brentq inside each eval
            bounds = (np.array([
                self.params['psi_bounds'][0],
                self.params['phi_bounds'][0],
                self.params['R_bounds'][0],
                self.params['Re_bounds'][0],
                self.params['M_1st_bounds'][0],
            ]), np.array([
                self.params['psi_bounds'][1],
                self.params['phi_bounds'][1],
                self.params['R_bounds'][1],
                self.params['Re_bounds'][1],
                self.params['M_1st_bounds'][1],
            ]))
            dimensions = 5
        else:
            # Keep your current 6D particle (includes r_m)
            bounds = (np.array([
                self.params['psi_bounds'][0],
                self.params['phi_bounds'][0],
                self.params['R_bounds'][0],
                self.params['Re_bounds'][0],
                self.params['r_m_bounds'][0],
                self.params['M_1st_bounds'][0],
            ]), np.array([
                self.params['psi_bounds'][1],
                self.params['phi_bounds'][1],
                self.params['R_bounds'][1],
                self.params['Re_bounds'][1],
                self.params['r_m_bounds'][1],
                self.params['M_1st_bounds'][1],
            ]))
            dimensions = 6
    
        # snapshot (unchanged)
        inp = dict(self.inputs)
        def pick(*names, default=None):
            for n in names:
                if n in inp and inp[n] is not None:
                    return inp[n]
            return default
    
        inputs_snapshot = {
            "p0_su": pick("p0_su", "P0_su", "P_su", "p_su"),
            "T0_su": pick("T0_su", "t0_su", "T_su", "t_su"),
            "p_ex" : pick("p_ex", "P_ex"),
            "mdot" : pick("mdot", "m_dot"),
            "W_dot": pick("W_dot", "W"),
        }
        missing = [k for k,v in inputs_snapshot.items() if v is None]
        if missing:
            raise ValueError(f"Missing required inputs for design: {missing}. Found keys: {sorted(self.inputs.keys())}")
    
        snapshot = {
            "cls": type(self),
            "fluid": self.fluid,
            "params": dict(self.params),
            "stage_params": getattr(self, "stage_params", None),
            "inputs": inputs_snapshot,
        }
    
        def objective_wrapper(X):
            X = np.asarray(X, dtype=float)
            with tqdm_joblib(tqdm(total=len(X), desc="Particles", unit="pt")):
                results = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
                    delayed(_eval_particle)(
                        xi, snapshot["cls"], snapshot["fluid"],
                        snapshot["params"], snapshot["stage_params"],
                        snapshot["inputs"], match_p_ex=match_p_ex
                    ) for xi in X
                )
            return results
    
        optimizer = ps.single.GlobalBestPSO(
            n_particles=20, dimensions=dimensions,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        patience, tol, max_iter = 5, 1e-3, 40
        no_improve, best_cost = 0, float("inf")
        for i in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            cur = optimizer.swarm.best_cost
            if cur < best_cost - tol:
                best_cost, no_improve = cur, 0
            else:
                no_improve += 1
            print(f"[{i+1}] Best cost: {best_cost:.6f}")
            if no_improve >= patience:
                print("Stopping early due to stagnation.")
                break
    
        best_pos = optimizer.swarm.best_pos
    
        # Finalize on best position
        if match_p_ex:
            self.design_system_match_p(best_pos)
        else:
            self.design_system(best_pos)
        self.cost_estimation()
    
        # prints...
        print(f"P_in : {self.stages[0].static_states['P'][1]} [Pa]")
        print(f"P_out: {self.stages[-1].static_states['P'][2]} [Pa]")
        print(f"r_m  : {self.r_m} [m]")
        print(f"Omega: {self.params['Omega']} [RPM]")
        print(f"eta_is: {self.eta_is}")
        print(f"W_dot : {self.W_dot} [W]")
        return best_pos


#%%

case_study = "TCO2_ORC"

Turb = AxialTurbineMeanLineDesign('Cyclopentane')

if case_study == 'Cuerva':

    Turb = AxialTurbineMeanLineDesign('Cyclopentane')
    
    Turb.set_inputs(
        mdot = 46.18, # kg/s
        W_dot = 4500*1e3, # W
        p0_su = 1230*1e3, # Pa
        T0_su = 273.15 + 158, # K
        p_ex = 78300, # Pa
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        AR_min = 0.8, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        Re_bounds = [1*1e5,1*1e6], # [-]
        psi_bounds = [1,2.5], # [-]
        phi_bounds = [0.4,0.8], # [-]
        R_bounds = [0.4,0.6], # [-]
        r_m_bounds = [0.1, 0.6], # [m]
        M_1_st = 0.3, # [-]
        damping = 0.2, # [-]
        delta_tip = 0.4*1e-3, # [m] : tip clearance
        N_lw = 0, # [-] : Number of lashing wires
        D_lw = 0, # [m] : Diameter of lashing wires
        e_blade = 0.002*1e-3, # [m] : blade roughness
        t_TE_o = 0.05, # [-] : trailing edge to throat opening ratio
        t_TE_min = 5*1e-4, # [m]
        )

elif case_study == 'Zorlu':
    
    Turb = AxialTurbineMeanLineDesign('Cyclopentane')

    Turb.set_inputs(
        mdot = 34.51, # kg/s
        W_dot = 2506000, # W
        p0_su = 767800, # Pa
        T0_su = 273.15 + 131, # K
        p_ex = 82000, # Pa
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        AR_min = 0.8, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        Re_bounds = [1*1e5,1*1e6], # [-]
        psi_bounds = [1,2.5], # [-]
        phi_bounds = [0.4,0.7], # [-]
        R_bounds = [0.45,0.55], # [-]
        r_m_bounds = [0.15, 0.5], # [m]
        M_1_st = 0.3, # [-]
        damping = 0.2, # [-]
        delta_tip = 0.4*1e-3, # [m] : tip clearance
        N_lw = 0, # [-] : Number of lashing wires
        D_lw = 0, # [m] : Diameter of lashing wires
        e_blade = 0.002*1e-3, # [m] : blade roughness
        t_TE_o = 0.05, # [-] : trailing edge to throat opening ratio
        t_TE_min = 5*1e-4, # [m]
        )
    
elif case_study == 'TCO2_ORC':

    Turb = AxialTurbineMeanLineDesign('CO2')

    Turb.set_inputs(
        mdot = 5*100, # kg/s
        W_dot = 5*4.69*1e6, # W
        p0_su = 140*1e5, # Pa
        T0_su = 273.15 + 121, # K
        p_ex = 39.8*1e5, # Pa
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        AR_min = 0.8, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.5, # [-]
        Re_bounds = [3*1e6,8*1e6], # [-]
        psi_bounds = [1.3,2.5], # [-]
        phi_bounds = [0.5,0.8], # [-]
        R_bounds = [0.45,0.55], # [-]
        M_1st_bounds = [0.2, 0.5], # [-]
        r_m_bounds = [0.1, 0.5], # [m]
        Omega_choices = [500,750,1000,1500,3000], # [RPM] : [500,750,1000,1500,3000]
        damping = 0.3, # [-]
        delta_tip = 0.4*1e-3, # [m] : tip clearance
        N_lw = 0, # [-] : Number of lashing wires
        D_lw = 0, # [m] : Diameter of lashing wires
        e_blade = 0.002*1e-3, # [m] : blade roughness
        t_TE_o = 0.05, # [-] : trailing edge to throat opening ratio
        t_TE_min = 5*1e-4, # [m]
        )

elif case_study == 'Salah_Case':

    Turb = AxialTurbineMeanLineDesign('CO2')

    Turb.set_inputs(
        mdot = 655.18, # kg/s
        W_dot = 100*1e6, # W
        p0_su = 250*1e5, # Pa
        T0_su = 923, # K
        p_ex = 100*1e5, # Pa
        )
    
    Turb.set_parameters(
        Zweifel = 0.8, # [-]
        AR_min = 0.8, # [-]
        r_hub_tip_max = 0.95, # [-]
        r_hub_tip_min = 0.6, # [-]
        Re_bounds = [1*1e6,7*1e6], # [-]
        psi_bounds = [1.5,2.5], # [-]
        phi_bounds = [0.6,0.9], # [-]
        R_bounds = [0.45,0.55], # [-]
        r_m_bounds = [0.15, 0.5], # [m]
        M_1_st = 0.5, #0.3, # [-]
        damping = 0.2, # [-]
        delta_tip = 0.4*1e-3, # [m] : tip clearance
        N_lw = 0, # [-] : Number of lashing wires
        D_lw = 0, # [m] : Diameter of lashing wires
        e_blade = 0.002*1e-3, # [m] : blade roughness
        t_TE_o = 0.05, # [-] : trailing edge to throat opening ratio
        t_TE_min = 5*1e-4, # [m]
        )
    
best_pos = Turb.design_parallel()

Turb.plot_geometry()
Turb.plot_n_blade()
Turb.plot_radius_verif()
Turb.plot_Mollier()

# (1.4528688807250987, 0.4939933148314975, 0.42406128096733176, 6364578.748939947, 0.4557043394023795)
