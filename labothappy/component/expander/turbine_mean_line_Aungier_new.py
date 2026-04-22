
from connector.mass_connector import MassConnector
from correlations.turbomachinery.aungier_axial_turbine import aungier_loss_model
from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
 
from toolbox.turbomachinery.mean_line_axial_turbine_mapping import map_plot, map_plot_clean, plot_power_eta_vs_mdot, filter_sparse_by_proximity
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, differential_evolution
import pyswarms as ps
 
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
 
import warnings
warnings.filterwarnings("ignore")
 
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
            self.total_states = {
                'H': np.empty(3, dtype=np.float64),
                'S': np.empty(3, dtype=np.float64),
                'P': np.empty(3, dtype=np.float64),
                'D': np.empty(3, dtype=np.float64),
                'A': np.empty(3, dtype=np.float64),
                'V': np.empty(3, dtype=np.float64),
            }
            self.static_states = {
                k: np.empty(3, dtype=np.float64) for k in self.total_states
            }
            
            self._last_update = {'tot': None, 'stat': None}
            
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
           


        def _maybe_update(self, kind, CP_INPUTS, a, b):
            key = (CP_INPUTS, float(a), float(b))
            if self._last_update[kind] != key:
                self.AS.update(CP_INPUTS, a, b)
                self._last_update[kind] = key
                return True   # did update
            return False      # skipped identical state
    
        def update_total_AS(self, CP_INPUTS, a, b, position):
            i = int(position-1)
            self._maybe_update('tot', CP_INPUTS, a, b)
            AS = self.AS
            T = self.total_states
            T['H'][i] = AS.hmass()
            T['S'][i] = AS.smass()
            T['P'][i] = AS.p()
            T['D'][i] = AS.rhomass()
            try:   T['A'][i] = AS.speed_sound()
            except: T['A'][i] = -1.0
            T['V'][i] = AS.viscosity()
    
        def update_static_AS(self, CP_INPUTS, a, b, position):
            i = int(position-1)
            self._maybe_update('stat', CP_INPUTS, a, b)
            AS = self.AS
            S = self.static_states
            S['H'][i] = AS.hmass()
            S['S'][i] = AS.smass()
            S['P'][i] = AS.p()
            S['D'][i] = AS.rhomass()
            try:   S['A'][i] = AS.speed_sound()
            except: S['A'][i] = -1.0
            S['V'][i] = AS.viscosity()

        def get_total_prop(self, prop, position):
            value = self.total_states[prop][position-1]
            return value
        
        def get_static_prop(self, prop, position):
            value = self.static_states[prop][position-1]
            return value
 
    # --- Anderson acceleration helper -------------------------------------------
    
    class AndersonMixer:
        """
        Anderson acceleration (type-I) for fixed-point problems x = g(x).
        Designed for small vectors (here: 2D). Stable and fast.
        """
        def __init__(self, m=2, beta=0.5, damping=0.3, reg=1e-10):
            self.m = int(m)            # history size
            self.beta = float(beta)    # relaxation on AA step
            self.damping = float(damping)  # fallback Picard damping
            self.reg = float(reg)      # Tikhonov regularization
            self.reset()
    
        def reset(self):
            self.x_hist = []
            self.f_hist = []
            self.k = 0
    
        def step(self, xk, gk):
            """Return next iterate given x_k and g(x_k)."""
            fk = gk - xk
    
            # first step: simple damped Picard
            if self.k == 0:
                x_next = (1.0 - self.damping) * xk + self.damping * gk
                self.x_hist = [xk.copy()]
                self.f_hist = [fk.copy()]
                self.k = 1
                return x_next
    
            # build ΔF and ΔX columns from last m steps
            x_stack = self.x_hist[-self.m:] + [xk]
            f_stack = self.f_hist[-self.m:] + [fk]
            dF_cols, dX_cols = [], []
            for j in range(1, len(f_stack)):
                dF_cols.append(f_stack[j] - f_stack[j-1])
                dX_cols.append(x_stack[j] - x_stack[j-1])
    
            if not dF_cols:  # safety
                return (1.0 - self.damping) * xk + self.damping * gk
    
            dF = np.column_stack(dF_cols)
            dX = np.column_stack(dX_cols)
    
            # solve (dF^T dF + reg I) alpha = dF^T f_k
            AtA = dF.T @ dF
            rhs = dF.T @ fk
            AtA.flat[::AtA.shape[0] + 1] += self.reg  # regularize
    
            try:
                alpha = np.linalg.solve(AtA, rhs)
                x_aa = gk - dX @ alpha               # AA-I formula
                x_next = (1.0 - self.beta) * gk + self.beta * x_aa
            except np.linalg.LinAlgError:
                x_next = (1.0 - self.damping) * xk + self.damping * gk
    
            # update history
            self.x_hist.append(xk.copy())
            self.f_hist.append(fk.copy())
            if len(self.x_hist) > self.m + 1:
                self.x_hist = self.x_hist[-(self.m+1):]
                self.f_hist = self.f_hist[-(self.m+1):]
            self.k += 1
            return x_next
    
    
    def solve_fixed_point(self, func, x0, mixer, tol=1e-6, max_iter=100):
        """
        Generic fixed-point driver. func(x) must return g(x).
        Uses the residual sum(|(x-g)/g|) like your profiler loop.
        """
        x = np.asarray(x0, dtype=float)
        for it in range(max_iter):
            print(f"x : {x}")
            
            g = np.asarray(func(x), dtype=float)
            res = np.sum(np.abs((x - g) / np.maximum(1e-12, np.abs(g))))
            if res <= tol:
                return g, it + 1
            x = mixer.step(x, g)
        return x, max_iter   
 
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
               
    # ---------------- Loss Models ------------------------------------------------------------------------
    
    def stator_blade_row_system(self, x):
                
        stage = self.stages[self.curr_stage_index]
 
        # 1) Guess outlet state
        h_static_out = x[0]*1e5
        p_static_out = x[1]*1e5
               
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
       
        stage.Vel_Tri_S['u'] = u = self.u
       
        A_flow = stage.h_blade_S*(2*np.pi*self.params['r_m'])
        stage.Vel_Tri_S['vm'] = vm = self.inputs['m_dot']/(stage.static_states['D'][1]*A_flow)
       
        if self.curr_stage_index == 0:
            stage.Vel_Tri_S['alpha1'] = alpha1 = stage.xhi_S1
            stage.Vel_Tri_S['vu1'] = vu1 = vm*np.tan(alpha1)
        else:
            stage.Vel_Tri_S['wu1'] = wu1 = np.tan(stage.Vel_Tri_S['beta1'])*vm
            stage.Vel_Tri_S['vu1'] = vu1 = wu1 + u
            stage.Vel_Tri_S['alpha1'] = alpha1 = np.arctan(vu1/vm)
       
       
        stage.Vel_Tri_S['v1'] = v1 = np.sqrt(stage.Vel_Tri_S['vm']**2 + stage.Vel_Tri_S['vu1']**2)
        stage.M1_S = v1/stage.static_states['A'][0]
 
        # 2) Compute total inlet state
        hin = stage.static_states['H'][0]
        h0in = hin + (vu1**2 + vm**2)/2  
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.static_states['S'][0], 1)            
       
        h02 = h0in
       
        stage.Vel_Tri_S['v2'] = v2 = np.sqrt(2*(h02 - h_static_out))
        stage.M2_S = v2/stage.static_states['A'][1]
        stage.Vel_Tri_S['alpha2'] = alpha2 = np.arctan2(np.sqrt(v2**2 - vm**2), vm)
 
        # 4) Compute cord, aspect ratio, blade pitch and blade number
       
        stage.Re_s = stage.chord_S*(stage.static_states['D'][1]*vm)/(stage.static_states['V'][1])
        stage.AR_S = stage.h_blade_S/stage.chord_S
               
        # 5) Loss model
       
        stage.beta_g_S = np.arcsin(stage.o_S/stage.pitch_S)
       
        self.compute_deviation_stator(stage)
        alpha2_calc = stage.xhi_S2 + stage.delta_S
       
        v2_new = vm/np.cos(alpha2_calc) 
       
        stage.Y_vec_S = aungier_loss_model(stage.Vel_Tri_S['alpha1'], alpha2_calc, stage.beta_g_S*180/np.pi, stage.xhi_S1, stage.chord_S,
                               0, self.params['D_lw'], self.params['e_blade'], stage.h_blade_S, stage.static_states['V'][1],
                               stage.M1_S, stage.M2_S, self.params['N_lw'], stage.R_c_S, stage.static_states['D'][1], stage.pitch_S, stage.t_blade_S, stage.t_TE_S,
                               vm, v2_new,1)
       
        Y = stage.Y_vec_S['Y_tot']
        p0_out = (stage.total_states['P'][0] + Y * p_static_out)/(1+Y)
 
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.total_states['S'][1]
       
        hout = stage.total_states['H'][1]-v2_new**2/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
               
        pout_calc = stage.static_states['P'][1]
 
        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.static_states['S'][0])
        hout_s = self.AS.hmass()
 
        stage.eta_is_S = (stage.static_states['H'][0]-stage.static_states['H'][1])/(stage.static_states['H'][0]-hout_s)
 
        return np.array([hout, pout_calc])*1e-5 # (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def rotor_blade_row_system(self, x):
        
        stage = self.stages[self.curr_stage_index]
        
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x*1e5
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 3)
        
        v = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['vu3']**2)
        w = np.sqrt(self.Vel_Tri['vm']**2 + self.Vel_Tri['wu3']**2)
        stage.M_R = max(v,w)/stage.get_static_prop('A',3)
        
        # 2) Compute total inlet state
        hin = stage.get_static_prop('H',2)
        h0in = hin + (self.Vel_Tri['wu2']**2 + self.Vel_Tri['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.get_static_prop('S',2), 2)            
        
        # 3) Compute A_flow and h_blade based on r_m 
        stage.A_flow_R = self.inputs['mdot']/(stage.get_static_prop('D',3)*self.Vel_Tri['vm'])
        stage.h_blade_R = stage.A_flow_R/(2*np.pi*self.r_m)
        
        if stage.h_blade_R > self.r_m*0.8:
            raise ValueError('')
        
        # 4) Compute cord, aspect ratio, pitch and blade number
        stage.chord_R = (self.params['Re_min']*stage.get_static_prop('V',3))/(stage.get_static_prop('D',3)*self.Vel_Tri['vm'])            
        stage.AR_R = stage.h_blade_R/stage.chord_R
        
        stage.pitch_R = stage.chord_R/self.solidityRotor
        stage.n_blade_R = round(2*np.pi*self.r_m/stage.pitch_R)
        
        self.compute_rotor_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha3'] - self.Vel_Tri['alpha2']
        
        stage.R_c_R = 2*np.sin(abs(camber)/2)*stage.chord_R # Geometrical estimation of blade curvature radius
                
        stage.M2_R = self.Vel_Tri['w2']/stage.get_static_prop('A',2)
        stage.M3_R = self.Vel_Tri['w3']/stage.get_static_prop('A',3)
        
        stage.beta_g_R = 0.5*(self.Vel_Tri['beta2'] + self.Vel_Tri['beta3'])  # mid-passage metal angle
        stage.o_R = abs(np.sin(stage.beta_g_R)*stage.pitch_R)
        
        
        stage.t_TE_R = max(self.params['t_TE_o']*stage.pitch_R * np.cos(self.Vel_Tri['alpha2'])/(1+self.params['t_TE_o']), self.params['t_TE_min'])
        
        stage.Y_vec_R = aungier_loss_model(-self.Vel_Tri['beta2'], -self.Vel_Tri['beta3'], stage.beta_g_R*180/np.pi, -self.Vel_Tri['beta2'], stage.chord_R, 
                               self.params['delta_tip'], self.params['D_lw'], self.params['e_blade'], stage.h_blade_R, stage.get_static_prop('V',3), 
                               stage.M2_R, stage.M3_R, self.params['N_lw'], stage.R_c_R, stage.get_static_prop('D',3), stage.pitch_R, stage.t_blade_R, stage.t_TE_R,
                               self.Vel_Tri['vm'], self.Vel_Tri['w3'],1)

        Y = stage.Y_vec_R['Y_tot']
                
        p0_out = (stage.get_total_prop('P',2) + Y * p_static_out)/(1+Y)

        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 3)
        sout = stage.get_total_prop('S',3)
        
        hout = h0in-(self.Vel_Tri['wu3']**2 + self.Vel_Tri['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 3)
        
        pout_calc = stage.get_static_prop('P',3)

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.get_static_prop('S',2))
        hout_s = self.AS.hmass()

        stage.eta_is_R = (stage.get_static_prop('H',2)-stage.get_static_prop('H',3))/(stage.get_static_prop('H',2)-hout_s)

        return np.array([hout, pout_calc])*1e-5 # (p_static_out - pout_calc)**2 + (h_static_out - hout)**2

    def last_blade_row_system(self, x):
        # 1) Guess outlet state
        [h_static_out, p_static_out] = x*1e5
        
        stage = self.stages[-1]
        
        stage.update_static_AS(CP.HmassP_INPUTS, h_static_out, p_static_out, 2)
        
        v = np.sqrt(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu2']**2)
        w = np.sqrt(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['wu2']**2)
        stage.M_S = max(v,w)/stage.get_static_prop('A',2)
        
        # 2) Compute total inlet state
        hin = stage.get_static_prop('H',1)
        h0in = hin + (self.Vel_Tri_Last_Stage['vu1']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2  
        
        stage.update_total_AS(CP.HmassSmass_INPUTS, h0in, stage.get_static_prop('S',1), 1)            
        
        # 3) Compute A_flow and h_blade based on r_m guess
        stage.A_flow_S = self.inputs['mdot']/(stage.get_static_prop('D',2)*self.Vel_Tri['vm'])
        stage.h_blade_S = stage.A_flow_S/(2*np.pi*self.r_m)

        if stage.h_blade_S > self.r_m*0.8:
            raise ValueError('')

        # 4) Compute cord, aspect ratio, blade pitch and blade number
                
        stage.chord_S = (self.params['Re_min']*stage.get_static_prop('V',2))/(stage.get_static_prop('D',2)*self.Vel_Tri['vm'])            
        stage.AR_S = stage.h_blade_S/stage.chord_S
        
        stage.pitch_S = stage.chord_S/self.solidityStator
        stage.n_blade_S = round(2*np.pi*self.r_m/stage.pitch_S)
        
        self.compute_stator_t_max(stage)
        
        # 5) Loss model
        
        camber = self.Vel_Tri['alpha2'] - self.Vel_Tri['alpha1']
        
        stage.R_c_S = 2*np.sin(abs(camber)/2)*stage.chord_S # Geometrical estimation of blade curvature radius
                
        stage.M1_S = self.Vel_Tri['v1']/stage.get_static_prop('A',1)
        stage.M2_S = self.Vel_Tri['v2']/stage.get_static_prop('A',2)
        
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
        DP_loss = (Yp+Yew)*(self.Vel_Tri_Last_Stage['vm']**2 + self.Vel_Tri_Last_Stage['vu1']**2)*stage.get_static_prop('D',1)/2
        p0_out = stage.get_total_prop('P',1)-DP_loss
                
        # Computation of static outlet pressure
        stage.update_total_AS(CP.HmassP_INPUTS, h0in, p0_out, 2)
        sout = stage.get_total_prop('S',2)
        
        hout = h0in-(self.Vel_Tri_Last_Stage['vu2']**2 + self.Vel_Tri_Last_Stage['vm']**2)/2
        stage.update_static_AS(CP.HmassSmass_INPUTS, hout, sout, 2)
        
        pout_calc = stage.get_static_prop('P',2)

        # Isentropic efficiency of the blade
        self.AS.update(CP.PSmass_INPUTS, pout_calc, stage.get_static_prop('S',1))
        hout_s = self.AS.hmass()

        stage.eta_is_S = (stage.get_static_prop('H',1)-stage.get_static_prop('H',2))/(stage.get_static_prop('H',1)-hout_s)

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
   
    # def computeBladeRow(self,stage_index,row_type):
    #     stage = self.stages[stage_index]
       
    #     self.curr_stage_index = stage_index
               
    #     if row_type == 'S': # Stator
           
    #         print("Stator")
       
    #         if 'P_ex' in self.inputs:
    #             RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
    #         else:
    #             RP_1_row = 5**(1/(2*self.nStages))      
            
    #         if self.Dh0_stage_guess !=0:
    #             h_out_guess = stage.static_states['H'][1] - self.Dh0_stage_guess/2    
    #         else:
    #             h_out_guess = stage.static_states['H'][1]*0.99
               
    #         pout_guess = stage.static_states['P'][1]/RP_1_row
    #         # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])        
           
    #         # Initial guess vector
    #         x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
            
    #         print(f"h_out_guess : {h_out_guess}")
    #         print(f"pout_guess : {pout_guess}")
            
    #         res = 1
    #         x_in = x0_disc
           
    #         c = 0
           
    #         while res > 1e-6:
                
    #             print(f"x_in : {x_in}")
                
    #             if c > 1000:
    #                 raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")
               
    #             # print(f"x_in : {x_in}")
               
    #             x_out = self.stator_blade_row_system(x_in)
 
    #             # print(f"x_out : {x_out}")
               
    #             res_vec = abs((x_in - x_out)/x_out)
    #             res = sum(res_vec)
               
    #             x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out
                             
    #             # print(f"new x_in : {x_in}")
 
    #             c += 1
               
    #         self.stator_blade_row_system(x_out)
                       
    #         # print(f'Y_S : {stage.Y_vec_S}')
 
    #     else: # Rotor
 
    #         # print("Rotor")
 
    #         if 'P_ex' in self.inputs:
    #             RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
    #         else:
    #             RP_1_row = 5**(1/(2*self.nStages))  
           
    #         if self.Dh0_stage_guess !=0:
    #             h_out_guess = stage.static_states['H'][2] - self.Dh0_stage_guess/2    
    #         else:
    #             h_out_guess = stage.static_states['H'][2]*0.99
           
    #         pout_guess = stage.static_states['P'][2]/RP_1_row
    #         # sol = minimize(self.rotor_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])    
           
    #         # Initial guess vector
    #         x0_disc = np.concatenate(([h_out_guess], [pout_guess]))*1e-5
 
    #         res = 1
    #         x_in = x0_disc
           
    #         c = 0
           
    #         while res > 1e-6:
 
    #             if c > 1000:
    #                 raise RuntimeError("Max iterations exceeded in computeBladeRow (stator/rotor/last stage).")
 
    #             # print(f"x_in : {x_in}")
 
    #             x_out = self.rotor_blade_row_system(x_in)
 
    #             # print(f"x_out : {x_out}")
 
    #             res_vec = abs((x_in - x_out)/x_out)
    #             res = sum(res_vec)
               
    #             x_in = (1-self.params['damping'])*x_in + self.params['damping'] * x_out
                       
    #             # print(f"new x_in : {x_in}")
 
    #             c += 1
           
    #         self.rotor_blade_row_system(x_out)
    #         self.compute_deviation_rotor(stage)
                       
    #         # print(f'Y_R : {stage.Y_vec_R}')
 
    #     return
    
    def computeBladeRow(self,stage_index,row_type):
        stage = self.stages[stage_index]

        self.curr_stage_index = stage_index
               
        if row_type == 'S': # Stator
            
            print("Stator")
        
            RP_1_row = (self.inputs['P_su']/self.inputs['P_ex'])**(1/(2*self.nStages))
            h_out_guess = stage.get_static_prop('H',1) *0.99
            pout_guess = stage.get_static_prop('P',1)/RP_1_row
            # sol = minimize(self.stator_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])         
            
            # Initial guess vector
            x0_disc = np.array([h_out_guess, pout_guess]) * 1e-5
            
            # Anderson-accelerated fixed point
            mixer = self.AndersonMixer(
                m=2,
                beta=0.5,
                damping=self.params.get('damping', 0.3),
                reg=1e-10
            )
            
            x_out, iters = self.solve_fixed_point(self.stator_blade_row_system, x0_disc, mixer, max_iter=100)
            
            # finalize at the converged point
            self.stator_blade_row_system(x_out)
            
            stage.xhi_S1 = self.Vel_Tri['alpha1']
            self.compute_deviation_stator(stage)
            stage.xhi_S2 = self.Vel_Tri['alpha2'] - stage.delta_S
            
            # print(f'Y_S : {stage.Y_vec_S}')

        else: # Rotor

            # print("Rotor")

            RP_1_row = (self.inputs['p0_su']/self.inputs['p_ex'])**(1/(2*self.nStages))
            h_out_guess = stage.get_static_prop('H',2) - self.Dh0Stage/2
            pout_guess = stage.get_static_prop('P',2)/RP_1_row
            # sol = minimize(self.rotor_blade_row_system, x0=(h_out_guess,pout_guess), args=(stage), bounds=[(stage.static_states['H'][1]-2*self.Dh0Stage, stage.static_states['H'][1]), (self.inputs['p_ex']*0.8, stage.static_states['P'][1])])    
            
            # Initial guess vector
            x0_disc = np.array([h_out_guess, pout_guess]) * 1e-5
            
            mixer = self.AndersonMixer(
                m=2,
                beta=0.5,
                damping=self.params.get('damping', 0.3),
                reg=1e-10
            )
            x_out, iters = self.solve_fixed_point(self.rotor_blade_row_system, x0_disc, mixer,
                                             tol=1e-6, max_iter=100)
            
            self.rotor_blade_row_system(x_out)
                        
            stage.xhi_R1 = self.Vel_Tri['beta2']
            self.compute_deviation_rotor(stage)
            stage.xhi_R2 = self.Vel_Tri['beta3'] - stage.delta_R
                        
            # print(f'Y_R : {stage.Y_vec_R}')
        return
    
    
    def computeRepeatingStages(self):
        
        print(f"{self.params['nStages']} stages.")

        self.nStages = self.params['nStages']
       
        for i in range(int(self.nStages)):
            
            print(f"Stage {i+1}")
            
            if i == 0:
                self.computeBladeRow(i, 'S')
                self.computeBladeRow(i, 'R')

            else:
                for prop in self.stages[i].static_states.keys():
                    self.stages[i].static_states[prop][0] = self.stages[i-1].static_states[prop][2]
                
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
        for prop in self.stages[-1].static_states.keys():
            stage.static_states[prop][0] = self.stages[-2].static_states[prop][2]
        
        h_out_guess = self.stages[-2].get_total_prop('H',3)
        pout_guess = self.stages[-2].get_total_prop('P',3)
        # sol = minimize(self.last_blade_row_system, x0=(h_out_guess,pout_guess), bounds=[(self.stages[-1].static_states['H'][1], h_out_guess), (self.stages[-1].static_states['P'][1], pout_guess)])         
        
        # Initial guess vector
        x0_disc = np.array([h_out_guess, pout_guess]) * 1e-5

        mixer = self.AndersonMixer(
            m=2,
            beta=0.5,
            damping=self.params.get('damping', 0.3),
            reg=1e-10
        )
        x_out, iters = self.solve_fixed_point(self.last_blade_row_system, x0_disc, mixer,
                                         tol=1e-6, max_iter=100)
        
        self.last_blade_row_system(x_out)
        
        stage.xhi_S1 = self.Vel_Tri['alpha1']
        self.compute_deviation_stator(stage)
        stage.xhi_S2 = self.Vel_Tri['alpha2'] - stage.delta_S
        
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
       
        self.check_calculable()
        self.check_parametrized()
       
        if not self.calculable:
            print("Component is not calculable. Check inputs.")
            return
       
        if not self.parametrized:
            print("Component is not parametrized. Check parameters.")
            return
       
        self.omega_rads = 2*np.pi*self.inputs['N_rot']/60
        self.u = self.omega_rads*self.params['r_m']*2
       
        self.stages[0].update_static_AS(CP.PT_INPUTS, self.su.p, self.su.T, 1)
        # self.stages.append(self.stage(self.fluid))
        
        print("OH")
        
        self.computeRepeatingStages()
       
        self.computeLastStage()
       
        hin = self.stages[0].total_states['H'][1]
        hout = self.stages[-1].static_states['H'][2]
       
        self.AS.update(CP.PSmass_INPUTS, self.stages[-1].static_states['P'][2], self.stages[0].static_states['S'][1])
 
        hout_s = self.AS.hmass()
       
        self.W_dot = self.inputs['m_dot']*(hin-hout)
               
        self.eta_is = (hin - hout)/(hin - hout_s)
       
        self.update_connectors()
        self.solved = True
        return
 
    def update_connectors(self):
 
        self.ex.set_fluid(self.su.fluid)        
        self.ex.set_p(self.stages[-1].total_states['P'][2])
        self.ex.set_h(self.stages[-1].total_states['H'][2])
        self.ex.set_m_dot(self.su.m_dot)
       
        return
 
 