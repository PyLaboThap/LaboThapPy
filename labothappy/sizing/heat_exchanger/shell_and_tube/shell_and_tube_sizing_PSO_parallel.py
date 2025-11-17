"""
INDEPENDENT VARIABLES
---------------------

D_o possible values : [1/2, 3/4, 1, 1+1/4, 1+1/2]*25.4*1e-3 # [m]

Shell_ID possible values : 
[8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27, 29, 31, 33, 35,
37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120]*25.4*1e-3 [m]

Central Spacing Limited Values : 
[Shell_ID/5; 74*D_o**0.75] # To put back in meters and computed with D_o in inches. 

L_shell has free values

"""

"""
COULD BE VARIABLES BUT FIXED
----------------------------

Tube_pass = 2 # (or 1)

Tube_layout_angle = 45 # [°] (or 60, 90) : 45 and 90 mean square / 60 means triangular

Pitch_ratio constrainted to values depending on D_o, could be varied from 1.2 to 1.5 on its own
Square:
-------
D_o = 1/2     [in] => Pitch_ratio = 1.25
D_o = 5/8     [in] => Pitch_ratio = (7/8)/D_o
D_o = 3/4     [in] => Pitch_ratio = (1)/D_o
D_o = 1       [in] => Pitch_ratio = (1+1/4)/D_o
D_o = 1 + 1/4 [in] => Pitch_ratio = (1+9/16)/D_o
D_o = 1 + 1/2 [in] => Pitch_ratio = (1+7/8)/D_o

Triangular:
-----------
D_o = 1/2     [in] => Pitch_ratio = 1.25
D_o = 5/8     [in] => Pitch_ratio = (25/32)/D_o
D_o = 3/4     [in] => Pitch_ratio = (15/16)/D_o
D_o = 1       [in] => Pitch_ratio = (1+1/4)/D_o
D_o = 1 + 1/4 [in] => Pitch_ratio = (1+9/16)/D_o
D_o = 1 + 1/2 [in] => Pitch_ratio = (1+7/8)/D_o

Baffle_cut = 0.25 # Could be varied from 0.15 to 0.4 but 0.25 is usual value for liquid flow
"""

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")


import __init__

# Connector import
from connector.mass_connector import MassConnector

# Component import
from component.base_component import BaseComponent
from component.heat_exchanger.hex_MB_charge_sensitive import HeatExchangerMB

# Cost model import
from correlations.heat_exchanger.STHE_cost_estimation import HeatExchangerCost

# Shell and tube related toolbox
from toolbox.heat_exchangers.shell_and_tubes.pitch_ratio_shell_and_tube import pitch_ratio_fun
from toolbox.heat_exchangers.shell_and_tubes.estimate_tube_in_shell import estimate_number_of_tubes
from toolbox.heat_exchangers.shell_and_tubes.shell_toolbox import shell_thickness
from toolbox.heat_exchangers.shell_and_tubes.tubesheet_toolbox import tube_sheet_thickness
from toolbox.heat_exchangers.shell_and_tubes.baffle_toolbox import baffle_thickness, find_divisors_between_bounds

# Piping toolbox
from toolbox.piping.pipe_thickness import carbon_steel_pipe_thickness_mm

# External imports
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import pandas as pd
import random
import numpy as np
import copy
import sys

# Ignore convergence warnings
import warnings
warnings.filterwarnings('ignore')

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

#%%

class ShellAndTubeSizingOpt(BaseComponent):
    class Particle(BaseComponent):
        def __init__(self, params = {}, su_S = None, ex_S = None, su_T = None, ex_T = None, choice_vectors = None, P_max_cycle = None, T_max_cycle = None, H_htc_Corr = None, C_htc_Corr = None, H_DP_Corr = None, C_DP_Corr = None):
            super().__init__()

            self.params = copy.deepcopy(params)
            
            self.choice_vectors = choice_vectors
            
            # For tube thickness study
            self.P_max_cycle = P_max_cycle
            self.T_max_cycle = T_max_cycle

            self.position = None
            self.velocity = None
            self.unmoved = {}

            self.score = None
            self.Q = None
            self.DP_h = None
            self.DP_c = None

            self.personnal_best_position = None
            self.personnal_best_score = 1e20

            # Will be Mass connectors
            self.su_S = su_S
            self.ex_S = ex_S

            self.su_T = su_T
            self.ex_T = ex_T
        
            # Correlations 
            self.H_htc_Corr = H_htc_Corr
            self.C_htc_Corr = C_htc_Corr

            self.H_DP_Corr = H_DP_Corr
            self.C_DP_Corr = C_DP_Corr  
    
        def set_position(self, position):
            self.position = position            
            return
        
        def set_velocity(self, velocity):
            self.velocity = velocity
            return 
        
        def set_score(self, score):
            self.score = score
            
            if self.personnal_best_score - 0.1 > score:
                self.personnal_best_score = self.score
                self.personnal_best_position = self.position

                return 1 
            
            return 0

        def check_reinit(self):
            re_init_flag = True
            
            for opt_var in self.unmoved:
                if self.unmoved[opt_var] < 3:
                    re_init_flag = False
                    return re_init_flag
            return re_init_flag

        def compute_geom(self, opt_inputs):
            """
            Compute rest of geometry
            """

            # pitch_ratio
            pitch_ratio = pitch_ratio_fun(self.position['D_o_inch'], self.position['tube_layout'])

            # Pipe length
            L_tube = self.position['L_shell']
            
            # Cross Passes
            Cross_Passes = round(self.position['L_shell']/self.position['Central_spac']) - 1

            D_o = self.position['D_o_inch']*25.4*1e-3
            
            if self.params['Shell_Side'] == 'H':
                p_shell = opt_inputs['P_su_H']
            else:
                p_shell = opt_inputs['P_su_C']
            
            # Pipe Thickness
            Tube_t = carbon_steel_pipe_thickness_mm(self.position['D_o_inch']*25.4/1e3, self.T_max_cycle, p_shell, self.P_max_cycle)
                        
            # if self.position['D_o_inch'] == 1:
            #     Tube_t = pipe_thickness['1']
            # else:
            #     Tube_t = pipe_thickness[str(self.position['D_o_inch'])]
            self.Tube_t = Tube_t
            # Tube_t = D_o /10
            
            Shell_ID = self.position['Shell_ID_inch']*25.4*1e-3
            
            # Number of tubes 
            min_tubes_in_row = 4 # 8
            n_tubes = estimate_number_of_tubes(Shell_ID, D_o, pitch_ratio*D_o, self.position['tube_layout'], min_tubes_in_row)[0]

            # HT Area and HTX volumes
            A_eff = n_tubes*L_tube*np.pi*D_o
            
            T_V_tot = L_tube*n_tubes*np.pi*((D_o - 2*Tube_t)/2)**2

            T_V_out = np.pi*(D_o/2)**2*L_tube*n_tubes
            S_V_tot = self.position['L_shell']*np.pi*(Shell_ID/2)**2 - T_V_out
                
            self.set_parameters( 
                            A_eff = A_eff, S_V_tot = S_V_tot, Shell_ID = Shell_ID, T_V_tot = T_V_tot, Tube_L = L_tube, 
                            Tube_OD = D_o, Tube_t = Tube_t, central_spacing = self.position['Central_spac'], Tube_pass = self.position["Tube_pass"],
                            cross_passes = Cross_Passes, n_tubes = n_tubes, pitch_ratio = pitch_ratio, tube_layout = self.position["tube_layout"],
                            Baffle_cut = self.position["Baffle_cut"]
                            ) 
                        
            return

        def HeatTransferRate(self, opt_inputs, opt_params):
                        
            self.HX = HeatExchangerMB('Shell&Tube')

            # Set HX inputs
            self.HX.set_inputs(**opt_inputs)

            self.compute_geom(opt_inputs)

            # Set HX params
            self.HX.set_parameters(**opt_params)

            "Correlation Loading And Setting"

            Corr_H = self.H_htc_Corr
            Corr_C = self.C_htc_Corr
            
            Corr_H_DP = self.H_DP_Corr 
            Corr_C_DP = self.C_DP_Corr 

            self.HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

            "Parameters Setting"
                        
            self.HX.set_parameters(
                A_eff = self.params['A_eff'], Baffle_cut = self.position['Baffle_cut'], S_V_tot = self.params['S_V_tot'],
                Shell_ID = self.params['Shell_ID'], T_V_tot = self.params['T_V_tot'], Tube_L = self.params['Tube_L'], 
                Tube_OD = self.params['Tube_OD'], Tube_pass = self.params['Tube_pass'], Tube_t = self.params['Tube_t'],
                central_spacing = self.params['central_spacing'], cross_passes = self.params['cross_passes'], foul_s = self.params['foul_s'],
                foul_t = self.params['foul_t'], n_series = self.params['n_series'], n_tubes = self.params['n_tubes'], 
                pitch_ratio = self.params['pitch_ratio'], tube_cond = self.params['tube_cond'], tube_layout = self.params['tube_layout'],

                Shell_Side = self.params['Shell_Side'],

                Flow_Type = self.params['Flow_Type'], H_DP_ON = self.params['H_DP_ON'], C_DP_ON = self.params['C_DP_ON'], n_disc = self.params['n_disc']) 

            # HX.set_DP(DP_type="User-Defined", UD_H_DP=1e4, UD_C_DP=1e4)
            self.HX.set_DP(DP_type = "Correlation", Corr_H=Corr_H_DP, Corr_C=Corr_C_DP)    
            
            if self.HX.params['n_tubes'] == 0:
                self.Q = 0
                self.DP_h = 0
                self.DP_c = 0
                
                return 0, 0, 0
                        
            try:
                self.HX.solve()
                self.Q = self.HX.Q
                self.DP_h = self.HX.DP_h
                self.DP_c = self.HX.DP_c
    
                return self.HX.Q, self.HX.DP_h, self.HX.DP_c
    
            except:
                self.Q = 0
                self.DP_h = 0
                self.DP_c = 0
                return 0, 0, 0

    def __init__(self, seed = None):
        super().__init__()

        self.params = {}

        self.particles = None
        self.global_best_position = None
        self.global_best_score = None
        self.best_particle = None
        
        self.suitable_param_set = []

        # Optimization related parameters/variables
        self.opt_vars = {}
        self.bounds = {}
        self.choice_vectors = {}

        # For tube thickness study
        self.P_max_cycle = None
        self.T_max_cycle = None

        # Will be Mass connectors
        self.su_S = None
        self.ex_S = None

        self.su_T = None
        self.ex_T = None

        # Correlations 
        self.H_htc_Corr = None
        self.C_htc_Corr = None

        self.H_DP_Corr = None
        self.C_DP_Corr = None        
        
        self.global_patience_iter = 15     # stop/refresh if no improvement for N iterations
        self.global_min_delta = 1e-1       # improvement threshold to count as progress
        self._global_no_improve = 0        # internal counter
        
        # self.global_action = "refresh"     # "refresh" or "stop"
        # self.global_refresh_frac = 0.3     # when refreshing, fraction of worst particles to reinit

        if seed == None: # random seed // Input a chosen seed for replicability
            import os
            seed = int.from_bytes(os.urandom(4), "little")
            
        self.rng = np.random.default_rng(seed)   # one RNG for the whole optimizer
        
        self.seen = {}
        from threading import RLock
        self.seen_lock = RLock()
        
    def __getstate__(self):
        state = self.__dict__.copy()
        # locks and other unpicklables must not be serialized
        state['seen_lock'] = None
        return state
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        # recreate the lock in the child process
        from threading import RLock
        self.seen_lock = RLock()

    #%% 
    
    def pos_key(self, pos):
        # quantize helper
        def q(x, step):  # snap to a step (e.g., 1e-3 m)
            return float(round(float(x) / step) * step)
    
        # quantize continuous vars
        cs = q(pos['Central_spac'], 1e-2)   # 1 mm
        L  = q(pos['L_shell'],    1e-2)     # 1 mm
    
        # derived integer feature
        cross = int(round(L / cs)) - 1
    
        return (
            round(float(pos['D_o_inch']), 6),
            round(float(pos['Shell_ID_inch']), 6),
            int(round(float(pos['Tube_pass']))),
            int(round(float(pos['tube_layout']))),
            round(float(pos['Baffle_cut']), 1),  # 0.1% step
            cs,
            L,
            cross,
            # include knobs that change physics (optional but recommended):
            self.params.get('Shell_Side', 'H'),
            tuple(sorted((self.H_htc_Corr or {}).items())),
            tuple(sorted((self.C_htc_Corr or {}).items())),
            self.H_DP_Corr,
            self.C_DP_Corr,
        )


    #%% SETTERS
    
    def set_opt_vars(self, opt_vars):
        for opt_var in opt_vars:
            self.opt_vars[opt_var] = None 
        return

    def set_opt_vars_values(self, opt_vars_val):
        for key, value in opt_vars_val.items():
            if key in self.opt_vars:
                self.opt_vars[key] = value
            else:
                print(f"Key {key} is not an optimization variable.")
        return

    def set_bounds(self, bounds):
        for key, value in bounds.items():
            self.bounds[key] = value 
        return

    def set_choice_vectors(self, choice_vectors):
        for key, value in choice_vectors.items():
            self.choice_vectors[key] = value 
        return

    def set_max_cycle_prop(self, T_max_cycle = None, p_max_cycle = None):
        self.T_max_cycle = T_max_cycle
        self.P_max_cycle = p_max_cycle
        return
    
    def set_constraints(self, Q_dot = None, DP_h = None, DP_c = None):
        
        self.Q_dot_constr = Q_dot 
        self.DP_h_constr = DP_h
        self.DP_c_constr = DP_c
        
        return

    def set_corr(self, H_Corr, C_Corr, H_DP, C_DP):
        self.H_htc_Corr = H_Corr
        self.C_htc_Corr = C_Corr

        self.H_DP_Corr = H_DP
        self.C_DP_Corr = C_DP  
        return

    #%%
    
    def random_multiple(self, lower_bound, upper_bound, multiple):
        """
        Generate a random number that is a multiple of `multiple` within the range [lower_bound, upper_bound].

        Parameters:
        - lower_bound: The lower bound of the range.
        - upper_bound: The upper bound of the range.
        - multiple: The number which the generated number must be a multiple of.

        Returns:
        - A random number between lower_bound and upper_bound that is a multiple of `multiple`.
        """

        # L_shell shall be a multiple of the central spacing to get an integer value of cross_passes as a simplification 

        # Find the smallest multiple within the range
        start = np.ceil(lower_bound / multiple) * multiple
        end   = np.floor(upper_bound / multiple) * multiple
        
        if start > end:
            print(f"No multiples of {multiple} in the range [{lower_bound}, {upper_bound}]")
            return 0.0
        num = int(round((end - start) / multiple)) + 1
        k = self.rng.integers(0, num)  # [0, num-1]
        
        return float(start + k * multiple)

    #%% PARTICLE MANAGEMENT

    def init_particle(self, particle):
    
        particle_position = {}
        particle_velocity = {}
    
        for opt_var in self.opt_vars:
            particle_velocity[opt_var] = float(np.round(self.rng.uniform(-1, 1), 2))
            
        for key, vec in self.choice_vectors.items():
            vec = np.asarray(vec)
            if isinstance(vec[0], str):
                particle_position[key] = float(pd.to_numeric(self.rng.choice(vec)))
            else:
                particle_position[key] = float(self.rng.choice(vec))
        
        low_bound_central_spac = (particle_position['Shell_ID_inch']/5)*25.4*1e-3 # [m]
        high_bound_central_spac = (74*particle_position['D_o_inch']**(0.75))*25.4*1e-3 # [m]
        
        low_bound_L_shell = max(self.bounds['L_shell'][0], particle_position['Shell_ID_inch']*25.4*1e-3*3)
        high_bound_L_shell = min(self.bounds['L_shell'][-1], particle_position['Shell_ID_inch']*25.4*1e-3*15)
        
        particle_position['Central_spac'] = float(np.round(self.rng.uniform(low_bound_central_spac, high_bound_central_spac), 2))
        particle_position['L_shell'] = float(np.round(self.random_multiple(low_bound_L_shell, high_bound_L_shell, particle_position['Central_spac']), 2))
        particle_position['Baffle_cut'] = float(np.round(self.rng.uniform(self.bounds['Baffle_cut'][0], self.bounds['Baffle_cut'][1]), 2))
    
        # Put these positions in the particles 
    
        particle.set_position(particle_position)
        particle.set_velocity(particle_velocity)
    
        for opt_var in particle.position.keys():
            particle.unmoved[opt_var] = 0
    
        return 

    def clone_Particle(self, particle):
        # Create a new Particle instance
        new_particle = self.Particle(
            params=particle.params, su_S=particle.su_S, ex_S=particle.ex_S, su_T=particle.su_T, ex_T=particle.ex_T, choice_vectors=particle.choice_vectors,
            P_max_cycle=particle.P_max_cycle, T_max_cycle=particle.T_max_cycle, H_htc_Corr=particle.H_htc_Corr, C_htc_Corr = particle.C_htc_Corr,
            H_DP_Corr = particle.H_DP_Corr, C_DP_Corr = particle.C_DP_Corr
        )
        
        # Manually copy attributes (excluding AbstractState)
        new_particle.position = particle.position.copy()
        new_particle.velocity = particle.velocity.copy()
        new_particle.masses = particle.masses.copy()
        new_particle.Q = particle.Q
        new_particle.DP_h = particle.DP_h
        new_particle.DP_c = particle.DP_c
        new_particle.personnal_best_position = particle.personnal_best_position.copy()
        new_particle.personnal_best_score = particle.personnal_best_score

        if hasattr(particle, "HX"):
            new_particle.HX = copy.copy(particle.HX)
        else:
            new_particle.HX = HeatExchangerMB("Shell&Tube")
        
        # ⚠️ Do NOT copy the AbstractState object
        
        return new_particle

    #%% SCORE RELATED
    
    def HX_Mass(self, HX_params, HX):
        
        rho_carbon_steel = 7850 # kg/m^3
        
        T_shell_m = (HX.su_H.T + HX.su_C.T)/2
            
        "Shell Mass"
                
        shell_t = shell_thickness(HX_params['Shell_ID'], T_shell_m, self.P_max_cycle)        
        HX_params['t_S'] = shell_t
        
        Shell_OD = HX_params['Shell_ID'] + 2*shell_t       
        Shell_volume = np.pi*((Shell_OD/2)**2 - (HX_params['Shell_ID']/2)**2)*HX_params['Tube_L'] + shell_t*np.pi*Shell_OD**2/4 
        Shell_mass = Shell_volume*rho_carbon_steel
        
        "Tube Mass"
        
        T_mass = np.pi*((HX_params['Tube_OD']/2)**2 - ((HX_params['Tube_OD']-2*HX_params['Tube_t'])/2)**2)*HX_params['Tube_L']*HX_params['n_tubes']*rho_carbon_steel*HX_params['n_series']

        "Tube Sheet Mass"
        
        TS_t = tube_sheet_thickness(HX_params['Tube_OD'],HX_params['Tube_OD']*HX_params['pitch_ratio'], T_shell_m, self.P_max_cycle, HX_params["Shell_ID"]) # HX_params["Shell_ID"] assumed to be the gasket diameter
        Full_Tube_sheet_A = np.pi*(HX_params["Shell_ID"]/2)**2
        Tube_in_tube_sheet_A = HX_params["n_tubes"]*np.pi*(HX_params["Tube_OD"]/2)**2
        
        TS_mass = TS_t*(Full_Tube_sheet_A - Tube_in_tube_sheet_A)*rho_carbon_steel*2*HX_params['n_series']
        
        HX_params['t_TS'] = TS_t
        
        "Baffle Mass"
        if self.params['Shell_Side'] == 'H':
            rho_shell = HX.su_H.D
            T_shell = self.inputs['T_su_H']
        else:
            rho_shell = HX.su_C.D
            T_shell = self.inputs['T_su_C']
        
        B_t = baffle_thickness(HX_params["Shell_ID"], HX_params["Baffle_cut"]/100, rho_shell, T_shell)
        Full_Baffle_A = np.pi*(HX_params["Shell_ID"]/2)**2 * (1-HX_params["Baffle_cut"]/100)
        Tube_in_Baffle_A = HX_params["n_tubes"]*(1-HX_params["Baffle_cut"]/100)*np.pi*(HX_params["Tube_OD"]/2)**2

        B_mass = HX_params["cross_passes"] * B_t * (Full_Baffle_A - Tube_in_Baffle_A)*rho_carbon_steel*HX_params['n_series']
        
        HX_params['t_B'] = B_t
        
        return abs(T_mass + Shell_mass + TS_mass + B_mass), Shell_mass, T_mass, TS_mass, B_mass 

    def constraint_Q_dot(self, Q_particle):
        if self.Q_dot_constr == None:
            return 0
        else:
            return max(self.Q_dot_constr - Q_particle,0) # [W] 

    def constraint_DP_h(self, DP_h_particle):
        if self.DP_h_constr == None:
            return 0
        else:
            return max(DP_h_particle - self.DP_h_constr,0) # [Pa] 

    def constraint_DP_c(self, DP_c_particle):
        if self.DP_c_constr == None:
            return 0
        else:
            return max(DP_c_particle - self.DP_c_constr,0) # [Pa]

    #%%

    def _eval_particle_pure(self, particle, objective_function, penalty_factor):
        warnings.filterwarnings("ignore")
        
        """Compute (Q, DP_h, DP_c, total_score, masses) for one particle (no cache write)."""
        particle.HeatTransferRate(self.inputs, self.params)
    
        score, S_mass, T_mass, TS_mass, B_mass = objective_function(particle.HX.params, particle.HX)
        masses = {'Shell': S_mass, 'Tubes': T_mass, 'Tubesheet': TS_mass, 'Baffles': B_mass, 'Total': B_mass+TS_mass+T_mass+S_mass}
    
        pen  = max(self.Q_dot_constr - particle.Q, 0.0)
        pen += max(particle.DP_h - self.DP_h_constr, 0.0)
        pen += max(particle.DP_c - self.DP_c_constr, 0.0)
    
        total = score + penalty_factor * pen
        return (particle.Q, particle.DP_h, particle.DP_c, total, masses)
    
    def _apply_eval(self, particle, res):
        """Apply results to particle and update cache."""
        Q, DP_h, DP_c, total, masses = res
        particle.Q, particle.DP_h, particle.DP_c = Q, DP_h, DP_c
        particle.masses = masses
        particle.set_score(total)
        key = self.pos_key(particle.position)
        with self.seen_lock:
            self.seen[key] = res
    
    def evaluate_population_parallel(self, particles, objective_function, penalty_factor, n_jobs=-1, desc="eval"):
        """Evaluate a list of particles using threads + tqdm, honoring the cache."""
        results = [None] * len(particles)
        misses = []
    
        # cache hits
        for i, p in enumerate(particles):
            key = self.pos_key(p.position)
            with self.seen_lock:
                hit = self.seen.get(key)
            if hit is not None:
                results[i] = hit
            else:
                misses.append(i)
    
        # parallel compute cache misses
        if misses:
            outs = Parallel(n_jobs=n_jobs, backend="loky", batch_size="auto")(
                delayed(self._eval_particle_pure)(particles[i], objective_function, penalty_factor)
                for i in misses
            )
            for idx, out in zip(misses, outs):
                results[idx] = out
    
        # apply + store
        for i, p in enumerate(particles):
            self._apply_eval(p, results[i])
    
        # NEW: return number of new evaluations so we can display it
        return results, len(misses)


    def cost_estimation(self):
        """
        Technoeconomic optimization of superalloy supercritical CO2 microtube
        shell-and-tube-heat exchangers
        
        Akshay Bharadwaj Krishna, Kaiyuan Jin, Portonovo S. Ayyaswamy, Ivan Catton,
        Timothy S. Fisher
        
        Cost in $ of 2023
        """
        
        A = 1.2 # [$/kg] : For carbon steel pipes // 255 for superalloy piping
        B = 5 # [$ * m]
        C = 14 # [$]
        D = 2 # [$*m]
        E = 2 # [$]
        F = 4000 # [$]
        
        BP = self.best_particle
        HX_params = BP.HX.params
        
        n_U_tubes = HX_params['n_tubes']/HX_params['Tube_pass']
        
        A_term = A*BP.masses['Total']
        B_term = B*n_U_tubes/(HX_params['Tube_OD']*1000)
        C_term = C*n_U_tubes
        D_term = D*n_U_tubes/HX_params['Tube_L']
        E_term = E*(HX_params['Tube_L']/HX_params['central_spacing'])*n_U_tubes
        
        self.CAPEX = {'HX' : A_term + B_term + C_term + D_term + E_term + F}
        self.CAPEX['Install'] = self.CAPEX['HX']*0.35
        self.CAPEX['Total'] = self.CAPEX['HX'] + self.CAPEX['Install']

        return

    #%%


    
    def particle_swarm_optimization(self, objective_function, bounds, num_particles=30, num_dimensions=2, max_iterations=50, 
                                inertia_weight=0.4, cognitive_constant=1.5, social_constant=1.5, constraints=None,
                                penalty_factor=1000):

        # --- initialization (unchanged) ---
        self.particles = [self.Particle(params=self.params, su_S=self.su_S, ex_S=self.ex_S,
                                        su_T=self.su_T, ex_T=self.ex_T, choice_vectors=self.choice_vectors,
                                        P_max_cycle=self.P_max_cycle, T_max_cycle=self.T_max_cycle,
                                        H_htc_Corr=self.H_htc_Corr, C_htc_Corr=self.C_htc_Corr,
                                        H_DP_Corr=self.H_DP_Corr, C_DP_Corr=self.C_DP_Corr)
                          for _ in range(num_particles*3)]
    
        self.all_scores = np.zeros((num_particles, max_iterations+1))
        self.particles_all_pos = {opt_var: np.zeros((num_particles + 1, max_iterations)) for opt_var in self.opt_vars}
    
        for p in self.particles:
            self.init_particle(p)
    
        # initial evaluation
        _, n_new = self.evaluate_population_parallel(self.particles, objective_function, penalty_factor, desc="init")
    
        particle_scores = np.array([particle.personnal_best_score for particle in self.particles])
        best_particle_indices = np.argsort(particle_scores)[:num_particles]
        self.particles = [self.particles[i] for i in best_particle_indices]
    
        self.personal_best_positions = np.array([self.particles[i].personnal_best_position for i in range(len(self.particles))])
        self.personal_best_scores = np.array([self.particles[i].personnal_best_score for i in range(len(self.particles))])
    
        for i in range(len(self.particles)):
            self.all_scores[i][0] = self.personal_best_scores[i]
    
        self.global_best_score = min(self.personal_best_scores)
        self.best_particle = self.clone_Particle(self.particles[np.argmin(self.personal_best_scores)])
        self.global_best_position = self.best_particle.position
        self.global_best_Q  = self.best_particle.Q
        self.global_best_DP_h = self.best_particle.DP_h
        self.global_best_DP_c = self.best_particle.DP_c
        self.best_particle.compute_geom(self.inputs)
    
        cognitive_velocity = {}
        social_velocity = {}
    
        # --- ONE progress bar only ---
        pbar = tqdm(
            total=max_iterations,
            desc="PSO",
            unit="iter",
            dynamic_ncols=True,
            leave=False,              # bar disappears after close()
            file=sys.stdout,
            mininterval=0.1
        )
    
        total_new_evals = 0
    
        for iteration in range(max_iterations):
            for i in range(num_particles):
    
                flag = self.particles[i].check_reinit()
                if flag:
                    self.init_particle(self.particles[i])
    
                for opt_var in self.opt_vars:
                    self.particles_all_pos[opt_var][i][iteration] = self.particles[i].position[opt_var]
                    if iteration > 1:
                        pos = round(self.particles_all_pos[opt_var][i][iteration], 2)
                        pos_prev = round(self.particles_all_pos[opt_var][i][iteration-1], 2)
                        if opt_var in ('L_shell', 'Central_spac'):
                            self.particles[i].unmoved[opt_var] = (self.particles[i].unmoved[opt_var] + 1) if abs(pos - pos_prev) < 0.3 else 0
                        else:
                            self.particles[i].unmoved[opt_var] = (self.particles[i].unmoved[opt_var] + 1) if abs(pos - pos_prev) < 0.03 else 0
    
                    personal_best_position_val = self.particles[i].personnal_best_position[opt_var]
                    current_position_val = self.particles[i].position[opt_var]
                    global_best_position_val = self.global_best_position[opt_var]
    
                    self.particles_all_pos[opt_var][-1][iteration] = global_best_position_val
    
                    u1 = self.rng.random(); u2 = self.rng.random()
    
                    if opt_var in self.choice_vectors.keys():
                        scale = self.choice_vectors[opt_var][0]
                        r1 = u1 * scale; r2 = u2 * scale
                    elif opt_var in self.bounds.keys():
                        if opt_var == 'L_shell':
                            scale = self.bounds[opt_var][0] * 5
                        elif opt_var == 'Baffle_cut':
                            scale = self.bounds[opt_var][0] / 5
                        else:
                            scale = self.bounds[opt_var][0] * 3
                        r1 = u1 * scale; r2 = u2 * scale
    
                    cognitive_velocity[opt_var] = cognitive_constant * r1 * (personal_best_position_val - current_position_val)
                    social_velocity[opt_var]    = social_constant    * r2 * (global_best_position_val - current_position_val)
    
                    self.particles[i].velocity[opt_var] = (
                        inertia_weight * self.particles[i].velocity[opt_var] +
                        cognitive_velocity[opt_var] + social_velocity[opt_var]
                    )
                    self.particles[i].position[opt_var] += round(self.particles[i].velocity[opt_var], 3)
    
                    if opt_var == 'Central_spac':
                        low_bound_central_spac = (self.particles[i].position['Shell_ID_inch']/5)*25.4*1e-3
                        high_bound_central_spac = (74*self.particles[i].position['D_o_inch']**0.75)*25.4*1e-3
                        L_shell_divs = find_divisors_between_bounds(self.particles[i].position['L_shell'],
                                                                    low_bound_central_spac, high_bound_central_spac)
                        if len(L_shell_divs) == 0:
                            self.particles[i].position[opt_var] = round(low_bound_central_spac, 2)
                        else:
                            self.particles[i].position[opt_var] = min(L_shell_divs, key=lambda x: abs(x - self.particles[i].position[opt_var]))
    
                    # snap discrete variables
                    if self.choice_vectors and opt_var in self.choice_vectors.keys():
                        if isinstance(self.choice_vectors[opt_var][0], str):
                            allowed = pd.to_numeric(self.choice_vectors[opt_var])
                            new_val = min(allowed, key=lambda x: abs(x - pd.to_numeric(self.particles[i].position[opt_var])))
                            self.particles[i].position[opt_var] = str(new_val)
                        else:
                            allowed = self.choice_vectors[opt_var]
                            new_val = min(allowed, key=lambda x: abs(x - self.particles[i].position[opt_var]))
                            self.particles[i].position[opt_var] = new_val
    
                bound_flag = 0
                for bound_key in self.bounds:
                    if bound_key == 'L_shell':
                        lowL = max(self.bounds['L_shell'][0], self.particles[i].position['Shell_ID_inch']*25.4*1e-3*3)
                        highL = min(self.bounds['L_shell'][-1], self.particles[i].position['Shell_ID_inch']*25.4*1e-3*15)
                        if self.particles[i].position[bound_key] < lowL:
                            self.particles[i].position[bound_key] = lowL
                            if lowL == self.particles[i].position['Shell_ID_inch']*25.4*1e-3*3:
                                index = np.where(np.array(self.choice_vectors['Shell_ID_inch']) == self.particles[i].position['Shell_ID_inch'])[0]
                                self.particles[i].position['Shell_ID_inch'] = self.choice_vectors['Shell_ID_inch'][int(index-1)]
                            self.particles[i].velocity[bound_key] *= -0.5; bound_flag = 1
                        if self.particles[i].position[bound_key] > highL:
                            self.particles[i].position[bound_key] = highL
                            self.particles[i].velocity[bound_key] *= -0.5; bound_flag = 1
                    else:
                        if self.particles[i].position[bound_key] < self.bounds[bound_key][0]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][0]
                            self.particles[i].velocity[bound_key] *= -0.5; bound_flag = 1
                        if self.particles[i].position[bound_key] > self.bounds[bound_key][1]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][1]
                            self.particles[i].velocity[bound_key] *= -0.5; bound_flag = 1
    
            # mark bound violations (penalize after eval)
            if bound_flag == 1:
                self.particles[i].bound_violated = True
            else:
                self.particles[i].bound_violated = False
    
            # evaluate this iteration
            (_, n_new) = self.evaluate_population_parallel(
                self.particles[:num_particles],
                objective_function, penalty_factor, desc=f"it{iteration+1}"
            )
            total_new_evals += n_new

            for p in self.particles[:num_particles]:
                if getattr(p, "bound_violated", False):
                    p.score += 1e6
                    p.personnal_best_score = min(p.personnal_best_score, p.score)
    
            for i in range(num_particles):
                self.all_scores[i][iteration + 1] = self.particles[i].score
    
            self.personal_best_positions = np.array([self.particles[i].personnal_best_position for i in range(len(self.particles))])
            self.personal_best_scores = np.array([self.particles[i].personnal_best_score for i in range(len(self.particles))])
    
            new_best = min(self.personal_best_scores)
    
            if (self.global_best_score is None) or ((self.global_best_score - new_best) > self.global_min_delta):
                self.global_best_score = new_best
                self.best_particle = self.clone_Particle(self.particles[np.argmin(self.personal_best_scores)])
                self.global_best_position = self.best_particle.position
                self.global_best_Q  = self.best_particle.Q
                self.global_best_DP_h = self.best_particle.DP_h
                self.global_best_DP_c = self.best_particle.DP_c
                self.best_particle.compute_geom(self.inputs)
                self._global_no_improve = 0
            else:
                self._global_no_improve += 1
    
            # ---- update the single bar IN-PLACE (no prints) ----
            pbar.update(1)
            pbar.set_description(f"PSO {iteration+1}/{max_iterations}")
            pbar.set_postfix(
                best=f"{self.global_best_score:.3f}",
                Q=f"{self.global_best_Q:.1f}",
                DP_h=f"{self.global_best_DP_h:.0f}",
                DP_c=f"{self.global_best_DP_c:.0f}",
                new=n_new,
                total=total_new_evals
            )
    
            # early stop -> finish bar cleanly (no extra bars)
            if self._global_no_improve >= self.global_patience_iter:
                tqdm.write(
                    f"[PSO] Early stop at iteration {iteration+1}: "
                    f"no improvement > {self.global_min_delta} for "
                    f"{self.global_patience_iter} iterations."
                )
                break
    
        pbar.close()  # <- only once
    
        return self.global_best_position, self.global_best_score, self.best_particle

    def opt_size(self):
        
        import numexpr as ne
        import os, multiprocessing
        
        # Detect number of CPU cores
        num_cores = multiprocessing.cpu_count()
        
        # Choose a safe dynamic cap (e.g., min(physical cores, 16))
        # You can tune this — often 50–75% of total cores is optimal for large data
        num_threads = min(num_cores, 16)
        
        # Apply to NumExpr
        ne.set_num_threads(num_threads)
        os.environ["NUMEXPR_MAX_THREADS"] = str(num_threads)        

        self.particle_swarm_optimization(objective_function = self.HX_Mass , bounds = self.bounds, num_particles = 30, num_dimensions = len(self.opt_vars), max_iterations = 50, inertia_weight = 0.5,
                                          cognitive_constant = 0.5, social_constant = 0.5, constraints = [self.constraint_Q_dot, self.constraint_DP_h, self.constraint_DP_c], penalty_factor = 1)
        
        self._eval_particle_pure(self.best_particle, self.HX_Mass, 1)

        self.cost_calculator = HeatExchangerCost(
            D_S_i=self.best_particle.HX.params['Shell_ID'],  
            t_S=self.best_particle.HX.params['t_S'], 
            r=self.best_particle.HX.params['n_series'], 
            L_B=self.best_particle.HX.params['central_spacing'], 
            t_B=self.best_particle.HX.params['t_B'], 
            B_c = self.best_particle.HX.params['Baffle_cut']/100, 
            N_TS=2, 
            D_r=2*0.05/self.best_particle.HX.params['Shell_ID'], 
            L_T=self.best_particle.HX.params['Tube_L'], 
            D_T_o=self.best_particle.HX.params['Tube_OD'], 
            D_T_i=self.best_particle.HX.params['Tube_OD']-2*self.best_particle.Tube_t, 
            N_T=self.best_particle.HX.params['n_tubes'],
            pitch_r=self.best_particle.HX.params['pitch_ratio'], 
            L_CH=0.3, 
            N_CH=self.best_particle.HX.params['n_series']*2, 
            t_TS=self.best_particle.HX.params['t_TS'], 
            N_FL=self.best_particle.HX.params['n_series']*2, 
            t_FL=0.015, 
            t_RC=self.best_particle.HX.params['t_S']*2,
            D_SP_e=0.006*2, 
            D_SP_i=0.006, 
            N_TR=self.best_particle.HX.params['Tube_L']/0.72, 
            D_TR=0.006, 
            L_TR=self.best_particle.HX.params['Tube_L'], 
            N_Bt=100
        )

        self.costs = self.cost_calculator.calculate_total_cost()
        self.cost_estimation()
        
        print("\n")
        print(f"Best Position")
        print(f"-------------")
        print(f"D_tube_o : {round(self.best_particle.HX.params['Tube_OD']*1e3,2)} [mm]")
        print(f"D_shell : {round(self.best_particle.HX.params['Shell_ID'],3)} [m]") 
        print(f"Tube_pass : {self.best_particle.HX.params['Tube_pass']} [-]") 
        print(f"Tube_layout : {self.best_particle.HX.params['tube_layout']} [°]")
        print(f"Baffle Distance : {round(self.best_particle.HX.params['central_spacing'],3)} [m]")
        print(f"Shell/Tube length : {round(self.best_particle.HX.params['Tube_L'],3)} [m]")
        print(f"Baffle cut : {round(self.best_particle.HX.params['Baffle_cut'],1)} [%]")
        
        print("\n")
        print(f"Results")
        print(f"-------------")
        print(f"A_HX : {round(self.best_particle.HX.params['A_eff'],2)} [m^2]")
        print(f"Q_dot : {round(self.best_particle.HX.Q,1)} [W]") 
        print(f"DP_c : {round(self.best_particle.HX.DP_c,1)} [Pa]")
        print(f"DP_h : {round(self.best_particle.HX.DP_h,1)} [Pa]")
        print(f"m_HX : {round(self.best_particle.masses['Total'],1)} [kg]")
        print(f"Manufacturing costs est. : {round(self.costs,1)} [€ (2025)]")
        print(f"CAPEX : {round(self.CAPEX['Total'],1)} [$ (2025)]")
        
        return self.global_best_position, self.global_best_score, self.best_particle
    
if __name__ == "__main__":
    
    """
    Instanciate Optimizer and test case choice
    """
    
    HX_test = ShellAndTubeSizingOpt()
    test_case = "CO2_CD"
    
    if test_case == "Methanol":
    
        """
        Optimization related parameters/variables
        """
        
        HX_test.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'])
        
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2], # [1,2,4,6,8,10]
                            'tube_layout' : [0,45,60]}
        
        """
        'D_o_inch' : [0.5, 0.75, 1, 1.25, 1.5],
        'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78,
                                84, 90, 96, 108, 120]
        """
        
        HX_test.set_choice_vectors(choice_vectors)
        
        """
        Max T and P for pipe thickness computation
        """
        
        # Worst Case
        P_max_cycle = 10*1e5 # Pa
        T_max_cycle = 273.15+110 # K 
        
        HX_test.set_max_cycle_prop(T_max_cycle = T_max_cycle, p_max_cycle = P_max_cycle)
        
        """
        Thermodynamical parameters : Inlet and Outlet Design States
        """
    
        HX_test.set_inputs(
            # First fluid
            fluid_H = 'Methanol',
            T_su_H = 273.15 + 95, # K
            P_su_H = 10*1e5, # Pa
            m_dot_H = 27.8, # kg/s
    
            # Second fluid
            fluid_C = 'Water',
            T_su_C = 273.15 + 25, # K
            P_su_C = 5*1e5, # Pa
            m_dot_C = 68.9, # kg/s  # Make sure to include fluid information
            )
    
        "Constraints Values"
        Q_dot_cstr = 4.34*1e6
        DP_h_cstr = 13.2*1e3
        DP_c_cstr = 4.3*1e3
    
        """
        Parameters Setting
        """
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                # OPTI -> Oui (regarder le papier pour déterminer ça)
    
                                foul_t = 0.0002, # (m^2 * K/W)
                                foul_s = 0.00033, # (m^2 * K/W)
                                tube_cond = 50, # W/(m*K)
                                Overdesign = 0,
                                
                                Shell_Side = 'H',
    
                                Flow_Type = 'Shell&Tube',
                                H_DP_ON = True,
                                C_DP_ON = True,
                                n_disc = 1
                              )
    
        H_Corr = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
        C_Corr = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
        
        H_DP = "Shell_Kern_DP"
        C_DP = "Gnielinski_DP"
        
        HX_test.set_corr(H_Corr, C_Corr, H_DP, C_DP)
    
    elif test_case == "R134a":
    
        """
        Optimization related parameters/variables
        """
        
        HX_test.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'])
        
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2], # [1,2,4,6,8,10]
                            'tube_layout' : [60]} # [0,45,60]}
        
        """
        'D_o_inch' : [0.5, 0.75, 1, 1.25, 1.5],
        'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78,
                                84, 90, 96, 108, 120]
        """
        
        HX_test.set_choice_vectors(choice_vectors)
        
        """
        Max T and P for pipe thickness computation
        """
        
        # Worst Case
        P_max_cycle = 5*1e5 # Pa
        T_max_cycle = 273.15+110 # K 
        
        HX_test.set_max_cycle_prop(T_max_cycle = T_max_cycle, p_max_cycle = P_max_cycle)
        
        """
        Thermodynamical parameters : Inlet and Outlet Design States
        """
        
        HX_test.set_inputs(
            # First fluid
            fluid_H = 'Water',
            T_su_H = 26 + 273.15, # K
            P_su_H = 5*1e5, # Pa
            m_dot_H = 5.35, # kg/s
    
            # Second fluid
            fluid_C = 'R134a',
            h_su_C = PropsSI('H','T', 273.15+7,'Q',0,'R134a')+1, # J/kg
            P_su_C = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # Pa
            m_dot_C = 1.62, # kg/s  # Make sure to include fluid information
            )
        
        "Constraints Values"
        Q_dot_cstr = 0.313*1e6
        DP_h_cstr = 8.2*1e3
        DP_c_cstr = 21.7*1e3
    
        """
        Parameters Setting
        """
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                # OPTI -> Oui (regarder le papier pour déterminer ça)
    
                                foul_t = 0.000176, # (m^2 * K/W)
                                foul_s =  0.000176, # (m^2 * K/W)
                                tube_cond = 50, # W/(m*K)
                                Overdesign = 0,
                                
                                Shell_Side = 'H',
    
                                Flow_Type = 'Shell&Tube',
                                H_DP_ON = True,
                                C_DP_ON = True,
                                n_disc = 50
                              )
    
        H_Corr = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
        C_Corr = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
        
        H_DP = "Shell_Kern_DP"
        C_DP = "Muller_Steinhagen_Heck_DP"
        
        HX_test.set_corr(H_Corr, C_Corr, H_DP, C_DP)
    
    elif test_case == "CO2_CD":
        """
        Optimization related parameters/variables
        """
        
        HX_test.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'])
        
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2], # [1,2,4,6,8,10]
                            'tube_layout' : [0,45,60]} # [0,45,60]}
        
        """
        'D_o_inch' : [0.5, 0.75, 1, 1.25, 1.5],
        'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78,
                                84, 90, 96, 108, 120]
        """
        
        HX_test.set_choice_vectors(choice_vectors)
        
        """
        Max T and P for pipe thickness computation
        """
        
        # Worst Case
        P_max_cycle = 160*1e5 # Pa
        T_max_cycle = 273.15+140 # K 
        
        HX_test.set_max_cycle_prop(T_max_cycle = T_max_cycle, p_max_cycle = P_max_cycle)
        
        """
        Thermodynamical parameters : Inlet and Outlet Design States
        """
        
        HX_test.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 273.15 + 40, # K
            P_su_H = 5050000, # Pa
            m_dot_H = 30, # kg/s
    
            # Second fluid
            fluid_C = 'Water',
            T_su_C = 3 + 273.15, # K
            P_su_C = 5*1e5, # Pa
            m_dot_C = 150, # kg/s  # Make sure to include fluid information
            )
        
        "Constraints Values"
        Q_dot_cstr = 6295150
        DP_h_cstr = 3*1e5
        DP_c_cstr = 100*1e3
    
        """
        Parameters Setting
        """
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                # OPTI -> Oui (regarder le papier pour déterminer ça)
    
                                foul_t = 0.000176, # (m^2 * K/W)
                                foul_s =  0.000176, # (m^2 * K/W)
                                tube_cond = 50, # W/(m*K)
                                Overdesign = 0,
                                
                                Shell_Side = 'C',
    
                                Flow_Type = 'Shell&Tube',
                                H_DP_ON = True,
                                C_DP_ON = True,
                                n_disc = 30
                              )
        
        Corr_H = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Thome_Condensation"}
        # Corr_H = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Gnielinski"}
        Corr_C = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

        Corr_H_DP = "Choi_DP"
        Corr_C_DP = "Shell_Kern_DP"
        
        HX_test.set_corr(Corr_H, Corr_C, Corr_H_DP, Corr_C_DP)
    
    """
    Bounds and Constraints
    """
    
    bounds = {
                "L_shell" : [1,15], # 10],
                "D_o_inch" : [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
                "Shell_ID_inch" : [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]],
                "Tube_pass" : [choice_vectors['Tube_pass'][0], choice_vectors['Tube_pass'][-1]],
                "tube_layout" : [choice_vectors['tube_layout'][0], choice_vectors['tube_layout'][-1]],
                "Baffle_cut" : [15, 45]
                }
    
    HX_test.set_bounds(bounds)
    HX_test.set_constraints(Q_dot = Q_dot_cstr, DP_h = DP_h_cstr, DP_c = DP_c_cstr)
    
    # global_best_position, global_best_score, best_particle = HX_test.opt_size()
    
    import time
    time_1 = []
    
    # for i in range(10):
    t0 = time.perf_counter()
    
    global_best_position, global_best_score, best_particle = HX_test.opt_size()
    
    elapsed = time.perf_counter() - t0
    
    print(f"Optimization completed in {elapsed:.2f} s")
    
    time_1.append(elapsed)
    
    def HX_price(T_mass, Shell_mass, Baffle_mass, TS_mass):
        """
        Technoeconomic optimization of superalloy supercritical CO2 microtube
        shell-and-tube-heat exchangers (2023)
        
        Akshay Bharadwaj Krishna, Kaiyuan Jin, Portonovo S. Ayyaswamy, Ivan Catton,
        Timothy S. Fisher
        """
        
        A_pipes = 2 # €/kg : price of carbon steel seamless pipes : 
        A_plate = 0.5 # €/kg : price of carbon steel plate 
        
        material_costs = A_pipes*T_mass + A_plate*(Shell_mass + Baffle_mass + TS_mass)
        print(f"mat_cost : {material_costs}")
        
        # CAPEX = A*HX_Mass + B*(n_tubes/2)/tube_od_mm + C * (n_tubes/2) + D *(n_tubes/2)/tube_l_m + E*n_baffle*(n_tubes/2) + F
        
        return material_costs # CAPEX
    
    
