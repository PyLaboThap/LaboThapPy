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

# Connector import
from labothappy.connector.mass_connector import MassConnector

# Component import
from labothappy.component.base_component import BaseComponent
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive

# Cost model import
from labothappy.correlations.heat_exchanger.STHE_cost_estimation import HeatExchangerCost, total_STHE_cost, krishna_cost_correlation_STHE

# Shell and tube related toolbox
from labothappy.toolbox.heat_exchangers.shell_and_tubes.pitch_ratio_shell_and_tube import pitch_ratio_fun
from labothappy.toolbox.heat_exchangers.shell_and_tubes.estimate_tube_in_shell import estimate_number_of_tubes
from labothappy.toolbox.heat_exchangers.shell_and_tubes.shell_toolbox import shell_thickness
from labothappy.toolbox.heat_exchangers.shell_and_tubes.tubesheet_toolbox import tube_sheet_thickness
from labothappy.toolbox.heat_exchangers.shell_and_tubes.baffle_toolbox import baffle_thickness, find_divisors_between_bounds

# Piping toolbox
from labothappy.toolbox.piping.pipe_thickness import carbon_steel_pipe_thickness_mm

# Inflation
from labothappy.toolbox.economics.cpi_data import actualize_price

# External imports
from CoolProp.CoolProp import PropsSI
import pandas as pd
import random
import numpy as np
import copy

# Parallel evaluation (sizing_parallel)
from joblib import Parallel, delayed
from tqdm import tqdm

# Ignore convergence warnings
import warnings
warnings.filterwarnings('ignore')

def make_hashable(obj):
    if isinstance(obj, dict):
        return tuple(sorted((k, make_hashable(v)) for k, v in obj.items()))
    return obj

#%% Top-level standalone versions of HX_Mass / HX_total_cost
#    (must not depend on `self` so they can be called from a joblib worker
#     in a separate process, without pickling the whole optimizer instance)

def _hx_mass_standalone(HX, HX_params, P_max_cycle, inputs, params):
    """Version top-level de ShellAndTubeSizingOpt.HX_Mass, sans dépendance à self."""

    rho_carbon_steel = 7850  # kg/m^3

    T_shell_m = (HX.su_H.T + HX.su_C.T) / 2

    "Shell Mass"

    shell_t = shell_thickness(HX_params['Shell_ID'], T_shell_m, P_max_cycle)
    HX_params['t_S'] = shell_t

    Shell_OD = HX_params['Shell_ID'] + 2 * shell_t
    Shell_volume = np.pi * ((Shell_OD / 2) ** 2 - (HX_params['Shell_ID'] / 2) ** 2) * HX_params['Tube_L'] + shell_t * np.pi * Shell_OD ** 2 / 4
    Shell_mass = Shell_volume * rho_carbon_steel * HX_params['n_parallel']

    "Tube Mass"

    T_mass = np.pi * ((HX_params['Tube_OD'] / 2) ** 2 - ((HX_params['Tube_OD'] - 2 * HX_params['Tube_t']) / 2) ** 2) * HX_params['Tube_L'] * HX_params['n_tubes'] * rho_carbon_steel * HX_params['n_series'] * HX_params['n_parallel']

    "Tube Sheet Mass"

    TS_t = tube_sheet_thickness(HX_params['Tube_OD'], HX_params['Tube_OD'] * HX_params['pitch_ratio'], T_shell_m, P_max_cycle, HX_params["Shell_ID"])
    Full_Tube_sheet_A = np.pi * (HX_params["Shell_ID"] / 2) ** 2
    Tube_in_tube_sheet_A = HX_params["n_tubes"] * np.pi * (HX_params["Tube_OD"] / 2) ** 2

    TS_mass = TS_t * (Full_Tube_sheet_A - Tube_in_tube_sheet_A) * rho_carbon_steel * 2 * HX_params['n_series'] * HX_params['n_parallel']

    HX_params['t_TS'] = TS_t

    "Baffle Mass"
    if params['Shell_Side'] == 'H':
        rho_shell = HX.su_H.D
        T_shell = inputs['T_su_H']
    else:
        rho_shell = HX.su_C.D
        T_shell = inputs['T_su_C']

    B_t = baffle_thickness(HX_params["Shell_ID"], HX_params["Baffle_cut"] / 100, rho_shell, T_shell)
    Full_Baffle_A = np.pi * (HX_params["Shell_ID"] / 2) ** 2 * (1 - HX_params["Baffle_cut"] / 100)
    Tube_in_Baffle_A = HX_params["n_tubes"] * (1 - HX_params["Baffle_cut"] / 100) * np.pi * (HX_params["Tube_OD"] / 2) ** 2

    B_mass = HX_params["cross_passes"] * B_t * (Full_Baffle_A - Tube_in_Baffle_A) * rho_carbon_steel * HX_params['n_series'] * HX_params['n_parallel']

    HX_params['t_B'] = B_t

    return abs(T_mass + Shell_mass + TS_mass + B_mass), Shell_mass, T_mass, TS_mass, B_mass


def _hx_total_cost_standalone(HX, n_y=10, h_per_y=7000, C_e=0.12, i=0.1, eta_pp=0.8):
    """Version top-level de ShellAndTubeSizingOpt.HX_total_cost, sans dépendance à self."""

    # CAPEX : Hall Correlation
    a_1 = 8000
    a_2 = 259.2
    a_3 = 0.91

    CAPEX = a_1 + a_2 * (HX.params['A_eff']) ** a_3

    # OPEX of year n
    if HX.params['Shell_Side'] == 'H':
        DP_t = HX.DP_c
        DP_s = HX.DP_h

        mdot_t = HX.su_C.m_dot
        mdot_s = HX.su_H.m_dot

        rho_t = HX.su_C.D
        rho_s = HX.su_H.D
    else:
        DP_t = HX.DP_h
        DP_s = HX.DP_c

        mdot_t = HX.su_H.m_dot
        mdot_s = HX.su_C.m_dot

        rho_t = HX.su_H.D
        rho_s = HX.su_C.D

    P = 1 / eta_pp * (DP_t * mdot_t / rho_t + DP_s * mdot_s / rho_s)
    OPEX_n = P * C_e * h_per_y / 1000

    OPEX = 0
    for y in np.array(range(n_y)) + 1:
        OPEX += OPEX_n / ((1 + i) ** y)

    return OPEX + CAPEX

#%%

class ShellAndTubeSizingOpt(BaseComponent):

    _CORR_KEYS = {'H_Corr', 'C_Corr', 'H_DP', 'C_DP'}
    _CONSTRAINT_KEYS = {'Q_dot', 'DP_h', 'DP_c'}
    _MAX_CYCLE_KEYS = {'T_max_cycle', 'p_max_cycle'}
    _OPT_VAR_KEYS = {'opt_vars'}

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
    
            self.reject = 0
            self.Q_guess = None
            
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

        def compute_geom(self, opt_inputs, Tube_t_flag = True):
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
            if Tube_t_flag == True: # Pipe thickness is pressure dependent
                Tube_t = carbon_steel_pipe_thickness_mm(self.position['D_o_inch']*25.4/1e3, self.T_max_cycle, p_shell, self.P_max_cycle)
            else: # Reference assumption
                Tube_t = D_o/10
            
            self.Tube_t = Tube_t
            
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
                            Baffle_cut = self.position["Baffle_cut"], n_parallel = self.position['n_parallel']
                            ) 
                        
            return

        def HeatTransferRate(self, opt_inputs, opt_params):
                        
            self.HX = HexMBChargeSensitive('Shell&Tube')

            self.HX.set_inputs(**opt_inputs)
            
            if "Tube_t_flag" not in opt_params:
                opt_params['Tube_t_flag'] = True    
            
            self.compute_geom(opt_inputs, Tube_t_flag = opt_params['Tube_t_flag'])

            "Correlation Loading And Setting"

            Corr_H = self.H_htc_Corr
            Corr_C = self.C_htc_Corr
            
            Corr_H_DP = self.H_DP_Corr 
            Corr_C_DP = self.C_DP_Corr # "Gnielinski_DP"

            self.HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

            "Parameters Setting"
                        
            self.HX.set_parameters(
                A_eff = self.params['A_eff'], Baffle_cut = self.position['Baffle_cut'], S_V_tot = self.params['S_V_tot'],
                Shell_ID = self.params['Shell_ID'], T_V_tot = self.params['T_V_tot'], Tube_L = self.params['Tube_L'], 
                Tube_OD = self.params['Tube_OD'], Tube_pass = self.params['Tube_pass'], Tube_t = self.params['Tube_t'],
                central_spacing = self.params['central_spacing'], cross_passes = self.params['cross_passes'], foul_s = self.params['foul_s'],
                foul_t = self.params['foul_t'], n_series = self.params['n_series'], n_tubes = self.params['n_tubes'], n_parallel = self.params['n_parallel'],
                pitch_ratio = self.params['pitch_ratio'], tube_cond = self.params['tube_cond'], tube_layout = self.params['tube_layout'],

                Shell_Side = self.params['Shell_Side'],

                Flow_Type = self.params['Flow_Type'], H_DP_ON = self.params['H_DP_ON'], C_DP_ON = self.params['C_DP_ON'], n_disc = self.params['n_disc']) 

            if self.Q_guess is not None:
                self.HX.set_parameters(Q_guess = self.Q_guess)

            self.HX.set_DP(DP_type = "Correlation_Disc", Corr_H=Corr_H_DP, Corr_C=Corr_C_DP)    
            
            if self.HX.params['n_tubes'] == 0:
                self.Q = 0
                self.DP_h = 0
                self.DP_c = 0
                
                return 0, 0, 0
                        
            try:
                self.HX.solve()
                self.Q = self.HX.Q.Q_dot
                self.DP_h = self.HX.DP_h
                self.DP_c = self.HX.DP_c
    
                self.Q_guess = self.HX.Q.Q_dot
    
                return self.HX.Q.Q_dot, self.HX.DP_h, self.HX.DP_c
            
            except:
                
                self.reject+=1
                
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
        
        if seed == None: # random seed // Input a chosen seed for replicability
            import os
            seed = int.from_bytes(os.urandom(4), "little")
            
        self.rng = np.random.default_rng(seed)   # one RNG for the whole optimizer
        
        self.seen = {}       # cache of previous positions
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
            tuple(sorted((self.H_DP_Corr or {}).items())),
            tuple(sorted((self.C_DP_Corr or {}).items())),
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

    def set_bounds(self, bounds, choice_vectors=None):
        for key, value in bounds.items():
            self.bounds[key] = value
        if choice_vectors:
            self.set_choice_vectors(choice_vectors)
        return

    def set_choice_vectors(self, choice_vectors):
        for key, value in choice_vectors.items():
            self.choice_vectors[key] = value 
        return

    def set_max_cycle_prop(self, T_max_cycle = None, p_max_cycle = None):
        self.T_max_cycle = T_max_cycle
        self.P_max_cycle = p_max_cycle
        return

    def set_thermo_BC(self, su_S = None, ex_S = None, su_T = None, ex_T = None):
        self.su_S = su_S
        self.ex_S = ex_S

        self.su_T = su_T
        self.ex_T = ex_T
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
    
    def _apply_deferred_parameters(self):
        """
        set_parameters() reste inchangé (hérité de BaseComponent) et stocke
        tout dans self.params. Cette méthode en extrait les clés spéciales
        (corrélations, contraintes, propriétés max de cycle, opt_vars) et
        les dispatche vers les setters internes correspondants — appelée
        une fois, juste avant l'optimisation.
        """
        corr_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CORR_KEYS}
        constraint_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CONSTRAINT_KEYS}
        max_cycle_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._MAX_CYCLE_KEYS}
        opt_var_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._OPT_VAR_KEYS}
    
        if corr_kwargs:
            self.set_corr(**corr_kwargs)
        if constraint_kwargs:
            self.set_constraints(**constraint_kwargs)
        if max_cycle_kwargs:
            self.set_max_cycle_prop(**max_cycle_kwargs)
        if 'opt_vars' in opt_var_kwargs:
            self.set_opt_vars(opt_var_kwargs['opt_vars'])

    #%%
    def export_params_dict(self):
        
        return {
            "type": "Shell and Tube",
            "Flow_Type": self.params["Flow_Type"],
            
            "Q_dot": self.best_particle.Q,
            "DP_h": self.best_particle.DP_h,
            "DP_c": self.best_particle.DP_c,
            
            "fluid_H": self.inputs['fluid_H'],            
            "m_dot_H": self.inputs['m_dot_H'],
            "T_su_H": self.inputs['T_su_H'],
            "P_su_H": self.inputs['P_su_H'],
            
            "fluid_C": self.inputs['fluid_C'],            
            "m_dot_C": self.inputs['m_dot_C'],
            "T_su_C": self.inputs['T_su_C'],
            "P_su_C": self.inputs['P_su_C'],

            "n_series": self.best_particle.params["n_series"],
            "n_parallel": self.best_particle.params["n_parallel"],
            "foul_t": self.best_particle.params["foul_t"],
            "foul_s": self.best_particle.params["foul_s"],      
            
            "tube_cond": self.best_particle.params["tube_cond"],
            "Overdesign": self.best_particle.params["Overdesign"],
            "Shell_Side": self.best_particle.params["Shell_Side"],
            "S_V_tot": self.best_particle.params["S_V_tot"],
            "T_V_tot": self.best_particle.params["T_V_tot"],
            
            "A_eff": self.best_particle.params["A_eff"],
            "Shell_ID": self.best_particle.params["Shell_ID"],
            "Tube_L": self.best_particle.params["Tube_L"],
            "Tube_OD": self.best_particle.params["Tube_OD"],
            "Tube_t": self.best_particle.params["Tube_t"],
            "central_spacing": self.best_particle.params["central_spacing"],
            "Tube_pass": self.best_particle.params["Tube_pass"],
            "cross_passes": self.best_particle.params["cross_passes"],
            "n_tubes": self.best_particle.params["n_tubes"],
            "pitch_ratio": self.best_particle.params["pitch_ratio"],
            "tube_layout": self.best_particle.params["tube_layout"],
            "Baffle_cut": self.best_particle.params["Baffle_cut"],
            
            "CAPEX": self.CAPEX,
        } 
    
    #%%

    def random_multiple(self, lower_bound, upper_bound, multiple):
        """
        Generate a random number that is a multiple of `multiple` within the range [lower_bound, upper_bound].
        """

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
    
        particle.set_position(particle_position)
        particle.set_velocity(particle_velocity)
    
        for opt_var in particle.position.keys():
            particle.unmoved[opt_var] = 0
    
        return 

    def clone_Particle(self, particle):
        new_particle = self.Particle(
            params=particle.params, su_S=particle.su_S, ex_S=particle.ex_S, su_T=particle.su_T, ex_T=particle.ex_T, choice_vectors=particle.choice_vectors,
            P_max_cycle=particle.P_max_cycle, T_max_cycle=particle.T_max_cycle, H_htc_Corr=particle.H_htc_Corr, C_htc_Corr = particle.C_htc_Corr,
            H_DP_Corr = particle.H_DP_Corr, C_DP_Corr = particle.C_DP_Corr
        )
        
        new_particle.position = particle.position.copy()
        new_particle.velocity = particle.velocity.copy()
        
        if self.obj == 'masses':
            new_particle.masses = particle.masses.copy()
            
        new_particle.Q = particle.Q
        new_particle.penalty = particle.penalty
        new_particle.DP_h = particle.DP_h
        new_particle.DP_c = particle.DP_c
        new_particle.personnal_best_position = particle.personnal_best_position.copy()
        new_particle.personnal_best_score = particle.personnal_best_score
        new_particle.HX = copy.copy(particle.HX)
        
        return new_particle

    def _ensure_particle_HX(self, particle):
        """
        Après une évaluation parallèle, particle.HX n'existe pas sur l'objet du
        process principal (seul un dict de résultats a été renvoyé par le worker).
        Recalcule .HX localement, nécessaire avant compute_geom()/clone_Particle().
        """
        particle.HeatTransferRate(self.inputs, self.params)
        return particle

    #%% SCORE RELATED
    
    def HX_Mass(self, HX, HX_params):
        return _hx_mass_standalone(HX, HX_params, self.P_max_cycle, self.inputs, self.params)

    def HX_total_cost(self, HX, n_y=10, h_per_y=7000, C_e=0.12, i=0.1, eta_pp=0.8):
        return _hx_total_cost_standalone(HX, n_y, h_per_y, C_e, i, eta_pp)

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

    def evaluate_with_penalty(self, objective_function, particle, constraints, penalty_factor,i):
        """
        Evaluates the objective function with a penalty for constraint violations.
        """
        
        key = self.pos_key(particle.position)

        with self.seen_lock:
            hit = self.seen.get(key)
        if hit is not None:
            if self.obj == 'mass':
                Q, DP_h, DP_c, total, masses = hit
                particle.Q, particle.DP_h, particle.DP_c = Q, DP_h, DP_c
                particle.masses = masses
                particle.set_score(total)
                return total
            else:
                Q, DP_h, DP_c, total = hit
                particle.Q, particle.DP_h, particle.DP_c = Q, DP_h, DP_c
                particle.set_score(total)
                return total
                
        particle.HeatTransferRate(self.inputs, self.params)
        
        if self.obj == 'mass':
            score, S_mass, T_mass, TS_mass, B_mass = objective_function(particle.HX, particle.HX.params)
            
            particle.masses = {'Shell' : S_mass,
                               'Tubes' : T_mass,
                               'Tubesheet' : TS_mass,
                               'Baffles' : B_mass,
                               'Total' : score}
    
            particle.penalty = 0
        
            constraints = [
                self.constraint_Q_dot(particle.Q),
                self.constraint_DP_h(particle.DP_h),
                self.constraint_DP_c(particle.DP_c)
            ]
        
            particle.penalty = sum(penalty_factor * abs(value) for value in constraints)
        
        else:
            score = objective_function(particle.HX)
            
            particle.penalty = 0

            constraints = [
                self.constraint_Q_dot(particle.Q),
            ]
            
            particle.penalty = sum(penalty_factor * abs(value) for value in constraints)
            
        total_score = score + particle.penalty
        flag = particle.set_score(total_score)
    
        if particle.penalty == 0 and particle.params not in self.suitable_param_set:
            self.suitable_param_set.append(copy.deepcopy(particle.params))
            self.suitable_param_set[-1]['score'] = particle.score
        
        with self.seen_lock:
            if self.obj == 'mass':
                self.seen[key] = (particle.Q, particle.DP_h, particle.DP_c, total_score, particle.masses)
            else:
                self.seen[key] = (particle.Q, particle.DP_h, particle.DP_c, total_score)
                
        return total_score

    #%%

    def particle_swarm_optimization(self, objective_function, bounds, num_particles=30, num_dimensions=2, max_iterations=50, 
                                inertia_weight=0.4, cognitive_constant=1.5, social_constant=1.5, constraints = None,
                                penalty_factor=1):
        """
        Perform Particle Swarm Optimization (PSO) to minimize the given objective function with constraints.
        """

        self.particles = [self.Particle(params = self.params, su_S = self.su_S, ex_S = self.ex_S, su_T = self.su_T, ex_T = self.ex_T, choice_vectors = self.choice_vectors, 
                                        P_max_cycle = self.P_max_cycle, T_max_cycle = self.T_max_cycle, H_htc_Corr=self.H_htc_Corr, C_htc_Corr = self.C_htc_Corr,
                                        H_DP_Corr = self.H_DP_Corr, C_DP_Corr = self.C_DP_Corr) for _ in range(num_particles*3)]

        self.all_scores = np.zeros((num_particles, max_iterations+1))

        self.particles_all_pos = {}

        for opt_var in self.opt_vars:
            self.particles_all_pos[opt_var] = np.zeros((num_particles + 1, max_iterations))

        for i in range(len(self.particles)):
            self.init_particle(self.particles[i])
            score = self.evaluate_with_penalty(objective_function, self.particles[i], constraints, penalty_factor,i)

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
        self.global_best_Q = self.best_particle.Q
        self.global_best_DP_h = self.best_particle.DP_h
        self.global_best_DP_c = self.best_particle.DP_c     

        self.best_particle.compute_geom(self.inputs, Tube_t_flag = self.params['Tube_t_flag'])

        cognitive_velocity = {}
        social_velocity = {}

        for iteration in range(max_iterations):

            for i in range(num_particles):

                flag = self.particles[i].check_reinit()
                if flag:
                    self.init_particle(self.particles[i])

                for opt_var in self.opt_vars:
                    self.particles_all_pos[opt_var][i][iteration] = self.particles[i].position[opt_var]
                    if iteration > 1:
                        pos = round(self.particles_all_pos[opt_var][i][iteration],2)
                        pos_previous = round(self.particles_all_pos[opt_var][i][iteration-1],2)
                        
                        if opt_var == 'L_shell' or opt_var == 'Central_spac':
                            if abs(pos - pos_previous) < 0.3:
                                self.particles[i].unmoved[opt_var] += 1
                            else:
                                self.particles[i].unmoved[opt_var] = 0
                        else:
                            if abs(pos - pos_previous) < 0.03:
                                self.particles[i].unmoved[opt_var] += 1
                            else:
                                self.particles[i].unmoved[opt_var] = 0                            

                    personal_best_position_val = self.particles[i].personnal_best_position[opt_var]
                    current_position_val = self.particles[i].position[opt_var]
                    global_best_position_val = self.global_best_position[opt_var]

                    self.particles_all_pos[opt_var][-1][iteration] = global_best_position_val
    
                    u1 = self.rng.random()
                    u2 = self.rng.random()
                    
                    if opt_var in self.choice_vectors.keys():
                        scale = self.choice_vectors[opt_var][0]
                        r1 = u1 * scale
                        r2 = u2 * scale
                    elif opt_var in self.bounds.keys():
                        if opt_var == 'L_shell':
                            scale = self.bounds[opt_var][0] * 5
                        elif opt_var == 'Baffle_cut':
                            scale = 0.02
                        else:
                            scale = self.bounds[opt_var][0] * 3
                        
                        r1 = u1 * scale
                        r2 = u2 * scale
    
                    cognitive_fact = cognitive_constant
                    social_fact = social_constant

                    cognitive_velocity[opt_var] = cognitive_fact * r1 * (personal_best_position_val - current_position_val)
                    social_velocity[opt_var] = social_fact * r2 * (global_best_position_val - current_position_val)
                    
                    self.particles[i].velocity[opt_var] = (inertia_weight * self.particles[i].velocity[opt_var] + cognitive_velocity[opt_var] + social_velocity[opt_var])

                    self.particles[i].position[opt_var] += round(self.particles[i].velocity[opt_var],3)
                    
                    if opt_var == 'Central_spac':
                        low_bound_central_spac = (self.particles[i].position['Shell_ID_inch']/5)*25.4*1e-3 # [m]
                        high_bound_central_spac = (74*self.particles[i].position['D_o_inch']**(0.75))*25.4*1e-3 # [m]
                        
                        L_shell_divisors = find_divisors_between_bounds(self.particles[i].position['L_shell'],low_bound_central_spac,high_bound_central_spac)

                        if len(L_shell_divisors) == 0:
                            self.particles[i].position[opt_var] = round(low_bound_central_spac,2)
                        else:
                            self.particles[i].position[opt_var] = min(L_shell_divisors, key=lambda x: abs(x - self.particles[i].position[opt_var]))
                        
                    if self.choice_vectors and opt_var in self.choice_vectors.keys():

                        if isinstance(self.choice_vectors[opt_var][0],str):
                            allowed_values = pd.to_numeric(self.choice_vectors[opt_var])
                            new_position_value = min(allowed_values, key=lambda x: abs(x - pd.to_numeric(self.particles[i].position[opt_var])))
                            self.particles[i].position[opt_var] = str(new_position_value)
                        else: 
                            allowed_values = self.choice_vectors[opt_var]
                            new_position_value = min(allowed_values, key=lambda x: abs(x - self.particles[i].position[opt_var]))
                            self.particles[i].position[opt_var] = new_position_value

                bound_flag = 0

                for bound_key in self.bounds:
                    
                    if bound_key == 'L_shell':
                        low_bound_L_shell = max(self.bounds['L_shell'][0], self.particles[i].position['Shell_ID_inch']*25.4*1e-3*3)
                        high_bound_L_shell = min(self.bounds['L_shell'][-1], self.particles[i].position['Shell_ID_inch']*25.4*1e-3*15)
                    
                        if self.particles[i].position[bound_key] < low_bound_L_shell:
                            self.particles[i].position[bound_key] = low_bound_L_shell
                            if low_bound_L_shell == self.particles[i].position['Shell_ID_inch']*25.4*1e-3*3:
                                index = np.where(np.array(self.choice_vectors['Shell_ID_inch']) == self.particles[i].position['Shell_ID_inch'])[0]
                                self.particles[i].position['Shell_ID_inch'] = self.choice_vectors['Shell_ID_inch'][int(index-1)]
                                self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key]
                            else:
                                self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key]
                                
                            bound_flag = 1
                            
                        if self.particles[i].position[bound_key] > high_bound_L_shell:
                            self.particles[i].position[bound_key] = high_bound_L_shell
                            self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key]
                            bound_flag = 1
                    
                    else:
                        if self.particles[i].position[bound_key] < self.bounds[bound_key][0]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][0]
                            self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key]
                            bound_flag = 1
    
                        if self.particles[i].position[bound_key] > self.bounds[bound_key][1]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][1]
                            self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key]
                            bound_flag = 1

                if bound_flag == 1:
                    new_score = self.evaluate_with_penalty(objective_function, self.particles[i], constraints, penalty_factor,i) + 1e6
                else:
                    new_score = self.evaluate_with_penalty(objective_function, self.particles[i], constraints, penalty_factor,i)

            for i in range(num_particles):
                self.all_scores[i][iteration + 1] = self.particles[i].score

            self.personal_best_positions = np.array([self.particles[i].personnal_best_position for i in range(len(self.particles))])
            self.personal_best_scores = np.array([self.particles[i].personnal_best_score for i in range(len(self.particles))])

            new_pot_global_best_score = min(self.personal_best_scores)

            if new_pot_global_best_score + 0.1 < self.global_best_score:
                self.global_best_score = new_pot_global_best_score
                self.best_particle = self.clone_Particle(self.particles[np.argmin(self.personal_best_scores)])
                self.global_best_position = self.best_particle.position
                self.global_best_Q = self.best_particle.Q
                self.global_best_DP_h = self.best_particle.DP_h
                self.global_best_DP_c = self.best_particle.DP_c

                self.best_particle.compute_geom(self.inputs, Tube_t_flag = self.params['Tube_t_flag'])

            if self.print_flag: 
                print("===========================")
                print(f"Iteration {iteration+1}/{max_iterations}, Global Best Score: {self.global_best_score}, Related Q: {self.global_best_Q}")
                print(f"Related DP_h: {self.global_best_DP_h}, Related DP_c: {self.global_best_DP_c}")
                print(f"Best Position : {self.global_best_position}")
                print(f"Best Part Velocity : {self.best_particle.velocity}")
                
        return self.global_best_position, self.global_best_score, self.best_particle
    
    def cost_estimation(self):
        
        m_HX = self.best_particle.masses['Total']
        n_tubes = self.best_particle.HX.params['n_tubes']
        tube_L = self.best_particle.HX.params['Tube_L']
        tube_OD = self.best_particle.HX.params['Tube_OD']

        n_B = np.round(tube_L/self.best_particle.HX.params['central_spacing']) - 1

        capex = krishna_cost_correlation_STHE(m_HX, n_B, n_tubes, tube_L, tube_OD)
        
        self.CAPEX = {"HX" : actualize_price(capex, 2023, "USD"),
                      "Currency" : "USD"}
        
        self.CAPEX["Install"] = self.CAPEX["HX"]*0.35
        self.CAPEX["Total"] = self.CAPEX["HX"] + self.CAPEX["Install"]
        
        return
    
    def opt_size(self, n_particles = 50, max_iter = 50, obj = 'mass', print_flag = 0):
        self._apply_deferred_parameters()

        self.obj = obj
        self.print_flag = print_flag
        
        if obj == 'mass':
            self.particle_swarm_optimization(objective_function = self.HX_Mass , bounds = self.bounds, num_particles = n_particles, num_dimensions = len(self.opt_vars), max_iterations = max_iter, inertia_weight = 0.5,
                                          cognitive_constant = 0.5, social_constant = 0.5, constraints = [self.constraint_Q_dot, self.constraint_DP_h, self.constraint_DP_c], penalty_factor = 1e6)
        else:
            self.particle_swarm_optimization(objective_function = self.HX_total_cost , bounds = self.bounds, num_particles = n_particles, num_dimensions = len(self.opt_vars), max_iterations = max_iter, inertia_weight = 0.5,
                                          cognitive_constant = 0.5, social_constant = 0.5, constraints = [self.constraint_Q_dot, self.constraint_DP_h, self.constraint_DP_c], penalty_factor = 1e6)

        self.best_particle.compute_geom(self.inputs, Tube_t_flag = self.params['Tube_t_flag'])
        
        total, S_mass, T_mass, TS_mass, B_mass = self.HX_Mass(self.best_particle.HX, self.best_particle.HX.params)
    
        self.best_particle.masses = {'Shell' : S_mass,
                           'Tubes' : T_mass,
                           'Tubesheet' : TS_mass,
                           'Baffles' : B_mass,
                           'Total' : total}
        
        self.best_particle.total_cost = self.HX_total_cost(self.best_particle.HX)

        if obj == 'mass':
            
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
    
            self.manuf_cost = self.cost_calculator.calculate_total_cost()
            self.cost_estimation()
            
        self.reject = 0
        for part in self.particles:
            self.reject += part.reject
        
        return self.global_best_position, self.global_best_score, self.best_particle

    #%% PARALLEL SIZING

    def sizing_parallel(self, n_jobs=-1, backend="loky", n_particles=30, max_iterations=50,
                         inertia_weight=0.5, cognitive_constant=0.5, social_constant=0.5,
                         penalty_factor=1e6, obj='mass', print_flag=0, show_progress=True):
        """
        Version parallélisée de particle_swarm_optimization/opt_size.
        Même logique PSO custom (bornes dépendantes, snapping discret,
        réinitialisation de particules stagnantes) mais l'évaluation
        HeatTransferRate (coûteuse, CoolProp) est distribuée sur plusieurs
        process via joblib. Le cache self.seen n'est pas partagé entre
        process : chaque particule est réévaluée à chaque itération.
        """
        self._apply_deferred_parameters()

        self.obj = obj
        self.print_flag = print_flag

        snapshot = {
            'params': dict(self.params),
            'inputs': dict(self.inputs),
            'choice_vectors': dict(self.choice_vectors),
            'P_max_cycle': self.P_max_cycle,
            'T_max_cycle': self.T_max_cycle,
            'H_htc_Corr': self.H_htc_Corr, 'C_htc_Corr': self.C_htc_Corr,
            'H_DP_Corr': self.H_DP_Corr, 'C_DP_Corr': self.C_DP_Corr,
            'Q_dot_constr': self.Q_dot_constr, 'DP_h_constr': self.DP_h_constr, 'DP_c_constr': self.DP_c_constr,
            'obj': obj,
        }

        def batch_evaluate(particles):
            results = Parallel(n_jobs=n_jobs, backend=backend)(
                delayed(_eval_particle_sht)(p.position, snapshot) for p in particles
            )
            for p, r in zip(particles, results):
                p.Q, p.DP_h, p.DP_c = r['Q'], r['DP_h'], r['DP_c']
                p.penalty = r['penalty']
                if r['masses'] is not None:
                    p.masses = r['masses']
                p.set_score(r['score'])
            return results

        # --- Init (3x particles, garder les meilleures) ---
        self.particles = [self.Particle(
            params=self.params, su_S=self.su_S, ex_S=self.ex_S, su_T=self.su_T, ex_T=self.ex_T,
            choice_vectors=self.choice_vectors, P_max_cycle=self.P_max_cycle, T_max_cycle=self.T_max_cycle,
            H_htc_Corr=self.H_htc_Corr, C_htc_Corr=self.C_htc_Corr,
            H_DP_Corr=self.H_DP_Corr, C_DP_Corr=self.C_DP_Corr,
        ) for _ in range(n_particles * 3)]

        for p in self.particles:
            self.init_particle(p)

        batch_evaluate(self.particles)

        particle_scores = np.array([p.personnal_best_score for p in self.particles])
        best_idx = np.argsort(particle_scores)[:n_particles]
        self.particles = [self.particles[i] for i in best_idx]
        
        self.personal_best_positions = np.array([p.personnal_best_position for p in self.particles])
        self.personal_best_scores = np.array([p.personnal_best_score for p in self.particles])
        
        self.global_best_score = min(self.personal_best_scores)
        winner = self.particles[np.argmin(self.personal_best_scores)]
        self._ensure_particle_HX(winner)
        self.best_particle = self.clone_Particle(winner)
        self.global_best_position = self.best_particle.position
        self.global_best_Q, self.global_best_DP_h, self.global_best_DP_c = (
            self.best_particle.Q, self.best_particle.DP_h, self.best_particle.DP_c
        )
        self.best_particle.compute_geom(self.inputs, Tube_t_flag=self.params['Tube_t_flag'])

        cognitive_velocity, social_velocity = {}, {}

        pbar = tqdm(total=max_iterations, desc="Shell&Tube", unit="iter") if show_progress else None

        for iteration in range(max_iterations):

            # --- Passe 1 : update des positions (séquentiel, comme particle_swarm_optimization) ---
            for i in range(n_particles):
                if self.particles[i].check_reinit():
                    self.init_particle(self.particles[i])
                    continue  # réévaluée avec le reste au batch ci-dessous

                for opt_var in self.opt_vars:
                    personal_best_val = self.particles[i].personnal_best_position[opt_var]
                    current_val = self.particles[i].position[opt_var]
                    global_best_val = self.global_best_position[opt_var]

                    u1, u2 = self.rng.random(), self.rng.random()

                    if opt_var in self.choice_vectors:
                        scale = self.choice_vectors[opt_var][0]
                    elif opt_var == 'L_shell':
                        scale = self.bounds[opt_var][0] * 5
                    elif opt_var == 'Baffle_cut':
                        scale = 0.02
                    elif opt_var == 'Central_spac':
                        low_bound_central_spac = (self.particles[i].position['Shell_ID_inch']/5)*25.4*1e-3
                        scale = low_bound_central_spac * 3
                    else:
                        scale = self.bounds[opt_var][0] * 3

                    cognitive_velocity[opt_var] = cognitive_constant * u1 * scale * (personal_best_val - current_val)
                    social_velocity[opt_var] = social_constant * u2 * scale * (global_best_val - current_val)

                    self.particles[i].velocity[opt_var] = (
                        inertia_weight * self.particles[i].velocity[opt_var]
                        + cognitive_velocity[opt_var] + social_velocity[opt_var]
                    )
                    self.particles[i].position[opt_var] += round(self.particles[i].velocity[opt_var], 3)

                    if opt_var == 'Central_spac':
                        low = (self.particles[i].position['Shell_ID_inch']/5)*25.4*1e-3
                        high = (74*self.particles[i].position['D_o_inch']**0.75)*25.4*1e-3
                        divisors = find_divisors_between_bounds(self.particles[i].position['L_shell'], low, high)
                        self.particles[i].position[opt_var] = (
                            round(low, 2) if not divisors
                            else min(divisors, key=lambda x: abs(x - self.particles[i].position[opt_var]))
                        )

                    if opt_var in self.choice_vectors:
                        vec = self.choice_vectors[opt_var]
                        if isinstance(vec[0], str):
                            allowed = pd.to_numeric(vec)
                            val = min(allowed, key=lambda x: abs(x - pd.to_numeric(self.particles[i].position[opt_var])))
                            self.particles[i].position[opt_var] = str(val)
                        else:
                            self.particles[i].position[opt_var] = min(vec, key=lambda x: abs(x - self.particles[i].position[opt_var]))

                for bound_key in self.bounds:
                    if bound_key == 'L_shell':
                        low = max(self.bounds['L_shell'][0], self.particles[i].position['Shell_ID_inch']*25.4*1e-3*3)
                        high = min(self.bounds['L_shell'][-1], self.particles[i].position['Shell_ID_inch']*25.4*1e-3*15)
                        if self.particles[i].position[bound_key] < low:
                            self.particles[i].position[bound_key] = low
                            self.particles[i].velocity[bound_key] *= -0.5
                        if self.particles[i].position[bound_key] > high:
                            self.particles[i].position[bound_key] = high
                            self.particles[i].velocity[bound_key] *= -0.5
                    else:
                        if self.particles[i].position[bound_key] < self.bounds[bound_key][0]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][0]
                            self.particles[i].velocity[bound_key] *= -0.5
                        if self.particles[i].position[bound_key] > self.bounds[bound_key][1]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][1]
                            self.particles[i].velocity[bound_key] *= -0.5

            # --- Passe 2 : évaluation en parallèle ---
            batch_evaluate(self.particles)

            self.personal_best_positions = np.array([p.personnal_best_position for p in self.particles])
            self.personal_best_scores = np.array([p.personnal_best_score for p in self.particles])

            new_best = min(self.personal_best_scores)
            if new_best + 0.1 < self.global_best_score:
                self.global_best_score = new_best
                winner = self.particles[np.argmin(self.personal_best_scores)]
                self._ensure_particle_HX(winner)
                self.best_particle = self.clone_Particle(winner)
                self.global_best_position = self.best_particle.position
                self.global_best_Q, self.global_best_DP_h, self.global_best_DP_c = (
                    self.best_particle.Q, self.best_particle.DP_h, self.best_particle.DP_c
                )
                self.best_particle.compute_geom(self.inputs, Tube_t_flag=self.params['Tube_t_flag'])

            if pbar is not None:
                pbar.set_postfix_str(f"best={self.global_best_score:.4f}")
                pbar.update(1)
            elif self.print_flag:
                print(f"Iter {iteration+1}/{max_iterations} — Best: {self.global_best_score:.4f}")

        if pbar is not None:
            pbar.close()

        self.best_particle.compute_geom(self.inputs, Tube_t_flag=self.params['Tube_t_flag'])
        total, S_mass, T_mass, TS_mass, B_mass = self.HX_Mass(self.best_particle.HX, self.best_particle.HX.params)
        self.best_particle.masses = {'Shell': S_mass, 'Tubes': T_mass, 'Tubesheet': TS_mass, 'Baffles': B_mass, 'Total': total}
        self.best_particle.total_cost = self.HX_total_cost(self.best_particle.HX)
        
        if obj == 'mass':
        
            self.cost_calculator = HeatExchangerCost(
                D_S_i=self.best_particle.HX.params['Shell_ID'],
                t_S=self.best_particle.HX.params['t_S'],
                r=self.best_particle.HX.params['n_series'],
                L_B=self.best_particle.HX.params['central_spacing'],
                t_B=self.best_particle.HX.params['t_B'],
                B_c=self.best_particle.HX.params['Baffle_cut']/100,
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
        
            self.manuf_cost = self.cost_calculator.calculate_total_cost()
            self.cost_estimation()
        
        return self.global_best_position, self.global_best_score, self.best_particle

#%% Worker top-level (nécessite ShellAndTubeSizingOpt.Particle déjà défini)

def _eval_particle_sht(position, snapshot):
    """Évalue une position de particule Shell&Tube (top-level, picklable par joblib)."""
    warnings.filterwarnings('ignore')

    particle = ShellAndTubeSizingOpt.Particle(
        params=snapshot['params'],
        choice_vectors=snapshot['choice_vectors'],
        P_max_cycle=snapshot['P_max_cycle'],
        T_max_cycle=snapshot['T_max_cycle'],
        H_htc_Corr=snapshot['H_htc_Corr'], C_htc_Corr=snapshot['C_htc_Corr'],
        H_DP_Corr=snapshot['H_DP_Corr'], C_DP_Corr=snapshot['C_DP_Corr'],
    )
    particle.set_position(position)

    try:
        particle.HeatTransferRate(snapshot['inputs'], snapshot['params'])
    except Exception:
        return {'Q': 0.0, 'DP_h': 0.0, 'DP_c': 0.0, 'score': 1e12, 'penalty': 1e12, 'masses': None}

    def constraint(val, constr, direction):
        if constr is None:
            return 0
        return max((constr - val) if direction == 'min' else (val - constr), 0)

    c_Q   = constraint(particle.Q,    snapshot['Q_dot_constr'], 'min')
    c_DPh = constraint(particle.DP_h, snapshot['DP_h_constr'],  'max')
    c_DPc = constraint(particle.DP_c, snapshot['DP_c_constr'],  'max')
    penalty = 1e6 * (abs(c_Q) + abs(c_DPh) + abs(c_DPc))

    masses = None
    if snapshot['obj'] == 'mass':
        score_mass, S_mass, T_mass, TS_mass, B_mass = _hx_mass_standalone(
            particle.HX, particle.HX.params, snapshot['P_max_cycle'], snapshot['inputs'], snapshot['params']
        )
        masses = {'Shell': S_mass, 'Tubes': T_mass, 'Tubesheet': TS_mass, 'Baffles': B_mass, 'Total': score_mass}
        total_score = score_mass + penalty
    else:
        total_score = _hx_total_cost_standalone(particle.HX) + penalty

    return {
        'Q': particle.Q, 'DP_h': particle.DP_h, 'DP_c': particle.DP_c,
        'score': total_score, 'penalty': penalty, 'masses': masses,
    }


#%%
if __name__ == "__main__":

    HX_test = ShellAndTubeSizingOpt()

    test_case = "Methanol"

    n_disc = 5
    Tube_t_flag = True
    obj = 'mass'
    print_flag = 1

    n_part = 50
    max_iter = 50

    if test_case == "Methanol":
    
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2],
                            'tube_layout' : [0,45,60],
                            'n_parallel' : [1,2,3]}
    
        HX_test.set_inputs(
            fluid_H = 'Methanol',
            T_su_H = 273.15 + 95, # K
            P_su_H = 10*1e5, # Pa
            m_dot_H = 27.8, # kg/s
    
            fluid_C = 'Water',
            T_su_C = 273.15 + 25, # K
            P_su_C = 5*1e5, # Pa
            m_dot_C = 68.9, # kg/s
            )
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                foul_t = 0.0002, foul_s = 0.00033, tube_cond = 50, Overdesign = 0,
                                Shell_Side = 'H',
                                Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True,
                                n_disc = n_disc, Tube_t_flag = Tube_t_flag,
    
                                opt_vars = ['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'],
                                T_max_cycle = 273.15+110, # K
                                p_max_cycle = 10*1e5, # Pa
    
                                H_Corr = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"},
                                C_Corr = {"1P" : "Gnielinski", "2P" : "Flow_boiling"},
                                H_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"},
                                C_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Choi_DP"},
    
                                Q_dot = 4.34*1e6,
                                DP_h = 13.2*1e3,
                                DP_c = 4.3*1e3,
                              )

    elif test_case == "R134a":
    
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2],
                            'tube_layout' : [60],
                            'n_parallel' : [1]}
    
        HX_test.set_inputs(
            fluid_H = 'Water',
            T_su_H = 26 + 273.15, # K
            P_su_H = 5*1e5, # Pa
            m_dot_H = 5.35, # kg/s
    
            fluid_C = 'R134a',
            h_su_C = PropsSI('H','T', 273.15+7,'Q',0,'R134a')+1, # J/kg
            P_su_C = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # Pa
            m_dot_C = 1.62, # kg/s
            )
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                foul_t = 0.000176, foul_s = 0.000176, tube_cond = 50, Overdesign = 0,
                                Shell_Side = 'H',
                                Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True,
                                n_disc = n_disc, Tube_t_flag = Tube_t_flag,
    
                                opt_vars = ['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'],
                                T_max_cycle = 273.15+110, # K
                                p_max_cycle = 5*1e5, # Pa
    
                                H_Corr = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"},
                                C_Corr = {"1P" : "Gnielinski", "2P" : "Flow_boiling"},
                                H_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"},
                                C_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Muller_Steinhagen_Heck_DP"},
    
                                Q_dot = 0.313*1e6,
                                DP_h = 8.2*1e3,
                                DP_c = 21.7*1e3,
                              )

    elif test_case == "CO2_CD":
    
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [1,2,4],
                            'tube_layout' : [0,45,60],
                            'n_parallel' : [1,2,3]} 
    
        HX_test.set_inputs(
            fluid_H = 'CO2',
            T_su_H = 289.64, # K
            P_su_H = 4234404, # Pa
            m_dot_H = 335.3, # kg/s
    
            fluid_C = 'Water',
            T_su_C = 0.1 + 273.15, # K
            P_su_C = 5*1e5, # Pa
            m_dot_C = 2000, # kg/s
            )
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                foul_t = 0.000176, foul_s = 0.000176, tube_cond = 50, Overdesign = 0,
                                Shell_Side = 'C',
                                Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True,
                                n_disc = 30,
    
                                opt_vars = ['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'],
                                T_max_cycle = 273.15+140, # K
                                p_max_cycle = 160*1e5, # Pa
    
                                H_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Thome_Condensation"},
                                C_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"},
                                H_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Choi_DP"},
                                C_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"},
    
                                Q_dot = 75285944,
                                DP_h = 2*1e5,
                                DP_c = 1*1e5,
                              )

    elif test_case == "CO2_GH":
    
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                                29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [1],
                            'tube_layout' : [0,45,60],
                            'n_parallel' : [1,2,3]}
    
        HX_test.set_inputs(
            fluid_C = 'CO2',
            T_su_C = 316.5, # K
            P_su_C = 12822693, # Pa
            m_dot_C = 41.85, # kg/s
    
            fluid_H = 'Water',
            T_su_H = 403.15, # K
            P_su_H = 5*1e5, # Pa
            m_dot_H = 38.85, # kg/s
            )
    
        HX_test.set_parameters(
                                n_series = 1, # [-]
                                foul_t = 0.000176, foul_s = 0.000176, tube_cond = 50, Overdesign = 0,
                                Shell_Side = 'H',
                                Flow_Type = 'Shell&Tube', H_DP_ON = True, C_DP_ON = True,
                                n_disc = 30,
    
                                opt_vars = ['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'],
                                T_max_cycle = 273.15+140, # K
                                p_max_cycle = 160*1e5, # Pa
    
                                H_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"},
                                C_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Flow_boiling"},
                                H_DP = {"SC" : "Shell_Kern_DP", "1P" : "Shell_Kern_DP", "2P" : "Shell_Kern_DP"},
                                C_DP = {"SC" : "Gnielinski_DP", "1P" : "Gnielinski_DP", "2P" : "Choi_DP"},
    
                                Q_dot = 8417198,
                                DP_h = 112284,
                                DP_c = 205160.5,
                          )

    bounds = {
                "L_shell" : [1,15],
                "D_o_inch" : [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
                "Shell_ID_inch" : [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]],
                "Tube_pass" : [choice_vectors['Tube_pass'][0], choice_vectors['Tube_pass'][-1]],
                "tube_layout" : [choice_vectors['tube_layout'][0], choice_vectors['tube_layout'][-1]],
                "Baffle_cut" : [15, 45]
                }
    
    HX_test.set_bounds(bounds, choice_vectors=choice_vectors)


    import time
    t0 = time.perf_counter()

    global_best_position, global_best_score, best_particle = HX_test.sizing_parallel(
        n_particles=n_part,
        max_iterations=max_iter,
        obj=obj,
        print_flag=print_flag,
    )

    elapsed = time.perf_counter() - t0
    print(f"Optimization completed in {elapsed:.2f} s")