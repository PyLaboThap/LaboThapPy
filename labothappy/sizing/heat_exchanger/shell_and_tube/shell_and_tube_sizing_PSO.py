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
import pandas as pd
import random
import numpy as np
import copy

# Ignore convergence warnings
import warnings
warnings.filterwarnings('ignore')

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

        def compute_geom(self):
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
            
            # Pipe Thickness
            Tube_t = carbon_steel_pipe_thickness_mm(self.position['D_o_inch']*25.4/1e3, self.T_max_cycle, self.su_S.p, self.P_max_cycle)
                        
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

        def HeatTransferRate(self):
            
            self.compute_geom()
            
            self.HX = HeatExchangerMB('Shell&Tube')

            if self.params['Shell_Side'] == 'C':

                self.HX.set_inputs(
                    # First fluid
                    fluid_H = self.su_T.fluid,
                    m_dot_H = self.su_T.m_dot, # kg/s

                    # Second fluid
                    fluid_C = self.su_S.fluid,
                    m_dot_C = self.su_S.m_dot, # kg/s  # Make sure to include fluid information
                )
                
                HX_req_inputs = self.HX.get_required_inputs()

                for input_var in self.su_S.variables_input:
                    input_str = input_var[0] + '_su_C'
                    if input_str in HX_req_inputs:    
                        self.HX.set_inputs(**{input_str: input_var[1]})
                    else:
                        input_str = input_var[0].lower() + '_su_C'
                        self.HX.set_inputs(**{input_str: input_var[1]})

                for input_var in self.su_T.variables_input:
                    input_str = input_var[0] + '_su_H'
                    if input_str in HX_req_inputs:    
                        self.HX.set_inputs(**{input_str: input_var[1]})
                    else:
                        input_str = input_var[0].lower() + '_su_H'
                        self.HX.set_inputs(**{input_str: input_var[1]})

            else:
                self.HX.set_inputs(
                    # First fluid
                    fluid_H = self.su_S.fluid,
                    m_dot_H = self.su_S.m_dot, # kg/s

                    # Second fluid
                    fluid_C = self.su_T.fluid,
                    m_dot_C = self.su_T.m_dot, # kg/s  # Make sure to include fluid information
                )    
                
                HX_req_inputs = self.HX.get_required_inputs()

                for input_var in self.su_T.variables_input:
                    input_str = input_var[0] + '_su_C'
                    if input_str in HX_req_inputs:    
                        self.HX.set_inputs(**{input_str: input_var[1]})
                    else:
                        input_str = input_var[0].lower() + '_su_C'
                        self.HX.set_inputs(**{input_str: input_var[1]})

                for input_var in self.su_S.variables_input:
                    input_str = input_var[0] + '_su_H'
                    if input_str in HX_req_inputs:    
                        self.HX.set_inputs(**{input_str: input_var[1]})
                    else:
                        input_str = input_var[0].lower() + '_su_H'
                        self.HX.set_inputs(**{input_str: input_var[1]})

            "Correlation Loading And Setting"

            Corr_H = self.H_htc_Corr
            Corr_C = self.C_htc_Corr

            # Corr_H = {"1P" : "Gnielinski", "2P" : "Flow_boiling", "SC" : "Liu_sCO2"}
            # Corr_C = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
            
            Corr_H_DP = self.H_DP_Corr 
            Corr_C_DP = self.C_DP_Corr # "Gnielinski_DP"
            
            # Corr_H_DP = "Shell_Kern_DP"
            # Corr_C_DP = "Muller_Steinhagen_Heck_DP"

            # Corr_H_DP = "Cheng_CO2_DP"
            # Corr_C_DP = "Shell_Kern_DP"

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
        new_particle.HX = copy.copy(particle.HX)
        
        # ⚠️ Do NOT copy the AbstractState object
        
        return new_particle

    #%% SCORE RELATED
    
    def HX_Mass(self, HX_params):
        
        rho_carbon_steel = 7850 # kg/m^3
        
        "Shell Mass"
                
        shell_t = shell_thickness(HX_params['Shell_ID'], (self.su_S.T + self.ex_S.T)/2, self.P_max_cycle)        
        HX_params['t_S'] = shell_t
        
        Shell_OD = HX_params['Shell_ID'] + 2*shell_t       
        Shell_volume = np.pi*((Shell_OD/2)**2 - (HX_params['Shell_ID']/2)**2)*HX_params['Tube_L'] + shell_t*np.pi*Shell_OD**2/4 
        Shell_mass = Shell_volume*rho_carbon_steel
        
        "Tube Mass"
        
        T_mass = np.pi*((HX_params['Tube_OD']/2)**2 - ((HX_params['Tube_OD']-2*HX_params['Tube_t'])/2)**2)*HX_params['Tube_L']*HX_params['n_tubes']*rho_carbon_steel*HX_params['n_series']

        "Tube Sheet Mass"
        
        TS_t = tube_sheet_thickness(HX_params['Tube_OD'],HX_params['Tube_OD']*HX_params['pitch_ratio'], self.su_S.T, self.P_max_cycle, HX_params["Shell_ID"]) # HX_params["Shell_ID"] assumed to be the gasket diameter
        Full_Tube_sheet_A = np.pi*(HX_params["Shell_ID"]/2)**2
        Tube_in_tube_sheet_A = HX_params["n_tubes"]*np.pi*(HX_params["Tube_OD"]/2)**2
        
        TS_mass = TS_t*(Full_Tube_sheet_A - Tube_in_tube_sheet_A)*rho_carbon_steel*2*HX_params['n_series']
        
        HX_params['t_TS'] = TS_t
        
        "Baffle Mass"
        B_t = baffle_thickness(HX_params["Shell_ID"], HX_params["Baffle_cut"]/100, self.su_S.D, self.su_S.T)
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

    def evaluate_with_penalty(self, objective_function, particle, constraints, penalty_factor,i):
        """
        Evaluates the objective function with a penalty for constraint violations.
        """
        
        key = self.pos_key(particle.position)

        # fast path: seen already
        with self.seen_lock:
            hit = self.seen.get(key)
        if hit is not None:
            Q, DP_h, DP_c, total, masses = hit
            particle.Q, particle.DP_h, particle.DP_c = Q, DP_h, DP_c
            particle.masses = masses
            particle.set_score(total)
            return total
        
        particle.HeatTransferRate()
        
        score, S_mass, T_mass, TS_mass, B_mass = objective_function(particle.HX.params)
        particle.masses = {'Shell' : S_mass,
                           'Tubes' : T_mass,
                           'Tubesheet' : TS_mass,
                           'Baffles' : B_mass}
        
        penalty = 0
    
        # List of constraint checks
        constraints = [
            self.constraint_Q_dot(particle.Q),
            self.constraint_DP_h(particle.DP_h),
            self.constraint_DP_c(particle.DP_c)
        ]
    
        # Compute total penalty
        penalty = sum(penalty_factor * abs(value) for value in constraints)
        
        # print(f"penalty : {penalty}")
    
        # Update particle score
        total_score = score + penalty
        flag = particle.set_score(total_score)
    
        # if flag == 1:
        #     print(f"Particle {i} PB !!!")
    
        # If no penalty, add unique params to the suitable set
        if penalty == 0 and particle.params not in self.suitable_param_set:
            self.suitable_param_set.append(copy.deepcopy(particle.params))
            self.suitable_param_set[-1]['score'] = particle.score
        
        # store for reuse
        with self.seen_lock:
            self.seen[key] = (particle.Q, particle.DP_h, particle.DP_c, total_score, particle.masses)
                
        return total_score

    #%%

    def particle_swarm_optimization(self, objective_function, bounds, num_particles=30, num_dimensions=2, max_iterations=50, 
                                inertia_weight=0.4, cognitive_constant=1.5, social_constant=1.5, constraints = None,
                                penalty_factor=1000):
        """
        Perform Particle Swarm Optimization (PSO) to minimize the given objective function with constraints.
        
        Parameters:
        - objective_function: Function to be minimized. Should take a particle position (array-like) as input.
        - bounds: List of tuples specifying (min, max) for each dimension.
        - num_particles: Number of particles in the swarm (default 30).
        - num_dimensions: Number of dimensions of the search space (default 2).
        - max_iterations: Maximum number of iterations (default 100).
        - inertia_weight: Inertia weight to balance exploration and exploitation (default 0.7).
        - cognitive_constant: Constant to control personal best influence (default 1.5).
        - social_constant: Constant to control global best influence (default 1.5).
        - discrete_indices: List of indices of variables that should be discrete.
        - discrete_values: Dictionary specifying allowed discrete values for each index.
        - constraints: List of constraint functions. Each should return a value, where negative or zero indicates that the constraint is satisfied.
        - penalty_factor: Factor to penalize constraint violations (default 1000).
        
        Returns:
        - global_best_position: Position of the best solution found.
        - global_best_score: Value of the objective function at the best solution.
        """

        # Initialize particle positions

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

        # Get scores for sorting
        particle_scores = np.array([particle.personnal_best_score for particle in self.particles])
        
        # Get indices of 50 best particles (sorted in ascending order)
        best_particle_indices = np.argsort(particle_scores)[:num_particles]
        
        # Keep only the top 50 particles
        self.particles = [self.particles[i] for i in best_particle_indices]
        
        # Update personal best positions and scores accordingly
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

        self.best_particle.compute_geom()

        # Initialize velocities and positions as dictionaries
        cognitive_velocity = {}
        social_velocity = {}

        # PSO loop
        for iteration in range(max_iterations):

            # print("==============================")
            # print(f"Iteration {iteration + 1}")

            for i in range(num_particles):

                flag = self.particles[i].check_reinit()
                if flag:
                    # print("Particle Reinitialized")
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

                    # Extract the positions from dictionaries
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
                    
                    # if L_flag == 1:
                    # Update velocity with inertia term
                    self.particles[i].velocity[opt_var] = (inertia_weight * self.particles[i].velocity[opt_var] + cognitive_velocity[opt_var] + social_velocity[opt_var])

                    # Update the position of the particle
                    self.particles[i].position[opt_var] += round(self.particles[i].velocity[opt_var],3)
                    
                        
                    if opt_var == 'Central_spac':
                        low_bound_central_spac = (self.particles[i].position['Shell_ID_inch']/5)*25.4*1e-3 # [m]
                        high_bound_central_spac = (74*self.particles[i].position['D_o_inch']**(0.75))*25.4*1e-3 # [m]
                        
                        L_shell_divisors = find_divisors_between_bounds(self.particles[i].position['L_shell'],low_bound_central_spac,high_bound_central_spac)

                        if len(L_shell_divisors) == 0:
                            self.particles[i].position[opt_var] = round(low_bound_central_spac,2)
                        else:
                            self.particles[i].position[opt_var] = min(L_shell_divisors, key=lambda x: abs(x - self.particles[i].position[opt_var]))
                        
                    # Handle discrete variables if needed (apply rounding or snapping to discrete values)
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

                # Enforce boundary conditions and handle discrete variables
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
                        # Bound constraints
                        if self.particles[i].position[bound_key] < self.bounds[bound_key][0]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][0]
                            self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key] #self.particles[i].velocity[bound_key]
                            bound_flag = 1
    
                        if self.particles[i].position[bound_key] > self.bounds[bound_key][1]:
                            self.particles[i].position[bound_key] = self.bounds[bound_key][1]
                            self.particles[i].velocity[bound_key] = -0.5*self.particles[i].velocity[bound_key] # self.particles[i].velocity[bound_key]
                            bound_flag = 1

                # Evaluate the new position with penalty for constraint violation
                
                # self.particles[i].compute_geom()
                
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
                # print("BEST SCORE BEATEN")
                self.global_best_score = new_pot_global_best_score
                self.best_particle = self.clone_Particle(self.particles[np.argmin(self.personal_best_scores)])
                self.global_best_position = self.best_particle.position
                self.global_best_Q = self.best_particle.Q
                self.global_best_DP_h = self.best_particle.DP_h
                self.global_best_DP_c = self.best_particle.DP_c

                self.best_particle.compute_geom()

            # Optionally, print progress
            print("===========================")
            print(f"Iteration {iteration+1}/{max_iterations}, Global Best Score: {self.global_best_score}, Related Q: {self.global_best_Q}")
            print(f"Related DP_h: {self.global_best_DP_h}, Related DP_c: {self.global_best_DP_c}")
            print(f"Best Position : {self.global_best_position}")
            print(f"Best Part Velocity : {self.best_particle.velocity}")
            
        return self.global_best_position, self.global_best_score, self.best_particle
    
    def opt_size(self):

        self.particle_swarm_optimization(objective_function = self.HX_Mass , bounds = self.bounds, num_particles = 50, num_dimensions = len(self.opt_vars), max_iterations = 50, inertia_weight = 0.5,
                                          cognitive_constant = 0.5, social_constant = 0.5, constraints = [self.constraint_Q_dot, self.constraint_DP_h, self.constraint_DP_c], penalty_factor = 1)

        self.best_particle.compute_geom()


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

        self.CAPEX,_,_ = self.cost_calculator.calculate_total_cost()
        
        return self.global_best_position, self.global_best_score, self.best_particle
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

    su_S = MassConnector()
    su_S.set_properties(T = 273.15 + 95, # K
                        P = 10*1e5, # 5*1e5, # Pa
                        m_dot = 27.8, # kg/s
                        fluid = 'Methanol'
                        )

    ex_S = MassConnector()
    ex_S.set_properties(T = 273.15 + 40, # K
                        P = 10*1e5, # 4.5*1e5, # Pa
                        m_dot = 27.8, # kg/s
                        fluid = 'Methanol'
                        )

    su_T = MassConnector()
    su_T.set_properties(T = 273.15 + 25, # K
                        P = 5*1e5, # 51.75*1e3, # Pa
                        m_dot = 68.9, # kg/s
                        fluid = 'Water'
                        )

    ex_T = MassConnector()
    ex_T.set_properties(T = 273.15 + 40, # K
                        P = 5*1e5, # Pa
                        m_dot = 68.9, # kg/s
                        fluid = 'Water'
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

    su_S = MassConnector()
    su_S.set_properties(T = 273.15 + 26, # K
                        P = 5*1e5, # 51.75*1e3, # Pa
                        m_dot = 5.35, # kg/s
                        fluid = 'Water'
                        )

    ex_S = MassConnector()
    ex_S.set_properties(T = 273.15 + 11.7, # K
                        P = 5*1e5, # Pa
                        m_dot = 5.35, # kg/s
                        fluid = 'Water'
                        )

    su_T = MassConnector()
    su_T.set_properties(P = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # K
                        H = PropsSI('H','T', 273.15+7,'Q',0,'R134a')+1, # Pa
                        m_dot = 1.62, # kg/s
                        fluid = 'R134a'
                        )

    ex_T = MassConnector()
    ex_T.set_properties(P = PropsSI('P','T', 273.15+7,'Q',1,'R134a') - 15*1e3, # K
                        H = PropsSI('H','T', 273.15+7,'Q',1,'R134a'), # Pa
                        m_dot = 1.62, # kg/s
                        fluid = 'R134a'
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
    P_max_cycle = 160*1e5 # Pa
    T_max_cycle = 273.15+140 # K 
    
    HX_test.set_max_cycle_prop(T_max_cycle = T_max_cycle, p_max_cycle = P_max_cycle)
    
    """
    Thermodynamical parameters : Inlet and Outlet Design States
    """
    
    HX_test.set_inputs(
        # First fluid
        fluid_H = 'CO2',
        T_su_H = 306.92, # K
        P_su_H = 5945330, # Pa
        m_dot_H = 33.36, # kg/s

        # Second fluid
        fluid_C = 'Water',
        T_su_C = 15 + 273.15, # K
        P_su_C = 5*1e5, # Pa
        m_dot_C = 100, # kg/s  # Make sure to include fluid information
        )
    
    "Constraints Values"
    Q_dot_cstr = 6295150
    DP_h_cstr = 2*1e5
    DP_c_cstr = 50*1e3

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
                            n_disc = 50
                          )

    H_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}    
    C_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    
    H_DP = "Choi_DP"
    C_DP = "Shell_Kern_DP"
    
    HX_test.set_corr(H_Corr, C_Corr, H_DP, C_DP)


HX_test.set_thermo_BC(su_S = su_S, ex_S = ex_S, su_T = su_T, ex_T = ex_T)

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


