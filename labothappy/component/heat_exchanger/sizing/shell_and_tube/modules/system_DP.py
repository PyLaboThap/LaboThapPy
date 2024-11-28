
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

n_passes = 2 # (or 1)

Tube_layout_angle = 45 # [°] (or 60, 90) : 45 and 90 mean square / 60 means triangular

Pitch_ratio constrainted to values depending on D_o, could be varied from 1.2 to 1.5 on its own
Square:
-------
D_o = 1/2     [in] => Pitch_ratio = 1.25
D_o = 3/4     [in] => Pitch_ratio = (1)/D_o
D_o = 1       [in] => Pitch_ratio = (1+1/4)/D_o
D_o = 1 + 1/4 [in] => Pitch_ratio = (1+9/16)/D_o
D_o = 1 + 1/2 [in] => Pitch_ratio = (1+7/8)/D_o
Triangular:
-----------
D_o = 1/2     [in] => Pitch_ratio = 1.25
D_o = 3/4     [in] => Pitch_ratio = (15/16)/D_o
D_o = 1       [in] => Pitch_ratio = (1+1/4)/D_o
D_o = 1 + 1/4 [in] => Pitch_ratio = (1+9/16)/D_o
D_o = 1 + 1/2 [in] => Pitch_ratio = (1+7/8)/D_o

Baffle_cut = 0.25 # Could be varied from 0.15 to 0.4 but 0.25 is usual value for liquid flow
"""

import __init__
import copy
from library.component.sizing.heat_exchanger.shell_and_tube.modules.tubes_toolbox import estimate_number_of_tubes, carbon_steel_pipe_thickness, pitch_ratio_fun
from scipy.interpolate import interp1d
from CoolProp.CoolProp import PropsSI
from library.component.sizing.heat_exchanger.basic_sizing_UA import find_UA
from central_spacing_comp import find_divisors_between_bounds
from library.connector.mass_connector import MassConnector
from component.base_component import BaseComponent

from library.component.steady_state.heat_exchanger.moving_boundary.charge_sensitive.simulation_model_AS_DP import HeatExchangerMB

import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings('ignore')

#%%

class ShellAndTubeSizingOpt(BaseComponent):
    class Particle(BaseComponent):
        def __init__(self, params = None, su_S = None, ex_S = None, su_T = None, ex_T = None, choice_vectors = None, P_max_cycle = None, T_max_cycle = None):
            super().__init__()

            self.params = params
            
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
            self.personnal_best_score = None

            # Will be Mass connectors
            self.su_S = su_S
            self.ex_S = ex_S

            self.su_T = su_T
            self.ex_T = ex_T
    
        def set_position(self, position):
            self.position = position            
            return
        
        def set_velocity(self, velocity):
            self.velocity = velocity
            return 
        
        def set_score(self, score):
            self.score = score
            
            if self.personnal_best_score is None or self.personnal_best_score > score:
                self.personnal_best_score = self.score
                self.personnal_best_position = self.position

            return

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

            # Pipe Thickness
            pipe_wall_T = (self.su_S.T + self.su_T.T)/2
            pipe_thickness = carbon_steel_pipe_thickness(self.choice_vectors['D_o_inch'], self.T_max_cycle, self.su_S.p, self.P_max_cycle)

            Tube_t = pipe_thickness[str(self.position['D_o_inch'])]

            # pitch_ratio
            pitch_ratio = pitch_ratio_fun(self.position['D_o_inch'], self.params['tube_layout'])

            # Pipe length
            L_tube = self.position['L_shell']*self.params['Tube_pass']

            # Cross Passes
            Cross_Passes = round(self.position['L_shell']/self.position['Central_spac']) - 1

            D_o = self.position['D_o_inch']*25.4*1e-3
            Shell_ID = self.position['Shell_ID_inch']*25.4*1e-3

            # Number of tubes 
            min_tubes_in_row = 8
            n_tubes = estimate_number_of_tubes(Shell_ID, D_o, pitch_ratio*D_o, self.params['tube_layout'], min_tubes_in_row)[0]

            # HT Area and HTX volumes
            A_eff = n_tubes*L_tube*np.pi*D_o
            
            T_V_tot = L_tube*n_tubes*np.pi*((D_o - 2*Tube_t)/2)**2

            T_V_out = np.pi*(D_o/2)**2*L_tube*n_tubes
            S_V_tot = self.position['L_shell']*np.pi*(Shell_ID/2)**2 - T_V_out
                
            self.set_parameters( 
                            A_eff = A_eff, S_V_tot = S_V_tot, Shell_ID = Shell_ID, T_V_tot = T_V_tot, Tube_L = L_tube, 
                            Tube_OD = D_o, Tube_t = Tube_t, central_spacing = self.position['Central_spac'], 
                            cross_passes = Cross_Passes, n_tubes = n_tubes, pitch_ratio = pitch_ratio
                            ) 
            return

        def HeatTransferRate(self):

            HX = HeatExchangerMB('Shell&Tube')

            if self.params['Shell_Side'] == 'C':
                HX.set_inputs(
                    # First fluid
                    Hsu_fluid = self.su_T.fluid,
                    Hsu_T = self.su_T.T, # K
                    Hsu_p = self.su_T.p, # Pa
                    Hsu_m_dot = self.su_T.m_dot, # kg/s

                    # Second fluid
                    Csu_fluid = self.su_S.fluid,
                    Csu_T = self.su_S.T, # K
                    Csu_p = self.su_S.p, # Pa
                    Csu_m_dot = self.su_S.m_dot, # kg/s  # Make sure to include fluid information
                )
            else:
                HX.set_inputs(
                    # First fluid
                    Hsu_fluid = self.su_S.fluid,
                    Hsu_T = self.su_S.T, # K
                    Hsu_p = self.su_S.p, # Pa
                    Hsu_m_dot = self.su_S.m_dot, # kg/s

                    # Second fluid
                    Csu_fluid = self.su_T.fluid,
                    Csu_T = self.su_T.T, # K
                    Csu_p = self.su_T.p, # Pa
                    Csu_m_dot = self.su_T.m_dot, # kg/s  # Make sure to include fluid information
                )          

            "Correlation Loading And Setting"

            Corr_H = {"1P" : "Gnielinski", "2P" : "Horizontal_Tube_Internal_Condensation"}
            Corr_C = {"1P" : "shell_htc_kern", "2P" : "shell_htc_kern"}

            HX.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C) # 'User-Defined' or 'Correlation' # 31

            "Parameters Setting"

            HX.set_parameters(
                A_eff = self.params['A_eff'], Baffle_cut = self.params['Baffle_cut'], S_V_tot = self.params['S_V_tot'],
                Shell_ID = self.params['Shell_ID'], T_V_tot = self.params['T_V_tot'], Tube_L = self.params['Tube_L'], 
                Tube_OD = self.params['Tube_OD'], Tube_pass = self.params['Tube_pass'], Tube_t = self.params['Tube_t'],
                central_spacing = self.params['central_spacing'], cross_passes = self.params['cross_passes'], foul_s = self.params['foul_s'],
                foul_t = self.params['foul_t'], n_series = self.params['n_series'], n_tubes = self.params['n_tubes'], 
                pitch_ratio = self.params['pitch_ratio'], tube_cond = self.params['tube_cond'], tube_layout = self.params['tube_layout'],

                Shell_Side = self.params['Shell_Side'],

                Flow_Type = self.params['Flow_Type'], H_DP_ON = self.params['H_DP_ON'], C_DP_ON = self.params['C_DP_ON'], n_disc = self.params['n_disc']) 

            Corr_H_DP = "pipe_internal_DP"
            Corr_C_DP = "shell_DP_kern"

            # HX.set_DP(DP_type="User-Defined", UD_H_DP=1e4, UD_C_DP=1e4)
            HX.set_DP(DP_type = "Correlation", Corr_H=Corr_H_DP, Corr_C=Corr_C_DP)    

            try: 
                HX.solve()
                self.Q = HX.Q
                self.DP_h = HX.DP_h
                self.DP_c = HX.DP_c
                return HX.Q
            except:
                self.Q = 0
                self.DP_h = 0
                self.DP_c = 0
                return 0

    def __init__(self):
        super().__init__()

        self.params = {}

        self.particles = None
        self.global_best_position = None
        self.global_best_score = None
        self.best_particle = None

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

    def Tube_Mass(self, particle_params):
        
        rho_carbon_steel = 7850 # kg/m^3
        T_mass = np.pi*((particle_params['Tube_OD']/2)**2 - ((particle_params['Tube_OD'] - particle_params['Tube_t'])/2)**2) * particle_params['Tube_L'] * rho_carbon_steel * particle_params['n_tubes'] * particle_params['n_series']

        return T_mass

    def constraint_Q_dot(self, Q_particle):
        Q_dot_val_cstr = 13300*1e3 # 0% Overdesign
        return Q_particle - Q_dot_val_cstr # [W] : Q_dot - 13300000 <= 0

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
        start = (lower_bound + multiple) // multiple * multiple

        # Find the largest multiple within the range
        end = upper_bound // multiple * multiple

        if start > end:
            print("No multiples of {} in the range [{}, {}]".format(multiple, lower_bound, upper_bound))
            return 0

        # Generate a random multiple within the range
        num_multiples = (end - start) // multiple + 1

        random_multiple = random.randint(0, int(num_multiples) - 1) * multiple + start

        return random_multiple

    def init_particle(self, particle):

        particle_position = {}
        particle_velocity = {}

        for opt_var in self.opt_vars:
            particle_position[opt_var] = None
            particle_velocity[opt_var] = round(random.uniform(-1, 1),2)

            # Choose values randomly from choice vectors

        for choice_vector_key in self.choice_vectors:
            if isinstance(self.choice_vectors[choice_vector_key][0], str):    
                particle_position[choice_vector_key] = pd.to_numeric(random.choice(self.choice_vectors[choice_vector_key]))
            else:
                particle_position[choice_vector_key] = random.choice(self.choice_vectors[choice_vector_key])

            # Compute values of central spacing, then of L_shell 

        low_bound_central_spac = (particle_position['Shell_ID_inch']/5)*25.4*1e-3 # [m]
        high_bound_central_spac = (74*particle_position['D_o_inch']**(0.75))*25.4*1e-3 # [m]
        
        particle_position['Central_spac'] = round(random.uniform(low_bound_central_spac, high_bound_central_spac),2)
        particle_position['L_shell'] = round(self.random_multiple(lower_bound = self.bounds['L_shell'][0], upper_bound = self.bounds['L_shell'][1], multiple = particle_position['Central_spac']),2)

        # Put these positions in the particles 

        particle.set_position(particle_position)
        particle.set_velocity(particle_velocity)

        for opt_var in particle.position.keys():
            particle.unmoved[opt_var] = 0

        return 

    def evaluate_with_penalty(self, objective_function, particle, constraints, penalty_factor):
        """
        Evaluates the objective function with a penalty for constraint violations.
        """
        score = objective_function(particle.params)
        penalty = 0

        # Q_dot constraint
        constraint_value = self.constraint_Q_dot(particle.Q)
        if constraint_value < 0:  # If constraint is violated, add a penalty
            penalty += penalty_factor * abs(constraint_value)

        particle.set_score(score + penalty)

        return score + penalty
        
    def particle_swarm_optimization(self, objective_function, bounds, num_particles=30, num_dimensions=2, max_iterations=100, 
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
                                        P_max_cycle = self.P_max_cycle, T_max_cycle = self.T_max_cycle) for _ in range(num_particles)]

        self.all_scores = np.zeros((num_particles, max_iterations+1))

        self.particles_all_pos = {}

        for opt_var in self.opt_vars:
            self.particles_all_pos[opt_var] = np.zeros((num_particles + 1, max_iterations))

        for i in range(num_particles):
            self.init_particle(self.particles[i])

        # Compute particle geometry
        for i in range(len(self.particles)):
            self.particles[i].compute_geom()       
            self.particles[i].HeatTransferRate()
            score = self.evaluate_with_penalty(objective_function, self.particles[i], constraints, penalty_factor)
            self.all_scores[i][0] = score

        self.personal_best_positions = np.array([self.particles[i].personnal_best_position for i in range(len(self.particles))])
        self.personal_best_scores = np.array([self.particles[i].personnal_best_score for i in range(len(self.particles))])

        self.global_best_score = min(self.personal_best_scores)
        self.global_best_position = self.personal_best_positions[np.argmin(self.personal_best_scores)]
        self.global_best_Q = self.particles[np.argmin(self.personal_best_scores)].Q
        self.global_best_DP_h = self.particles[np.argmin(self.personal_best_scores)].DP_h
        self.global_best_DP_c = self.particles[np.argmin(self.personal_best_scores)].DP_c
        self.best_particle = copy.copy(self.particles[np.argmin(self.personal_best_scores)])

        # for i in range(num_particles):
        #     print("==============================")
        #     print(f"Initial score ({i}) : {self.particles[i].score}")
        #     print(f"Initial position ({i}) : {self.particles[i].position}")
        #     print(f"Initial velocity ({i}) : {self.particles[i].velocity}")
        #     print("\n")

        # Initialize velocities and positions as dictionaries
        cognitive_velocity = {}
        social_velocity = {}

        # PSO loop
        for iteration in range(max_iterations):

            print("==============================")
            print(f"Iteration {iteration + 1}")

            for i in range(num_particles):

                print(f"New score ({i}) : {self.particles[i].score}")
                print(f"Related Q ({i}) : {self.particles[i].Q}")
                print(f"Related DP_h ({i}) : {self.particles[i].DP_h}")
                print(f"Related DP_c ({i}) : {self.particles[i].DP_c}")
                print(f"New position ({i}) : {self.particles[i].position}")
                print(f"New velocity ({i}) : {self.particles[i].velocity}")
                print(f"\n")

                flag = self.particles[i].check_reinit()
                if flag:
                    print("Particle Reinitialized")
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

                    # Update velocity for each optimization variable (opt_var)
                    if opt_var in self.choice_vectors.keys():
                        if opt_var == 'D_o_inch':
                            r1 = np.random.rand()*self.choice_vectors[opt_var][0]/3
                            r2 = np.random.rand()*self.choice_vectors[opt_var][0]/3
                        else:
                            r1 = np.random.rand()*self.choice_vectors[opt_var][0]/10
                            r2 = np.random.rand()*self.choice_vectors[opt_var][0]/10
                    elif opt_var in self.bounds.keys():
                        if opt_var == 'L_shell':
                            r1 = np.random.rand()*self.bounds[opt_var][0]/3
                            r2 = np.random.rand()*self.bounds[opt_var][0]/3
                        else:
                            r1 = np.random.rand()*self.bounds[opt_var][0]/10
                            r2 = np.random.rand()*self.bounds[opt_var][0]/10                            
                    

                    cognitive_cos_fact = np.cos((iteration / max_iterations) * 2 * np.pi)
                    social_cos_fact = 1 - cognitive_cos_fact

                    cognitive_fact = cognitive_constant*(0.5 + cognitive_cos_fact)
                    social_fact = social_constant*(0.5 + social_cos_fact)


                    cognitive_velocity[opt_var] = cognitive_fact * r1 * (personal_best_position_val - current_position_val)
                    social_velocity[opt_var] = social_fact * r2 * (global_best_position_val - current_position_val)
                    
                    # Update velocity with inertia term
                    self.particles[i].velocity[opt_var] = (inertia_weight * self.particles[i].velocity[opt_var] +
                                                        cognitive_velocity[opt_var] + 
                                                        social_velocity[opt_var])

                    # Update the position of the particle
                    self.particles[i].position[opt_var] += round(self.particles[i].velocity[opt_var],2)
                    
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
                    # Bound constraints
                    if self.particles[i].position[bound_key] < self.bounds[bound_key][0]:
                        self.particles[i].position[bound_key] = self.bounds[bound_key][0]
                        bound_flag = 1

                    if self.particles[i].position[bound_key] > self.bounds[bound_key][1]:
                        self.particles[i].position[bound_key] = self.bounds[bound_key][1]
                        bound_flag = 1

                # Evaluate the new position with penalty for constraint violation
                if bound_flag == 1:
                    self.particles[i].compute_geom()
                    self.particles[i].Q = self.particles[i].HeatTransferRate()
                    new_score = self.evaluate_with_penalty(objective_function, self.particles[i], constraints, penalty_factor) + 1e6

                else:
                    self.particles[i].compute_geom()
                    self.particles[i].HeatTransferRate()
                    new_score = self.evaluate_with_penalty(objective_function, self.particles[i], constraints, penalty_factor)

                self.particles[i].set_score(new_score)

            for i in range(num_particles):
                self.all_scores[i][iteration + 1] = self.particles[i].score

            self.personal_best_positions = np.array([self.particles[i].personnal_best_position for i in range(len(self.particles))])
            self.personal_best_scores = np.array([self.particles[i].personnal_best_score for i in range(len(self.particles))])

            new_pot_global_best_score = min(self.personal_best_scores)

            if new_pot_global_best_score < self.global_best_score:
                self.global_best_score = new_pot_global_best_score
                self.global_best_position = self.personal_best_positions[np.argmin(self.personal_best_scores)]
                self.global_best_Q = self.particles[np.argmin(self.personal_best_scores)].Q
                self.global_best_DP_h = self.particles[np.argmin(self.personal_best_scores)].DP_h
                self.global_best_DP_c = self.particles[np.argmin(self.personal_best_scores)].DP_c            
                self.best_particle = copy.copy(self.particles[np.argmin(self.personal_best_scores)])


            # Optionally, print progress
            print("===========================")
            print(f"Iteration {iteration+1}/{max_iterations}, Global Best Score: {self.global_best_score}, Related Q: {self.global_best_Q}")
            print(f"Related DP_h: {self.global_best_DP_h}, Related DP_c: {self.global_best_DP_c}")
            print(f"Best Position : {self.global_best_position}")
        


        return self.global_best_position, self.global_best_score, self.best_particle
    
    def opt_size(self):

        return self.particle_swarm_optimization(objective_function = self.Tube_Mass , bounds = self.bounds, num_particles = 50, num_dimensions = len(self.opt_vars), max_iterations = 15, inertia_weight = 0.6,
                                         cognitive_constant = 1, social_constant = 1, constraints = [self.constraint_Q_dot], penalty_factor = 1)

        

"""
Preliminary Sizing Method : Heat Exchangers (Kakac, Liu, Pramuanjaroenkij) - Section 9.3.1
"""

def ShellAndTube_PrelimSizing(Q_dot, T_c_i, T_c_o, T_h_i, T_h_o, params, h_coefs):

    UA_req = find_UA(Q_dot = Q_dot, T_c_i = T_c_i, T_c_o = T_c_o, T_h_i = T_h_i, T_h_o = T_h_o, flow = 'Shell&Tube', params = params)
    U = 0

    for key, value in h_coefs.items():
        U += 1/value

    U = 1/U

    # Required Area
    A_req = UA_req/U
    A_req_OD = A_req*1.2

    # Number of Tubes
    
    # Tube count constant
    if params['Tube_pass'] == 1:
        CTP = 0.93
    elif params['Tube_pass'] == 2:
        CTP = 0.9
    elif params['Tube_pass'] == 3:
        CTP = 0.85
    else:
        print("Tube_pass number not considered.")

    # Tube layout constant
    if params['tube_layout'] == 45 or params['tube_layout'] == 90:
        CL = 1
    elif params['tube_layout'] == 30 or params['tube_layout'] == 60:
        CL = 0.87

    D_s = 0.637*(CL/CTP)**(0.5)*((A_req_OD*((params['pitch_ratio'])**2)*params['Tube_OD'])/params['L_tube'])**(1/2)
    N_tubes = 0.785*(CTP/CL)*(D_s**2)/(params['pitch_ratio']**2*params['Tube_OD']**2)

    return A_req, A_req_OD, D_s, N_tubes

HX_test = ShellAndTubeSizingOpt()

"""
Optimization related parameters/variables
"""

HX_test.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac'])

choice_vectors = {
                    'D_o_inch' : [0.5, 0.75, 1, 1.25, 1.5],
                    'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,
                        29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78,
                        84, 90, 96, 108, 120]
}

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
P_max_cycle = 1048*1e3 # Pa
T_max_cycle = 273.15+147

HX_test.set_max_cycle_prop(T_max_cycle = T_max_cycle, p_max_cycle = P_max_cycle)

"""
Thermodynamical parameters : Inlet and Outlet Design States
"""

su_S = MassConnector()
su_S.set_properties(T = 273.15 + 24, # K
                    P = 101.3*1e3, # Pa
                    m_dot = 1000, # kg/s
                    fluid = 'Water'
                    )

ex_S = MassConnector()
ex_S.set_properties(T = 273.15 + 27.29, # K
                    P = 101.3*1e3, # Pa
                    m_dot = 1000, # kg/s
                    fluid = 'Water'
                    )

"From Aitor Code"

su_T = MassConnector()
su_T.set_properties(T = 273.15 + 38.43, # K
                    P = 51.75*1e3, # Pa
                    m_dot = 32.85, # kg/s
                    fluid = 'Cyclopentane'
                    )

ex_T = MassConnector()
ex_T.set_properties(T = 273.15 + 27.21, # K
                    P = 51.75*1e3, # Pa
                    m_dot = 32.85, # kg/s
                    fluid = 'Cyclopentane'
                    )

HX_test.set_thermo_BC(su_S = su_S, ex_S = ex_S, su_T = su_T, ex_T = ex_T)

"""
Parameters Setting
"""

HX_test.set_parameters(
                        tube_layout = 45, # [°]
                        Tube_pass = 2, # [-]
                        n_series = 1, # [-]
                        Baffle_cut = 25, # [%]
                        foul_t = 0,
                        foul_s = 0,
                        tube_cond = 50, # W/(m*K)
                        
                        Shell_Side = 'C',

                        Flow_Type = 'Shell&Tube', 
                        H_DP_ON = True, 
                        C_DP_ON = True, 
                        n_disc = 30
                      )

"""
Geometry Computation Test
"""

# HX_test.set_opt_vars_values(
#                             {'D_o_inch' : 0.75,
#                             'L_shell' : 5,
#                             'Shell_ID_inch' : 90,
#                             'Central_spac' : 0.5}
#                             )

# HX_test.compute_geom()

bounds = {
            "L_shell" : [3,10],
            "D_o_inch" : [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
            "Shell_ID_inch" : [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]]
            }

HX_test.set_bounds(bounds)

global_best_position, global_best_score, best_particle = HX_test.opt_size()

all_scores = HX_test.all_scores

for row in all_scores: 
    plt.plot(row)

plt.show()

for row in all_scores: 
    filtered_row = [x if x <= 10000 else None for x in row]
    plt.plot(filtered_row)
    plt.axis([0,50, 0,10000])

plt.show()

for opt_var in HX_test.particles_all_pos:
    opt_var_matrix = HX_test.particles_all_pos[opt_var]
    for particle in opt_var_matrix:
        if sum(particle) == sum(opt_var_matrix[-1]):
            plt.plot(particle, color = 'k')
        else:
            plt.plot(particle)
    
    plt.title(opt_var)
    plt.show()

print(f"Best global position : {global_best_position}")
print(f"Best global score : {global_best_score}")
print(f"Best particle parameters : {best_particle.params}")
print(f"Related Q : {HX_test.global_best_Q} [W]")