
import __init__
import numpy as np 

from labothappy.connector.mass_connector import MassConnector
from component.base_component import BaseComponent

class SizingAHTX(BaseComponent):

    def __init__(self):

        super().__init__()

        self.su_H = MassConnector()
        self.su_C = MassConnector()

        self.ex_H = MassConnector()
        self.ex_C = MassConnector()

        self.fouling = None
        self.HTX_Type = None
        self.geom_input = None

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return['fouling', 'HTX_Type', 'geom_input']
    
    def sync_inputs(self):
        
        if self.fouling is not None:
            self.inputs['fouling'] = self.fouling
        if self.HTX_Type is not None:
            self.inputs['HTX_Type'] = self.HTX_Type
        if self.fouling is not None:
            self.inputs['geom_input'] = self.geom_input

        return

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'fouling' in self.inputs:
            self.fouling = self.inputs['fouling']
        if 'HTX_Type' in self.inputs:
            self.HTX_Type = self.inputs['HTX_Type']
        if 'geom_input' in self.inputs:
            self.geom_input = self.inputs['geom_input']

        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_H_fluid', 'su_H_T', 'su_H_cp', 'su_H_m_dot']
    
    def get_required_parameters(self):
        return [
            'su_C_fluid', 'su_C_T', 'su_C_m_dot', 'su_C_p', 'su_H_fluid', 'su_H_T', 'su_H_p', 'su_H_m_dot'
        ]
    
    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")

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


    def particle_swarm_optimization(objective_function, bounds, num_particles=30, num_dimensions=2, max_iterations=100, 
                                    inertia_weight=0.7, cognitive_constant=1.5, social_constant=1.5, 
                                    discrete_indices=None, discrete_values=None, 
                                    constraints=None, penalty_factor=1000):
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
        
        def evaluate_with_penalty(position):
            """
            Evaluates the objective function with a penalty for constraint violations.
            """
            score = objective_function(position)
            penalty = 0
            
            if constraints:
                for constraint in constraints:
                    constraint_value = constraint(position)
                    if constraint_value > 0:  # If constraint is violated, add a penalty
                        penalty += penalty_factor * constraint_value ** 2

            return score + penalty

        # Initialize particle positions and velocities
        particle_positions = np.random.uniform(low=[b[0] for b in bounds], 
                                            high=[b[1] for b in bounds], 
                                            size=(num_particles, num_dimensions))
        particle_velocities = np.random.uniform(low=-1, high=1, size=(num_particles, num_dimensions))

        # Initialize personal best positions and global best position
        personal_best_positions = np.copy(particle_positions)
        personal_best_scores = np.array([evaluate_with_penalty(p) for p in personal_best_positions])

        global_best_position = personal_best_positions[np.argmin(personal_best_scores)]
        global_best_score = np.min(personal_best_scores)

        # PSO loop
        for iteration in range(max_iterations):
            for i in range(num_particles):
                # Update velocity
                r1, r2 = np.random.rand(), np.random.rand()  # Random coefficients
                cognitive_velocity = cognitive_constant * r1 * (personal_best_positions[i] - particle_positions[i])
                social_velocity = social_constant * r2 * (global_best_position - particle_positions[i])
                particle_velocities[i] = inertia_weight * particle_velocities[i] + cognitive_velocity + social_velocity

                # Update position
                particle_positions[i] += particle_velocities[i]

                # Enforce boundary conditions and handle discrete variables
                for d in range(num_dimensions):
                    # Bound constraints
                    if particle_positions[i, d] < bounds[d][0]:
                        particle_positions[i, d] = bounds[d][0]
                    if particle_positions[i, d] > bounds[d][1]:
                        particle_positions[i, d] = bounds[d][1]
                    
                    # Apply discrete constraints for specific variables
                    if discrete_indices and d in discrete_indices:
                        # Snap the position to the nearest allowed discrete value
                        allowed_values = discrete_values[d]
                        particle_positions[i, d] = min(allowed_values, key=lambda x: abs(x - particle_positions[i, d]))

                # Evaluate the new position with penalty for constraint violation
                new_score = evaluate_with_penalty(particle_positions[i])

                # Update personal best
                if new_score < personal_best_scores[i]:
                    personal_best_positions[i] = particle_positions[i]
                    personal_best_scores[i] = new_score

            # Update global best
            min_personal_best_score = np.min(personal_best_scores)
            if min_personal_best_score < global_best_score:
                global_best_position = personal_best_positions[np.argmin(personal_best_scores)]
                global_best_score = min_personal_best_score

            # Optionally, print progress
            print(f"Iteration {iteration+1}/{max_iterations}, Global Best Score: {global_best_score}")

        return global_best_position, global_best_score
    
    # def opt_size(self):
        
    #     self.particle_swarm_optimization(bounds = )

    #     return

#%%

HTX = SizingAHTX()

HTX.set_parameters(
                    su_C_fluid = 'Water',
                    su_C_T = 273.15 + 24, # K
                    su_C_p = 101.3*1e3, # Pa
                    su_C_m_dot = 1000, # kg/s

                    su_H_fluid = 'Cyclopentane',
                    su_H_T = 273.15 + 38.43, # K
                    su_H_p = 51.75*1e3, # Pa
                    su_H_m_dot = 32.85 # kg/s                    
                    )



HTX.set_inputs(
                fouling = 0, HTX_Type = 'Shell&Tube', geom_input = ['Tube_OD','L_shell','Central_Spacing','Shell_ID']
                )

HTX.opt_size()
