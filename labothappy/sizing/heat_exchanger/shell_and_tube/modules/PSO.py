# import numpy as np

# def particle_swarm_optimization(objective_function, bounds, num_particles=30, num_dimensions=2, max_iterations=100, 
#                                 inertia_weight=0.7, cognitive_constant=1.5, social_constant=1.5, discrete_indices=None, discrete_values=None):
#     """
#     Perform Particle Swarm Optimization (PSO) to minimize the given objective function.
    
#     Parameters:
#     - objective_function: Function to be minimized. Should take a particle position (array-like) as input.
#     - bounds: List of tuples specifying (min, max) for each dimension.
#     - num_particles: Number of particles in the swarm (default 30).
#     - num_dimensions: Number of dimensions of the search space (default 2).
#     - max_iterations: Maximum number of iterations (default 100).
#     - inertia_weight: Inertia weight to balance exploration and exploitation (default 0.7).
#     - cognitive_constant: Constant to control personal best influence (default 1.5).
#     - social_constant: Constant to control global best influence (default 1.5).
#     - discrete_indices: List of indices of variables that should be discrete.
#     - discrete_values: Dictionary specifying allowed discrete values for each index.
    
#     Returns:
#     - global_best_position: Position of the best solution found.
#     - global_best_score: Value of the objective function at the best solution.
#     """
    
#     # Initialize particle positions and velocities
#     particle_positions = np.random.uniform(low=[b[0] for b in bounds], 
#                                            high=[b[1] for b in bounds], 
#                                            size=(num_particles, num_dimensions))
#     particle_velocities = np.random.uniform(low=-1, high=1, size=(num_particles, num_dimensions))

#     # Initialize personal best positions and global best position
#     personal_best_positions = np.copy(particle_positions)
#     personal_best_scores = np.array([objective_function(p) for p in personal_best_positions])

#     global_best_position = personal_best_positions[np.argmin(personal_best_scores)]
#     global_best_score = np.min(personal_best_scores)

#     # PSO loop
#     for iteration in range(max_iterations):
#         for i in range(num_particles):
#             # Update velocity
#             r1, r2 = np.random.rand(), np.random.rand()  # Random coefficients
#             cognitive_velocity = cognitive_constant * r1 * (personal_best_positions[i] - particle_positions[i])
#             social_velocity = social_constant * r2 * (global_best_position - particle_positions[i])
#             particle_velocities[i] = inertia_weight * particle_velocities[i] + cognitive_velocity + social_velocity

#             # Update position
#             particle_positions[i] += particle_velocities[i]

#             # Enforce boundary conditions and handle discrete variables
#             for d in range(num_dimensions):
#                 # Bound constraints
#                 if particle_positions[i, d] < bounds[d][0]:
#                     particle_positions[i, d] = bounds[d][0]
#                 if particle_positions[i, d] > bounds[d][1]:
#                     particle_positions[i, d] = bounds[d][1]
                
#                 # Apply discrete constraints for specific variables
#                 if discrete_indices and d in discrete_indices:
#                     # Snap the position to the nearest allowed discrete value
#                     allowed_values = discrete_values[d]
#                     particle_positions[i, d] = min(allowed_values, key=lambda x: abs(x - particle_positions[i, d]))

#             # Evaluate the new position
#             new_score = objective_function(particle_positions[i])

#             # Update personal best
#             if new_score < personal_best_scores[i]:
#                 personal_best_positions[i] = particle_positions[i]
#                 personal_best_scores[i] = new_score

#         # Update global best
#         min_personal_best_score = np.min(personal_best_scores)
#         if min_personal_best_score < global_best_score:
#             global_best_position = personal_best_positions[np.argmin(personal_best_scores)]
#             global_best_score = min_personal_best_score

#         # Optionally, print progress
#         print(f"Iteration {iteration+1}/{max_iterations}, Global Best Score: {global_best_score}")

#     return global_best_position, global_best_score


# # Example usage of the general PSO function with discrete variables
# if __name__ == "__main__":
#     # Define the objective function (Example: minimize f(x, y) = x^2 + y^2)
#     def objective_function(x):
#         return x[0]**2 + x[1]**2

#     # Define the bounds for each dimension (Example: for x and y between -10 and 10)
#     bounds = [(-10, 10), (-10, 10)]

#     # Define the discrete constraints
#     discrete_indices = [1]  # The second dimension (y) will have discrete values
#     discrete_values = {
#         1: [0, 1, 2]  # Allowed discrete values for the second dimension (y)
#     }

#     # Call the PSO function with discrete constraints
#     best_position, best_score = particle_swarm_optimization(objective_function, bounds, num_particles=50, num_dimensions=2, 
#                                                             max_iterations=200, discrete_indices=discrete_indices, discrete_values=discrete_values)

#     # Output the results
#     print("Best position found:", best_position)
#     print("Best score found:", best_score)

import numpy as np

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


# Example usage with constraints
if __name__ == "__main__":
    # Define the objective function (Example: minimize f(x, y) = x^2 + y^2)
    def objective_function(x):
        return x[0]**2 + x[1]**2

    # Define the bounds for each dimension (Example: for x and y between -10 and 10)
    bounds = [(-10, 10), (-10, 10)]

    # Define constraints
    # Example: constraint 1: x + y <= 5 (inequality constraint)
    # Example: constraint 2: x - y = 1 (equality constraint can be treated as inequality with tolerance)
    def constraint1(x):
        return x[0] + x[1] - 5  # x + y <= 5 -> x + y - 5 <= 0

    def constraint2(x):
        return abs(x[0] - x[1] - 1) - 0.001  # x - y = 1, treat as |x - y - 1| <= 0.001

    constraints = [constraint1, constraint2]

    # Call the PSO function with constraints
    best_position, best_score = particle_swarm_optimization(objective_function, bounds, num_particles=50, num_dimensions=2, 
                                                            max_iterations=200, constraints=constraints)

    # Output the results
    print("Best position found:", best_position)
    print("Best score found:", best_score)

