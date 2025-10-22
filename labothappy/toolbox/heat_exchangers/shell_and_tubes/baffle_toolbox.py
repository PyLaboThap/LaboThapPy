# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 17:42:12 2025

@author: basil
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

def baffle_thickness(D_s, BC, rho_S, T_S):
    """
    Inputs
    ------
    
    D_s : Shell diameter [m]
    BC : Baffle cut [%]
    rho_S : Shell fluid density [kg/m^3]
    T_S : Shell fluid temperature [K]

    Outputs
    -------
    t_B : Baffle thickness [m]
    
    Reference
    ---------
    Pressure Vessel Manual (2004)
    """
    
    "Maximum allowable stress - Carbon Steel"
        
    T_S_interp = np.array([0, 93.33, 204.444, 315.556, 371.111,
                          398.889, 426.667]) + 273.15  # [K] : Temperature vector
    # [MPa] : Max stress with respect to temperature vector
    S_interp = np.array([158.57942, 158.57942, 158.57942,
                        134.44777, 131.00039, 103.42136, 82.737088])
    
    S_fun = interp1d(T_S_interp, S_interp, kind='linear')
    
    S_calc = S_fun(T_S)*1e6  # [Pa] : Young Modulus
    
    "Baffle cord (b) and height (a) computation"
    
    # Define the center angle of the cord equation using a lambda function
    equation = lambda theta: theta - np.sin(theta) - BC * np.pi
    
    # Solve for theta with an initial guess
    theta = fsolve(equation, np.pi / 2)[0]
    
    a = (D_s/2)*(1-np.cos(theta/2)) + D_s/2
    b = D_s*np.sin(theta/2)

    """beta, gamma coefficients
    One short edge free, three edges simply supported, uniformly decreasing load to the free edge """

    ratio_val =   [0.25 ,   0.5,  0.75,    1,  1.5,     2,   2.5,     3,   3.5,    4]
    beta_1_val =  [0.5  ,  0.11,  0.16,  0.2, 0.28,  0.32,  0.35,  0.36,  0.37, 0.37]
    gamma_1_val = [0.013, 0.026, 0.033, 0.04, 0.05, 0.058, 0.064, 0.067, 0.069, 0.07]

    beta_1_fun = interp1d(ratio_val, beta_1_val, kind='linear')
    gamma_1_fun = interp1d(ratio_val, gamma_1_val, kind='linear')

    beta_1 = beta_1_fun(a/b)
    gamma_1 = gamma_1_fun(a/b)

    "Maximum Load"
    # Relative gravity
    rho_w = 1000 # kg/m^3
    S_g = rho_S/rho_w

    p = 62.4*a*S_g/144 # psi
    p = 0.0689475729*p*1e5 # Pa

    "Proposed Baffle Thickness"    
    t_b = np.sqrt(beta_1*p*b**2 / S_calc)

    "Verification of the thickness"
    E = 200*1e9 # Steel Young Modulus

    # Deflection
    delta = (p*gamma_1*b**4)/(E*t_b**3)

    if delta >= t_b/2 and delta >= b/360:
        
        def equation(t_b):
            delta = (p * gamma_1 * b**4) / (E * t_b**3)
            return max(delta - t_b/2, delta - b/360)  # Both constraints should be satisfied
        
        # Initial guess for t_b (should be positive)
        t_b_initial_guess = t_b
        
        # Solve for t_b
        t_b_solution = fsolve(equation, t_b_initial_guess)
        
        # Print result
        t_b = t_b_solution[0] if t_b_solution[0] > 0 else None
        
    t_b_min = 25.4*1e-3*1/16
    
    t_b = max(t_b, t_b_min)

    corrosion_allowance = 1.6*1e-3 # m

    return t_b + corrosion_allowance

#%% CENTRAL SPAC RELATED FUN

import numpy as np
import random

# Function to find divisors between bounds for floating-point numbers
def find_divisors_between_bounds(num_to_divide, lower_bound, upper_bound, step=0.001, tolerance=1e-5):
    """Vectorized, fast version — finds floating-point divisors within bounds."""
    # Generate possible divisors
    potential_divisors = np.arange(lower_bound, upper_bound, step)

    # Compute remainders all at once
    remainders = np.mod(num_to_divide, potential_divisors)

    # Check near-zero or near-divisor remainders (within tolerance)
    mask = (remainders < tolerance) | (np.abs(remainders - potential_divisors) < tolerance)

    # Optionally round for readability when returning
    valid_divisors = np.round(potential_divisors[mask], 3)

    return valid_divisors.tolist()

# Function to generate a random divisor between two bounds
def random_divisor_between_bounds(num_to_divide, lower_bound, upper_bound):
    valid_divisors = find_divisors_between_bounds(num_to_divide, lower_bound, upper_bound)
    
    if valid_divisors:
        # Choose a random divisor from the valid divisors list
        return random.choice(valid_divisors)
    else:
        return None  # No divisors found in the given range

if __name__ == "__main__":
    # Test case
    test = 0
    
    if test == 1: 
        low = 0.15747999999999998
        high = 2.547615488977596
        To_divide = 13.86
    
        divisors = find_divisors_between_bounds(To_divide, low, high)
    
        print(f"N Divisors: {len(divisors)}")
        print(f"N Divisors: {divisors}")
