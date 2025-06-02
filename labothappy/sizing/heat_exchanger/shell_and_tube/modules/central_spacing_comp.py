#%% CENTRAL SPAC RELATED FUN

import numpy as np
import random

# Function to find divisors between bounds for floating-point numbers
def find_divisors_between_bounds(num_to_divide, lower_bound, upper_bound, tolerance=1e-5):
    """ Finds floating-point divisors of num_to_divide within a specific range. """
    valid_divisors = []
    
    # Loop through potential divisors in small increments between lower and upper bounds
    potential_divisors = np.arange(round(lower_bound, 2), round(upper_bound, 2), 0.001)  # Adjust step size as needed

    for divisor in potential_divisors:
        divisor_rounded = round(divisor, 3)  # Use sufficient precision
        if abs(np.mod(num_to_divide, divisor_rounded)) < tolerance or abs(np.mod(num_to_divide, divisor_rounded) - divisor_rounded) < tolerance:  # Check if the remainder is within the tolerance
            valid_divisors.append(divisor_rounded)
    
    return valid_divisors

# Function to generate a random divisor between two bounds
def random_divisor_between_bounds(num_to_divide, lower_bound, upper_bound):
    valid_divisors = find_divisors_between_bounds(num_to_divide, lower_bound, upper_bound)
    
    if valid_divisors:
        # Choose a random divisor from the valid divisors list
        return random.choice(valid_divisors)
    else:
        return None  # No divisors found in the given range

# Test case
test = 0

if test == 1: 
    low = 0.15747999999999998
    high = 2.547615488977596
    To_divide = 13.86

    divisors = find_divisors_between_bounds(To_divide, low, high)

    print(f"N Divisors: {len(divisors)}")
    print(f"N Divisors: {divisors}")
