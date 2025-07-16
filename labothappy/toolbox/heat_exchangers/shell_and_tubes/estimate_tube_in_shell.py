
import numpy as np

# Function to calculate the number of tubes in a row
def tubes_in_row(D_main, P_h, y_pos, tube_radius):
    """ Calculate how many tubes can fit in a row at a given y position """
    # Effective available diameter at this y position
    effective_diameter = 2 * np.sqrt((D_main / 2) ** 2 - y_pos ** 2)

    # Calculate how many tubes fit within this effective diameter
    if y_pos < D_main/2:
        tubes_per_row = int(effective_diameter // P_h)
    
        return tubes_per_row
    
    else:
        return 0

# Function to estimate number of tubes in a staggered tube bank arrangement
def estimate_number_of_tubes(D_main, d_tube, P, config_angle, min_tube_row):
    """
    Estimate how many tubes can fit within a given diameter based on the pitch for a staggered tube bank.
    
    Args:
        D_main (float): Main diameter of the containing circle.
        d_tube (float): Tube diameter.
        P (float): Center-to-center distance (pitch).
        config_angle (float): Configuration angle (in degrees), typically 60 for triangular or 45 for staggered square.
        min_tube_row (int) : Minimum number of tubes in a row
        
    Returns:
        int: Estimated number of tubes.
        list: List containing number of tubes in each row.
    """
    
    # Convert configuration angle to radians
    config_angle_rad = np.pi * config_angle / 180
    
    # Vertical and horizontal pitches based on the configuration angle
    if config_angle == 45 or config_angle == 60:
        P_v = 2 * P * np.sin(config_angle_rad)  # Vertical pitch
        P_h = 2 * P * np.cos(config_angle_rad)  # Horizontal pitch
    elif config_angle == 0 or config_angle == 90:
        P_v = P   # Vertical pitch
        P_h = P   # Horizontal pitch
    else:
        print("Config angle not set to 0, 45, 60, 90 degrees")
        return None

    # Tube radius
    tube_radius = d_tube / 2
    
    # Total number of tubes
    tubes_total = 0
    
    # Store number of tubes per row
    tubes_per_row = []

    # Start at the top (y_pos = 0) and move row by row vertically down
    y_pos = P
    
    # Iterate row by row while staying within the half diameter of the circle
    while y_pos < D_main / 2:
        # Calculate number of tubes in the current row
        tubes_in_current_row = tubes_in_row(D_main, P_h, y_pos, tube_radius)
        
        if tubes_in_current_row >= min_tube_row:  # Only add rows with at least one tube
            tubes_per_row.append(tubes_in_current_row)
            tubes_total += tubes_in_current_row
            
            # If staggered (triangular), add a staggered row below
            if config_angle != 90 and config_angle != 0:
                tubes_in_next_row = tubes_in_row(D_main, P_h, y_pos + P_v / 2, tube_radius)
                if tubes_in_next_row >= min_tube_row and y_pos < D_main/2:
                    tubes_per_row.append(tubes_in_next_row)
                    tubes_total += tubes_in_next_row

        # Move to the next row down by vertical pitch
        y_pos += P_v
    
    return tubes_total*2, tubes_per_row[::-1] + tubes_per_row# -*- coding: utf-8 -*-


if __name__ == "__main__":

    test = 0
    
    if test == 1:
        # Example inputs
        clearance = 0.017
        D_main =  - clearance + 0.743  # Main diameter
        d_tube = 0.015875    # Tube diameter
        P = d_tube*1.4        # Pitch
        config_angle = 45  # Triangular (60) / Square (45) arrangement
        min_tube_row = 8
    
        # Estimate how many tubes fit in half of the shell (consider 2 passes)
        num_tubes, tubes_per_row = estimate_number_of_tubes(D_main, d_tube, P, config_angle,min_tube_row)
    
        print(f"Estimated number of tubes: {num_tubes}")
        print(f"Tubes per row : {tubes_per_row}")
    
    if test == 2:
        # Example inputs
        clearance = 0.017
        D_main =  - clearance + 0.81  # Main diameter
        d_tube = 0.015    # Tube diameter
        P = 0.0187       # Pitch
        config_angle = 60  # Triangular (60) / Square (45) arrangement
        min_tube_row = 1
    
        # Estimate how many tubes fit in half of the shell (consider 2 passes)
        num_tubes, tubes_per_row = estimate_number_of_tubes(D_main, d_tube, P, config_angle,min_tube_row)
    
        print(f"Estimated number of tubes: {num_tubes}")
        print(f"Tubes per row : {tubes_per_row}")


