import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

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
    
    return tubes_total*2, tubes_per_row[::-1] + tubes_per_row

#%% PIPE THICKNESS RELATED FUN

def Internal_Max_P_carbon_steel(D_o,t,T_tube):
    """
    Inputs
    ----------
        - D_o : Input outer diameter [m]
        - t : Tube Thickness [m]
        - T_tube : Input tube temperature [K]s
    
    Outputs
    -------
        - P_max_calc : Maximum allowable pressure [Pa]
        
    Reference
    ---------
    2007 ASME BPV Code + ASME B31
    
    """
    # Interpolation data from engineering Tooblox
    # For steel
    
    T_S_interp = np.array([0, 93.33, 204.444, 315.556, 371.111,
                          398.889, 426.667]) + 273.15  # [K] : Temperature vector
    # [MPa] : Max stress with respect to temperature vector
    S_interp = np.array([158.57942, 158.57942, 158.57942,
                        134.44777, 131.00039, 103.42136, 82.737088])
    
    S_fun = interp1d(T_S_interp, S_interp, kind='linear')

    "Compute P_max for inputs"
    
    S_tube_calc = S_fun(T_tube)*1e6  # [Pa]

    P_max = S_tube_calc*((2*t - 0.01*D_o)/(D_o - (t-0.005*D_o)))
    
    if D_o/t > 20: # Thin wall : Barlow
        P_max = (2 * S_tube_calc * t) / D_o
    else: # Thick wall : Lame
        r_o = D_o/2
        r_i = D_o/2 - t
        P_max = S_tube_calc * ((r_o**2 - r_i**2) / (r_o**2 + r_i**2))
    
    return P_max

def External_Max_P_carbon_steel(D_o,t,T_tube):
    """
    Inputs
    ----------
        - D_o : Input outer diameter [m]
        - t : Tube Thickness [m]
        - T_tube : Input tube temperature [K]
    
    Outputs
    -------
        - P_max_calc : Maximum allowable pressure [Pa]
        
    Reference
    ---------
    2007 ASME BPV Code 
    
    """
    
    sigma = 250*1e6 # Pa : Yield strength for steel
    E = 200*1e9 # Pa : Yield strength for steel
        
    """
    
    Max allowable pressure depending on pipe outside diameter and thickness
    If under critical pressure, associated saturation temperature
    
    """

    "Compute P_max for inputs"
    
    mu = 0.3 # Poisson Ratio

    # from in to m
    D_o = D_o

    r = D_o/2

    P_max_el = 2*E/(1-mu**2) * (t/r)**3
    P_max_pl = 2*sigma/3 * t/r    

    P_max = min(P_max_el, P_max_pl)

    return P_max

def carbon_steel_pipe_thickness(D_o_vect, tube_T, ext_p, int_p):    
    # ASTM A312 standard
    
    schedules = ['5S','10S','40S','80S']

    thickness = np.array([
                # [0.035, 0.049, 0.068, 0.095], # 1/8
                # [0.049, 0.065, 0.088, 0.119], # 1/4
                [0.049, 0.065, 0.091, 0.126], # 3/8
                [0.065, 0.083, 0.109, 0.147], # 1/2
                [0.065, 0.083, 0.111, 0.150], # 5/8
                [0.065, 0.083, 0.113, 0.154], # 3/4
                # [0.065, 0.109, 0.133, 0.179], # 1
                # [0.065, 0.109, 0.140, 0.191], # 1 + 1/4
                # [0.065, 0.109, 0.145, 0.2  ]  # 1 + 1/2
                ])*25.4*1e-3 # m

    thickness_df = pd.DataFrame(index = D_o_vect, columns = schedules, data = thickness)

    for D_o in thickness_df.index: 
        for standard in thickness_df.columns:
            P_max_int = Internal_Max_P_carbon_steel(pd.to_numeric(D_o, errors='coerce')*25.4*1e-3, thickness_df[standard][D_o], tube_T)
            P_max_ext = External_Max_P_carbon_steel(pd.to_numeric(D_o, errors='coerce')*25.4*1e-3, thickness_df[standard][D_o], tube_T)

            # print("-----------------------")

            # print(f"P_max_int : {P_max_int}")
            # print(f"P_max_ext : {P_max_ext}")

            # print(f"ext_p : {ext_p}")
            # print(f"int_p : {int_p}")

            if int_p > P_max_int or ext_p > P_max_ext:
                thickness_df[standard][D_o] = 10000     

    thickness_dic = {}

    for i in range(len(D_o_vect)):
        D_o = D_o_vect[i]
        thickness_dic[str(D_o)] = min(thickness_df.loc[D_o].values)

    return thickness_dic

#%%

def pitch_ratio_fun(D_o, Layout_angle_deg):

    if Layout_angle_deg == 45 or Layout_angle_deg == 0: # Square arrangement 
        if D_o == 1/8:
            Pitch_ratio = 1.25  
        elif D_o == 1/4:
            Pitch_ratio = 1.25  
        elif D_o == 3/8:
            Pitch_ratio = 1.33    
        elif D_o == 1/2:
            Pitch_ratio = 1.25
        elif D_o == 5/8:
            Pitch_ratio = 1.25
        elif D_o == 3/4:
            Pitch_ratio = (1)/D_o
        elif D_o == 1:
            Pitch_ratio = (1+1/4)/D_o
        elif D_o == 1+1/4:
            Pitch_ratio = (1+9/16)/D_o
        elif D_o == 1+1/2:
            Pitch_ratio = (1+7/8)/D_o
        else:
            print("This outer diameter is not considered.")
            return -1

    elif Layout_angle_deg == 30 or Layout_angle_deg == 60: # Square arrangement 
        if D_o == 1/8:
            Pitch_ratio = 1.25  
        elif D_o == 1/4:
            Pitch_ratio = 1.25  
        elif D_o == 3/8:
            Pitch_ratio = 1.33     
        elif D_o == 1/2:
            Pitch_ratio = 1.25
        elif D_o == 5/8:
            Pitch_ratio = 1.25
        elif D_o == 3/4:
            Pitch_ratio = (15/16)/D_o
        elif D_o == 1:
            Pitch_ratio = (1+1/4)/D_o
        elif D_o == 1+1/4:
            Pitch_ratio = (1+9/16)/D_o
        elif D_o == 1+1/2:
            Pitch_ratio = (1+7/8)/D_o
        else:
            print("This outer diameter is not considered.")
            return -1
    else:
        print("This tube arrangement is not considered.")
        return -1

    return Pitch_ratio

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

