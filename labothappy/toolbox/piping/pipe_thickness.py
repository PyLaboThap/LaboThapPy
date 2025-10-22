import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

#%% PIPE THICKNESS RELATED FUN

def Internal_Max_P_carbon_steel(D_o,t,T_tube):
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
    # ASTM A179 standard (BWG for thickness)
        
    # Define BWG and corresponding thickness (in mm)
    schedules_BWG = [str(int(x)) for x in np.linspace(0, 26, 27)]
    thickness = [0.34, 0.3, 0.284, 0.259, 0.238, 0.22, 0.203, 0.18, 0.165, 0.148, 0.134,
                 0.12, 0.109, 0.095, 0.083, 0.072, 0.065, 0.058, 0.049, 0.042, 0.035,
                 0.032, 0.028, 0.025, 0.022, 0.02, 0.018]
    
    BWG_dict = dict(zip(schedules_BWG, thickness))
    
    # Define pipe size thickness lists (in mm)
    # list_1_4 = [BWG_dict['22'], BWG_dict['24'], BWG_dict['26']] 
    list_3_8 = [BWG_dict['18'], BWG_dict['20'], BWG_dict['22'], BWG_dict['24']]
    list_1_2 = [BWG_dict['16'], BWG_dict['18'], BWG_dict['20'], BWG_dict['22']]
    list_5_8 = [BWG_dict[x] for x in schedules_BWG[12:21]]     # 0.109 to 0.035
    list_3_4 = [BWG_dict[x] for x in schedules_BWG[12:21]]
    list_1   = [BWG_dict[x] for x in schedules_BWG[8:21]]      # 0.165 to 0.035
    list_5_4 = [BWG_dict[x] for x in schedules_BWG[7:21]]      # 0.18 to 0.035
    list_3_2 = [BWG_dict[x] for x in schedules_BWG[10:17]]     # 0.134 to 0.058
    
    # Pad lists to match length of full thickness list
    def pad_list(lst, target_len):
        return lst + [max(lst)] * (target_len - len(lst))

    # list_1_4 = pad_list(list_1_4, len(thickness))    
    list_3_8 = pad_list(list_3_8, len(thickness))
    list_1_2 = pad_list(list_1_2, len(thickness))
    list_5_8 = pad_list(list_5_8, len(thickness))
    list_3_4 = pad_list(list_3_4, len(thickness))
    list_1   = pad_list(list_1, len(thickness))
    list_5_4 = pad_list(list_5_4, len(thickness))
    list_3_2 = pad_list(list_3_2, len(thickness))
    
    # Create final thickness array (in meters)
    thickness_array = np.array([
        # list_1_4,
        list_3_8,
        list_1_2,
        list_5_8,
        list_3_4,
        list_1,
        list_5_4,
        list_3_2
    ]) * 25.4e-3  # Convert mm to meters

    thickness_df = pd.DataFrame(index = D_o_vect, columns = schedules_BWG, data = thickness_array)

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

def carbon_steel_pipe_thickness_mm(D_o, tube_T, ext_p, int_p):    
    # ASTM A179 standard (BWG for thickness)
    D_o_choices = np.array([0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5])*25.4*1e-3
    
    # Define BWG and corresponding thickness (in mm)
    schedules_BWG = [str(int(x)) for x in np.linspace(0, 26, 27)]
    thickness = [0.34, 0.3, 0.284, 0.259, 0.238, 0.22, 0.203, 0.18, 0.165, 0.148, 0.134,
                 0.12, 0.109, 0.095, 0.083, 0.072, 0.065, 0.058, 0.049, 0.042, 0.035,
                 0.032, 0.028, 0.025, 0.022, 0.02, 0.018]
    
    BWG_dict = dict(zip(schedules_BWG, thickness))
    
    # Define pipe size thickness lists (in mm)
    list_1_4 = [BWG_dict['22'], BWG_dict['24'], BWG_dict['26']] 
    list_3_8 = [BWG_dict['18'], BWG_dict['20'], BWG_dict['22'], BWG_dict['24']]
    list_1_2 = [BWG_dict['16'], BWG_dict['18'], BWG_dict['20'], BWG_dict['22']]
    list_5_8 = [BWG_dict[x] for x in schedules_BWG[12:21]]     # 0.109 to 0.035
    list_3_4 = [BWG_dict[x] for x in schedules_BWG[12:21]]
    list_1   = [BWG_dict[x] for x in schedules_BWG[8:21]]      # 0.165 to 0.035
    list_5_4 = [BWG_dict[x] for x in schedules_BWG[7:21]]      # 0.18 to 0.035
    list_3_2 = [BWG_dict[x] for x in schedules_BWG[10:17]]     # 0.134 to 0.058
    
    # Pad lists to match length of full thickness list
    def pad_list(lst, target_len):
        return lst + [max(lst)] * (target_len - len(lst))

    list_1_4 = pad_list(list_1_4, len(thickness))    
    list_3_8 = pad_list(list_3_8, len(thickness))
    list_1_2 = pad_list(list_1_2, len(thickness))
    list_5_8 = pad_list(list_5_8, len(thickness))
    list_3_4 = pad_list(list_3_4, len(thickness))
    list_1   = pad_list(list_1, len(thickness))
    list_5_4 = pad_list(list_5_4, len(thickness))
    list_3_2 = pad_list(list_3_2, len(thickness))
    
    # Create final thickness array (in meters)
    thickness_array = np.array([
        list_1_4,
        list_3_8,
        list_1_2,
        list_5_8,
        list_3_4,
        list_1,
        list_5_4,
        list_3_2
    ]) * 25.4e-3  # Convert mm to meters

    idx = np.argmin(np.abs(D_o_choices - D_o))
    
    t_list = np.sort(thickness_array[idx])
        
    for t in t_list:
        
        P_max_int = Internal_Max_P_carbon_steel(D_o, t, tube_T)
        P_max_ext = External_Max_P_carbon_steel(D_o, t, tube_T)
        
        if P_max_int >= int_p and P_max_ext >= ext_p:
            t_min = t
            return t_min
                
    return 1000

if __name__ == "__main__":
    
    t_test = carbon_steel_pipe_thickness_mm(1.5*25.4*1e-3, 273.15+26, 100*1e5, 100*1e5)

