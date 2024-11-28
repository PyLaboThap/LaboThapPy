# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:44:25 2023

@author: Basile
"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os

def nozzle(m_dot, T_in, cp_in, v_in, P_in, A1, A2):
    """
    Inputs :
        - m_dot : air mass flowrate [kg/s]
        - T_in : inlet air temperature [K]
        - cp_in : inlet air heat capacity [J/(kg*K)]
        - v_in : air inlet speed [m/s]
        - P_in : air inlet pressure [Pa] 
        - A1 : inlet area [m^2]
        - A2 : outlet area [m^2]
   -------------------------------------------------------
    Outputs :    
        - u_2 : air outlet speed [m/s]
        - T_2 : outlet air temperature [K]
        - P_2 : air outlet pressure [Pa] 
    """
    
    
    "1) Compressible flow functions"
    
    # # Define the relative path to the Excel file
    # file_name = 'compr_flow.xlsx'
        
    # # Use pandas to read the Excel file
    # data_frame = pd.read_excel(file_name)

    # Get the directory of the Python script
    script_directory = os.path.dirname(os.path.abspath(__file__))

    # Construct the full path to the Excel file
    file_name = os.path.join(script_directory, 'compr_flow.xlsx')

    # Check if the file exists
    if os.path.exists(file_name):
        # Read the Excel file using pandas
        data_frame = pd.read_excel(file_name)
        # Convert the DataFrame to a NumPy array
        data = data_frame.to_numpy()
    else:
        print(f"Error: {file_name} not found.")

    "2) Data"
    
    gamma = 1.4 # [/]
    R = 285.71429 # J/(kg*K)
    
    "2.1) Inlet"
            
    rho_in = P_in/(T_in*R)
    a_in = (gamma*R*T_in)**(1/2)
    
    "2.2) Total Properties"
    
    T_tot = T_in + (v_in**2)/(2*cp_in)
    rho_tot = rho_in/((T_in/T_tot)**(1/(gamma-1)))
    P_tot = rho_tot*R*T_tot
    
    # Critical area
    A_crit = m_dot/(rho_in*a_in)
    
    "2.3) Outlet properties"
    
    # Data Interpolation
    f_M = interp1d(data[:,3], data[:,0])
    f_T = interp1d(data[:,3], data[:,2])
    f_P = interp1d(data[:,3], data[:,1])
    
    # Relative area
    A_r_2 = A2/A_crit
    
    # Outlet State
    P_2 = P_tot*f_P(A_r_2)
    T_2 = T_tot*f_T(A_r_2)
    
    # Outlet Speed
    a_2 = (gamma*R*T_2)**(1/2)
    M_2 = f_M(A_r_2)
    u_2 = M_2*a_2
    
    return (u_2, T_2, P_2)
