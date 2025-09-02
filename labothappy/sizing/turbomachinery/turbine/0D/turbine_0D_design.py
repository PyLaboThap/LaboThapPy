
import __init__
import numpy as np
from scipy.interpolate import interp1d
from CoolProp.CoolProp import PropsSI

def cordier_line():

    spec_speed = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                            2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30])

    spec_diam = np.array([9.66, 6.98, 5.65, 4.7, 4.05, 3.57, 3.25,
                          3.01, 1.89, 1.57, 1.40, 1.26, 1.19, 1.11,
                          1.04, 0.97, 0.93, 0.64, 0.48])

    line = interp1d(spec_speed,spec_diam, fill_value='extrapolate')

    return line

def turbine_0D_design(Omega, Q, h_1, h_2, rho_1):

    g = 9.81 # m/s
    cordier = cordier_line()

    # Enthalpy difference
    Dh = (h_1 - h_2)

    # Specific speed
    Omega_s = (Omega*np.sqrt(Q))/((Dh)**(3/4))

    Omega_s_imp = 2733*Omega_s

    if Omega_s_imp < 1500:
        geom = "Radial-Vane Area"
    elif Omega_s_imp < 4000:
        geom = "Francis-Vane Area"
    elif Omega_s_imp < 7500:
        geom = "Mixed-Flow Area"
    else:
        geom = "Axial-Flow Area"
    
    # Specific Diameter using Cordier Line
    D_s = cordier(Omega_s)

    # Real diameter 
    D = D_s*((Q)**(1/2))/(Dh)**(1/4) # [m]

    return Omega_s, D_s, D, geom

if __name__ == "__main__":
    case_study = 'TCO2'
    
    if case_study == "TCO2":
        # Working fluid
        fluid = 'CO2'
    
        # Rotating speed
        Omega_pp = 1500 # RPM
        Omega_pp_Hz = Omega_pp/60 # Hz
        Omega_pp_rads = Omega_pp_Hz*(2*np.pi) # rad/s
    
        # Static pressures
        p_1 = 140*1e5 # Pa
        p_2 = 40*1e5 # Pa
    
        # Temperatures and densities
        T_1 = 121 + 273.15 # K
        T_2 = 40 + 273.15 # K
    
        rho_1 = PropsSI('D', 'P', p_1, 'T', T_1, fluid) # kg/m^3
        rho_2 = PropsSI('D', 'P', p_2, 'T', T_2, fluid) # kg/m^3
    
        h_1 = PropsSI('H', 'P', p_1, 'T', T_1, fluid) # kg/m^3
        h_2 = PropsSI('H', 'P', p_2, 'T', T_2, fluid) # kg/m^3
    
        # Flowrate
        m_dot = 500 # kg/s
        Q = m_dot/rho_1 # m^3/s
    
        Omega_s,D_s,D,geom = turbine_0D_design(Omega_pp_rads, Q, h_1, h_2, rho_1)
        
      
print(f"Omega_s: {Omega_s}")
print(f"D_s: {D_s}")
print(f"D: {D}")
print(f"geom: {geom}")
