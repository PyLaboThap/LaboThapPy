
import __init__
import numpy as np
from scipy.interpolate import interp1d
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

def cordier_line():

    spec_speed = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                            2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30])

    spec_diam = np.array([9.66, 6.98, 5.65, 4.7, 4.05, 3.57, 3.25,
                          3.01, 1.89, 1.57, 1.40, 1.26, 1.19, 1.11,
                          1.04, 0.97, 0.93, 0.64, 0.48])

    line = interp1d(spec_speed,spec_diam, fill_value='extrapolate')

    return line

def pump_0D_design(Omega, Q, h_1, h_2, v_1, v_2, p_1, p_2, rho_1):

    g = 9.81 # m/s
    cordier = cordier_line()

    # Height difference
    DH = (p_2 - p_1)/(rho_1*g) + (h_2 - h_1) + (v_2**2 - v_1**2)/(2*g)

    # Specific speed
    Omega_s = (Omega*np.sqrt(Q))/((g*DH)**(3/4))

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
    D = D_s*((Q)**(1/2))/(g*DH)**(1/4) # [m]

    return Omega_s, D_s, D, geom

case = 'Torrecid'

if case == 'Cuerva':
    # Working fluid
    fluid = 'Cyclopentane'

    # Rotating speed
    Omega_pp = 3000 # RPM
    Omega_pp_Hz = Omega_pp/60 # Hz
    Omega_pp_rads = Omega_pp_Hz*(2*np.pi) # rad/s

    # ORC pump assumptions for heights and velocities
    h_1 = 0 # m
    h_2 = 0 # m

    v_1 = 0 # m/s
    v_2 = 0 # m/s

    # Static pressures
    p_1 = 53.2*1e3 # Pa
    p_2 = 1330*1e3 # Pa

    # Temperatures and densities
    T_1 = 30 + 273.15 # K
    T_2 = 30.8 + 273.15 # K

    rho_1 = PropsSI('D', 'P', p_1, 'T', T_1, fluid) # kg/m^3
    rho_2 = PropsSI('D', 'P', p_2, 'T', T_2, fluid) # kg/m^3

    # Flowrate
    m_dot = 46.18 # kg/s
    Q = m_dot/rho_1 # m^3/s

    Omega_s,D_s,D,geom = pump_0D_design(Omega_pp_rads, Q, h_1, h_2, v_1, v_2, p_1, p_2, rho_1)

elif case == "Zorlu":
    # Working fluid
    fluid = 'Cyclopentane'

    # Rotating speed
    Omega_pp = 3000 # RPM
    Omega_pp_Hz = Omega_pp/60 # Hz
    Omega_pp_rads = Omega_pp_Hz*(2*np.pi) # rad/s

    # ORC pump assumptions for heights and velocities
    h_1 = 0 # m
    h_2 = 0 # m

    v_1 = 0 # m/s
    v_2 = 0 # m/s

    # Static pressures
    p_1 = 56.82*1e3 # Pa
    p_2 = 917.8*1e3 # Pa

    # Temperatures and densities
    T_1 = 31.7 + 273.15 # K
    T_2 = 32.24 + 273.15 # K

    rho_1 = PropsSI('D', 'P', p_1, 'T', T_1, fluid) # kg/m^3
    rho_2 = PropsSI('D', 'P', p_2, 'T', T_2, fluid) # kg/m^3

    # Flowrate
    m_dot = 34.51 # kg/s
    Q = m_dot/rho_1 # m^3/s

    Omega_s,D_s,D,geom = pump_0D_design(Omega_pp_rads, Q, h_1, h_2, v_1, v_2, p_1, p_2, rho_1)
  
elif case == "Torrecid":
    # Working fluid
    fluid = 'Cyclopentane'

    # ORC pump assumptions for heights and velocities
    h_1 = 0 # m
    h_2 = 0 # m

    v_1 = 0 # m/s
    v_2 = 0 # m/s

    # Static pressures
    p_1 = 53.63*1e3 # Pa
    p_2 = 3390*1e3 # Pa

    # Temperatures and densities
    T_1 = 30.8 + 273.15 # K
    T_2 = 32.6 + 273.15 # K

    rho_1 = PropsSI('D', 'P', p_1, 'T', T_1, fluid) # kg/m^3
    rho_2 = PropsSI('D', 'P', p_2, 'T', T_2, fluid) # kg/m^3

    # Flowrate
    m_dot = 3.58 # kg/s
    Q = m_dot/rho_1 # m^3/s

    Omega_pp_vec = np.linspace(3000, 3000, 100)
    Omega_s_vec = np.zeros(len(Omega_pp_vec))
    D_s_vec = np.zeros(len(Omega_pp_vec))
    D_vec = np.zeros(len(Omega_pp_vec))

    for i in range(len(Omega_pp_vec)):
        # Rotating speed
        Omega_pp = Omega_pp_vec[i] # RPM
        Omega_pp_Hz = Omega_pp/60 # Hz
        Omega_pp_rads = Omega_pp_Hz*(2*np.pi) # rad/s

        Omega_s_vec[i],D_s_vec[i],D_vec[i],geom = pump_0D_design(Omega_pp_rads, Q, h_1, h_2, v_1, v_2, p_1, p_2, rho_1)

    plt.plot(Omega_pp_vec, Omega_s_vec)
    plt.show()

    plt.plot(Omega_pp_vec, D_s_vec)
    plt.show()

    plt.plot(Omega_pp_vec, D_vec)
    plt.show()

# print(f"Omega_s: {Omega_s}")
# print(f"D_s: {D_s}")
# print(f"D: {D}")
# print(f"geom: {geom}")
