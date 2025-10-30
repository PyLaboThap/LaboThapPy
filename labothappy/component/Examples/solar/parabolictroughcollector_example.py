
import __init__
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

#import component.Examples.solar.parabolictroughcollector_example as parabolictroughcollector_example
# Muted by Titouan, do not work without this line on Elise computer 
from toolbox.geometries.solar.parabolictrough_geometry import PT_Collector_Geom
from component.solar.parabolictroughcollector import PT_collector

PT_geom = PT_Collector_Geom()
PT_geom.set_parameters("Soponova_MicroCSP")

case = "plant_sizing"

if case == "study_disc":
    n_disc_vec = np.linspace(1,50,50)
    Q_dot_vec = np.zeros(len(n_disc_vec))

    for i in range(len(n_disc_vec)):

        PT = PT_collector()

        PT.set_parameters(coll_eff = PT_geom.coll_eff, L = 10*PT_geom.L, W = PT_geom.W, A = PT_geom.A, 
                            A_r = PT_geom.A_r, m = PT_geom.m, L_f = PT_geom.L_f,

                            Tube_OD = PT_geom.Tube_OD, Tube_V = PT_geom.Tube_V,

                            alpha_r = PT_geom.alpha_r, refl_m = PT_geom.refl_m, epsilon_r = PT_geom.epsilon_r, envel_tau = PT_geom.envel_tau, eta_opt = PT_geom.eta_opt,
                            
                            V_w_max_tr = PT_geom.V_w_max_tr, V_w_max_st = PT_geom.V_w_max_st, Vdot_min = PT_geom.Vdot_min,
                            Vdot_max = PT_geom.Vdot_max, T_f_min = PT_geom.T_f_min, T_f_max = PT_geom.T_f_max,

                            a = PT_geom.a, n_disc = int(n_disc_vec[i])
        )

        PT.set_inputs(fluid = 'Water',
                    m_dot = 0.6, # kg/s
                    P_su = 5*1e5, # Pa
                    T_su = 100 + 273.15, # K
                    T_amb = 25 + 273.15, # K
                    DNI = 900, # W/m^2
                    Theta = 10*np.pi/180, # rad
                    v_wind = 5  # m/s
                    )

        PT.solve()

        Q_dot_vec[i] = PT.Q_dot

    plt.plot(n_disc_vec, Q_dot_vec*1e-3)
    plt.grid()
    plt.xlabel("Discretizations")
    plt.ylabel("Absorbed Power [kW]")

if case == "plant_sizing":
    
    Q_dot_need = 14*1e6 # W
    DT_need = 10 # K
    p_circ = 5*1e5 # Pa
    fluid = "Water"

    DT_HX = 10 # K
    T_out_HX = 135 + 273.15 # K

    T_in_HX = T_out_HX + DT_HX # K

    h_in = PropsSI('H', 'T', T_in_HX, 'P', p_circ, fluid)
    h_out = PropsSI('H', 'T', T_out_HX, 'P', p_circ, fluid)

    m_dot_needed = Q_dot_need/(h_in - h_out)
    print(f"m_dot_needed : {m_dot_needed}")

    D_in = PropsSI('D', 'T', T_out_HX, 'P', p_circ, fluid)
    m_dot_max_coll = PT_geom.Vdot_max*D_in

    nb_coll_parallel = m_dot_needed/m_dot_max_coll

    print(f"nb_coll_parallel : {nb_coll_parallel}")

    PT = PT_collector()

    PT.set_parameters(coll_eff = PT_geom.coll_eff, L = 10*PT_geom.L, W = PT_geom.W, A = PT_geom.A, 
                        A_r = PT_geom.A_r, m = PT_geom.m, L_f = PT_geom.L_f,

                        Tube_OD = PT_geom.Tube_OD, Tube_V = PT_geom.Tube_V,

                        alpha_r = PT_geom.alpha_r, refl_m = PT_geom.refl_m, epsilon_r = PT_geom.epsilon_r, envel_tau = PT_geom.envel_tau, eta_opt = PT_geom.eta_opt,
                        
                        V_w_max_tr = PT_geom.V_w_max_tr, V_w_max_st = PT_geom.V_w_max_st, Vdot_min = PT_geom.Vdot_min,
                        Vdot_max = PT_geom.Vdot_max, T_f_min = PT_geom.T_f_min, T_f_max = PT_geom.T_f_max,

                        a = PT_geom.a, n_disc = 20
    )

    T_ex = T_out_HX
    n_series = 0

    while T_ex < T_in_HX:
        T_in_coll = T_ex
        n_series += 1

        PT.set_inputs(
                fluid = 'Water',
                m_dot = m_dot_max_coll, # kg/s
                P_su = p_circ, # Pa
                T_su = T_in_coll, # K
                T_amb = 20 + 273.15, # K
                DNI = 600, # W/m^2
                Theta = 10*np.pi/180, # rad
                v_wind = 5 # m/s
                )

        PT.solve()
        T_ex = PT.ex.T
        print(f"Q_dot {PT.Q_dot}")
        print(f"T_ex {T_ex}")

    print(f"n_series {n_series}")

    Filling_factor = 0.6
    A = np.ceil(nb_coll_parallel)*(np.ceil(n_series))*PT.params['W']*PT.params['L']/Filling_factor
    print(f"Required solar field : {A} [m^2]")
    print(f"Required solar field : {A*1e-6} [km^2]")


