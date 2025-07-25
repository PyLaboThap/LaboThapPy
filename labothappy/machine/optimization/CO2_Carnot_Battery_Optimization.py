# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

@author: Basile
"""

#%% Imports

from machine.examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import IHX_CO2_HP, IHX_EXP_CO2_HP
from connector.mass_connector import MassConnector

import numpy as np
from CoolProp.CoolProp import PropsSI
from pyswarm import pso
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm

#%% Optimization Parmeters

study_case = "NoExp"

m_dot = 0.04

def opt_CO2_HP(x):
    
    study_case = "NoExp"
    
    T_hot_set_point = 130+273.15
    
    T_cold_source = 0.1+273.15
    T_hot_source = 15+273.15

    eta_is_cp = 0.7
    eta_gc = 0.9
    eta_IHX = 0.7
    eta_exp = 0.8

    PPTD_ev = 5
    SH_ev = 0.1
    
    P_high = x[0]
    m_dot_HS = x[1]
    
    P_sat_T_CSource = PropsSI('P', 'T', T_cold_source, 'Q', 0.5, 'CO2')
    P_crit_CO2 = PropsSI('PCRIT', 'CO2')

    P_low_guess = 0.5*min(P_sat_T_CSource, P_crit_CO2)

    HSource = MassConnector()
    HSource.set_properties(fluid='Water', T=T_hot_source,
                           p=5e5, m_dot=m_dot_HS)  # 0.1 # 62.5
    
    try:
        if study_case == "Exp":
            CO2_HP = IHX_EXP_CO2_HP(HSource, eta_is_cp, eta_gc, eta_IHX, eta_exp, PPTD_ev, SH_ev, P_low_guess, P_high, m_dot) # 0.16 # 100
        else:
            CO2_HP = IHX_CO2_HP(HSource, eta_is_cp, eta_gc, eta_IHX, PPTD_ev, SH_ev, P_low_guess, P_high, m_dot, print_flag = 0) # 0.16 # 100
    
        CO2_HP.mute_print()
    
        CO2_HP.solve()
    
        if study_case == "Exp":
            COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/(CO2_HP.components['Compressor'].model.W.W_dot - CO2_HP.components['Expander'].model.W_exp.W_dot)
        else:
            COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/CO2_HP.components['Compressor'].model.W.W_dot
    except:
        return 100      
      
    if CO2_HP.components['GasCooler'].model.ex_C.T - T_hot_set_point < 0:
        penalty = abs(CO2_HP.components['GasCooler'].model.ex_C.T - T_hot_set_point)
    else:
        penalty = 0
      
    return -COP + penalty, CO2_HP

# Define bounds
P_high_min = 100 * 1e5
P_high_max = 200 * 1e5

m_dot_HS_min = m_dot*0.5
m_dot_HS_max = m_dot*1

bounds = (np.array([P_high_min, m_dot_HS_min]), np.array([P_high_max, m_dot_HS_max]))

# Objective wrapper for 1D input (reshape required by pyswarms)
def objective_wrapper(x):
    return np.array([opt_CO2_HP(xi) for xi in x])

# Initialize the optimizer
optimizer = GlobalBestPSO(
    n_particles=5,
    dimensions=2,
    options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
    bounds=bounds
)

# Custom stopping logic
patience = 5
tol = 1e-3
max_iter = 10
no_improve_counter = 0
best_cost = np.inf
convergence_curve = []

# Optimization loop with progress
for i in tqdm(range(max_iter), desc="Optimizing", ncols=80):
    optimizer.optimize(objective_wrapper, iters=1, verbose=False)
    current_best, CO2_HP = optimizer.swarm.best_cost

    convergence_curve.append(current_best)

    if current_best < best_cost - tol:
        best_cost = current_best
        no_improve_counter = 0
    else:
        no_improve_counter += 1

    print(f"[{i+1:03}] Best cost: {best_cost:.6f}")

    if no_improve_counter >= patience:
        print("Stopping early due to stagnation.")
        break

# Final result
best_P_high = optimizer.swarm.best_pos[0]
best_m_dot_HS = optimizer.swarm.best_pos[1]

print(f"\nOptimal P_high: {best_P_high:.2f} Pa")
print(f"\nOptimal m_dot_HS: {best_m_dot_HS:.4f} kg/s")
print(f"Best score: {-best_cost:.4f}")


