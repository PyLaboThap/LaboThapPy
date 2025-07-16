# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

@author: Basile
"""

#%% Imports

from machine.examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import IHX_CO2_HP, IHX_EXP_CO2_HP
from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from pyswarm import pso

#%% Optimization Parmeters

study_case = "NoExp"

def opt_CO2_HP(x):
    
    study_case = "NoExp"
    
    T_cold_source = 0.1+273.15
    T_hot_source = 15+273.15

    eta_is_cp = 0.7
    eta_gc = 0.9
    eta_IHX = 0.7
    eta_exp = 0.8
    
    m_dot = 0.04

    PPTD_ev = 5
    SH_ev = 0.1
    
    P_high = x[0]
    
    P_sat_T_CSource = PropsSI('P', 'T', T_cold_source, 'Q', 0.5, 'CO2')
    P_crit_CO2 = PropsSI('PCRIT', 'CO2')

    P_low_guess = 0.8*min(P_sat_T_CSource, P_crit_CO2)

    HSource = MassConnector()
    HSource.set_properties(fluid='Water', T=T_hot_source,
                           p=5e5, m_dot=0.625*m_dot)  # 0.1 # 62.5

    CSource = MassConnector()
    CSource.set_properties(fluid='Water', T=T_cold_source,
                           p=5e5, m_dot=1000*m_dot)  # 1000 # 625000
    
    try:
        if study_case == "Exp":
            CO2_HP = IHX_EXP_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, eta_exp, PPTD_ev, SH_ev, P_low_guess, P_high, m_dot) # 0.16 # 100
        else:
            CO2_HP = IHX_CO2_HP(HSource, CSource, eta_is_cp, eta_gc, eta_IHX, PPTD_ev, SH_ev, P_low_guess, P_high, m_dot, print_flag =0) # 0.16 # 100
    
        CO2_HP.mute_print()
    
        CO2_HP.solve()
    
        if study_case == "Exp":
            COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/(CO2_HP.components['Compressor'].model.W.W_dot - CO2_HP.components['Expander'].model.W_exp.W_dot)
        else:
            COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/CO2_HP.components['Compressor'].model.W.W_dot
    except:
        return 100    
      
    return -COP


# Set bounds for P_high (in Pa)
P_high_min = 100*1e5
P_high_max = 200*1e5

# Run PSO
best_x, best_neg_COP = pso(opt_CO2_HP, [P_high_min], [P_high_max], swarmsize=5, maxiter=20)

# Display result
best_COP = -best_neg_COP
print(f"Optimal P_high: {best_x[0]:.2f} Pa")
print(f"Maximum COP: {best_COP:.4f}")

