# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

@author: Basile
"""

#%% Imports

from machine.examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import IHX_CO2_HP, IHX_EXP_CO2_HP
from machine.examples.CO2_Transcritical_Circuits.CO2_Transcritical_circuit import REC_CO2_TC
from connector.mass_connector import MassConnector

import numpy as np
from CoolProp.CoolProp import PropsSI
from pyswarm import pso
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm

#%%

class CO2RCOptimizer(object):
    
    def __init__(self, fluid):
        
        self.fluid = fluid
        
        self.HP = None
        self.RC = None
        
        self.inputs = {}
        self.params = {}
        self.it_var = {}
        self.obj = {}
        
        self.HSource = MassConnector()
        self.CSource = MassConnector()

#%% Set Methods

    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
        for key, value in parameters.items():
            self.params[key] = value

    def set_it_var(self, **parameters):
        for key, value in parameters.items():
            self.it_var[key] = value    

    def set_obj(self, **parameters):
        for key, value in parameters.items():
            self.obj[key] = value   

#%% Cycle Set

    def set_RC(self):
        
        self.HSource.set_properties(m_dot = self.it_var['mdot_HS'])  # 0.1 # 62.5
                        
        P_sat_T_CSource = PropsSI('P', 'T', self.CSource.T, 'Q', 0.5, self.fluid)
        P_crit_CO2 = PropsSI('PCRIT', self.fluid)

        P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)  

        self.RC = REC_CO2_TC(self.HSource, self.CSource, self.params['PP_gh'], self.params['PP_rec'], self.params['eta_pp'], self.params['eta_exp'], self.params['eta_gh'], 
                              self.params['eta_rec'], self.params['PP_cd'], self.params['SC_cd'], P_low_guess, self.it_var['P_high'], self.it_var['mdot'], mute_print_flag = 1)

        return # PP_rec
        
#%% Cycle Optimization
    
    def system_RC(self, x):
        
        self.it_var['P_high'] = x[0]
        self.it_var['mdot'] = x[1]
        self.it_var['mdot_HS'] = x[1]*x[2]
        
        # print(self.it_var['P_high'])
        # print(self.it_var['mdot'])
        # print(self.it_var['mdot_HS'])
        
        self.set_RC()
        
        try:       
            CO2_RC = self.RC
    
            CO2_RC.solve()

            # Pump

            DP = 50*1e3
            rho = CO2_RC.components['GasHeater'].model.su_H.D
            mdot = CO2_RC.components['GasHeater'].model.su_H.m_dot
            h_ex = CO2_RC.components['GasHeater'].model.ex_H.h
            
            # Cooling Tower
            
            T_amb = 15+273.15
            h_ex_req = PropsSI('H', 'T', T_amb, 'P', 101325, 'Water')
            
            Q_dot_req = mdot*(h_ex - h_ex_req)
            self.W_dot_fan = 0.2*Q_dot_req
            
            eta_pp = 0.8
            self.pp_power = DP*mdot/(rho*eta_pp)
            
            self.W_dot_net = CO2_RC.components['Expander'].model.W.W_dot - CO2_RC.components['Pump'].model.W.W_dot - self.pp_power
    
            self.eta = (self.W_dot_net) / CO2_RC.components['GasHeater'].model.Q_dot.Q_dot 
            

        except:
            return 1000
              
        if not CO2_RC.converged:
            return 100
              
        self.penalty_1 = 0
        # self.penalty_2 = 0
        
        # if abs((CO2_RC.components['GasCooler'].model.ex_C.T - self.obj['T_high'])/self.obj['T_high']) > 1e-2:
        #     self.penalty_1 = abs((CO2_RC.components['GasCooler'].model.ex_C.T - self.obj['T_high'])/self.obj['T_high'])*1000
            
        if abs((self.W_dot_net - self.obj['W_dot'])/self.obj['W_dot']) > 1e-2:
            self.penalty_1 = abs((self.W_dot_net - self.obj['W_dot'])/self.obj['W_dot'])*10        
        
        self.penalty = self.penalty_1 # + self.penalty_2
          
        return -self.eta + self.penalty

    def opt_RC(self):
        
        self.bounds = (np.array([self.params['P_high_min'], self.params['m_dot_min'], self.params['m_dot_HS_fact_min']]), 
                        np.array([self.params['P_high_max'], self.params['m_dot_max'], self.params['m_dot_HS_fact_max']]))

        # Objective wrapper for 1D input (reshape required by pyswarms)
        def objective_wrapper(x):
            return np.array([self.system_RC(xi) for xi in x])

        # Initialize the optimizer
        optimizer = GlobalBestPSO(
            n_particles=20,
            dimensions=3,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=self.bounds
        )

        # Custom stopping logic
        patience = 5
        tol = 1e-3
        max_iter = 30
        no_improve_counter = 0
        best_cost = np.inf
        convergence_curve = []

        # Optimization loop with progress
        for i in tqdm(range(max_iter), desc="Optimizing", ncols=80):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = optimizer.swarm.best_cost

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
        best_m_dot = optimizer.swarm.best_pos[1]
        best_m_dot_HS = optimizer.swarm.best_pos[2]

        self.system_RC([best_P_high, best_m_dot, best_m_dot_HS])

        print(f"\nOptimal P_high: {best_P_high:.2f} Pa")
        print(f"\nOptimal m_dot: {best_m_dot:.4f} kg/s")
        print(f"Best score: {-best_cost:.4f}")
        
        return

#%% Parameter Setting

Optimizer = CO2RCOptimizer('CO2')

T = 95

if T == 95:
    Optimizer.set_parameters(
        eta_pp = 0.95, # [-]
        eta_gh = 0.95, # [-]
        eta_rec = 0.9, # [-]
        eta_exp = 0.9, # [-]
        PP_cd = 5, # [K]
        PP_gh = 5, # [K]
        PP_rec = 0, # [K]
        SC_cd = 0.1, # [K]
        
        P_high_min = 80*1e5, # [Pa]
        P_high_max = 200*1e5, # [Pa]
        
        m_dot_min = 70, # [kg/s]
        m_dot_max = 90, # [kg/s]
        
        m_dot_HS_fact_min = 0.5, # [-]
        m_dot_HS_fact_max = 1, # [-]
        )
    
    Optimizer.set_it_var(
        P_high = 100*1e5, # [Pa]
        mdot_HS = 70, # [kg/s]
        mdot = 70 # [kg/s]
        )
    
    Optimizer.set_obj(
        W_dot = 1e6, # [W]
        )
    
    Optimizer.CSource.set_properties(
        T = 15 + 273.15, # [K]
        P = 5*1e5, # [Pa]
        fluid = 'Water' # []
        )
    
    Optimizer.HSource.set_properties(
        T = 95 + 273.15, # [K]
        P = 5*1e5, # [Pa]
        fluid = 'Water', # []
        )
    
    Optimizer.set_RC()
    Optimizer.opt_RC()
    
else:
    
    T_vec = np.linspace(100, 150, 6) + 273.15
    eta_vec = []
    P_high_vec = []
    m_dot_vec = []
    m_dot_HS_vec = []
    T_h_ex_vec = []
    
    for T in T_vec:
        Optimizer.set_parameters(
            eta_pp = 0.8, # [-]
            eta_gh = 0.95, # [-]
            eta_rec = 0.9, # [-]
            eta_exp = 0.9, # [-]
            PP_cd = 5, # [K]
            PP_gh = 5, # [K]
            PP_rec = 0, # [K]
            SC_cd = 0.1, # [K]
            
            P_high_min = 80*1e5, # [Pa]
            P_high_max = 200*1e5, # [Pa]
            
            m_dot_HS_fact_min = 0.6, # [-]
            m_dot_HS_fact_max = 1, # [-]
            
            m_dot_min = 10, # [kg/s]
            m_dot_max = 100 # [kg/s]
            )
        
        Optimizer.set_it_var(
            P_high = 100*1e5, # [Pa]
            mdot_HS = 20, # [kg/s]
            mdot = 17 # [kg/s]
            )
        
        Optimizer.set_obj(
            W_dot = 1e6, # [W]
            )
        
        Optimizer.CSource.set_properties(
            T = 15 + 273.15, # [K]
            P = 5*1e5, # [Pa]
            fluid = 'Water' # []
            )
        
        Optimizer.HSource.set_properties(
            T = T, # 130 + 273.15, # [K]
            P = 10*1e5, # [Pa]
            fluid = 'Water', # []
            )
    
        Optimizer.set_RC()
        Optimizer.opt_RC()

        eta_vec.append(Optimizer.eta)
        P_high_vec.append(Optimizer.it_var['P_high'])
        m_dot_vec.append(Optimizer.it_var['mdot'])
        m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
        T_h_ex_vec.append(Optimizer.RC.components['GasHeater'].model.ex_H.T)
        
    import matplotlib.pyplot as plt
    
    # Plot 1: Efficiency vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, eta_vec, linewidth=2)
    plt.title("Efficiency vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Efficiency")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Plot 2: P_high vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, P_high_vec, linewidth=2, label='P_high')
    plt.title("High Pressure vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("P_high (units?)")  # Replace with actual units
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot 3: Mass Flow Rate (mdot) vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, m_dot_vec, linewidth=2, label='mdot')
    plt.title("Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Mass Flow Rate (kg/s)")  # Adjust unit if needed
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot 4: Heat Source Mass Flow Rate (mdot_HS) vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, m_dot_HS_vec, linewidth=2, label='mdot_HS')
    plt.title("Heat Source Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Heat Source Mass Flow Rate (kg/s)")  # Adjust unit if needed
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# eta_vec = [0.06760463125362191, 0.10175411487372459, 0.11253675580506668, 
#            0.12287728838101378, 0.13310575584372045, 0.14309211472631098]

# P_high_vec = [9574334.736041067, 13254944.76096478, 13747434.495208107, 
#            13977935.909633664, 15020414.816827059, 15947783.739003953]

# mdot_vec = [78.4014749812898, 47.20629954763076, 41.82720799121871, 
#             37.89369424554172, 33.90218368813903, 30.971937420801787]

# mdot_HS_vec = [55.95542253387586, 41.61756307479707, 28.580367647727293, 
#            25.719270520004418, 22.25176029471058, 19.404497699132122]
