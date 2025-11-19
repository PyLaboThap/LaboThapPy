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
import CoolProp.CoolProp as CP

from joblib import Parallel, delayed

#%%

class CO2HPOptimizer(object):
    
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

        self.AS = CP.AbstractState('HEOS', fluid)

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

    def set_HP(self):
        
        self.HSource.set_properties(m_dot = self.it_var['mdot_HS'])  # 0.1 # 62.5
                        
        self.AS.update(CP.QT_INPUTS, 0.5,  self.CSource.T)
        
        P_sat_T_CSource = self.AS.p()
        P_crit_CO2 = self.AS.p_critical()

        P_low_guess = 0.8*min(P_sat_T_CSource, P_crit_CO2)
        
        if self.params['HP_ARCH'] == 'IHX_EXP':
            self.HP = IHX_EXP_CO2_HP(self.HSource, self.CSource.T, self.params['eta_cp'], self.params['eta_gc'], self.params['eta_IHX'], self.params['eta_exp'], self.params['PP_ev'], self.params['SH_ev'], P_low_guess, self.it_var['P_high'], self.it_var['mdot'], mute_print_flag = 1)
            
        elif self.params['HP_ARCH'] == 'IHX':
            self.HP = IHX_CO2_HP(self.HSource, self.CSource.T, self.params['eta_cp'], self.params['eta_gc'], self.params['eta_IHX'], self.params['PP_ev'], self.params['SH_ev'], P_low_guess, self.it_var['P_high'], self.it_var['mdot'], mute_print_flag = 1)
        
        return
        
#%% Cycle Optimization
    
    def system_HP(self, x):
        
        self.it_var['P_high'] = x[0]
        self.it_var['mdot'] = x[1]
        self.it_var['mdot_HS'] = x[1]*x[2]
                
        self.set_HP()
        
        try:       
            CO2_HP = self.HP
            
            # CO2_HP.mute_print()
            CO2_HP.solve()
    
            if self.params['HP_ARCH'] == 'IHX_EXP':     
                self.COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/(CO2_HP.components['Compressor'].model.W.W_dot - CO2_HP.components['Expander'].model.W_exp.W_dot)
                
            if self.params['HP_ARCH'] == 'IHX':     
                self.COP = CO2_HP.components['GasCooler'].model.Q_dot.Q_dot/CO2_HP.components['Compressor'].model.W.W_dot
            
        except:
            return 100 
              
        if not CO2_HP.converged:
            print("OUCH")
            return 100  
              
        self.penalty_1 = 0
        self.penalty_2 = 0
        
        if abs((CO2_HP.components['GasCooler'].model.ex_C.T - self.obj['T_high'])/self.obj['T_high']) > 1e-2:
            self.penalty_1 = abs((CO2_HP.components['GasCooler'].model.ex_C.T - self.obj['T_high'])/self.obj['T_high'])*300
            
        if abs((CO2_HP.components['Compressor'].model.W.W_dot - self.obj['W_dot'])/self.obj['W_dot']) > 1e-2:
            self.penalty_2 = abs((CO2_HP.components['Compressor'].model.W.W_dot - self.obj['W_dot'])/self.obj['W_dot'])*100          
        
        self.penalty = self.penalty_1 + self.penalty_2
        self.score = -self.COP + self.penalty
                  
        return self.score

    def opt_HP(self):
        
        self.bounds = (np.array([self.params['P_high_min'], self.params['m_dot_min'], self.params['m_dot_HS_fact_min']]), 
                        np.array([self.params['P_high_max'], self.params['m_dot_max'], self.params['m_dot_HS_fact_max']]))

        # Objective wrapper for 1D input (reshape required by pyswarms)
        def objective_wrapper(x):
            return np.array([self.system_HP(xi) for xi in x])

        # def objective_wrapper(x):
        #     with parallel_backend('threading', n_jobs=-1):
        #         return np.array([self.system_HP(xi) for xi in x])


        # Initialize the optimizer
        self.optimizer = GlobalBestPSO(
            n_particles=50,
            dimensions=3,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=self.bounds
        )

        # Custom stopping logic
        patience = 10
        tol = 1e-3
        max_iter = 30
        no_improve_counter = 0
        best_cost = np.inf
        convergence_curve = []

        # Optimization loop with progress
        for i in tqdm(range(max_iter), desc="Optimizing", ncols=80):
            self.optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = self.optimizer.swarm.best_cost

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
        best_P_high = self.optimizer.swarm.best_pos[0]
        best_m_dot = self.optimizer.swarm.best_pos[1]
        best_m_dot_HS_fact = self.optimizer.swarm.best_pos[2]

        self.system_HP([best_P_high, best_m_dot, best_m_dot_HS_fact])

        print(f"\nOptimal P_high: {best_P_high:.2f} Pa")
        print(f"\nOptimal m_dot: {best_m_dot:.4f} kg/s")
        print(f"\nOptimal m_dot_HS: {best_m_dot_HS_fact*best_m_dot:.4f} kg/s")
        print(f"Best score: {-best_cost:.4f}")
        
        return        

# %% Parameter Setting

Optimizer = CO2HPOptimizer('CO2')

T = 100

if T == 110:
    Optimizer.set_parameters(
        HP_ARCH = 'IHX_EXP', # IHX and IHX_EXP
        eta_cp = 0.8, # [-]
        eta_gc = 0.95, # [-]
        eta_IHX = 0.7, # [-]
        PP_ev = 5, # [K]
        SH_ev = 0.1, # [K]
        eta_exp = 0.7, # [-]
        
        P_high_min = 100*1e5, # [Pa]
        P_high_max = 150*1e5, # [Pa]
        
        m_dot_HS_fact_min = 0.3, # [-]
        m_dot_HS_fact_max = 1, # [-]
            
        m_dot_min = 5, # [kg/s]
        m_dot_max = 25 # [kg/s]
        )
    
    Optimizer.set_it_var(
        P_high = 140*1e5, # [Pa]
        mdot_HS = 20, # [kg/s]
        mdot = 17 # [kg/s]
        )
    
    Optimizer.set_obj(
        T_high = 110 + 273.15, # [K]
        W_dot = 1e6, # [W]
        )
    
    Optimizer.CSource.set_properties(
        T = 15 + 273.15,
        P = 1*1e5,
        fluid = 'Water'
        )
    
    Optimizer.HSource.set_properties(
        T = 15 + 273.15,
        P = 10*1e5,
        fluid = 'Water',
        )
    
    Optimizer.set_HP()
    Optimizer.opt_HP()

else:
    
    T_vec = np.linspace(100, 150, 6) + 273.15
    COP_vec = []
    P_high_vec = []
    m_dot_vec = []
    m_dot_HS_vec = []
    T_act_vec = []
    
    m_dot_HS_fact_min = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2])
    m_dot_HS_fact_max = np.array([1, 1, 1, 1, 1, 1])
    P_high_min = np.array([ 90,  90,  90,  90,  90,  90])*1e5
    P_high_max = np.array([130, 180, 180, 180, 200, 200])*1e5
    
    for i in range(len(T_vec)):
        
        Optimizer.set_parameters(
            HP_ARCH = 'IHX_EXP', # IHX and IHX_EXP
            eta_cp = 0.8, # [-]
            eta_gc = 0.95, # [-]
            eta_IHX = 0.8, # [-]
            PP_ev = 5, # [K]
            SH_ev = 0.1, # [K]
            eta_exp = 0.7, # [-]
            
            P_high_min = P_high_min[i], # [Pa]
            P_high_max = P_high_max[i], # [Pa]
            
            m_dot_HS_fact_min = m_dot_HS_fact_min[i], # [-]
            m_dot_HS_fact_max = m_dot_HS_fact_max[i], # [-]
            
            m_dot_min = 5, # [kg/s]
            m_dot_max = 25 # [kg/s]
            )
        
        Optimizer.set_it_var(
            P_high = 140*1e5, # [Pa]
            mdot_HS = 20, # [kg/s]
            mdot = 17 # [kg/s]
            )
        
        Optimizer.set_obj(
            T_high = T_vec[i], # [K]
            W_dot = 1e6, # [W]
            )
        
        Optimizer.CSource.set_properties(
            T = 15 + 273.15,
            P = 1*1e5,
            fluid = 'Water'
            )
        
        Optimizer.HSource.set_properties(
            T = 15 + 273.15,
            P = 10*1e5,
            fluid = 'Water',
            )
    
        Optimizer.set_HP()
        Optimizer.opt_HP()

        COP_vec.append(Optimizer.COP)
        P_high_vec.append(Optimizer.it_var['P_high'])
        m_dot_vec.append(Optimizer.it_var['mdot'])
        m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
        T_act_vec.append(Optimizer.HP.components['GasCooler'].model.ex_C.T - 273.15)

    import matplotlib.pyplot as plt
    
    # Plot 1: Efficiency vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, COP_vec, linewidth=2)
    plt.title("Efficiency vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("COP")
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

