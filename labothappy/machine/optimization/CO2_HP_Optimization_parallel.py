# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

Adapted for multiprocessing-based PSO with optimizer return
"""

#%% Imports

from machine.examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import IHX_CO2_HP, IHX_EXP_CO2_HP
from connector.mass_connector import MassConnector

import numpy as np
import CoolProp.CoolProp as CP
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm
from multiprocessing import Pool
from joblib import Parallel, delayed

#%% Externalized cost function for multiprocessing

def system_HP_parallel(x, input_data):
    fluid = input_data['fluid']
    params = input_data['params']
    obj = input_data['obj']

    HSource = MassConnector()
    CSource = MassConnector()
    HSource.set_properties(**input_data['HSource'])
    CSource.set_properties(**input_data['CSource'])

    AS = CP.AbstractState('HEOS', fluid)
    AS.update(CP.QT_INPUTS, 0.5, CSource.T)
    P_sat_T_CSource = AS.p()
    P_crit_CO2 = AS.p_critical()
    P_low_guess = 0.8 * min(P_sat_T_CSource, P_crit_CO2)

    P_high, mdot, mdot_HS_fact = x
    mdot_HS = mdot * mdot_HS_fact
    HSource.set_properties(m_dot=mdot_HS)

    try:
        if params['HP_ARCH'] == 'IHX_EXP':
            HP = IHX_EXP_CO2_HP(HSource, CSource.T, params['eta_cp'], params['eta_gc'],
                                params['eta_IHX'], params['eta_exp'],
                                params['PP_ev'], params['SH_ev'],
                                P_low_guess, P_high, mdot, mute_print_flag=1)
        else:
            HP = IHX_CO2_HP(HSource, CSource.T, params['eta_cp'], params['eta_gc'],
                            params['eta_IHX'], params['PP_ev'], params['SH_ev'],
                            P_low_guess, P_high, mdot, mute_print_flag=1)

        HP.solve()

        if not HP.converged:
            return 100

        if params['HP_ARCH'] == 'IHX_EXP':
            COP = HP.components['GasCooler'].model.Q_dot.Q_dot / (
                HP.components['Compressor'].model.W.W_dot - HP.components['Expander'].model.W_exp.W_dot)
        else:
            COP = HP.components['GasCooler'].model.Q_dot.Q_dot / HP.components['Compressor'].model.W.W_dot

        T_err = abs((HP.components['GasCooler'].model.ex_C.T - obj['T_high']) / obj['T_high'])
        W_err = abs((HP.components['Compressor'].model.W.W_dot - obj['W_dot']) / obj['W_dot'])

        penalty = 0
        if T_err > 1e-2:
            penalty += 300 * T_err
        if W_err > 1e-2:
            penalty += 100 * W_err

        return -COP + penalty

    except Exception:
        return 100

#%% Optimizer Class

class CO2HPOptimizer:
    def __init__(self, fluid):
        self.fluid = fluid
        self.params = {}
        self.it_var = {}
        self.obj = {}
        self.HSource = MassConnector()
        self.CSource = MassConnector()
        self.COP = None
        self.optimizer = None

    def set_parameters(self, **kwargs):
        self.params.update(kwargs)

    def set_it_var(self, **kwargs):
        self.it_var.update(kwargs)

    def set_obj(self, **kwargs):
        self.obj.update(kwargs)

    def opt_HP(self):
        import multiprocessing
        n_cores = multiprocessing.cpu_count()
        
        bounds = (
            np.array([self.params['P_high_min'], self.params['m_dot_min'], self.params['m_dot_HS_fact_min']]),
            np.array([self.params['P_high_max'], self.params['m_dot_max'], self.params['m_dot_HS_fact_max']])
        )
    
        input_data = {
            'fluid': self.fluid,
            'params': self.params,
            'obj': self.obj,
            'HSource': {
                'T': self.HSource.T,
                'P': self.HSource.p,
                'fluid': self.HSource.fluid,
                'm_dot': self.it_var['mdot_HS']
            },
            'CSource': {
                'T': self.CSource.T,
                'P': self.CSource.p,
                'fluid': self.CSource.fluid
            }
        }
    
        # ✅ Open pool once and keep it open during the loop
        def objective_wrapper(X):
            return np.array(Parallel(n_jobs=n_cores - 1, backend='loky', prefer="processes")(
                delayed(system_HP_parallel)(x, input_data) for x in X
            ))
            
        self.optimizer = GlobalBestPSO(
            n_particles=40,
            dimensions=3,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )

        best_cost = np.inf
        no_improve_counter = 0
        patience = 5
        tol = 1e-3
        max_iter = 30

        for i in tqdm(range(max_iter), desc="Optimizing", ncols=80):
            self.optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = self.optimizer.swarm.best_cost

            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1

            print(f"[{i+1:03}] Best cost: {best_cost:.6f}")
            if no_improve_counter >= patience:
                print("Stopping early due to stagnation.")
                break
    
        # After optimization loop is done
        best_P_high = self.optimizer.swarm.best_pos[0]
        best_mdot = self.optimizer.swarm.best_pos[1]
        best_mdot_HS_fact = self.optimizer.swarm.best_pos[2]
        best_mdot_HS = best_mdot * best_mdot_HS_fact
        
        # Store best variables
        self.it_var['P_high'] = best_P_high
        self.it_var['mdot'] = best_mdot
        self.it_var['mdot_HS'] = best_mdot_HS
        
        # 🔁 Rebuild the HP circuit in the main process
        self.HSource.set_properties(m_dot=best_mdot_HS)
        AS = CP.AbstractState('HEOS', self.fluid)
        AS.update(CP.QT_INPUTS, 0.5, self.CSource.T)
        P_sat_T_CSource = AS.p()
        P_crit_CO2 = AS.p_critical()
        P_low_guess = 0.8 * min(P_sat_T_CSource, P_crit_CO2)
        
        if self.params['HP_ARCH'] == 'IHX_EXP':
            self.HP = IHX_EXP_CO2_HP(self.HSource, self.CSource.T, self.params['eta_cp'],
                                      self.params['eta_gc'], self.params['eta_IHX'], self.params['eta_exp'],
                                      self.params['PP_ev'], self.params['SH_ev'],
                                      P_low_guess, best_P_high, best_mdot, mute_print_flag=1)
        else:
            self.HP = IHX_CO2_HP(self.HSource, self.CSource.T, self.params['eta_cp'],
                                 self.params['eta_gc'], self.params['eta_IHX'],
                                 self.params['PP_ev'], self.params['SH_ev'],
                                 P_low_guess, best_P_high, best_mdot, mute_print_flag=1)
        
        # 🔁 Solve in main process and store COP
        try:
            self.HP.solve()
            if self.params['HP_ARCH'] == 'IHX_EXP':
                self.COP = self.HP.components['GasCooler'].model.Q_dot.Q_dot / (
                    self.HP.components['Compressor'].model.W.W_dot - self.HP.components['Expander'].model.W_exp.W_dot)
            else:
                self.COP = self.HP.components['GasCooler'].model.Q_dot.Q_dot / \
                           self.HP.components['Compressor'].model.W.W_dot
        except Exception as e:
            print(f"⚠️ Failed to solve final HP circuit: {e}")
            self.COP = None
            self.HP = None
    
        return self.optimizer

#%% Optimizer call

if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    # Initialize optimizer
    Optimizer = CO2HPOptimizer('CO2')
    
    # Define temperatures for sweep
    T_vec = np.linspace(100, 150, 6) + 273.15
    COP_vec = []
    P_high_vec = []
    m_dot_vec = []
    m_dot_HS_vec = []
    T_act_vec = []
    
    # Sweep parameters
    m_dot_HS_fact_min = np.array([0.2] * 6)
    m_dot_HS_fact_max = np.array([0.8] * 6)
    P_high_min = np.array([100, 100, 120, 130, 140, 150]) * 1e5
    P_high_max = np.array([130, 150, 170, 180, 190, 200]) * 1e5
    
    for i in range(len(T_vec)):
        Optimizer.set_parameters(
            HP_ARCH='IHX_EXP',
            eta_cp=0.8,
            eta_gc=0.95,
            eta_IHX=0.8,
            eta_exp=0.7,
            PP_ev=5,
            SH_ev=0.1,
            P_high_min=P_high_min[i],
            P_high_max=P_high_max[i],
            m_dot_min=10,
            m_dot_max=20,
            m_dot_HS_fact_min=m_dot_HS_fact_min[i],
            m_dot_HS_fact_max=m_dot_HS_fact_max[i]
        )
    
        Optimizer.set_it_var(
            P_high=140e5,
            mdot=17,
            mdot_HS=20
        )
    
        Optimizer.set_obj(
            T_high=T_vec[i],
            W_dot=1e6
        )
    
        Optimizer.CSource.set_properties(
            T=15 + 273.15,
            P=1e5,
            fluid='Water'
        )
    
        Optimizer.HSource.set_properties(
            T=15 + 273.15,
            P=10e5,
            fluid='Water'
        )
    
        # Run optimization and get optimizer instance
        opt = Optimizer.opt_HP()
    
        # Store results
        COP_vec.append(Optimizer.COP)
        P_high_vec.append(opt.swarm.best_pos[0])
        m_dot_vec.append(opt.swarm.best_pos[1])
        m_dot_HS_vec.append(opt.swarm.best_pos[1] * opt.swarm.best_pos[2])
        T_act_vec.append(T_vec[i] - 273.15)
    
    # Plotting
    plt.figure(figsize=(8, 5))
    plt.plot(T_act_vec, COP_vec, linewidth=2)
    plt.title("Efficiency vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("COP")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8, 5))
    plt.plot(T_act_vec, P_high_vec, linewidth=2, label='P_high')
    plt.title("High Pressure vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("P_high (Pa)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8, 5))
    plt.plot(T_act_vec, m_dot_vec, linewidth=2, label='mdot')
    plt.title("Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Mass Flow Rate (kg/s)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8, 5))
    plt.plot(T_act_vec, m_dot_HS_vec, linewidth=2, label='mdot_HS')
    plt.title("Heat Source Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Heat Source Mass Flow Rate (kg/s)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
# COP_vec = [4.697249244448596, 4.3895555255356085, 4.050117334016632, 
#            4.026683891241342, 3.7880968812034888, 3.5983675455471418]

# P_high_vec = [12331249.06403394, 13392550.320759114, 16319829.434295718, 
#               17739417.348679177, 17709543.51810449, 19081338.907870803]

# mdot_vec = [17.30013398744241, 15.670235783606566, 13.076689048124113, 
#                12.787677902485882, 11.890533504136116, 11.045665619017555]

# mdot_HS_vec = [12.392235492258568, 10.226879105272955, 7.796114170455033, 
#                7.481034247657304, 6.4890594730472175, 5.70901076559547]
