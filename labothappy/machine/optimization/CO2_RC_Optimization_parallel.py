# # -*- coding: utf-8 -*-
# """
# Created on Wed Jul 16 11:12:02 2025

# @author: Basile
# """

# #%% Imports

# from machine.examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import IHX_CO2_HP, IHX_EXP_CO2_HP
# from machine.examples.CO2_Transcritical_Circuits.CO2_Transcritical_circuit import REC_CO2_TC, basic_CO2_TC
# from connector.mass_connector import MassConnector

# import numpy as np
# from CoolProp.CoolProp import PropsSI
# from pyswarms.single import GlobalBestPSO
# from tqdm import tqdm
# from joblib import Parallel, delayed

# import warnings
# warnings.filterwarnings('ignore')

# #%% Define parallel system evaluation outside the class

# def system_RC_parallel(x, input_data):
#     warnings.filterwarnings('ignore')

#     fluid = input_data['fluid']
#     params = input_data['params']
#     obj = input_data['obj']
#     hs_props = input_data['HSource']
#     cs_props = input_data['CSource']

#     HSource = MassConnector()
#     CSource = MassConnector()

#     HSource.set_properties(
#         T=hs_props['T'], P=hs_props['P'], fluid=hs_props['fluid'], m_dot=x[1]*x[2])
#     CSource.set_properties(
#         T=cs_props['T'], P=cs_props['P'], fluid=cs_props['fluid'])

#     P_sat_T_CSource = PropsSI('P', 'T', CSource.T, 'Q', 0.5, fluid)
#     P_crit_CO2 = PropsSI('PCRIT', fluid)
#     P_low_guess = min(1.3*P_sat_T_CSource, 0.8*P_crit_CO2)

    
#     if input_data['RC_ARCH'] == 'REC':
#         RC = REC_CO2_TC(HSource, CSource.T, params['PP_gh'], params['PP_rec'], params['eta_pp'],
#                         params['eta_exp'], params['eta_gh'], params['eta_rec'], params['PP_cd'], params['SC_cd'], 
#                         P_low_guess, x[0], x[1], DP_h_rec = params['DP_h_rec'], DP_c_rec = params['DP_c_rec'], 
#                         DP_h_gh = params['DP_h_gh'], DP_c_gh = params['DP_c_gh'], DP_cond = params['DP_cond'], mute_print_flag=1)

#     elif input_data['RC_ARCH'] == 'basic':
#         RC = basic_CO2_TC(HSource, CSource.T, params['PP_gh'], params['PP_rec'], params['eta_pp'],
#                         params['eta_exp'], params['eta_gh'], params['PP_cd'], params['SC_cd'], P_low_guess, 
#                         x[0], x[1], mute_print_flag=1)

#     try:
#         RC.solve()

#         DP = 50e3
#         rho = RC.components['GasHeater'].model.su_H.D
#         mdot = RC.components['GasHeater'].model.su_H.m_dot
#         h_ex = RC.components['GasHeater'].model.ex_H.h

#         T_amb = 15 + 273.15
#         h_ex_req = PropsSI('H', 'T', T_amb, 'P', 101325, 'Water')
#         Q_dot_req = mdot * (h_ex - h_ex_req)
#         W_dot_fan = 0.2 * Q_dot_req

#         eta_pp = 0.8
#         pp_power = DP * mdot / (rho * eta_pp)

#         W_dot_net = RC.components['Expander'].model.W_exp.W_dot*0.95 - RC.components['Pump'].model.W_pp.W_dot/0.95 - pp_power/0.95
        
#         eta = W_dot_net / RC.components['GasHeater'].model.Q_dot.Q_dot

#         Th_out = RC.components['GasHeater'].model.ex_H.T
#         Th_out_obj = 15 + 273.15

#         penalty = 0

#         if abs((W_dot_net - obj['W_dot'])/obj['W_dot']) > 1e-2:
#             penalty = abs((W_dot_net - obj['W_dot'])/obj['W_dot']) * 100

#         # if abs((Th_out-Th_out_obj)/Th_out_obj) > 1e-2:
#         #     penalty = abs((Th_out-Th_out_obj)/Th_out_obj) * 1

#         return -eta + penalty

#     except:
#         return 1000

# #%% Optimizer Class

# class CO2RCOptimizer:

#     def __init__(self, fluid):
#         self.fluid = fluid
#         self.inputs = {}
#         self.params = {}
#         self.it_var = {}
#         self.obj = {}
#         self.HSource = MassConnector()
#         self.CSource = MassConnector()

#     def set_inputs(self, **parameters):
#         self.inputs.update(parameters)

#     def set_parameters(self, **parameters):
#         self.params.update(parameters)

#     def set_it_var(self, **parameters):
#         self.it_var.update(parameters)

#     def set_obj(self, **parameters):
#         self.obj.update(parameters)

#     def opt_RC(self):
#         import multiprocessing
#         n_cores = multiprocessing.cpu_count()
        
#         bounds = (
#             np.array([self.params['P_high_bounds'][0], self.params['m_dot_bounds'][0], self.params['m_dot_HS_fact_bounds'][0]]),
#             np.array([self.params['P_high_bounds'][1], self.params['m_dot_bounds'][1], self.params['m_dot_HS_fact_bounds'][1]])
#         )
    
#         input_data = {
#             'fluid': self.fluid,
#             'params': self.params,
#             'obj': self.obj,
#             'HSource': {
#                 'T': self.HSource.T,
#                 'P': self.HSource.p,
#                 'fluid': self.HSource.fluid
#             },
#             'CSource': {
#                 'T': self.CSource.T,
#                 'P': self.CSource.p,
#                 'fluid': self.CSource.fluid
#             },
#             'RC_ARCH' : self.params['RC_ARCH']
#         }
    
#         def objective_wrapper(X):
#             return np.array(Parallel(n_jobs=n_cores - 1, backend='loky')(
#                 delayed(system_RC_parallel)(x, input_data) for x in X
#             ))
    
#         self.optimizer = GlobalBestPSO(
#             n_particles=20,
#             dimensions=3,
#             options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
#             bounds=bounds
#         )
    
#         best_cost = np.inf
#         no_improve_counter = 0
#         patience = 5
#         tol = 1e-3
#         max_iter = 10
    
#         pbar = tqdm(total=max_iter, desc="Optimizing", ncols=80)
        
#         for i in range(max_iter):
#             self.optimizer.optimize(objective_wrapper, iters=1, verbose=False)
#             current_best = self.optimizer.swarm.best_cost
        
#             if current_best < best_cost - tol:
#                 best_cost = current_best
#                 no_improve_counter = 0
#             else:
#                 no_improve_counter += 1
        
#             pbar.set_postfix(best_cost=f"{best_cost:.6f}")
#             pbar.update(1)
        
#             if no_improve_counter >= patience:
#                 pbar.set_description("Stopped (no improvement)")
#                 break
        
#         pbar.close()
    
#         # Recompute system with best parameters
#         best_P_high, best_m_dot, best_m_dot_HS_fact = self.optimizer.swarm.best_pos
    
#         self.it_var['P_high'] = best_P_high
#         self.it_var['mdot'] = best_m_dot
#         self.it_var['mdot_HS'] = best_m_dot_HS = best_m_dot * best_m_dot_HS_fact
    
#         self.HSource.set_properties(m_dot=best_m_dot_HS)
    
#         # Estimate low pressure for initialization
#         P_sat_T_CSource = PropsSI('P', 'T', self.CSource.T, 'Q', 0.5, self.fluid)
#         P_crit_CO2 = PropsSI('PCRIT', self.fluid)
#         P_low_guess = min(1.3 * P_sat_T_CSource, 0.8 * P_crit_CO2)
    
#         try:
#             self.RC = REC_CO2_TC(
#                 self.HSource, self.CSource.T,
#                 self.params['PP_gh'], self.params['PP_rec'],
#                 self.params['eta_pp'], self.params['eta_exp'],
#                 self.params['eta_gh'], self.params['eta_rec'],
#                 self.params['PP_cd'], self.params['SC_cd'],
#                 P_low_guess, best_P_high, best_m_dot,
#                 mute_print_flag=1
#             )
    
#             self.RC.solve()
    
#             # Fan and pump calculations
#             DP = 50e3
#             rho = self.RC.components['GasHeater'].model.su_H.D
#             mdot = self.RC.components['GasHeater'].model.su_H.m_dot
#             h_ex = self.RC.components['GasHeater'].model.ex_H.h
    
#             h_ex_req = PropsSI('H', 'T', 15 + 273.15, 'P', 101325, 'Water')
#             Q_dot_req = mdot * (h_ex - h_ex_req)
#             W_dot_fan = 0.2 * Q_dot_req
#             W_dot_pp = DP * mdot / (rho * 0.8)
    
#             W_net = self.RC.components['Expander'].model.W_exp.W_dot*0.95 - \
#                     self.RC.components['Pump'].model.W_pp.W_dot/0.95 - W_dot_pp/0.95
    
#             eta_final = W_net / self.RC.components['GasHeater'].model.Q_dot.Q_dot
    
#             self.eta = eta_final
#             self.W_dot_net = W_net
#             self.pp_power = W_dot_pp
#             self.W_dot_fan = W_dot_fan
#             self.final_RC = self.RC

#             self.Q_dot_waste = self.RC.components['GasHeater'].model.ex_H.m_dot*(self.RC.components['GasHeater'].model.ex_H.h- PropsSI('H', 'T', 273.15 +15, 'P', self.RC.components['GasHeater'].model.ex_H.p, self.RC.components['GasHeater'].model.ex_H.fluid))
    
#             print(f"\nOptimal P_high: {best_P_high:.2f} Pa")
#             print(f"Optimal m_dot: {best_m_dot:.4f} kg/s")
#             print(f"Optimal m_dot_HS: {best_m_dot_HS:.4f} kg/s")
#             print(f"Final net W: {W_net:.2f} W")
#             print(f"Final η (thermal efficiency): {eta_final:.4f}")
#             print(f"Best score: {-best_cost:.4f}")
    
#         except Exception as e:
#             print(f"⚠️ Failed to solve final RC circuit: {e}")
#             self.RC = None
#             self.eta = None

#         return self.optimizer

# #%% Optimizer call

# if __name__ == "__main__":
    
#     import matplotlib.pyplot as plt
    
#     # Define temperature sweep (°C to K)
#     # T_vec = np.linspace(100, 150, 6) + 273.15
#     T_vec = np.array([130])+273.15
    
#     # Output vectors
#     eta_vec = []
#     P_high_vec = []
#     m_dot_vec = []
#     m_dot_HS_vec = []
#     T_h_ex_vec = []
#     Q_dot_waste = []
        
#     n_MW = 1 # W
#     W_dot_test = n_MW*1e6 # W
    
#     # Create optimizer instance
#     Optimizer = CO2RCOptimizer('CO2')
    
#     # Sweep parameters
#     m_dot_HS_fact_bounds = [0.5,1]
#     P_high_bounds = np.array([80, 180]) * 1e5
#     m_dot_bounds = np.array([30,80])*n_MW
    
#     # Sweep loop
#     for i in range(len(T_vec)):
#         # Set model parameters
#         Optimizer.set_parameters(
#             RC_ARCH= 'REC', # 'REC'
            
#             # Pump
#             eta_pp=0.8,
            
#             # GasHeater
#             eta_gh=0.95,
#             PP_gh=5,
#             DP_h_gh = 50*1e3,
#             DP_c_gh = 50*1e3,

#             # DP_h_gh = 50*1e3,
#             # DP_c_gh = 2*1e5,
    
#             # Recuperator
#             eta_rec=0.8,
#             PP_rec=0,
#             DP_h_rec = 50*1e3,
#             DP_c_rec = 50*1e3,

#             # DP_h_rec = 1*1e5,
#             # DP_c_rec = 2*1e5,
            
#             # Expander
#             eta_exp=0.9,
            
#             # Condenser
#             PP_cd=5,
#             SC_cd=0.1,
#             DP_cond = 50*1e3,

#             # DP_cond = 1*1e5,
            
#             # Bounds
#             P_high_bounds=P_high_bounds,
#             m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
#             m_dot_bounds = m_dot_bounds
#         )
    
#         # Initial guess
#         Optimizer.set_it_var(
#             P_high=100e5,
#             mdot=27,
#             mdot_HS=20
#         )
    
#         # Objective
#         Optimizer.set_obj(W_dot=1e6)
    
#         # Source definitions
#         Optimizer.CSource.set_properties(
#             T=15 + 273.15,
#             P=5e5,
#             fluid='Water'
#         )
    
#         Optimizer.HSource.set_properties(
#             T=T_vec[i],
#             P=10e5,
#             fluid='Water'
#         )
    
#         # Prepare model and optimize
#         Optimizer.opt_RC()
    
#         # Collect results
#         eta_vec.append(Optimizer.eta)
#         P_high_vec.append(Optimizer.it_var['P_high'])
#         m_dot_vec.append(Optimizer.it_var['mdot'])
#         m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
#         T_h_ex_vec.append(Optimizer.RC.components['GasHeater'].model.ex_H.T)
#         Q_dot_waste.append(Optimizer.Q_dot_waste)

#     #%% Plotting Results
    
#     T_C = T_vec - 273.15  # convert to °C
    
#     # Plot 1: Efficiency
#     plt.figure(figsize=(8, 5))
#     plt.plot(T_C, eta_vec, linewidth=2)
#     plt.title("Efficiency vs Temperature")
#     plt.xlabel("Temperature (°C)")
#     plt.ylabel("Efficiency")
#     plt.grid(True)
#     plt.tight_layout()
#     plt.show()
    
#     # Plot 2: P_high
#     plt.figure(figsize=(8, 5))
#     plt.plot(T_C, P_high_vec, linewidth=2, label='P_high')
#     plt.title("High Pressure vs Temperature")
#     plt.xlabel("Temperature (°C)")
#     plt.ylabel("P_high (Pa)")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
    
#     # Plot 3: mdot
#     plt.figure(figsize=(8, 5))
#     plt.plot(T_C, m_dot_vec, linewidth=2, label='mdot')
#     plt.title("Mass Flow Rate vs Temperature")
#     plt.xlabel("Temperature (°C)")
#     plt.ylabel("Mass Flow Rate (kg/s)")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
    
#     # Plot 4: mdot_HS
#     plt.figure(figsize=(8, 5))
#     plt.plot(T_C, m_dot_HS_vec, linewidth=2, label='mdot_HS')
#     plt.title("Heat Source Mass Flow Rate vs Temperature")
#     plt.xlabel("Temperature (°C)")
#     plt.ylabel("Heat Source Mass Flow Rate (kg/s)")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

# # eta_vec = [0.09086802588180352, 0.10128587354579623, 0.11146548353387467, 
# #            0.1214949488715974, 0.13104562249758014, 0.140919652539438]

# # P_high_vec = [12240858.880586991, 12984383.432874506, 13896082.028326932, 
# #               14910851.41300104, 14942289.443446713, 16668855.969638687]

# # mdot_vec = [54.192112753242355, 47.38773716711351, 42.00913879832734,
# #             37.49690624913382, 33.69443605681327, 30.290940206479696]

# # mdot_HS_vec = [43.24829905355492, 33.858992654942654, 35.71657513820893, 
# #                29.40146975511797, 32.027319612867196, 27.070350435856742]

# # 4266224.080146256, 4892408.350373052, 4262821.838781694, 3254882.4953648252, 2746835.699992969, 5196392.521202499
# # 2518221.691336776, 1963451.9608898002, 1891513.8494420673, 2062909.7198840023, 2717955.854383429, 1760254.759218207

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

@author: Basile
"""

#%% Imports

from machine.examples.CO2_Heat_Pumps.CO2_HeatPump_circuit import IHX_CO2_HP, IHX_EXP_CO2_HP
from machine.examples.CO2_Transcritical_Circuits.CO2_Transcritical_circuit import REC_CO2_TC, basic_CO2_TC
from connector.mass_connector import MassConnector

import numpy as np
from CoolProp.CoolProp import PropsSI
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm
from joblib import Parallel, delayed

import warnings
warnings.filterwarnings('ignore')

#%% Define parallel system evaluation outside the class

def system_RC_parallel(x, input_data):
    warnings.filterwarnings('ignore')

    fluid = input_data['fluid']
    params = input_data['params']
    obj = input_data['obj']
    hs_props = input_data['HSource']
    cs_props = input_data['CSource']

    HSource = MassConnector()
    CSource = MassConnector()

    HSource.set_properties(
        T=hs_props['T'], P=hs_props['P'], fluid=hs_props['fluid'], m_dot=x[1]*x[2])
    CSource.set_properties(
        T=cs_props['T'], P=cs_props['P'], fluid=cs_props['fluid'])

    P_sat_T_CSource = PropsSI('P', 'T', CSource.T, 'Q', 0.5, fluid)
    P_crit_CO2 = PropsSI('PCRIT', fluid)
    P_low_guess = min(1.3*P_sat_T_CSource, 0.8*P_crit_CO2)

    
    if input_data['RC_ARCH'] == 'REC':
        RC = REC_CO2_TC(HSource, CSource, params['PP_gh'], params['PP_rec'], params['eta_pp'],
                        params['eta_exp'], params['eta_gh'], params['eta_rec'], params['PP_cd'], params['SC_cd'], 
                        P_low_guess, x[0], x[1], DP_h_rec = params['DP_h_rec'], DP_c_rec = params['DP_c_rec'], 
                        DP_h_gh = params['DP_h_gh'], DP_c_gh = params['DP_c_gh'], DP_h_cond = params['DP_cond'], mute_print_flag=1)

    elif input_data['RC_ARCH'] == 'basic':
        RC = basic_CO2_TC(HSource, CSource, params['PP_gh'], params['PP_rec'], params['eta_pp'],
                        params['eta_exp'], params['eta_gh'], params['PP_cd'], params['SC_cd'], P_low_guess, 
                        x[0], x[1], mute_print_flag=1)

    try:
        RC.solve()

        DP = 50e3
        rho = RC.components['GasHeater'].model.su_H.D
        mdot = RC.components['GasHeater'].model.su_H.m_dot
        h_ex = RC.components['GasHeater'].model.ex_H.h

        T_amb = 15 + 273.15
        h_ex_req = PropsSI('H', 'T', T_amb, 'P', 101325, 'Water')
        Q_dot_req = mdot * (h_ex - h_ex_req)
        W_dot_fan = 0.2 * Q_dot_req

        eta_pp = 0.8
        pp_power = DP * mdot / (rho * eta_pp)

        W_dot_net = RC.components['Expander'].model.W.W_dot*0.95 - RC.components['Pump'].model.W.W_dot/0.95 - pp_power/0.95
        
        eta = W_dot_net / RC.components['GasHeater'].model.Q_dot.Q_dot

        Th_out = RC.components['GasHeater'].model.ex_H.T
        Th_out_obj = 15 + 273.15

        penalty = 0

        if abs((W_dot_net - obj['W_dot'])/obj['W_dot']) > 1e-2:
            penalty = abs((W_dot_net - obj['W_dot'])/obj['W_dot']) * 100

        # if abs((Th_out-Th_out_obj)/Th_out_obj) > 1e-2:
        #     penalty = abs((Th_out-Th_out_obj)/Th_out_obj) * 1

        return -eta + penalty

    except:
        return 1000

#%% Optimizer Class

class CO2RCOptimizer:

    def __init__(self, fluid):
        self.fluid = fluid
        self.inputs = {}
        self.params = {}
        self.it_var = {}
        self.obj = {}
        self.HSource = MassConnector()
        self.CSource = MassConnector()

    def set_inputs(self, **parameters):
        self.inputs.update(parameters)

    def set_parameters(self, **parameters):
        self.params.update(parameters)

    def set_it_var(self, **parameters):
        self.it_var.update(parameters)

    def set_obj(self, **parameters):
        self.obj.update(parameters)

    def opt_RC(self):
        import multiprocessing
        n_cores = multiprocessing.cpu_count()
        
        bounds = (
            np.array([self.params['P_high_bounds'][0], self.params['m_dot_bounds'][0], self.params['m_dot_HS_fact_bounds'][0]]),
            np.array([self.params['P_high_bounds'][1], self.params['m_dot_bounds'][1], self.params['m_dot_HS_fact_bounds'][1]])
        )
    
        input_data = {
            'fluid': self.fluid,
            'params': self.params,
            'obj': self.obj,
            'HSource': {
                'T': self.HSource.T,
                'P': self.HSource.p,
                'fluid': self.HSource.fluid
            },
            'CSource': {
                'T': self.CSource.T,
                'P': self.CSource.p,
                'fluid': self.CSource.fluid
            },
            'RC_ARCH' : self.params['RC_ARCH']
        }
    
        def objective_wrapper(X):
            return np.array(Parallel(n_jobs=n_cores - 1, backend='loky')(
                delayed(system_RC_parallel)(x, input_data) for x in X
            ))
    
        self.optimizer = GlobalBestPSO(
            n_particles=20,
            dimensions=3,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        best_cost = np.inf
        no_improve_counter = 0
        patience = 5
        tol = 1e-3
        max_iter = 10
    
        pbar = tqdm(total=max_iter, desc="Optimizing", ncols=80)
        
        for i in range(max_iter):
            self.optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = self.optimizer.swarm.best_cost
        
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
        
            pbar.set_postfix(best_cost=f"{best_cost:.6f}")
            pbar.update(1)
        
            if no_improve_counter >= patience:
                pbar.set_description("Stopped (no improvement)")
                break
        
        pbar.close()
    
        # Recompute system with best parameters
        best_P_high, best_m_dot, best_m_dot_HS_fact = self.optimizer.swarm.best_pos
    
        self.it_var['P_high'] = best_P_high
        self.it_var['mdot'] = best_m_dot
        self.it_var['mdot_HS'] = best_m_dot_HS = best_m_dot * best_m_dot_HS_fact
    
        self.HSource.set_properties(m_dot=best_m_dot_HS)
    
        # Estimate low pressure for initialization
        P_sat_T_CSource = PropsSI('P', 'T', self.CSource.T, 'Q', 0.5, self.fluid)
        P_crit_CO2 = PropsSI('PCRIT', self.fluid)
        P_low_guess = min(1.3 * P_sat_T_CSource, 0.8 * P_crit_CO2)
    
    
        try:
            self.RC = REC_CO2_TC(
                self.HSource, self.CSource,
                self.params['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'],
                self.params['eta_gh'], self.params['eta_rec'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, best_P_high, best_m_dot,
                mute_print_flag=1
            )
    
            self.RC.solve()
    
            # Fan and pump calculations
            DP = 50e3
            rho = self.RC.components['GasHeater'].model.su_H.D
            mdot = self.RC.components['GasHeater'].model.su_H.m_dot
            h_ex = self.RC.components['GasHeater'].model.ex_H.h
    
            h_ex_req = PropsSI('H', 'T', 15 + 273.15, 'P', 101325, 'Water')
            Q_dot_req = mdot * (h_ex - h_ex_req)
            W_dot_fan = 0.2 * Q_dot_req
            W_dot_pp = DP * mdot / (rho * 0.8)
    
            W_net = self.RC.components['Expander'].model.W.W_dot*0.95 - \
                    self.RC.components['Pump'].model.W.W_dot/0.95 - W_dot_pp/0.95
    
            eta_final = W_net / self.RC.components['GasHeater'].model.Q_dot.Q_dot
    
            self.eta = eta_final
            self.W_dot_net = W_net
            self.pp_power = W_dot_pp
            self.W_dot_fan = W_dot_fan
            self.final_RC = self.RC
    
            self.Q_dot_waste = self.RC.components['GasHeater'].model.ex_H.m_dot*(self.RC.components['GasHeater'].model.ex_H.h- PropsSI('H', 'T', 273.15 +15, 'P', self.RC.components['GasHeater'].model.ex_H.p, self.RC.components['GasHeater'].model.ex_H.fluid))
    
            print(f"\nOptimal P_high: {best_P_high:.2f} Pa")
            print(f"Optimal m_dot: {best_m_dot:.4f} kg/s")
            print(f"Optimal m_dot_HS: {best_m_dot_HS:.4f} kg/s")
            print(f"Final net W: {W_net:.2f} W")
            print(f"Final η (thermal efficiency): {eta_final:.4f}")
            print(f"Best score: {-best_cost:.4f}")
        
        except Exception as e:
            print(f"⚠️ Failed to solve final RC circuit: {e}")
            self.RC = None
            self.eta = None

        return self.optimizer

#%% Optimizer call

if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    # Define temperature sweep (°C to K)
    # T_vec = np.linspace(100, 150, 6) + 273.15
    T_vec = np.array([130])+273.15
    
    # Output vectors
    eta_vec = []
    P_high_vec = []
    m_dot_vec = []
    m_dot_HS_vec = []
    T_h_ex_vec = []
    Q_dot_waste = []
        
    n_MW = 1 # W
    W_dot_test = n_MW*1e6 # W
    
    # Create optimizer instance
    Optimizer = CO2RCOptimizer('CO2')
    
    # Sweep parameters
    m_dot_HS_fact_bounds = [0.5,1]
    P_high_bounds = np.array([80, 180]) * 1e5
    m_dot_bounds = np.array([30,80])*n_MW
    
    
    # Sweep loop
    for i in range(len(T_vec)):
        # Set model parameters
        Optimizer.set_parameters(
            RC_ARCH= 'REC', # 'REC'
            
            # Pump
            eta_pp=0.8,
            
            # GasHeater
            eta_gh=0.95,
            PP_gh=5,
            DP_h_gh = 50*1e3,
            DP_c_gh = 50*1e3,

            # DP_h_gh = 50*1e3,
            # DP_c_gh = 2*1e5,
    
            # Recuperator
            eta_rec=0.8,
            PP_rec=0,
            DP_h_rec = 50*1e3,
            DP_c_rec = 50*1e3,

            # DP_h_rec = 1*1e5,
            # DP_c_rec = 2*1e5,
            
            # Expander
            eta_exp=0.9,
            
            # Condenser
            PP_cd=5,
            SC_cd=0.1,
            DP_cond = 50*1e3,

            # DP_cond = 1*1e5,
            
            # Bounds
            P_high_bounds=P_high_bounds,
            m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
            m_dot_bounds = m_dot_bounds
        )
    
        # Initial guess
        Optimizer.set_it_var(
            P_high=100e5,
            mdot=27,
            mdot_HS=20
        )
    
        # Objective
        Optimizer.set_obj(W_dot=1e6)
    
        # Source definitions
        Optimizer.CSource.set_properties(
            T=15 + 273.15,
            P=5e5,
            fluid='Water'
        )
    
        Optimizer.HSource.set_properties(
            T=T_vec[i],
            P=10e5,
            fluid='Water'
        )
    
        # Prepare model and optimize
        Optimizer.opt_RC()
    
        # Collect results
        eta_vec.append(Optimizer.eta)
        P_high_vec.append(Optimizer.it_var['P_high'])
        m_dot_vec.append(Optimizer.it_var['mdot'])
        m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
        T_h_ex_vec.append(Optimizer.RC.components['GasHeater'].model.ex_H.T)
        Q_dot_waste.append(Optimizer.Q_dot_waste)

    #%% Plotting Results
    
    T_C = T_vec - 273.15  # convert to °C
    
    # Plot 1: Efficiency
    plt.figure(figsize=(8, 5))
    plt.plot(T_C, eta_vec, linewidth=2)
    plt.title("Efficiency vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Efficiency")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Plot 2: P_high
    plt.figure(figsize=(8, 5))
    plt.plot(T_C, P_high_vec, linewidth=2, label='P_high')
    plt.title("High Pressure vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("P_high (Pa)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot 3: mdot
    plt.figure(figsize=(8, 5))
    plt.plot(T_C, m_dot_vec, linewidth=2, label='mdot')
    plt.title("Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Mass Flow Rate (kg/s)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot 4: mdot_HS
    plt.figure(figsize=(8, 5))
    plt.plot(T_C, m_dot_HS_vec, linewidth=2, label='mdot_HS')
    plt.title("Heat Source Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Heat Source Mass Flow Rate (kg/s)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# eta_vec = [0.09086802588180352, 0.10128587354579623, 0.11146548353387467, 
#            0.1214949488715974, 0.13104562249758014, 0.140919652539438]

# P_high_vec = [12240858.880586991, 12984383.432874506, 13896082.028326932, 
#               14910851.41300104, 14942289.443446713, 16668855.969638687]

# mdot_vec = [54.192112753242355, 47.38773716711351, 42.00913879832734,
#             37.49690624913382, 33.69443605681327, 30.290940206479696]

# mdot_HS_vec = [43.24829905355492, 33.858992654942654, 35.71657513820893, 
#                29.40146975511797, 32.027319612867196, 27.070350435856742]

# 4266224.080146256, 4892408.350373052, 4262821.838781694, 3254882.4953648252, 2746835.699992969, 5196392.521202499
# 2518221.691336776, 1963451.9608898002, 1891513.8494420673, 2062909.7198840023, 2717955.854383429, 1760254.759218207