# # -*- coding: utf-8 -*-
# """
# Created on Wed Jul 16 11:12:02 2025

# Adapted for multiprocessing-based PSO with optimizer return
# + optimisation du débit de fluide froid (CSource)
# """

# #%% Imports

# from labothappy.machine.examples.Heat_Pumps.fpi_TC_HeatPump_example import IHX_CO2_HP, IHX_EXP_CO2_HP
# from labothappy.connector.mass_connector import MassConnector

# import numpy as np
# import CoolProp.CoolProp as CP
# from pyswarms.single import GlobalBestPSO
# from tqdm import tqdm
# from joblib import Parallel, delayed

# import warnings
# warnings.filterwarnings("ignore")

# #%% Helper : décodage du vecteur de variables

# def unpack_x(x, hp_arch):
#     """
#     Ordre du vecteur de variables :
#       [P_high, mdot, mdot_HS_fact, mdot_CS_fact]            (IHX, Ejector_IHX)
#       [P_high, mdot, mdot_HS_fact, mdot_CS_fact, eta_IHX]   (EXP_IHX)
#     """
#     if hp_arch == 'EXP_IHX':
#         P_high, mdot, fact_HS, fact_CS, eta_IHX = x
#     else:
#         P_high, mdot, fact_HS, fact_CS = x
#         eta_IHX = None
#     return P_high, mdot, fact_HS, fact_CS, eta_IHX

# #%% Externalized cost function for multiprocessing

# def system_HP_parallel(x, input_data):
#     fluid = input_data['fluid']
#     params = input_data['params']
#     obj = input_data['obj']

#     HSource = MassConnector()
#     CSource = MassConnector()
#     HSource.set_properties(**input_data['HSource'])
#     CSource.set_properties(**input_data['CSource'])

#     AS = CP.AbstractState('HEOS', fluid)
#     AS.update(CP.QT_INPUTS, 0.5, CSource.T)
#     P_sat_T_CSource = AS.p()
#     P_crit_CO2 = AS.p_critical()
#     P_low_guess = 0.8 * min(P_sat_T_CSource, P_crit_CO2)

#     P_high, mdot, fact_HS, fact_CS, eta_IHX = unpack_x(x, params['HP_ARCH'])

#     mdot_HS = mdot * fact_HS
#     mdot_CS = mdot * fact_CS
#     HSource.set_properties(m_dot=mdot_HS)
#     CSource.set_properties(m_dot=mdot_CS)   # <-- débit froid optimisé

#     if params['HP_ARCH'] == 'EXP_IHX':
#         HP = IHX_EXP_CO2_HP(HSource, CSource, params['eta_cp'], params['eta_gc'],
#                             eta_IHX, params['eta_exp'],
#                             params['PP_ev'], params['SH_ev'],
#                             P_low_guess, P_high, mdot, mute_print_flag=1)

#     elif params['HP_ARCH'] == 'IHX':
#         HP = IHX_CO2_HP(HSource, CSource, params['eta_cp'], params['eta_gc'],
#                         params['eta_IHX'], params['PP_ev'], params['SH_ev'],
#                         P_low_guess, P_high, mdot, mute_print_flag=1)

#     elif params['HP_ARCH'] == 'Ejector_IHX':
#         HP = IHX_CO2_HP(HSource, CSource.T, params['eta_cp'], params['eta_gc'],
#                         params['eta_IHX'], params['PP_ev'], params['SH_ev'],
#                         P_low_guess, P_high, mdot, mute_print_flag=1)

#     HP.solve()
    
#     # if not HP.converged:
#     #     return 1000
    
#     if params['HP_ARCH'] == 'EXP_IHX':
#         W_dot_net = HP.components['Compressor'].model.W.W_dot - HP.components['Expander'].model.W.W_dot
#     else:
#         W_dot_net = HP.components['Compressor'].model.W.W_dot

#     COP = HP.components['GasCooler'].model.Q.Q_dot / W_dot_net

#     T_err = abs((HP.components['GasCooler'].model.ex_C.T - obj['T_high']) / obj['T_high'])
#     W_err = abs((W_dot_net - obj['W_dot']) / obj['W_dot'])

#     penalty = 0
#     if T_err > 1e-2:
#         penalty += 300 * T_err
#     if W_err > 1e-2:
#         penalty += 100 * W_err

#     return -COP + penalty

# #%% Optimizer Class

# class CO2HPOptimizer:
#     def __init__(self, fluid):
#         self.fluid = fluid
#         self.params = {}
#         self.it_var = {}
#         self.obj = {}
#         self.HSource = MassConnector()
#         self.CSource = MassConnector()
#         self.COP = None
#         self.optimizer = None

#     def set_parameters(self, **kwargs):
#         self.params.update(kwargs)

#     def set_it_var(self, **kwargs):
#         self.it_var.update(kwargs)

#     def set_obj(self, **kwargs):
#         self.obj.update(kwargs)

#     def opt_HP(self, n_jobs=1, patience=10, max_iter=50, n_part=100):

#         p = self.params

#         lb = [p['P_high_min'], p['m_dot_min'], p['m_dot_HS_fact_min'], p['m_dot_CS_fact_min']]
#         ub = [p['P_high_max'], p['m_dot_max'], p['m_dot_HS_fact_max'], p['m_dot_CS_fact_max']]
#         if p['HP_ARCH'] == 'EXP_IHX':
#             lb.append(p['eta_IHX_min'])
#             ub.append(p['eta_IHX_max'])
#         bounds = (np.array(lb), np.array(ub))

#         input_data = {
#             'fluid': self.fluid,
#             'params': self.params,
#             'obj': self.obj,
#             'HSource': {
#                 'T': self.HSource.T,
#                 'P': self.HSource.p,
#                 'fluid': self.HSource.fluid,
#                 'm_dot': self.it_var['mdot_HS']
#             },
#             'CSource': {
#                 'T': self.CSource.T,
#                 'P': self.CSource.p,
#                 'fluid': self.CSource.fluid,
#                 'm_dot': self.it_var['mdot_CS']
#             }
#         }

#         def objective_wrapper(X):
#             return np.array(Parallel(n_jobs=n_jobs, backend='loky', prefer="processes")(
#                 delayed(system_HP_parallel)(x, input_data) for x in X
#             ))

#         self.optimizer = GlobalBestPSO(
#             n_particles=n_part,
#             dimensions=len(bounds[0]),
#             options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
#             bounds=bounds
#         )

#         best_cost = np.inf
#         no_improve_counter = 0
#         tol = 1e-3

#         for i in tqdm(range(max_iter), desc="Optimizing", ncols=80):
#             self.optimizer.optimize(objective_wrapper, iters=1, verbose=False)
#             current_best = self.optimizer.swarm.best_cost

#             if current_best < best_cost - tol:
#                 best_cost = current_best
#                 no_improve_counter = 0
#             else:
#                 no_improve_counter += 1

#             print(f"[{i+1:03}] Best cost: {best_cost:.6f}")
#             if no_improve_counter >= patience:
#                 print("Stopping early due to stagnation.")
#                 break

#         # Store best variables
#         best_P_high, best_mdot, best_fact_HS, best_fact_CS, best_eta_IHX = unpack_x(
#             self.optimizer.swarm.best_pos, p['HP_ARCH'])

#         best_mdot_HS = best_mdot * best_fact_HS
#         best_mdot_CS = best_mdot * best_fact_CS

#         self.it_var['P_high'] = best_P_high
#         self.it_var['mdot'] = best_mdot
#         self.it_var['mdot_HS'] = best_mdot_HS
#         self.it_var['mdot_CS'] = best_mdot_CS
#         if best_eta_IHX is not None:
#             self.it_var['eta_IHX'] = best_eta_IHX

#         # Rebuild the HP circuit in the main process
#         self.HSource.set_properties(m_dot=best_mdot_HS)
#         self.CSource.set_properties(m_dot=best_mdot_CS)

#         AS = CP.AbstractState('HEOS', self.fluid)
#         AS.update(CP.QT_INPUTS, 0.5, self.CSource.T)
#         P_sat_T_CSource = AS.p()
#         P_crit_CO2 = AS.p_critical()
#         P_low_guess = 0.8 * min(P_sat_T_CSource, P_crit_CO2)

#         if p['HP_ARCH'] == 'EXP_IHX':
#             self.HP = IHX_EXP_CO2_HP(self.HSource, self.CSource, p['eta_cp'],
#                                      p['eta_gc'], self.it_var['eta_IHX'], p['eta_exp'],
#                                      p['PP_ev'], p['SH_ev'],
#                                      P_low_guess, best_P_high, best_mdot, mute_print_flag=1)
#         else:
#             self.HP = IHX_CO2_HP(self.HSource, self.CSource, p['eta_cp'],
#                                  p['eta_gc'], p['eta_IHX'],
#                                  p['PP_ev'], p['SH_ev'],
#                                  P_low_guess, best_P_high, best_mdot, mute_print_flag=1)

#         # Solve in main process and store COP
#         try:
#             self.HP.solve()
#             if p['HP_ARCH'] == 'EXP_IHX':
#                 self.COP = self.HP.components['GasCooler'].model.Q.Q_dot / (
#                     self.HP.components['Compressor'].model.W.W_dot - self.HP.components['Expander'].model.W.W_dot)
#             else:
#                 self.COP = self.HP.components['GasCooler'].model.Q.Q_dot / \
#                            self.HP.components['Compressor'].model.W.W_dot
#         except Exception as e:
#             print(f"⚠️ Failed to solve final HP circuit: {e}")
#             self.COP = None
#             self.HP = None

#         return self.optimizer

# #%% Optimizer call

# if __name__ == "__main__":

#     import matplotlib.pyplot as plt
#     import multiprocessing

#     Optimizer = CO2HPOptimizer('CO2')

#     n_MW = 1 # 4*1e-3

#     T_vec = np.linspace(100, 100, 1) + 273.15
#     COP_vec = []
#     P_high_vec = []
#     m_dot_vec = []
#     m_dot_HS_vec = []
#     m_dot_CS_vec = []
#     T_act_vec = []
#     eta_IHX_vec = []

#     m_dot_HS_fact_min = np.array([0.2])
#     m_dot_HS_fact_max = np.array([0.9])
#     m_dot_CS_fact_min = np.array([10])   # <-- bornes du facteur de débit froid
#     m_dot_CS_fact_max = np.array([10])
#     P_high_min = np.array([80]) * 1e5
#     P_high_max = np.array([130]) * 1e5
#     eta_IHX_min = np.array([0.4])
#     eta_IHX_max = np.array([0.9])

#     for i in range(len(T_vec)):
#         Optimizer.set_parameters(
#             HP_ARCH='IHX',
#             eta_cp=0.7,
#             eta_gc=0.95,
#             eta_IHX=0.8,
#             eta_exp=0.8,
#             PP_ev=5,
#             SH_ev=3,
#             P_high_min=P_high_min[i],
#             P_high_max=P_high_max[i],
#             m_dot_min=10*n_MW,
#             m_dot_max=20*n_MW,
#             m_dot_HS_fact_min=m_dot_HS_fact_min[i],
#             m_dot_HS_fact_max=m_dot_HS_fact_max[i],
#             m_dot_CS_fact_min=m_dot_CS_fact_min[i],
#             m_dot_CS_fact_max=m_dot_CS_fact_max[i],
#             eta_IHX_min=eta_IHX_min[i],
#             eta_IHX_max=eta_IHX_max[i],
#         )

#         Optimizer.set_it_var(
#             P_high=100e5,
#             mdot=20*n_MW,
#             mdot_HS=20*n_MW,
#             mdot_CS=20*n_MW      # valeur initiale du débit froid
#         )

#         Optimizer.set_obj(
#             T_high=T_vec[i],
#             W_dot=n_MW*1e6
#         )

#         Optimizer.CSource.set_properties(
#             T=15 + 273.15,
#             P=5e5,
#             fluid='Water',
#         )

#         Optimizer.HSource.set_properties(
#             T=15 + 273.15,
#             P=10e5,
#             fluid='INCOMP::T66',
#         )

#         n_cores = multiprocessing.cpu_count()

#         opt = Optimizer.opt_HP(n_jobs=n_cores - 1)

#         # Store results
#         COP_vec.append(Optimizer.COP)
#         P_high_vec.append(Optimizer.it_var['P_high'])
#         m_dot_vec.append(Optimizer.it_var['mdot'])
#         m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
#         m_dot_CS_vec.append(Optimizer.it_var['mdot_CS'])
#         T_act_vec.append(T_vec[i] - 273.15)

#         if 'eta_IHX' in Optimizer.it_var:
#             eta_IHX_vec.append(Optimizer.it_var['eta_IHX'])

#         Optimizer.HP.plot_cycle_Ts()

#     # Plotting
#     def make_plot(y, title, ylabel, label=None):
#         plt.figure(figsize=(8, 5))
#         plt.plot(T_act_vec, y, linewidth=2, marker='o', label=label)
#         plt.title(title)
#         plt.xlabel("Temperature (°C)")
#         plt.ylabel(ylabel)
#         plt.grid(True)
#         if label:
#             plt.legend()
#         plt.tight_layout()
#         plt.show()

#     make_plot(COP_vec, "Efficiency vs Temperature", "COP")
#     make_plot(P_high_vec, "High Pressure vs Temperature", "P_high (Pa)", 'P_high')
#     make_plot(m_dot_vec, "Mass Flow Rate vs Temperature", "Mass Flow Rate (kg/s)", 'mdot')
#     make_plot(m_dot_HS_vec, "Heat Source Mass Flow Rate vs Temperature",
#               "Heat Source Mass Flow Rate (kg/s)", 'mdot_HS')
#     make_plot(m_dot_CS_vec, "Cold Source Mass Flow Rate vs Temperature",
#               "Cold Source Mass Flow Rate (kg/s)", 'mdot_CS')


# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

Adapted for multiprocessing-based PSO with optimizer return
"""

#%% Imports

from labothappy.machine.examples.Heat_Pumps.fpi_TC_HeatPump_example import IHX_CO2_HP, IHX_EXP_CO2_HP
from labothappy.connector.mass_connector import MassConnector

import numpy as np
import CoolProp.CoolProp as CP
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm
from multiprocessing import Pool
from joblib import Parallel, delayed

import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")

#%% Externalized cost function for multiprocessing

def system_HP_parallel(x, input_data):
    fluid = input_data['fluid']
    params = input_data['params']
    obj = input_data['obj']

    HSource = MassConnector()
    CSource = MassConnector()
    HSource.set_properties(**input_data['HSource'])
    CSource.set_properties(**input_data['CSource'])
    CSource.set_m_dot(100)

    AS = CP.AbstractState('HEOS', fluid)
    AS.update(CP.QT_INPUTS, 0.5, CSource.T)
    P_sat_T_CSource = AS.p()
    P_crit_CO2 = AS.p_critical()
    P_low_guess = 0.8 * min(P_sat_T_CSource, P_crit_CO2)

    if params['HP_ARCH'] == 'EXP_IHX':
        P_high, mdot, mdot_HS_fact, eta_IHX = x
        mdot_HS = mdot * mdot_HS_fact
        HSource.set_properties(m_dot=mdot_HS)
    else:
        P_high, mdot, mdot_HS_fact, mdot_CS_fact = x
        mdot_HS = mdot * mdot_HS_fact
        mdot_CS = mdot * mdot_CS_fact
        HSource.set_properties(m_dot=mdot_HS)
        CSource.set_properties(m_dot=mdot_CS)

    # try:
    if params['HP_ARCH'] == 'EXP_IHX':
        HP = IHX_EXP_CO2_HP(HSource, CSource, params['eta_cp'], params['eta_gc'],
                            eta_IHX, params['eta_exp'],
                            params['PP_ev'], params['SH_ev'],
                            P_low_guess, P_high, mdot, mute_print_flag=1)
        
    elif params['HP_ARCH'] == 'IHX':
        HP = IHX_CO2_HP(HSource, CSource, params['eta_cp'], params['eta_gc'],
                        params['eta_IHX'], params['PP_ev'], params['SH_ev'],
                        P_low_guess, P_high, mdot, mute_print_flag=1)

    elif params['HP_ARCH'] == 'Ejector_IHX':
        HP = IHX_CO2_HP(HSource, CSource.T, params['eta_cp'], params['eta_gc'],
                        params['eta_IHX'], params['PP_ev'], params['SH_ev'],
                        P_low_guess, P_high, mdot, mute_print_flag=1)


    HP.solve()

    if not HP.converged:
        return 100

    if params['HP_ARCH'] == 'EXP_IHX':
        COP = HP.components['GasCooler'].model.Q.Q_dot / (
            HP.components['Compressor'].model.W.W_dot - HP.components['Expander'].model.W.W_dot)
    else:
        COP = HP.components['GasCooler'].model.Q.Q_dot / HP.components['Compressor'].model.W.W_dot

    T_err = abs((HP.components['GasCooler'].model.ex_C.T - obj['T_high']) / obj['T_high'])

    if params['HP_ARCH'] == 'EXP_IHX':
        W_dot_net = HP.components['Compressor'].model.W.W_dot - HP.components['Expander'].model.W.W_dot
    else:
        W_dot_net = HP.components['Compressor'].model.W.W_dot
        
    W_err = abs((W_dot_net - obj['W_dot']) / obj['W_dot'])

    penalty = 0
    if T_err > 1e-2:
        penalty += 300 * T_err
    if W_err > 1e-2:
        penalty += 100 * W_err

    return -COP + penalty

    # except Exception:
    #     return 100

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

    def opt_HP(self, n_jobs = 1, patience = 10, max_iter = 50):

        if self.params['HP_ARCH'] == 'EXP_IHX':
            bounds = (
                np.array([self.params['P_high_min'], self.params['m_dot_min'], self.params['m_dot_HS_fact_min'], self.params['eta_IHX_min']]),
                np.array([self.params['P_high_max'], self.params['m_dot_max'], self.params['m_dot_HS_fact_max'], self.params['eta_IHX_max']])
            )
        else:
            bounds = (
                np.array([self.params['P_high_min'], self.params['m_dot_min'], self.params['m_dot_HS_fact_min'], self.params['m_dot_CS_fact_min']]), 
                np.array([self.params['P_high_max'], self.params['m_dot_max'], self.params['m_dot_HS_fact_max'], self.params['m_dot_CS_fact_max']])
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
                'fluid': self.CSource.fluid,
                'm_dot': self.it_var['mdot_CS']
            }
        }
    
        # ✅ Open pool once and keep it open during the loop
        def objective_wrapper(X):
            return np.array(Parallel(n_jobs=n_jobs, backend='loky', prefer="processes")(
                delayed(system_HP_parallel)(x, input_data) for x in X
            ))
            
        self.optimizer = GlobalBestPSO(
            n_particles=100,
            dimensions=len(bounds[0]),
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )

        best_cost = np.inf
        no_improve_counter = 0
        tol = 1e-3

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
        
        # Store best variables
        if self.params['HP_ARCH'] == 'EXP_IHX':
            best_P_high = self.optimizer.swarm.best_pos[0]
            best_mdot = self.optimizer.swarm.best_pos[1]
            best_mdot_HS_fact = self.optimizer.swarm.best_pos[2]
            best_mdot_HS = best_mdot * best_mdot_HS_fact
            best_eta_IHX = self.optimizer.swarm.best_pos[3]
            
            best_mdot_CS_fact = self.optimizer.swarm.best_pos[4]
            best_mdot_CS = best_mdot * best_mdot_CS_fact
            
            self.it_var['P_high'] = best_P_high
            self.it_var['mdot'] = best_mdot
            self.it_var['mdot_HS'] = best_mdot_HS
            self.it_var['eta_IHX'] = best_eta_IHX
        else:
            best_P_high = self.optimizer.swarm.best_pos[0]
            best_mdot = self.optimizer.swarm.best_pos[1]
            best_mdot_HS_fact = self.optimizer.swarm.best_pos[2]
            best_mdot_HS = best_mdot * best_mdot_HS_fact
            
            best_mdot_CS_fact = self.optimizer.swarm.best_pos[3]
            best_mdot_CS = best_mdot * best_mdot_CS_fact
            
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
        
        self.CSource.set_properties(m_dot=best_mdot_CS)
        
        if self.params['HP_ARCH'] == 'EXP_IHX':
            self.HP = IHX_EXP_CO2_HP(self.HSource, self.CSource, self.params['eta_cp'],
                                      self.params['eta_gc'], self.it_var['eta_IHX'], self.params['eta_exp'],
                                      self.params['PP_ev'], self.params['SH_ev'],
                                      P_low_guess, best_P_high, best_mdot, mute_print_flag=1)
        else:
            self.HP = IHX_CO2_HP(self.HSource, self.CSource, self.params['eta_cp'],
                                 self.params['eta_gc'], self.params['eta_IHX'],
                                 self.params['PP_ev'], self.params['SH_ev'],
                                 P_low_guess, best_P_high, best_mdot, mute_print_flag=1)
        
        # 🔁 Solve in main process and store COP
        try:
            self.HP.solve()
            if self.params['HP_ARCH'] == 'EXP_IHX':
                self.COP = self.HP.components['GasCooler'].model.Q.Q_dot / (
                    self.HP.components['Compressor'].model.W.W_dot - self.HP.components['Expander'].model.W.W_dot)
            else:
                self.COP = self.HP.components['GasCooler'].model.Q.Q_dot / \
                           self.HP.components['Compressor'].model.W.W_dot
        except Exception as e:
            print(f"⚠️ Failed to solve final HP circuit: {e}")
            self.COP = None
            self.HP = None
        
        bp = self.optimizer.swarm.best_pos
        print("\n" + "="*55)
        print("  OPTIMAL RESULT")
        print("="*55)
        print(f"  P_high            : {bp[0]/1e5:.2f}  bar")
        print(f"  m_dot (CO2)       : {bp[1]:.4f}  kg/s")
        print(f"  m_dot_HS_factor   : {bp[2]:.4f}  [-]")
        print(f"  m_dot_HS          : {bp[1]*bp[2]:.4f}  kg/s")
        print(f"  m_dot_CS_factor   : {bp[3]:.4f}  [-]")
        print(f"  m_dot_CS          : {bp[1]*bp[3]:.4f}  kg/s")
        
        return self.optimizer

#%% Optimizer call

if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    # Initialize optimizer
    Optimizer = CO2HPOptimizer('CO2')
    
    n_MW = 4.4*1e-3
    
    # Define temperatures for sweep
    T_vec = np.linspace(80, 80, 1) + 273.15
    COP_vec = []
    P_high_vec = []
    m_dot_vec = []
    m_dot_HS_vec = []
    T_act_vec = []
    eta_IHX_vec = []
    
    # Sweep parameters
    m_dot_HS_fact_min = np.array([2])
    m_dot_HS_fact_max = np.array([3])
    m_dot_CS_fact_min = np.array([20])
    m_dot_CS_fact_max = np.array([30])
    P_high_min = np.array([90]) * 1e5
    P_high_max = np.array([130]) * 1e5
    eta_IHX_min = np.array([0.8] * 6)
    eta_IHX_max = np.array([0.9] * 6)
    
    for i in range(len(T_vec)):
        Optimizer.set_parameters(
            HP_ARCH='IHX',
            eta_cp=0.7,
            eta_gc=0.8,
            eta_IHX=0.8,
            eta_exp=0.8,
            PP_ev=3,
            SH_ev=0.1,
            P_high_min=P_high_min[i],
            P_high_max=P_high_max[i],
            m_dot_min=8*n_MW,
            m_dot_max=30*n_MW,
            m_dot_HS_fact_min=m_dot_HS_fact_min[i],
            m_dot_HS_fact_max=m_dot_HS_fact_max[i],
            m_dot_CS_fact_min=m_dot_CS_fact_min[i],
            m_dot_CS_fact_max=m_dot_CS_fact_max[i],
            eta_IHX_min = eta_IHX_min[i],
            eta_IHX_max = eta_IHX_max[i],
        )
    
        Optimizer.set_it_var(
            P_high=100e5,
            mdot=17*n_MW,
            mdot_HS=20*n_MW,
            mdot_CS=20*n_MW
        )
    
        Optimizer.set_obj(
            T_high=T_vec[i],
            W_dot=n_MW*1e6
        )
    
        Optimizer.CSource.set_properties(
            T=15 + 273.15,
            P=1e5,
            fluid='Water'
        )
    
        Optimizer.HSource.set_properties(
            T=20 + 273.15,
            P=10e5,
            fluid='INCOMP::T66'
        )
        
        import multiprocessing
        n_cores = multiprocessing.cpu_count()
        
        # Run optimization and get optimizer instance
        opt = Optimizer.opt_HP(n_jobs = n_cores - 1, patience = 5)
    
        # Store results
        COP_vec.append(Optimizer.COP)
        P_high_vec.append(opt.swarm.best_pos[0])
        m_dot_vec.append(opt.swarm.best_pos[1])
        m_dot_HS_vec.append(opt.swarm.best_pos[1] * opt.swarm.best_pos[2])
        T_act_vec.append(T_vec[i] - 273.15)
        
        if 'eta_IHX' in Optimizer.it_var:
            eta_IHX_vec.append(Optimizer.it_var['eta_IHX'])
    
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
