# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 13:31:47 2025

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize, root_scalar

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyswarms as ps
import warnings
warnings.filterwarnings("ignore")


#%%

# # ---- joblib worker (process-based) ----
# import os, numpy as np
# from joblib import Parallel, delayed

# from contextlib import contextmanager
# from tqdm import tqdm
# import joblib
# from joblib.parallel import BatchCompletionCallBack

# @contextmanager
# def tqdm_joblib(tqdm_object):
#     """Context manager to patch joblib to report into tqdm progress bar."""
#     class TqdmBatchCompletionCallback(BatchCompletionCallBack):
#         def __call__(self, *args, **kwargs):
#             tqdm_object.update(n=self.batch_size)
#             return super().__call__(*args, **kwargs)

#     old_cb = joblib.parallel.BatchCompletionCallBack
#     joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
#     try:
#         yield tqdm_object
#     finally:
#         joblib.parallel.BatchCompletionCallBack = old_cb
#         tqdm_object.close()

# _SOLVER = None  # cached per-process

# # --- worker for joblib ---
# def _eval_particle(x, cls, fluid, params, stage_params, inputs):
#     """Evaluate one particle using a per-process cached solver."""
#     os.environ.setdefault("OMP_NUM_THREADS", "1")
#     os.environ.setdefault("MKL_NUM_THREADS", "1")
#     os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
#     os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")

#     global _SOLVER
#     if _SOLVER is None:
#         s = cls(fluid)
#         s.set_parameters(**params)
#         if stage_params:
#             s.set_stage_parameters(**stage_params)
#         s.set_inputs(**inputs)
#         _SOLVER = s

#     # Re-apply inputs every call to avoid cross-particle contamination
#     _SOLVER.set_inputs(**inputs)
#     _SOLVER.W_dot = 0  # start clean for this evaluation

#     x = np.asarray(x, dtype=float)
#     cost = float(_SOLVER.design_system(x))
#     return cost, float(_SOLVER.W_dot)

# #%%
#     def design_parallel(self, n_jobs=-1, backend="loky", chunksize="auto"):
#         import numpy as np
#         import pyswarms as ps
    
#         # --- always 6D swarm (psi, phi, R, Re_min, r_m, M_1_st) ---
#         bounds = (np.array([
#             self.params['psi_bounds'][0],
#             self.params['phi_bounds'][0],
#             self.params['R_bounds'][0],
#             self.params['Re_bounds'][0],
#             self.params['r_m_bounds'][0],
#             self.params['M_1st_bounds'][0],
#         ]), np.array([
#             self.params['psi_bounds'][1],
#             self.params['phi_bounds'][1],
#             self.params['R_bounds'][1],
#             self.params['Re_bounds'][1],
#             self.params['r_m_bounds'][1],
#             self.params['M_1st_bounds'][1],
#         ]))
#         dimensions = 6
    
#         # snapshot of class + parameters
#         inp = dict(self.inputs)
    
#         def pick(*names, default=None):
#             for n in names:
#                 if n in inp and inp[n] is not None:
#                     return inp[n]
#             return default
    
#         # this dict is updated each iteration to carry the latest W_dot target
#         inputs_snapshot = {
#             "p0_su": pick("p0_su", "P0_su", "P_su", "p_su"),
#             "T0_su": pick("T0_su", "t0_su", "T_su", "t_su"),
#             "p_ex" : pick("p_ex", "P_ex"),
#             "mdot" : pick("mdot", "m_dot"),
#             "W_dot": pick("W_dot", "W"),
#         }
    
#         snapshot = {
#             "cls": type(self),
#             "fluid": self.fluid,
#             "params": dict(self.params),
#             "stage_params": getattr(self, "stage_params", None),
#             "inputs": inputs_snapshot,
#         }
    
#         def objective_wrapper(X):
#             X = np.asarray(X, dtype=float)
#             with tqdm_joblib(tqdm(total=len(X), desc="Particles", unit="pt")):
#                 results = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
#                     delayed(_eval_particle)(
#                         xi, snapshot["cls"], snapshot["fluid"],
#                         snapshot["params"], snapshot["stage_params"],
#                         snapshot["inputs"]
#                     ) for xi in X
#                 )
#             # results: list of (cost, W_dot)
#             costs = [r[0] for r in results]
#             wdots = [r[1] for r in results]
#             self._last_batch_max_wdot = max(wdots) if wdots else self.inputs.get("W_dot", 0.0)
#             return np.asarray(costs, dtype=float)
    
#         # --- PSO optimizer ---
#         optimizer = ps.single.GlobalBestPSO(
#             n_particles=40, dimensions=dimensions,
#             options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
#             bounds=bounds
#         )
    
#         patience, tol, max_iter = 5, 1e-3, 40
#         no_improve, best_cost = 0, float("inf")
    
#         for i in range(max_iter):
#             # ensure workers see the current target THIS iteration
#             snapshot["inputs"]["W_dot"] = self.inputs.get("W_dot", snapshot["inputs"]["W_dot"])
    
#             optimizer.optimize(objective_wrapper, iters=1, verbose=False)
#             cur = optimizer.swarm.best_cost
    
#             # --- between-iteration W_dot raise ---
#             batch_best = getattr(self, "_last_batch_max_wdot", self.inputs.get("W_dot", 0.0))
#             if batch_best > self.inputs.get("W_dot", 0.0):
#                 self.inputs["W_dot"] = batch_best
#                 # optional trace:
#                 # print(f"[iter {i+1}] raised target W_dot to {self.inputs['W_dot']:.3f} W")
    
#             if cur < best_cost - tol:
#                 best_cost, no_improve = cur, 0
#             else:
#                 no_improve += 1
#             print(f"[{i+1}] Best cost: {best_cost:.6f}")
#             if no_improve >= patience:
#                 print("Stopping early due to stagnation.")
#                 break
    
#         best_pos = optimizer.swarm.best_pos
    
#         # Finalize
#         self.design_system(best_pos)
#         self.cost_estimation()

#         print(f"Work Coef : {self.psi}")
#         print(f"Flow Coef : {self.phi}")
#         print(f"Reaction : {self.R}")
#         print(f"Re : {self.params['Re_min']}")
#         print(f"r_m  : {self.r_m} [m]")
#         print(f"M_1st  : {self.params['M_1_st']}")
    
#         print(f"P_in : {self.stages[0].static_states['P'][1]} [Pa]")
#         print(f"P_out: {self.stages[-1].static_states['P'][2]} [Pa]")
#         print(f"Omega: {self.params['Omega']} [RPM]")
#         print(f"eta_is: {self.eta_is}")
#         print(f"W_dot : {self.W_dot} [W]")
#         return best_pos

#%%

class RadialTurbineMeanLineDesign(object):

    def __init__(self, fluid):
        # Inputs
        self.inputs = {}
        
        # Params
        self.params = {}  

        # Abstract State 
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Blade Dictionnary
        self.stages = []

        # Velocity Triangle Data
        self.Vel_Tri_R = {}
        self.Vel_Tri_S = {}
        
        # Blade Row Efficiency
        self.eta_blade_row = None
        
        self.total_states  = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3,4,5])
        self.static_states = pd.DataFrame(columns=['H','S','P','D','A','V'], index = [1,2,3,4,5])
        self.AS = CP.AbstractState('HEOS', fluid)
            
    def update_total_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.total_states['H'][position] = self.AS.hmass()            
        self.total_states['S'][position] = self.AS.smass()            
        self.total_states['P'][position] = self.AS.p()            
        self.total_states['D'][position] = self.AS.rhomass()            

        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1  
            
        self.total_states['V'][position] = self.AS.viscosity()            
        
        return
    
    def update_static_AS(self, CP_INPUTS, input_1, input_2, position):
        self.AS.update(CP_INPUTS, input_1, input_2)
        
        self.static_states['H'][position] = self.AS.hmass()            
        self.static_states['S'][position] = self.AS.smass()            
        self.static_states['P'][position] = self.AS.p()            
        self.static_states['D'][position] = self.AS.rhomass()    
        
        try:        
            self.static_states['A'][position] = self.AS.speed_sound()            
        except:
            self.static_states['A'][position] = -1            
            
        self.static_states['V'][position] = self.AS.viscosity()            

        return
        
    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value
                
    # ---------------- Blade row ----------------------------------------------------------------------
    
    def computeStator(self):
   
        p0loss_volute = 0
        
        p0_2 = self.total_states['P'][1] - p0loss_volute
        T0_2 = self.inputs['T0_su']

        self.update_total_AS(CP.PT_INPUTS, p0_2, T0_2, 2)
  
        h02 = self.total_states['H'][2]
        h2 = h02 - self.Vel_Tri_S['v2']**2 / 2
        s2 = self.total_states['S'][2]
        
        self.update_static_AS(CP.HmassSmass_INPUTS, h2, s2, 2)
        
        h03 = h02
        
        h3 = h03 - (self.Vel_Tri_S['v3']**2)/2            
        h3_s = h2 - (h2-h3)/self.eta_blade_row_stator
        
        self.AS.update(CP.HmassSmass_INPUTS, h3_s, self.total_states['S'][2])
        p3 = self.AS.p()
        
        self.update_static_AS(CP.HmassP_INPUTS, h3, p3, 3)    
        self.update_total_AS(CP.HmassSmass_INPUTS, h03, self.static_states['S'][3], 3)
        
        return
     
    def computeRotor(self):

        self.update_total_AS(CP.HmassP_INPUTS, self.total_states['H'][3], self.total_states['P'][3], 4)
        self.update_static_AS(CP.HmassP_INPUTS, self.static_states['H'][3], self.static_states['P'][3], 4)
        
        h4 = self.static_states['H'][4]
        
        h05 = self.total_states['H'][4] - self.Dh0
        
        h5 = h05 - self.Vel_Tri_R['v5']**2 /2             
        h5_s = h4 - (h4-h5)/self.eta_blade_row_rotor
        
        self.AS.update(CP.HmassSmass_INPUTS, h5_s, self.total_states['S'][4])
        p5 = self.AS.p()
        
        self.update_static_AS(CP.HmassP_INPUTS, h5, p5, 5)        
        self.update_total_AS(CP.HmassSmass_INPUTS, h05, self.static_states['S'][5], 5)        
        
        return          

    def computeVelTriangles(self):
        
        "1) -------- Rotor ------------------------"
        self.Vel_Tri_R['u4'] = np.sqrt(self.Dh0/self.inputs['psi'])
        self.Vel_Tri_R['vm5'] = self.Vel_Tri_R['u4']*self.inputs['phi']
        self.Vel_Tri_R['vm4'] = self.Vel_Tri_R['vm5']/self.inputs['xhi']

        # Rotor Geometry
        self.params['r4'] = self.Vel_Tri_R['u4']/self.omega_rads
        self.params['r5'] = self.params['r4']*self.params['r5_r4_ratio']

        self.Vel_Tri_R['alpha5'] = 0 # No outlet swirl
        self.Vel_Tri_R['vu5'] = self.Vel_Tri_R['vm5']*np.tan(self.Vel_Tri_R['alpha5'])
        self.Vel_Tri_R['v5'] = np.sqrt(self.Vel_Tri_R['vm5']**2 + self.Vel_Tri_R['vu5']**2)

        # Rotor Outlet
        self.Vel_Tri_R['u5'] = self.omega_rads*self.params['r5']
        self.Vel_Tri_R['wu5'] = self.Vel_Tri_R['vu5'] - self.Vel_Tri_R['u5']
        self.Vel_Tri_R['w5'] = np.sqrt(self.Vel_Tri_R['wu5']**2 + self.Vel_Tri_R['vm5']**2)
        self.Vel_Tri_R['beta5'] = np.arctan(self.Vel_Tri_R['wu5']/self.Vel_Tri_R['vm5'])
        
        # Euler Equation determines vt4
        self.Vel_Tri_R['vu4'] = (self.Dh0 + self.Vel_Tri_R['u5']*self.Vel_Tri_R['vm5'])/self.Vel_Tri_R['u4']
        self.Vel_Tri_R['alpha4'] = np.arctan(self.Vel_Tri_R['vu4']/self.Vel_Tri_R['vm4'])
        self.Vel_Tri_R['v4'] = np.sqrt(self.Vel_Tri_R['vu4']**2 + self.Vel_Tri_R['vm4']**2)

        self.Vel_Tri_R['wu4'] = self.Vel_Tri_R['vu4'] - self.Vel_Tri_R['u4']
        self.Vel_Tri_R['w4'] = np.sqrt(self.Vel_Tri_R['wu4']**2 + self.Vel_Tri_R['vm4']**2)
        self.Vel_Tri_R['beta4'] = np.arctan(self.Vel_Tri_R['wu4']/self.Vel_Tri_R['vm4'])

        "2) -------- Stator ------------------------"
        # Outlet
        self.Vel_Tri_S['alpha3'] = self.Vel_Tri_R['alpha4']
        self.Vel_Tri_S['beta3'] = self.Vel_Tri_R['beta4']
        self.Vel_Tri_S['vm3'] = self.Vel_Tri_R['vm4']
        self.Vel_Tri_S['u3'] = self.Vel_Tri_R['u4']
        self.Vel_Tri_S['vu3'] = self.Vel_Tri_R['vu4']
        self.Vel_Tri_S['wu3'] = self.Vel_Tri_R['wu4']
        self.Vel_Tri_S['v3'] = self.Vel_Tri_R['v4']
        self.Vel_Tri_S['w3'] = self.Vel_Tri_R['w4']
        
        # Inlet
        self.Vel_Tri_S['alpha2'] = 0 # design choice
        self.Vel_Tri_S['vm2'] = self.Vel_Tri_S['vm3']
        self.Vel_Tri_S['vu2'] = self.Vel_Tri_S['vm2']*np.tan(self.Vel_Tri_S['alpha2'])
        self.Vel_Tri_S['v2'] = np.sqrt(self.Vel_Tri_S['vm2']**2 + self.Vel_Tri_S['vu2']**2)
        
        self.Vel_Tri_S['wu2'] = self.Vel_Tri_S['vu2'] - self.Vel_Tri_S['u3']
        self.Vel_Tri_R['beta2'] = np.arctan(self.Vel_Tri_S['vu2']/self.Vel_Tri_S['vm2'])

        self.Vel_Tri_S['w2'] = np.sqrt(self.Vel_Tri_S['vm2']**2 + self.Vel_Tri_S['wu2']**2)
        
        return
     
    def design_system(self, x):    
                
        self.inputs['psi'] = x[0]
        self.inputs['phi'] = x[1]
        self.inputs['xhi'] = x[2]
        self.params['r5_r4_ratio'] = x[3]
        
        self.update_total_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        self.update_static_AS(CP.PT_INPUTS, self.inputs['p0_su'], self.inputs['T0_su'], 1)
        
        "------------- 1) Isentropic Expansion Calculation -----------------------------------------------" 
        s_in = self.total_states['S'][1]
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], s_in)
        
        h_is_ex = self.AS.hmass()
        Dh0s = self.total_states['H'][1] - h_is_ex
                
        self.Dh0 = self.inputs['W_dot']/self.inputs['mdot']
        self.eta_is = self.Dh0/Dh0s

        hout = self.total_states['H'][1] - Dh0s * self.eta_is
        
        self.AS.update(CP.HmassP_INPUTS, hout, self.inputs['p_ex'])
        rho_out = self.AS.rhomass()
        
        "------------- 2) Preliminary design parameters -----------------------------------------------" 
        
        if self.eta_is >= 0.87:
            self.omega_s = 0.55
        else:   
            def eta_is_func(omega_s):
                return 0.87 - 1.07*(omega_s - 0.55)**2 - 0.5*(omega_s - 0.55)**3
            
            def find_omega_s(eta_target, bounds=(0.1, 0.55)): # bounds can be inverted to 0.1 to 0.55 
                func = lambda omega: eta_is_func(omega) - eta_target
                sol = root_scalar(func, bracket=bounds, method='brentq')
                if sol.converged:
                    return sol.root
                else:
                    raise ValueError("No solution found in given bounds")
        
            self.omega_s = find_omega_s(self.eta_is)
        
        V_dot = self.inputs['mdot']/rho_out

        self.omega_rads = self.omega_s * (Dh0s)**0.75 /np.sqrt(V_dot)    
        self.Omega = self.omega_rads*60/(2*np.pi)
        
        self.v_s = 0.737*self.omega_s**0.2
        
        "------------- 3) Rotor Sizing -------------------------------------"
        
        # 3.1 Rotor Tip Sizing
        
        v_0s = np.sqrt(2*Dh0s) # discharge spouting velocity
        self.Vel_Tri_R['u4'] = self.v_s*v_0s
        self.Vel_Tri_R['vu4'] = self.Vel_Tri_R['u4']*self.eta_is/(2*self.v_s**2)
        
        h5 = self.static_states['H'][1] - self.Dh0
        self.params['r4'] = self.Vel_Tri_R['u4']/self.omega_rads
        
        Pt4_guess = self.total_states['P'][1]-self.static_states['D'][1]*Dh0s*(1-self.eta_is)/4

        # 3.2 Rotor Design Specifications
        self.Vel_Tri_R['alpha4'] = np.pi*(10.8+14.2*self.omega_s**2)/180
        
        self.params['tb4'] = 0.04*self.params['r4']
        self.params['tb5'] = 0.02*self.params['r4']
        self.params['rh5'] = 0.185*self.params['r4']
        
        self.Vel_Tri_R['vu5'] = 0 # assumed

        self.params['n_blade_R'] = np.round(12+0.03*(33 - self.Vel_Tri_R['alpha4']*180/np.pi),0)       


        "------------- 4) Velocity Triangle Computation -------------------------------------"         

        self.computeVelTriangles()
        
        "------------- 5) Find eta_blade_row and stator sizing ------------------------------------"         
        
        def find_eta_blade(x):
            self.eta_blade_row_stator = x[0]
            self.eta_blade_row_rotor = x[1]
            
            self.computeStator()
    
            self.computeStator()
            self.computeRotor()
    
            pn_comp = self.static_states['P'][5]

            return (self.inputs["p_ex"] - pn_comp)**2
        
        try:
            self.sol = minimize(find_eta_blade, [1,1], bounds=[(self.eta_is-0.1, 1), (self.eta_is-0.1, 1)], tol = 1e-4)
        except:
            obj = 1000
            print(f"-------------------------")
            print(f"eta_is : {obj}")
            return obj
    
        self.exit_loss = self.inputs['mdot']*(self.Vel_Tri_R['v5']**2)/2        
        
        self.params['b4'] = self.inputs['mdot']/(2*np.pi*self.static_states['D'][4]*self.params['r4']*self.Vel_Tri_R['vm4'])
        self.params['b5'] = self.inputs['mdot']/(2*np.pi*self.static_states['D'][4]*self.params['r5']*self.Vel_Tri_R['vm5'])
        self.params['b3'] = self.params['b4']
        
        self.AS.update(CP.PSmass_INPUTS, self.inputs['p_ex'], self.static_states['S'][1])
        h5s = self.AS.hmass()
        
        self.eta_is = (self.total_states['H'][1] - self.static_states['H'][5])/(self.total_states['H'][1] - h5s)
        
        if self.sol.success:
            obj = -self.eta_is
        else:
            obj = 1000
        
        print(f"-------------------------")
        print(f"eta_is : {obj}")
        print(f"eta_is : {self.eta_blade_row_stator}")
        print(f"eta_is : {self.eta_blade_row_rotor}")
        
        return obj

#%%

    def design(self):
        bounds = (np.array([
            self.params['psi_bounds'][0],
            self.params['phi_bounds'][0],
            self.params['xhi_bounds'][0],
            self.params['r5_r4_bounds'][0],
        ]),
        np.array([
            self.params['psi_bounds'][1],
            self.params['phi_bounds'][1],
            self.params['xhi_bounds'][1],
            self.params['r5_r4_bounds'][1],
        ]))
    
        def objective_wrapper(x):
            rounded_x = np.copy(x)
            costs, wdots = [], []
            for xi in rounded_x:
                c = self.design_system(xi)
                costs.append(c)
            return np.asarray(costs, dtype=float)
    
        optimizer = ps.single.GlobalBestPSO(
            n_particles=10,
            dimensions=4,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        patience = 5
        tol = 1e-3
        max_iter = 40
        no_improve_counter = 0
        best_cost = np.inf
    
        for i in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = optimizer.swarm.best_cost
    
            # between-iteration W_dot raise
            batch_best = getattr(self, "_last_batch_max_wdot", self.inputs.get("W_dot", 0.0))
            if batch_best > self.inputs.get("W_dot", 0.0):
                self.inputs["W_dot"] = batch_best
                # print(f"[iter {i+1}] raised target W_dot to {self.inputs['W_dot']:.3f} W")
    
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
            # print(f"[{i+1}] Best cost: {best_cost:.6f}")
            if no_improve_counter >= patience:
                print("Stopping early due to stagnation.")
                break
    
        best_pos = optimizer.swarm.best_pos
        self.design_system(best_pos)
    
        print(f"Parameters : {self.inputs['psi'], self.inputs['phi'], self.inputs['xhi'], self.params['r5_r4_bounds']}")
        print(f"Turbine rotation speed: {self.Omega} [RPM]")
        print(f"Turbine total-to-static efficiency : {self.eta_is} [-]")
        print(f"Turbine Generation : {self.W_dot} [W]")
        return best_pos

Turb = RadialTurbineMeanLineDesign('CO2')

Turb.set_inputs(
    mdot = 100, # kg/s
    W_dot = 4.69*1e6, # W
    p0_su = 140*1e5, # Pa
    T0_su = 273.15 + 121, # K
    p_ex = 39.8*1e5, # Pa
    )

Turb.set_parameters(
    r5_r4_bounds = [0.3,0.7], # [-] : r5/r4 ratio
    psi_bounds = [0.5, 1.5],
    phi_bounds = [0.3, 0.6],
    xhi_bounds = [0.3, 0.6]
    )
    
Turb.design()

