"""
Author : Basile Chaudoir
"""

from component.base_component import BaseComponent
from component.heat_exchanger.hex_MB_charge_sensitive import HeatExchangerMB

from toolbox.heat_exchangers.PCHE.thicknesses import PCHE_thicknesses

import pyswarms as ps
import numpy as np
import warnings

#%%

# ---- joblib worker (process-based) ----
import os, numpy as np
from joblib import Parallel, delayed

from contextlib import contextmanager
from tqdm import tqdm
import joblib
from joblib.parallel import BatchCompletionCallBack

import warnings
warnings.filterwarnings('ignore')


from contextlib import contextmanager
from joblib.parallel import BatchCompletionCallBack
import joblib

@contextmanager
def patch_joblib_tqdm(tqdm_object, close_on_exit=False):
    class _TqdmBatchCompletionCallback(BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)
    old_cb = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = _TqdmBatchCompletionCallback
    try:
        yield
    finally:
        joblib.parallel.BatchCompletionCallBack = old_cb
        if close_on_exit:
            tqdm_object.close()

# @contextmanager
# def patch_joblib_tqdm(tqdm_object, close_on_exit=False):
#     """Route joblib batch completions into a *persistent* tqdm bar."""
#     class _TqdmBatchCompletionCallback(BatchCompletionCallBack):
#         def __call__(self, *args, **kwargs):
#             tqdm_object.update(n=self.batch_size)
#             return super().__call__(*args, **kwargs)

#     old_cb = joblib.parallel.BatchCompletionCallBack
#     joblib.parallel.BatchCompletionCallBack = _TqdmBatchCompletionCallback
#     try:
#         yield
#     finally:
#         joblib.parallel.BatchCompletionCallBack = old_cb
#         if close_on_exit:
#             tqdm_object.close()

# def tqdm_joblib(tqdm_object):
#     warnings.filterwarnings('ignore')

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

#%%

class PCHESizingOpt(BaseComponent):
    
    def __init__(self):
        super().__init__()
        
        self.HX = HeatExchangerMB('PCHE')
        self.bounds = {}
        
        return

    #%%

    def set_corr(self, H_Corr=None, C_Corr=None, H_DP=None, C_DP=None, htc_type = "Correlation", DP_type = "Correlation", UD_H_HTC = None, UD_C_HTC = None, UD_H_DP = None, UD_C_DP = None):

        # Set HTC
        if htc_type == "Correlation":
            self.HX.set_htc(Corr_H = H_Corr, Corr_C = C_Corr) 
        else:
            self.HX.set_htc(UD_H_HTC = UD_H_HTC, UD_C_HTC = UD_C_HTC)
            
        self.HX.params['H_DP_ON'] = self.params['H_DP_ON']
        self.HX.params['C_DP_ON'] = self.params['C_DP_ON']
        
        # Set DP
        if DP_type == "Correlation":
            self.HX.set_DP(DP_type = 'Correlation', Corr_H = H_DP, Corr_C = C_DP)
        else:
            self.HX.set_DP(DP_type = 'User-Defined', UD_H_DP = UD_H_DP, UD_C_DP = UD_C_DP)
        
        return

    def set_bounds(self, **bounds):
        for key, value in bounds.items():
            self.bounds[key] = value

    def set_constraints(self, Q_dot = None, DP_h = None, DP_c = None):
        
        self.Q_dot_constr = Q_dot 
        self.DP_h_constr = DP_h
        self.DP_c_constr = DP_c
        
        return

    #%%

    def compute_geom(self):
        
        P_max = max(self.inputs['P_su_H'], self.inputs['P_su_C'])*1.5
        T_max = max(self.inputs['T_su_H'], self.inputs['T_su_C'])*1.2
        
        self.params['t_2'], self.params['t_3'] = t_2, t_3 = PCHE_thicknesses(self.params['D_c'], P_max, T_max)
        self.params['L_c'] = L_c = self.params['L_x']/np.cos(np.pi*self.params['alpha']/180)
        
        self.params['N_p'] = N_p = np.floor((self.params['L_y']-t_3)/(self.params['D_c']/2 + t_3))
        self.params['N_c'] = N_c = np.floor((self.params['L_z']-t_2)/(self.params['D_c'] + t_2))
        
        V_channel = np.pi*self.params['D_c']**2 / 8 * L_c
        
        self.params['C_V_tot'] = 1/(1+self.params['R_p']) * N_p * N_c * V_channel
        self.params['H_V_tot'] = self.params['R_p']/(1+self.params['R_p']) * N_p * N_c * V_channel
        
        return

    #%%
    
    def compute_score(self):
        
        PF = 0.1
        
        # Objective Function : HX Mass
        rho_mat = 7850 # kg/m^3
        self.m_HX = rho_mat*(self.params['L_x'] * self.params['L_y'] * self.params['L_z'] - (self.params['C_V_tot'] + self.params['H_V_tot'])) 
        
        # Penalties 
        if self.Q_dot_constr:
            self.pen_Q = max(self.Q_dot_constr - self.HX.Q,0)
        else:
            self.pen_Q = 0
        
        if self.DP_h_constr:
            self.pen_DP_h = max(self.HX.DP_h - self.DP_h_constr,0)
        else:
            self.pen_DP_h = 0
            
        if self.DP_c_constr:
            self.pen_DP_c = max(self.HX.DP_c - self.DP_c_constr,0)
        else:
            self.pen_DP_c = 0

        self.penalty = PF*(self.pen_DP_c + self.pen_DP_h + self.pen_Q)
        
        score = self.m_HX + self.penalty
        
        # self.cost_estimation()
        
        # return self.CAPEX + penalty
        return score
    
    def cost_estimation(self):
        """
        Analysis of Supercritical CO2 Brayton Cycle Recuperative Heat Exchanger Size and Capital Cost
        with Variation of Layout Design (2018)
        
        Kyle R. Zada, Ryan Kim, Aaron Wildberger, Carl P. Schalansky
        """
        # C_m = 4 # €/kg : https://mepsinternational.com/gb/en/products/europe-stainless-steel-prices for 316L plates
        
        # C_UA = 1.77 # $/UA : 
        
        self.U = sum((self.HX.Qvec_h/self.HX.LMTD)*self.HX.w)
        self.UA = self.U*(1/(1/self.HX.A_h + 1/self.HX.A_c))
        
        """
        Evaluation of thermal-hydraulic performance and economics of Printed Circuit Heat Exchanger (PCHE)
        for recuperators of Sodiumcooled Fast Reactors (SFRs) using CO2 and N2 as working fluids (2022)
        
        Su Won Lee, Seong Min Shin, SungKun Chung, HangJin Jo
        """
        
        self.CAPEX = 40*self.m_HX # $ (2022) # C_UA*self.UA
        
        return
    
    #%%
    
    def simulate_HX(self, x):
        warnings.filterwarnings('ignore')
        
        self.params['alpha'] = x[0]
        self.params['D_c'] = x[1]
        self.params['L_x'] = x[2]
        self.params['L_y'] = x[3]
        self.params['L_z'] = x[4]        
        
        # Set HX inputs
        self.HX.set_inputs(**self.inputs)
        
        # Compute_geometry
        self.compute_geom()
        
        # Set HX params
        self.HX.set_parameters(**self.params)
        
        # Compute HX 
        try:
            self.HX.solve()
        except Exception:
            self.score = 1_000_000.0
            self.m_HX = np.nan
            self.penalty = np.nan
            self.pen_Q = self.pen_DP_h = self.pen_DP_c = np.nan
            return self.score
    
        return self.compute_score()
        
    #%%
    
    def design(self):
        
        # Choose a fixed order for variables
        ORDER = ['alpha', 'D_c', 'L_x', 'L_y', 'L_z']
        
        def bounds_dict_to_arrays(bounds_dict, order=ORDER):
            lb = np.array([bounds_dict[k][0] for k in order], dtype=float)
            ub = np.array([bounds_dict[k][1] for k in order], dtype=float)
            if np.any(lb > ub):
                raise ValueError("Lower bound > upper bound for at least one variable.")
            return lb, ub
        
        lb, ub = bounds_dict_to_arrays(self.bounds, ORDER)
        D = lb.size
                
        def objective_wrapper(X):
        # X shape: (n_particles, D)
            costs = np.empty(X.shape[0], dtype=float)
            for i, xi in enumerate(X):
                # clip to bounds (safety), then snap alpha to discrete set
                xi = np.clip(xi, lb, ub)
                # if your simulateHX expects a dict, convert here:
                # params = vector_to_params(xi, ORDER)
                # c = self.simulateHX(params)
                c = self.simulate_HX(xi)  # if your simulateHX already accepts vector
                costs[i] = float(c)
            return costs

        optimizer = ps.single.GlobalBestPSO(
            n_particles=40,
            dimensions=D,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=(lb, ub)
        )
    
        patience, tol, max_iter = 5, 1e-3, 40
        no_improve, best_cost = 0, np.inf
    
        for _ in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = optimizer.swarm.best_cost
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve = 0
            else:
                no_improve += 1
            if no_improve >= patience:
                print("Stopping early due to stagnation.")
                break
    
        best_pos = optimizer.swarm.best_pos
    
        # Final evaluation
        self.simulate_HX(best_pos)  # or best_params if simulateHX expects dict
        self.cost_estimation()

        return best_pos  # or return best_params
    
    #%%
    
    # def design_parallel(self, n_jobs=-1, backend="loky", chunksize="auto"):
    #     ORDER = ['alpha', 'D_c', 'L_x', 'L_y', 'L_z']
    
    #     def bounds_dict_to_arrays(bounds_dict, order=ORDER):
    #         lb = np.array([bounds_dict[k][0] for k in order], dtype=float)
    #         ub = np.array([bounds_dict[k][1] for k in order], dtype=float)
    #         if np.any(lb > ub):
    #             raise ValueError("Lower bound > upper bound for at least one variable.")
    #         return lb, ub
    
    #     lb, ub = bounds_dict_to_arrays(self.bounds, ORDER)
    #     D = lb.size
    
    #     optimizer = ps.single.GlobalBestPSO(
    #         n_particles=50, dimensions=D,
    #         options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
    #         bounds=(lb, ub)
    #     )
    
    #     patience, tol = 20, 1e-3
    #     max_iter = patience*10
    #     no_improve, best_cost = 0, np.inf
    
    #     # --- single progress bar for all function evaluations ---
    #     total_evals = max_iter * optimizer.swarm.n_particles
    #     pbar = tqdm(total=total_evals, desc="Function evaluations", unit="eval")
    
    #     def objective_wrapper(X):
    #         Xp = [np.clip(xi, lb, ub) for xi in np.asarray(X)]
    #         # joblib will call our patched callback and update *the same* bar
    #         costs = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
    #             delayed(self.simulate_HX)(xi) for xi in Xp
    #         )
    #         return np.asarray(costs, dtype=float)
    
    #     with patch_joblib_tqdm(pbar, close_on_exit=False):
    #         for i in range(max_iter):
    #             optimizer.optimize(objective_wrapper, iters=1, verbose=False)
    #             curr = optimizer.swarm.best_cost
    #             if curr < best_cost - tol:
    #                 best_cost, no_improve = curr, 0
    #             else:
    #                 no_improve += 1
    #             print(f"[{i+1}] Best cost: {best_cost:.6f}")
    #             if no_improve >= patience:
    #                 print("Stopping early due to stagnation.")
    #                 break
    
    #     pbar.close()
    
    #     best_pos = np.clip(optimizer.swarm.best_pos, lb, ub)
    #     self.simulate_HX(best_pos)
    #     self.cost_estimation()
    
    #     if __name__ == "__main__":
    #         self.HX.plot_cells()
    
    #     print("\nBest Position\n-------------")
    #     print(f"alpha : {round(self.params['alpha'],2)} [°]")
    #     print(f"D_c   : {round(self.params['D_c']*1e3,2)} [mm]")
    #     print(f"L_x   : {round(self.params['L_x'],3)} [m]")
    #     print(f"L_y   : {round(self.params['L_y'],3)} [m]")
    #     print(f"L_z   : {round(self.params['L_z'],3)} [m]")
    
    #     print("\nResults\n-------------")
    #     print(f"A_h   : {round(self.HX.A_h,2)} [m^2]")
    #     print(f"A_c   : {round(self.HX.A_c,2)} [m^2]")
    #     print(f"Q_dot : {round(self.HX.Q,1)} [W]")
    #     print(f"DP_c  : {round(self.HX.DP_c,1)} [Pa]")
    #     print(f"DP_h  : {round(self.HX.DP_h,1)} [Pa]")
    #     print(f"m_HX  : {round(self.m_HX,1)} [kg]")
    #     print(f"CAPEX : {round(self.CAPEX,1)} [$ (2022)]")
    
    #     return best_pos
   
    def design_parallel(self, n_jobs=-1, backend="loky", chunksize="auto", verbose=True):
        ORDER = ['alpha', 'D_c', 'L_x', 'L_y', 'L_z']
    
        def bounds_dict_to_arrays(bounds_dict, order=ORDER):
            lb = np.array([bounds_dict[k][0] for k in order], dtype=float)
            ub = np.array([bounds_dict[k][1] for k in order], dtype=float)
            if np.any(lb > ub):
                raise ValueError("Lower bound > upper bound for at least one variable.")
            return lb, ub
    
        lb, ub = bounds_dict_to_arrays(self.bounds, ORDER)
        D = lb.size
    
        optimizer = ps.single.GlobalBestPSO(
            n_particles=50, dimensions=D,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=(lb, ub)
        )
    
        patience, tol = 20, 1e-3
        max_iter = patience * 10
        no_improve, best_cost = 0, np.inf
    
        total_evals = max_iter * optimizer.swarm.n_particles
        pbar = tqdm(
            total=total_evals,
            desc="Function evaluations",
            unit="eval",
            disable=not verbose,
            position=0,
            leave=False,
            dynamic_ncols=True,
            mininterval=0.2,   # throttle redraws
            maxinterval=1.0,
            smoothing=0.3,
        )
    
        def objective_wrapper(X):
            Xp = [np.clip(xi, lb, ub) for xi in np.asarray(X)]
            costs = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
                delayed(self.simulate_HX)(xi) for xi in Xp
            )
            return np.asarray(costs, dtype=float)
    
        # only touch the postfix occasionally
        POSTFIX_EVERY = 3  # iterations
        with patch_joblib_tqdm(pbar, close_on_exit=False):
            for i in range(max_iter):
                optimizer.optimize(objective_wrapper, iters=1, verbose=False)
                curr = optimizer.swarm.best_cost
    
                improved = curr < best_cost - tol
                if improved:
                    best_cost, no_improve = curr, 0
                else:
                    no_improve += 1
    
                if verbose and (improved or (i+1) % POSTFIX_EVERY == 0):
                    pbar.set_postfix_str(f"iter={i+1} best={best_cost:.3f}", refresh=False)
    
                if no_improve >= patience:
                    if verbose:
                        # finish the current count cleanly and close
                        pbar.total = pbar.n
                        pbar.refresh()
                        pbar.write("Stopping early due to stagnation.")
                    break
    
        pbar.close()
    
        best_pos = np.clip(optimizer.swarm.best_pos, lb, ub)
        self.simulate_HX(best_pos)
        self.cost_estimation()
    
        if verbose and __name__ == "__main__":
            self.HX.plot_cells()
    
        if verbose:
            print("\nBest Position\n-------------")
            print(f"alpha : {round(self.params['alpha'],2)} [°]")
            print(f"D_c   : {round(self.params['D_c']*1e3,2)} [mm]")
            print(f"L_x   : {round(self.params['L_x'],3)} [m]")
            print(f"L_y   : {round(self.params['L_y'],3)} [m]")
            print(f"L_z   : {round(self.params['L_z'],3)} [m]")
    
            print("\nResults\n-------------")
            print(f"A_h   : {round(self.HX.A_h,2)} [m^2]")
            print(f"A_c   : {round(self.HX.A_c,2)} [m^2]")
            print(f"Q_dot : {round(self.HX.Q,1)} [W]")
            print(f"DP_c  : {round(self.HX.DP_c,1)} [Pa]")
            print(f"DP_h  : {round(self.HX.DP_h,1)} [Pa]")
            print(f"m_HX  : {round(self.m_HX,1)} [kg]")
            print(f"CAPEX : {round(self.CAPEX,1)} [$ (2022)]")
    
        return best_pos
    
#%%

if __name__ == "__main__":
    
    HX_opt = PCHESizingOpt()
    
    case_study = "TCO2_RC"
    
    if case_study == "CO2_Test":
        
        HX_opt.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 249 + 273.15, # K
            P_su_H = 96.4*1e5, # Pa
            m_dot_H = 5.35, # kg/s
    
            # Second fluid
            fluid_C = 'CO2',
            T_su_C = 52.77 + 273.15, # K
            P_su_C = 165.4*1e5, # Pa
            m_dot_C = 5.35, # kg/s  # Make sure to include fluid information
            )
    
        HX_opt.set_parameters(
            k_cond = 60, # plate conductivity
            R_p = 1, # n_hot_channel_row / n_cold_channel_row
            
            n_disc = 50,
            
            Flow_Type = 'CounterFlow', 
            H_DP_ON = True, 
            C_DP_ON = True,
            )
    
        H_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
        C_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
        
        H_DP = "Gnielinski_DP"
        C_DP = "Gnielinski_DP"
        
        HX_opt.set_corr(H_Corr, C_Corr, H_DP, C_DP)
        
        # Source for bounds
        # Design and Dynamic Modeling of Printed Circuit Heat Exchangers for Supercritical Carbon Dioxide Brayton Power Cycles
        # Yuan Jiang, Eric Liese, Stephen E. Zitney, Debangsu Bhattacharyya
        
        HX_opt.set_bounds(
            alpha = [10,40], # [°]
            D_c = [1*1e-3, 3*1e-3], # [m]
            L_x = [0.1, 1.5], # [m] : 1.5 limit fixed by Heatric (PCHE manufacturer) : Fluid direction
            L_y = [0.1, 2.3], # [m] : 2.3 limit for shipping requirements : Vertical direction
            L_z = [0.1, 0.6], # [m] : 0.6 limit fixed by Heatric (PCHE manufacturer) : Width
            )
        
        Q_dot_cstr = 1.4*1e6
        DP_c_cstr = 50*1e3
        DP_h_cstr = 50*1e3
        
        HX_opt.set_constraints(Q_dot = Q_dot_cstr, DP_h = DP_h_cstr, DP_c = DP_c_cstr)
    
    if case_study == "TCO2_RC":
        
        HX_opt.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 321.88, # K
            P_su_H = 5742510, # Pa
            m_dot_H = 415.93, # kg/s
    
            # Second fluid
            fluid_C = 'CO2',
            T_su_C = 305.61, # K
            P_su_C = 14153425, # Pa
            m_dot_C = 415.93, # kg/s  # Make sure to include fluid information
            )
    
        HX_opt.set_parameters(
            k_cond = 60, # plate conductivity
            R_p = 1, # n_hot_channel_row / n_cold_channel_row
            
            n_disc = 30,
            
            Flow_Type = 'CounterFlow', 
            H_DP_ON = True, 
            C_DP_ON = True,
            )
    
        H_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
        C_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"}
        
        H_DP = "Darcy_Weisbach"
        C_DP = "Darcy_Weisbach"
        
        HX_opt.set_corr(H_Corr, C_Corr, H_DP, C_DP)
    
        L_x_bounds = np.array([0.1, 1.5])*2
        L_y_bounds = np.array([0.1, 2.3])*2
        L_z_bounds = np.array([0.1, 0.6])*2
    
        HX_opt.set_bounds(
            alpha = [10,40], # [°]
            D_c = [1*1e-3, 3*1e-3], # [m]
            L_x = L_x_bounds, # [m] : 1.5 limit fixed by Heatric (PCHE manufacturer) : Fluid direction
            L_y = L_y_bounds, # [m] : 2.3 limit for shipping requirements : Vertical direction
            L_z = L_z_bounds, # [m] : 0.6 limit fixed by Heatric (PCHE manufacturer) : Width
            )
        
        Q_dot_cstr = 8*1e6 # 7*1e6
        DP_c_cstr = 100*1e3
        DP_h_cstr = 200*1e3
        
        HX_opt.set_constraints(Q_dot = Q_dot_cstr, DP_h = DP_h_cstr, DP_c = DP_c_cstr)
    
    best_pos = HX_opt.design_parallel()
    
    