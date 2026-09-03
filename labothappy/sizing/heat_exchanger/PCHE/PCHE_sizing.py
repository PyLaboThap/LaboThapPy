"""
Author : Basile Chaudoir
"""

import os
# os.environ.setdefault("OMP_NUM_THREADS", "1")
# os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
# os.environ.setdefault("MKL_NUM_THREADS", "1")
# os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
# os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

from labothappy.component.base_component import BaseComponent
from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive

from labothappy.toolbox.heat_exchangers.PCHE.thicknesses import PCHE_thicknesses
from labothappy.toolbox.economics.cpi_data import actualize_price

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

warnings.filterwarnings('ignore')

@contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar."""
    class TqdmBatchCompletionCallback(BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_cb = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_cb
        tqdm_object.close()

_SOLVER = None  # cached per-process

#%%

class PCHESizingOpt(BaseComponent):
    
    _CORR_KEYS = {
        'H_Corr', 'C_Corr', 'H_DP', 'C_DP',
        'htc_type', 'DP_type',
        'UD_H_HTC', 'UD_C_HTC', 'UD_H_DP', 'UD_C_DP',
    }
    _CONSTRAINT_KEYS = {'Q_dot', 'DP_h', 'DP_c'}
    
    def __init__(self):
        super().__init__()
        
        self.HX = HexMBChargeSensitive('PCHE')
        self.bounds = {}
        
        return

    #%%

    def set_corr(self, H_Corr=None, C_Corr=None, H_DP=None, C_DP=None, htc_type="Correlation", DP_type="Correlation", UD_H_HTC=None, UD_C_HTC=None, UD_H_DP=None, UD_C_DP=None):
    
        # Sauvegarder pour réutilisation dans simulate_HX
        self.corr_params = {
            'H_Corr': H_Corr, 'C_Corr': C_Corr,
            'H_DP': H_DP, 'C_DP': C_DP,
            'htc_type': htc_type, 'DP_type': DP_type,
            'UD_H_HTC': UD_H_HTC, 'UD_C_HTC': UD_C_HTC,
            'UD_H_DP': UD_H_DP, 'UD_C_DP': UD_C_DP,
        }
    
        self._apply_corr(self.HX)
    
    def _apply_corr(self, HX):
        """Applique les corrélations sauvegardées sur n'importe quel objet HX."""
        c = self.corr_params
    
        if c['htc_type'] == "Correlation":
            HX.set_htc(Corr_H=c['H_Corr'], Corr_C=c['C_Corr'])
        else:
            HX.set_htc(UD_H_HTC=c['UD_H_HTC'], UD_C_HTC=c['UD_C_HTC'])
    
        HX.params['H_DP_ON'] = self.params['H_DP_ON']
        HX.params['C_DP_ON'] = self.params['C_DP_ON']
    
        if c['DP_type'] == "Correlation":
            HX.set_DP(DP_type='Correlation_Disc', Corr_H=c['H_DP'], Corr_C=c['C_DP'])
        else:
            HX.set_DP(DP_type='User-Defined', UD_H_DP=c['UD_H_DP'], UD_C_DP=c['UD_C_DP'])


    def set_bounds(self, **bounds):
        for key, value in bounds.items():
            self.bounds[key] = value
    
    def _apply_deferred_parameters(self):
        """
        set_parameters() reste inchangé (hérité de BaseComponent) et stocke
        tout dans self.params, y compris les clés de corrélation et de
        contraintes. Cette méthode les en extrait et les dispatche vers
        set_corr()/set_constraints() — appelée juste avant l'optimisation.
        """
        corr_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CORR_KEYS}
        constraint_kwargs = {k: self.params.pop(k) for k in list(self.params) if k in self._CONSTRAINT_KEYS}
    
        if corr_kwargs:
            self.set_corr(**corr_kwargs)
        if constraint_kwargs:
            self.set_constraints(**constraint_kwargs)

    def set_constraints(self, Q_dot = None, DP_h = None, DP_c = None):
        
        self.Q_dot_constr = Q_dot 
        self.DP_h_constr = DP_h
        self.DP_c_constr = DP_c
        
        return

    #%%
    def export_params_dict(self):
        
        return {
            "type": "PCHE",
            "type_channel": "Zigzag",
            "Flow_Type": self.params["Flow_Type"],
            
            "Q_dot": self.HX.Q.Q_dot,
            "DP_h": self.HX.DP_h,
            "DP_c": self.HX.DP_c,
            
            "fluid_H": self.inputs['fluid_H'],            
            "m_dot_H": self.inputs['m_dot_H'],
            "T_su_H": self.inputs['T_su_H'],
            "P_su_H": self.inputs['P_su_H'],
            
            "fluid_C": self.inputs['fluid_C'],            
            "m_dot_C": self.inputs['m_dot_C'],
            "T_su_C": self.inputs['T_su_C'],
            "P_su_C": self.inputs['P_su_C'],

            "alpha": self.params["alpha"],
            "D_c": self.params["D_c"],
            "k_cond": self.params["k_cond"],
            "L_c": self.params["L_c"],            
            "N_c": self.params["N_c"],
            "N_p": self.params["N_p"],
            "R_p": self.params["R_p"],
            "t_2": self.params["t_2"],
            "t_3": self.params["t_3"],
            "L_x": self.params["L_x"],
            "L_y": self.params["L_y"],
            "L_z": self.params["L_z"],
            
            "CAPEX": self.CAPEX,
        } 

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
        
        PF = 1
        
        # Objective Function : HX Mass
        rho_mat = 7850 # kg/m^3
        self.m_HX = self.params['n_parallel']*rho_mat*(self.params['L_x'] * self.params['L_y'] * self.params['L_z'] - (self.params['C_V_tot'] + self.params['H_V_tot']))
        
        # Penalties 
        if self.Q_dot_constr:
            self.pen_Q = pen_Q = max(self.Q_dot_constr - self.HX.Q.Q_dot,0)
        else:
            pen_Q = 0
        
        if self.DP_h_constr:
            self.pen_DP_h = pen_DP_h = max(self.HX.DP_h - self.DP_h_constr,0)
        else:
            pen_DP_h = 0
            
        if self.DP_c_constr:
            self.pen_DP_c = pen_DP_c = max(self.HX.DP_c - self.DP_c_constr,0)
        else:
            pen_DP_c = 0
        
        self.penalty = PF*(abs(pen_DP_c) + abs(pen_DP_h) + abs(pen_Q))
        
        self.score = self.m_HX + self.penalty
        
        return self.score
    
    def cost_estimation(self):
        """
        Analysis of Supercritical CO2 Brayton Cycle Recuperative Heat Exchanger Size and Capital Cost
        with Variation of Layout Design (2018)
        
        Kyle R. Zada, Ryan Kim, Aaron Wildberger, Carl P. Schalansky
        """
        
        C_m = 40 # €/kg : Cost per kg 
        # Value from Evaluation of thermal-hydraulic performance and economics of
        # Printed Circuit Heat Exchanger (PCHE) for recuperators of Sodiumcooled Fast Reactors 
        # (SFRs) using CO2 and N2 as working fluids (2022)
        # Su Won Lee, Seong Min Shin, SungKun Chung, HangJin Jo
        
        C_UA = 1.77 # $/UA : 
        
        self.U = sum((self.HX.Qvec_h/self.HX.LMTD)*self.HX.w)
        self.UA = self.U*(1/(1/self.HX.A_h + 1/self.HX.A_c))
        
        self.CAPEX_alt = 49.45*self.UA**0.7544
        
        self.CAPEX = {"HX" : actualize_price(C_m*self.m_HX, 2022, "USD"),
                      "Currency" : "USD"}
        
        self.CAPEX["Install"] = self.CAPEX["HX"]*0.35
        self.CAPEX["Total"] = self.CAPEX["HX"] + self.CAPEX["Install"]
        
        return
    
    #%%
    
    def simulate_HX(self, x):
        warnings.filterwarnings('ignore')
    
        alpha      = x[0]
        D_c        = x[1]
        L_x        = x[2]
        L_y        = x[3]
        L_z        = x[4]
        n_parallel = np.round(x[5])
        n_series   = max(1, np.round(x[6]))
    
        # ---- local geometry ----
        P_max = max(self.inputs['P_su_H'], self.inputs['P_su_C']) * 1.5
        T_max = max(self.inputs['T_su_H'], self.inputs['T_su_C']) * 1.2
    
        t_2, t_3 = PCHE_thicknesses(D_c, P_max, T_max)
        L_c = L_x / np.cos(np.pi * alpha / 180)
        N_p = int(np.floor((L_y - t_3) / (D_c / 2 + t_3)))
        N_c = int(np.floor((L_z - t_2) / (D_c + t_2)))
    
        V_channel = np.pi * D_c**2 / 8 * L_c
        R_p = self.params['R_p']
        C_V_tot = 1 / (1 + R_p) * N_p * N_c * V_channel
        H_V_tot = R_p / (1 + R_p) * N_p * N_c * V_channel
    
        # ---- feasibility ----
        V_block = L_x * L_y * L_z
        if (C_V_tot + H_V_tot) >= V_block:
            return 1e6
    
        # ---- build local param dict ----
        local_params = {
            **self.params,
            'alpha': alpha, 'D_c': D_c, 'L_x': L_x, 'L_y': L_y, 'L_z': L_z,
            'n_parallel': n_parallel, 'n_series': n_series,
            't_2': t_2, 't_3': t_3, 'L_c': L_c,
            'N_p': N_p, 'N_c': N_c,
            'C_V_tot': C_V_tot, 'H_V_tot': H_V_tot,
        }
    
        # ---- simulate ----
        HX = HexMBChargeSensitive('PCHE')
        HX.set_inputs(**self.inputs)
        self._apply_corr(HX)
        HX.set_parameters(**local_params)
    
        try:
            HX.solve()
        except Exception:
            return 1e6
    
        # ---- score (fully local, no self writes) ----
        rho_mat = 7850
        m_HX = n_parallel * n_series * rho_mat * (V_block - C_V_tot - H_V_tot)
    
        PF = 1
        pen_Q    = max(self.Q_dot_constr - HX.Q.Q_dot, 0) if self.Q_dot_constr else 0
        pen_DP_h = max(HX.DP_h - self.DP_h_constr,    0) if self.DP_h_constr  else 0
        pen_DP_c = max(HX.DP_c - self.DP_c_constr,    0) if self.DP_c_constr  else 0
    
        score = m_HX + PF * (pen_Q + pen_DP_h + pen_DP_c)
    
        return score

    #%%
    
    def sizing(self):
        self._apply_deferred_parameters()
        
        # Choose a fixed order for variables
        ORDER = ['alpha', 'D_c', 'L_x', 'L_y', 'L_z', 'n_parallel', 'n_series']
        
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
    
        self.best_pos = optimizer.swarm.best_pos
    
        # Final evaluation
        self.simulate_HX(best_pos)  # or best_params if simulateHX expects dict
    
        self.cost_estimation()
    
        return best_pos  # or return best_params
    
    #%%
    
    def sizing_parallel(self, n_jobs=-1, backend="threading", chunksize="auto", n_particles = 30, max_iter = None, patience = 10):
        self._apply_deferred_parameters()
        
        # ---- fixed order + bounds ----
        ORDER = ['alpha', 'D_c', 'L_x', 'L_y', 'L_z', 'n_parallel', 'n_series']
        def bounds_dict_to_arrays(bounds_dict, order=ORDER):
            lb = np.array([bounds_dict[k][0] for k in order], dtype=float)
            ub = np.array([bounds_dict[k][1] for k in order], dtype=float)
            if np.any(lb > ub):
                raise ValueError("Lower bound > upper bound for at least one variable.")
            return lb, ub
    
        lb, ub = bounds_dict_to_arrays(self.bounds, ORDER)   # <-- pass ORDER
        D = lb.size
        
        optimizer = ps.single.GlobalBestPSO(
            n_particles=n_particles, dimensions=D,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=(lb, ub)
        )
        
        tol = 1e-3
        
        if max_iter is None:
            max_iter = patience * 10

        # after creating `optimizer`
        try:
            n_particles = int(optimizer.swarm.n_particles)
        except AttributeError:
            # fallback if API differs
            n_particles = int(getattr(optimizer, "n_particles", len(optimizer.swarm.position)))
        
        total_evals = n_particles * int(max_iter)
        # pbar = tqdm(total=total_evals, desc="Function evaluations", unit="eval")
    
        def objective_wrapper(X):
            Xp = [np.clip(xi, lb, ub) for xi in np.asarray(X)]
            costs = Parallel(n_jobs=n_jobs, backend=backend, batch_size=chunksize)(
                delayed(self.simulate_HX)(xi) for xi in Xp
            )
            # one update per PSO iteration
            pbar.update(len(Xp))
            return np.asarray(costs, dtype=float)
    
        no_improve, best_cost = 0, np.inf
    
        pbar = tqdm(total=total_evals, unit="eval", dynamic_ncols=True, desc='PCHE')

        for i in range(max_iter):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            curr = optimizer.swarm.best_cost
        
            pbar.clear()   # remove previous bar
            pbar.set_description(f"Iter {i+1} | Best: {curr:.3f}")
            pbar.refresh()
        
            if curr < best_cost - tol:
                best_cost = curr
                no_improve = 0
            else:
                no_improve += 1
        
            if no_improve >= patience:
                break
        
        best_pos = optimizer.swarm.best_pos
        self.score = optimizer.swarm.best_cost
        
        n_series_best = max(1, np.round(best_pos[6]))
        
        # Safe single-threaded final eval — write results to self
        self.params.update({
            'alpha': best_pos[0], 'D_c': best_pos[1],
            'L_x': best_pos[2],   'L_y': best_pos[3], 'L_z': best_pos[4],
            'n_parallel': np.round(best_pos[5]),
            'n_series': n_series_best,
        })
        self.compute_geom()
        
        self.HX = HexMBChargeSensitive('PCHE')
        self.HX.set_inputs(**self.inputs)
        self._apply_corr(self.HX)
        self.HX.set_parameters(**self.params)
        self.HX.solve()
        
        self.m_HX = np.round(best_pos[5]) * n_series_best * 7850 * (
            self.params['L_x'] * self.params['L_y'] * self.params['L_z']
            - self.params['C_V_tot'] - self.params['H_V_tot']
        )
        
        self.cost_estimation()

        self.cost_estimation()
        
        pbar.close()
            
        if __name__ == "__main__":
            print("\nBest Position")
            print("-------------")
            print(f"alpha : {round(self.params['alpha'],2)} [°]")
            print(f"D_c : {round(self.params['D_c']*1e3,2)} [mm]")
            print(f"L_x : {round(self.params['L_x'],3)} [m]")
            print(f"L_y : {round(self.params['L_y'],3)} [m]")
            print(f"L_z : {round(self.params['L_z'],3)} [m]")
            print(f"n_parallel : {self.params['n_parallel']} [-]")
            print(f"n_series : {self.params['n_series']} [-]")       
            
            print("\nResults")
            print("-------------")
            print(f"A_h : {round(self.HX.A_h,2)} [m^2]")
            print(f"A_c : {round(self.HX.A_c,2)} [m^2]")
            print(f"Q_dot : {round(self.HX.Q.Q_dot,1)} [W]")
            print(f"DP_c : {round(self.HX.DP_c,1)} [Pa]")
            print(f"DP_h : {round(self.HX.DP_h,1)} [Pa]")
            print(f"m_HX : {round(self.m_HX,1)} [kg]")
            print(f"CAPEX est. : {round(self.CAPEX['Total'],1)} [$ (2025)]")
    
        return best_pos

#%%

if __name__ == "__main__":
    HX_opt = PCHESizingOpt()

    case_study = 'REC'

    if case_study == "Reference": 
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
        
    elif case_study == "REC":
        HX_opt.set_inputs(
            # First fluid
            fluid_H = 'CO2',
            T_su_H = 298.26, # K
            P_su_H = 4292504, # Pa
            m_dot_H = 298.27, # kg/s
    
            # Second fluid
            fluid_C = 'CO2',
            T_su_C = 286.1, # K
            P_su_C = 15561277, # Pa
            m_dot_C = 298.27, # kg/s  # Make sure to include fluid information
         )
        
    HX_opt.set_parameters(
        k_cond = 60,
        R_p = 1,
        n_disc = 50,
        Flow_Type = 'CounterFlow', 
        H_DP_ON = True, 
        C_DP_ON = True,
        
        H_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"},
        C_Corr = {"1P" : "Gnielinski", "SC" : "Gnielinski"},
        H_DP   = {"1P" : "Gnielinski_DP", "SC" : "Gnielinski_DP"},
        C_DP   = {"1P" : "Gnielinski_DP", "SC" : "Gnielinski_DP"},
        
        Q_dot = 3*1e6,
        DP_h  = 4*1e5,
        DP_c  = 2*1e5,
        )
    
    HX_opt.set_bounds(
        alpha = [10,40],
        D_c = [1*1e-3, 3*1e-3],
        L_x = [0.2, 1.5],
        L_y = [0.3, 2.3],
        L_z = [0.2, 0.6],
        n_parallel = [1,10],
        n_series = [1,10]
        )

    best_pos = HX_opt.sizing_parallel()
    
    