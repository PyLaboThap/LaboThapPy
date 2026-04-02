# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

@author: Basile
"""

#%% Imports

from labothappy.machine.examples.CO2_Transcritical_Circuits.CO2_Transcritical_circuit import REC_CO2_TC, basic_CO2_TC
from connector.mass_connector import MassConnector

import numpy as np
from CoolProp.CoolProp import PropsSI
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm
from joblib import Parallel, delayed

from labothappy.sizing.turbomachinery.turbine.axial.design_1D.mean_line_axial_turbine_loss_model_design_aungier import AxialTurbineMeanLineDesign
from labothappy.sizing.heat_exchanger.shell_and_tube.shell_and_tube_sizing_PSO_parallel import ShellAndTubeSizingOpt
from labothappy.sizing.heat_exchanger.PCHE.PCHE_PSO import PCHESizingOpt
from labothappy.sizing.turbomachinery.pump.radial.radial_pump_0D_design import RadialPumpODDesign

import warnings
warnings.filterwarnings('ignore')

#%%

def system_RC_parallel(x, input_data):
    warnings.filterwarnings('ignore')

    x = np.array(x, dtype=float)

    # --- Discrétisation des variables 3,4,5,6 ---
    discrete_vars = input_data.get('discrete_vars', {})
    for idx, allowed_vals in discrete_vars.items():
        allowed_vals = np.array(allowed_vals, dtype=float)
        x[idx] = allowed_vals[np.argmin(np.abs(allowed_vals - x[idx]))]
        
    # --------------------------------------------

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
        T=cs_props['T'], P=cs_props['P'], fluid=cs_props['fluid'], m_dot=x[1]*x[7])

    P_sat_T_CSource = PropsSI('P', 'T', CSource.T, 'Q', 0.5, fluid)
    P_crit_CO2 = PropsSI('PCRIT', fluid)
    P_low_guess = min(1.3*P_sat_T_CSource, 0.8*P_crit_CO2)
    
    if input_data['RC_ARCH'] == 'REC':
        RC = REC_CO2_TC(HSource, CSource, x[4], params['PP_rec'], params['eta_pp'],
                        params['eta_exp'], x[3], x[5], x[6], params['SC_cd'], 
                        P_low_guess, x[0], x[1], DP_h_rec=params['DP_h_rec'], DP_c_rec=params['DP_c_rec'], 
                        DP_h_gh=params['DP_h_gh'], DP_c_gh=params['DP_c_gh'], DP_h_cond=params['DP_h_cond'],
                        DP_c_cond=params['DP_c_cond'], mute_print_flag=1)

    elif input_data['RC_ARCH'] == 'basic':
        RC = basic_CO2_TC(HSource, CSource, params['PP_gh'], params['PP_rec'], params['eta_pp'],
                          params['eta_exp'], params['eta_gh'], params['PP_cd'], params['SC_cd'],
                          P_low_guess, x[0], x[1], mute_print_flag=1)

    try:
        RC.solve(method='wegstein')
        
        if not RC.converged:
            # cost, penalty, eta
            return 1000, np.inf, np.nan
    
        DP = 50e3
        rho = RC.components['GasHeater'].model.su_H.D
        mdot = RC.components['GasHeater'].model.su_H.m_dot #+ RC.components['GasHeater'].model.su_C.m_dot
        
        eta_pp = 0.8
        pp_power = DP * mdot / (rho * eta_pp)
    
        W_dot_net = (RC.components['Expander'].model.W.W_dot * 0.95
                     - RC.components['Pump'].model.W.W_dot / 0.95
                     - pp_power / 0.95)
        
        eta = W_dot_net / RC.components['GasHeater'].model.Q.Q_dot
    
        penalty = 0.0
        
        RC.eta = eta
        RC.W_dot_net = W_dot_net
        
        Q_cond = RC.components['Condenser'].model.Q.Q_dot
        RC.components['Condenser'].model.equivalent_effectiveness()
        eta_cond = RC.components['Condenser'].model.epsilon
        
        Q_rec = RC.components['Recuperator'].model.Q.Q_dot
        eta_rec = RC.components['Recuperator'].model.epsilon
    
        Q_gh = RC.components['GasHeater'].model.Q.Q_dot
        eta_gh = RC.components['GasHeater'].model.epsilon
        
        eps = 1e-6
        eta_gh  = np.clip(eta_gh,  0.0, 1.0 - eps)
        eta_rec = np.clip(eta_rec, 0.0, 1.0 - eps)
        eta_cond= np.clip(eta_cond,0.0, 1.0 - eps)
        
        objective = (Q_gh*(-np.log(1-eta_gh)) + Q_rec*(-np.log(1-eta_rec)) + Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec + Q_gh)
        
        cost = objective + penalty
                
        # Return all three so we can use them in objective_wrapper
        return cost, penalty, eta

    except Exception:
        return 1000.0, np.inf, np.nan

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
        self.CAPEX = {}
        
        self.allowable_positions = []
        
    def set_inputs(self, **parameters):
        self.inputs.update(parameters)

    def set_parameters(self, **parameters):
        self.params.update(parameters)

    def set_it_var(self, **parameters):
        self.it_var.update(parameters)

    def set_obj(self, **parameters):
        self.obj.update(parameters)


    def opt_RC(self, n_jobs=None, n_particles=30, max_iter=30, patience=10, tol=1e-4, ntop = 10):
        
        import multiprocessing, numexpr as ne
        n_cores = multiprocessing.cpu_count()
        
        if n_jobs is None:
            n_jobs = n_cores - 1
    
        # keep NumExpr small-ish per process
        ne.set_num_threads(min(2, max(1, n_jobs)))
    
        # create the pool ONCE
        if n_jobs > 1:
            parallel = Parallel(n_jobs=n_jobs, backend='loky')
        else:
            parallel = None
        
        # Raccourcis locaux
        eta_gh_disc   = self.params['eta_gh_disc']
        PP_gh_disc    = self.params['PP_gh_disc']
        eta_rec_disc  = self.params['eta_rec_disc']
        PP_cd_disc    = self.params['PP_cd_disc']
    
        bounds = (
            np.array([
                self.params['P_high_bounds'][0],        # 0 P_high
                self.params['m_dot_bounds'][0],         # 1 m_dot
                self.params['m_dot_HS_fact_bounds'][0], # 2 m_dot_HS_fact (continu)
                eta_gh_disc[0],                         # 3 eta_gh (disc)
                PP_gh_disc[0],                          # 4 PP_gh (disc)
                eta_rec_disc[0],                        # 5 eta_rec (disc)
                PP_cd_disc[0],                          # 6 PP_cd (disc)
                self.params['m_dot_CS_fact_bounds'][0]  # 7 m_dot_CS_fact (continu ici)
            ]),
            np.array([
                self.params['P_high_bounds'][1],
                self.params['m_dot_bounds'][1],
                self.params['m_dot_HS_fact_bounds'][1],
                eta_gh_disc[-1],
                PP_gh_disc[-1],
                eta_rec_disc[-1],
                PP_cd_disc[-1],
                self.params['m_dot_CS_fact_bounds'][1]
            ])
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
                'fluid': self.CSource.fluid,
                'T': self.CSource.T,
                'P': self.CSource.p,
                # 'm_dot': self.CSource.m_dot # !!!
            },
            'RC_ARCH': self.params['RC_ARCH'],
            
            'discrete_vars': {
                3: eta_gh_disc,   # eta_gh
                4: PP_gh_disc,    # PP_gh
                5: eta_rec_disc,  # eta_rec
                6: PP_cd_disc,    # PP_cd
            }}
                      
        #%% 1) Optimize with pre-set params
        
        def objective_wrapper(X):
            discrete_vars = input_data.get('discrete_vars', {})
    
            def discretize(x):
                x = np.array(x, dtype=float)
                for idx, allowed_vals in discrete_vars.items():
                    allowed_vals = np.array(allowed_vals, dtype=float)
                    x[idx] = allowed_vals[np.argmin(np.abs(allowed_vals - x[idx]))]
                return x
    
            if parallel is None or n_jobs == 1:
                # pure serial evaluation, no joblib overhead
                results = [system_RC_parallel(x, input_data) for x in X]
            else:
                results = parallel(
                    delayed(system_RC_parallel)(x, input_data) for x in X
                )
    
            results = np.array(results)
            costs    = results[:, 0]
            penalties = results[:, 1]
            etas      = results[:, 2]
    
            eta_obj = self.obj['eta']
    
            for x_i, pen_i, cost in zip(X, penalties, costs):
                # print(f"x_i : {x_i}")
                # print(f"pen_i : {pen_i}")
                # print(f"cost : {cost}")
                
                if pen_i == 0 and np.isfinite(cost):
                    x_disc = discretize(x_i)
                    self.allowable_positions.append({
                        'x': x_disc.copy(),
                        'score': float(cost)
                    })
            
            print(costs)
            
            return costs

        self.optimizer = GlobalBestPSO(
            n_particles=n_particles,
            dimensions=len(bounds[0]),
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        best_cost = np.inf
        no_improve_counter = 0
        
        patience = patience
        tol = tol
        max_iter = max_iter
    
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
                pbar.set_description("Stopped (patience)")
                break
        
        pbar.close()
    
        # Recompute system with best parameters
        best_pos = self.optimizer.swarm.best_pos.copy()
        discrete_vars = input_data.get('discrete_vars', {})
        
        for idx, allowed_vals in discrete_vars.items():
            allowed_vals = np.array(allowed_vals, dtype=float)
            best_pos[idx] = allowed_vals[np.argmin(np.abs(allowed_vals - best_pos[idx]))]
        
        (best_P_high,
         best_m_dot,
         best_m_dot_HS_fact,
         best_eta_gh,
         best_PP_gh,
         best_eta_rec,
         best_PP_cd,
         best_m_dot_CS_fact) = best_pos
        
        # --- build unique positions with their best (lowest) score ---
        unique_positions = {}
        
        for entry in self.allowable_positions:
            x = np.array(entry['x'], dtype=float)
            score = float(entry['score'])
        
            # make a hashable key (round a bit to avoid float noise)
            key = tuple(np.round(x, 8))
        
            # keep the best score for each unique position
            if key not in unique_positions or score < unique_positions[key]['score']:
                unique_positions[key] = {
                    'x': x,
                    'score': score,
                }
        
        #%% 2) Get best positions and size the component (better perf estimate + CAPEX)
        
        # convert dict → list
        unique_list = list(unique_positions.values())
        
        # sort by score (ascending = best first)
        unique_list.sort(key=lambda e: e['score'])
        
        # keep only the 10 best unique positions
        self.top_positions = unique_list[:ntop]

    def cycle_design(self, n_jobs=None, n_particles=30, max_iter=30, patience=10, tol=1e-4, ntop = 5):
                
        self.criterion = 0
        n_it_max = 10
        
        it = 0
        
        #%% 1) Import RC
        
        # while self.criterion == 0 and it < n_it_max:
        
        self.opt_RC(n_jobs=n_jobs, n_particles=n_particles, max_iter=max_iter, patience=patience, tol=tol, ntop = ntop)
        
        return self.optimizer

#%% Optimizer call

if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    T_test = 130 + 273.15 # K
    
    n_MW = 1 # W
    W_dot_obj = n_MW*1e6 # W
    
    eta_obj = 0.12
    
    # Create optimizer instance
    Optimizer = CO2RCOptimizer('CO2')
    
    # Sweep parameters
    m_dot_HS_fact_bounds = [0.5,1]
    m_dot_CS_fact_bounds = [10,200]
    P_high_bounds = np.array([100, 180]) * 1e5
    m_dot_bounds = np.array([10,80])*n_MW
    
    # Discrete Variable choices
    eta_gh_disc = np.arange(0.8,0.99,0.01)
    PP_gh_disc = np.arange(1,10,1)
    eta_rec_disc= np.arange(0.7,0.98,0.02)
    PP_cd_disc = np.arange(1,15,1)
    
    # # GasHeater
    # eta_gh_disc = np.array([0.95])
    # PP_gh_disc = np.array([5])
    # eta_rec_disc = np.array([0.8])
    # PP_cd_disc = np.array([5])

    # Set model parameters
    Optimizer.set_parameters(
        RC_ARCH= 'REC', # 'REC'
        
        # Pump
        eta_pp=0.8,
        
        # GasHeater
        DP_h_gh = 50*1e3, # 100*1e3,
        DP_c_gh = 50*1e3, # 4*1e5,

        # Recuperator
        PP_rec=0,
        DP_h_rec = 50*1e3, # 4*1e5,
        DP_c_rec = 50*1e3, # 2*1e5,
        
        # Expander
        eta_exp=0.9,
        
        # Condenser
        SC_cd=0.1,
        DP_h_cond = 50*1e3, # 2*1e5,
        DP_c_cond = 50*1e3, # 100*1e3,
        
        # Bounds
        P_high_bounds=P_high_bounds,
        m_dot_HS_fact_bounds=m_dot_HS_fact_bounds,
        m_dot_CS_fact_bounds=m_dot_CS_fact_bounds,
        m_dot_bounds = m_dot_bounds,
        
        # Discrete Values
        eta_gh_disc=eta_gh_disc,
        PP_gh_disc=PP_gh_disc,
        eta_rec_disc=eta_rec_disc,
        PP_cd_disc=PP_cd_disc,
    )

    # Initial guess
    Optimizer.set_it_var(
        P_high=100e5,
        mdot=17,
        mdot_HS=20,
        
        # GasHeater
        eta_gh=0.95,
        PP_gh=5,
        
        # Recuperator
        eta_rec=0.8,
        
        # Condenser
        PP_cd=5,
    )

    # Objective
    Optimizer.set_obj(
        W_dot=W_dot_obj,
        eta=eta_obj
        )

    # Source definitions
    Optimizer.CSource.set_properties(
        T=15 + 273.15,
        P=5e5,
        fluid='Water',
    )

    Optimizer.HSource.set_properties(
        T=T_test,
        P=5e5,
        fluid='Water',
    )

    # Prepare model and optimize
    import time
    
    # --- First run (cold) ---
    # t0 = time.perf_counter()
    # Optimizer.opt_RC(n_particles=10, max_iter=5, patience=5)
    # t1 = time.perf_counter()
    # print(f"First run time (cold): {t1 - t0:.4f} s")
    
    # # --- Second run (warm) ---
    # t2 = time.perf_counter()
    # Optimizer.opt_RC(n_particles=10, max_iter=5, patience=5)
    # t3 = time.perf_counter()
    # print(f"Second run time (warm): {t3 - t2:.4f} s")

    Optimizer.cycle_design(ntop = 5, n_particles=100)

