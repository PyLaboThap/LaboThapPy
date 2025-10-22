"""
Author : Basile Chaudoir
"""

from component.base_component import BaseComponent
from component.heat_exchanger.hex_MB_charge_sensitive import HeatExchangerMB
from CoolProp.CoolProp import PropsSI

# Shell and tube related toolbox
from toolbox.heat_exchangers.shell_and_tubes.pitch_ratio_shell_and_tube import pitch_ratio_fun
from toolbox.heat_exchangers.shell_and_tubes.estimate_tube_in_shell import estimate_number_of_tubes
from toolbox.heat_exchangers.shell_and_tubes.shell_toolbox import shell_thickness
from toolbox.heat_exchangers.shell_and_tubes.tubesheet_toolbox import tube_sheet_thickness
from toolbox.heat_exchangers.shell_and_tubes.baffle_toolbox import baffle_thickness, find_divisors_between_bounds

# Piping toolbox
from toolbox.piping.pipe_thickness import carbon_steel_pipe_thickness, carbon_steel_pipe_thickness_mm

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

warnings.filterwarnings("ignore")

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

# --- worker for joblib ---
def _eval_particle(x, cls, fluid, params, stage_params, inputs):
    """Evaluate one particle using a per-process cached solver."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")

    warnings.filterwarnings("ignore")      # <- suppress in this child process
    np.seterr(all="ignore")                # <- silence NumPy runtime warnings (optional)

    global _SOLVER
    if _SOLVER is None:
        s = cls()
        s.set_parameters(**params)
        if stage_params:
            s.set_stage_parameters(**stage_params)
        s.set_inputs(**inputs)
        _SOLVER = s

    # Re-apply inputs every call to avoid cross-particle contamination
    _SOLVER.set_inputs(**inputs)
    _SOLVER.W_dot = 0  # start clean for this evaluation

    x = np.asarray(x, dtype=float)
    cost = float(_SOLVER.design_system(x))
    return cost

#%%

class ShellAndTubeSizingOpt(BaseComponent):
    
    def __init__(self):
        super().__init__()
        
        self.HX = HeatExchangerMB('Shell&Tube')
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

    def set_gene_space(self, gene_space):
        
        self.gene_space = gene_space
        
        return

    #%%

    def define_spacing(self):
        
        s_min = (self.params['Shell_ID']/5) # [m]
        s_max = (74*(self.params['Tube_OD']/(25.4*1e-3))**(0.75))*25.4*1e-3 # [m]
        
        # feasible integer range for N_b (because spacing = L_eff/(N_b+1) must be in [s_min, s_max])
        N_b_min = max(0, int(np.ceil(self.params['Tube_L']/s_max) - 1))
        N_b_max = max(N_b_min, int(np.floor(self.params['Tube_L']/s_min) - 1))
        
        N_b = int(round(N_b_min + self.u * (N_b_max - N_b_min)))
        s = self.params['Tube_L']/(N_b+1)
        
        return s

    def compute_geom(self):
        """
        Compute rest of geometry
        """

        # pitch_ratio
        self.params['pitch_ratio'] = pitch_ratio_fun(self.params['Tube_OD']*1e3/25.4, self.params['tube_layout'])
        
        # Cross Passes
        self.params['central_spacing'] = self.define_spacing()
        self.params['cross_passes'] = round(self.params['Tube_L']/self.params['central_spacing']) - 1

        # Pipe Thickness
        T_max = self.inputs['T_su_H']
        P_max = max(self.inputs['P_su_H'], self.inputs['P_su_C'])
        
        if self.params['Shell_Side'] == 'H':
            P_shell = self.inputs['P_su_H']
        else:
            P_shell = self.inputs['P_su_C']

        self.params['Tube_t'] = carbon_steel_pipe_thickness_mm(self.params['Tube_OD'], T_max, P_shell, P_max)
                
        # Number of tubes 
        min_tubes_in_row = 4 # 8
        self.params['n_tubes'] = estimate_number_of_tubes(self.params['Shell_ID'], self.params['Tube_OD'], 
                                 self.params['pitch_ratio']*self.params['Tube_OD'], self.params['tube_layout'], min_tubes_in_row)[0]

        # HT Area and HTX volumes
        self.params['A_eff'] = self.params['n_tubes']*self.params['Tube_L']*np.pi*self.params['Tube_OD']
        
        self.params['T_V_tot'] = self.params['Tube_L']*self.params['n_tubes']*np.pi*((self.params['Tube_OD'] - 2*self.params['Tube_t'])/2)**2

        T_V_out = np.pi*(self.params['Tube_OD']/2)**2*self.params['Tube_L']*self.params['n_tubes']
        self.params['S_V_tot'] = self.params['Tube_L']*np.pi*(self.params['Shell_ID']/2)**2 - T_V_out

    #%%
    
    # def compute_score(self):
        
    #     PF = 100
        
    #     # Objective Function : HX Mass
    #     rho_mat = 7850 # kg/m^3
    #     self.m_HX = rho_mat*(self.params['L_x'] * self.params['L_y'] * self.params['L_z'] - (self.params['C_V_tot'] + self.params['H_V_tot'])) 
        
    #     # Penalties 
    #     if self.Q_dot_constr:
    #         pen_Q = max(self.Q_dot_constr - self.HX.Q,0)
    #     else:
    #         pen_Q = 0
        
    #     if self.DP_h_constr:
    #         pen_DP_h = max(self.HX.DP_h - self.DP_h_constr,0)
    #     else:
    #         pen_DP_h = 0
            
    #     if self.DP_c_constr:
    #         pen_DP_c = max(self.HX.DP_c - self.DP_c_constr,0)
    #     else:
    #         pen_DP_c = 0

    #     penalty = PF*(pen_DP_c + pen_DP_h + pen_Q)
        
    #     score = self.m_HX + penalty
        
    #     return score
    
    # def cost_estimation(self):
    #     """
    #     Analysis of Supercritical CO2 Brayton Cycle Recuperative Heat Exchanger Size and Capital Cost
    #     with Variation of Layout Design (2018)
        
    #     Kyle R. Zada, Ryan Kim, Aaron Wildberger, Carl P. Schalansky
    #     """
    #     C_m = 4 # €/kg : https://mepsinternational.com/gb/en/products/europe-stainless-steel-prices for 316L plates
        
    #     C_UA = 1.77 # $/UA : 
        
    #     self.U = sum((self.HX.Qvec_h/self.HX.LMTD)*self.HX.w)
    #     self.UA = self.U*(1/(1/self.HX.A_h + 1/self.HX.A_c))
        
    #     self.CAPEX = C_UA*self.UA
        
    #     return
    
    #%%
    
    def simulate_HX(self, x):
                
        self.params['Tube_L'] = x[0]
        self.params['Tube_OD'] = x[1]
        self.params['Shell_ID'] = x[2]
        self.params['N_pass'] = x[3]
        self.params['tube_layout'] = x[4]        
        self.params['Baffle_Cut'] = x[5]        
        self.u = x[6]        
        
        # Set HX inputs
        self.HX.set_inputs(**self.inputs)
        
        # Compute_geometry
        try:
            self.compute_geom()
        except:
            print("Error in Geom Gen")
            return 10000
            
        # Set HX params
        self.HX.set_parameters(**self.params)
        
        # Compute HX 
        try:
            self.HX.solve()
        except:
            print("Error in Solving")
            return 10000
        
        # Score
        score = self.compute_score()
        self.cost_estimation()
        
        return score
    
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
    
        return best_pos  # or return best_params
    
    #%%
    
    def design_parallel(self, 
                  num_generations=80,
                  sol_per_pop=80,
                  num_parents_mating=20,
                  crossover_probability=0.9,
                  mutation_probability=0.1,
                  parallel_workers=-1,   # threads; set to None to disable
                  early_stop_patience=12,
                  verbose=True):
        import numpy as np, warnings
        import pygad
    
        ORDER = ['L_shell', 'Tube_OD', 'Shell_ID', 'N_pass', 'Tube_arrang', 'Baffle_Cut', 'u']
    
        gene_space = self.gene_space
        if gene_space is None:
            raise ValueError("self.params['gene_space'] must be set (list of choices per gene).")
        if len(gene_space) != len(ORDER):
            raise ValueError(f"gene_space length {len(gene_space)} != expected {len(ORDER)}")
    
        # ensure lists of choices
        gene_space = [np.asarray(g).ravel().tolist() for g in gene_space]
        num_genes = len(gene_space)
    
        best_fitness_so_far = -np.inf
        stagnation = 0
    
        # NEW SIGNATURE REQUIRED BY PyGAD 2.20.0+
        def fitness_func(ga_instance, solution, solution_idx):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sol = np.array(solution, dtype=float)
                # ensure integer where needed
                sol[3] = int(round(sol[3]))  # N_pass
                cost = float(self.simulate_HX(sol))
                return -cost  # PyGAD maximizes
    
        def on_generation(ga_instance):
            nonlocal best_fitness_so_far, stagnation
            bf = ga_instance.best_solution()[1]
            if bf > best_fitness_so_far + 1e-6:
                best_fitness_so_far, stagnation = bf, 0
            else:
                stagnation += 1
            if verbose:
                print(f"[Gen {ga_instance.generations_completed:03d}] best_fitness={bf:.6f}  (no_improve={stagnation}/{early_stop_patience})")
            if stagnation >= early_stop_patience:
                ga_instance.stop_generation = True
            
        n_workers = os.cpu_count()-1
    
        ga = pygad.GA(
            num_generations=num_generations,
            num_parents_mating=num_parents_mating,
            sol_per_pop=sol_per_pop,
            num_genes=num_genes,
            fitness_func=fitness_func,            # <- new signature
            gene_space=gene_space,                # <- your choice vectors
            parent_selection_type="tournament",
            crossover_probability=crossover_probability,
            mutation_probability=mutation_probability,
            keep_parents=2,
            allow_duplicate_genes=True,
            on_generation=on_generation,
            parallel_processing=("thread", n_workers)
        )
    
        ga.run()
    
        solution, fitness, _ = ga.best_solution()
        best_pos = np.array(solution, dtype=float)
        best_pos[3] = int(round(best_pos[3]))  # N_pass
    
        # final eval + optional plot
        self.simulate_HX(best_pos)
        if hasattr(self, "HX") and hasattr(self.HX, "plot_cells"):
            try: self.HX.plot_cells()
            except Exception: pass
    
        print("\nBest Position")
        print("-------------")
        print(f"L_shell     : {best_pos[0]:.3f} [m]")
        print(f"Tube_OD     : {best_pos[1]*1e3:.2f} [mm]")
        print(f"Shell_ID    : {best_pos[2]*1e3:.2f} [mm]")
        print(f"N_pass      : {int(best_pos[3])} [-]")
        print(f"Tube_arrang : {best_pos[4]:.0f} [°]")
        print(f"Baffle_Cut  : {best_pos[5]:.1f} [%]")
    
        if hasattr(self, "HX"):
            print("\nResults")
            print("-------------")
            for name in ["A_h", "A_c", "Q", "DP_c", "DP_h"]:
                if hasattr(self.HX, name):
                    val = getattr(self.HX, name)
                    if name in ("A_h", "A_c"):
                        print(f"{name} : {val:.2f} [m^2]")
                    elif name == "Q":
                        print(f"Q_dot : {val:.1f} [W]")
                    else:
                        print(f"{name} : {val:.1f} [Pa]")
    
        return best_pos

    
#%%

if __name__ == "__main__":
    HX_opt = ShellAndTubeSizingOpt()

    HX_opt.set_inputs(
                  # Hot Fluid
                  T_su_H = 273.15 + 26, # K
                  P_su_H = 10*1e5, # 51.75*1e3, # Pa
                  m_dot_H = 5.35, # kg/s
                  fluid_H = 'Water',
                  
                  # Cold Fluid
                  h_su_C = PropsSI('H','T', 273.15+7,'Q',0,'R134a')+1, # K
                  P_su_C = PropsSI('P','T', 273.15+7,'Q',0,'R134a'), # 51.75*1e3, # Pa
                  m_dot_C = 1.62, # kg/s
                  fluid_C = 'R134a'
                  )

    gene_space = [
        np.arange(1,15,0.01),  # L_shell [m]
        np.array([0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5])*25.4*1e-3, # Tube_OD [m]
        np.array([8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
            29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120])*25.4*1e-3, # Shell_ID [m]
        np.array([1,2,4,6,8,10]), # Tube_pass number [-]
        np.array([0,45,60]), # Tube arrang [°]
        np.arange(15,45,0.1), # Baffle Cut [%]
        np.arange(0,1,50) # u gene for central spacing
        ]

    HX_opt.set_parameters(
            n_series = 1, # [-]
            # OPTI -> Oui (regarder le papier pour déterminer ça)

            foul_t = 0.000176, # (m^2 * K/W)
            foul_s =  0.000176, # (m^2 * K/W)
            tube_cond = 50, # W/(m*K)
            Overdesign = 0,
            
            Shell_Side = 'H',

            Flow_Type = 'Shell&Tube',
            H_DP_ON = True,
            C_DP_ON = True,
            n_disc = 30,
            )

    HX_opt.set_gene_space(gene_space)

    H_Corr = {"1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
    C_Corr = {"1P" : "Gnielinski", "2P" : "Flow_boiling"}
    
    H_DP = "Shell_Kern_DP"
    C_DP = "Muller_Steinhagen_Heck_DP"
    
    HX_opt.set_corr(H_Corr, C_Corr, H_DP, C_DP)

    Q_dot_cstr = 0.313*1e6
    DP_h_cstr = 8.2*1e3
    DP_c_cstr = 21.7*1e3
    
    HX_opt.set_constraints(Q_dot = Q_dot_cstr, DP_h = DP_h_cstr, DP_c = DP_c_cstr)

    # best_pos = HX_opt.design_parallel()
    
    HX_opt.simulate_HX([4.18, 0.375*25.4*1e-3, 0.2032, 2, 60, 25, 0.5]) 
    