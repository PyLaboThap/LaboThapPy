# -*- coding: utf-8 -*-
"""
CO2 Transcritical Rankine Cycle Optimizer
- Parallel PSO evaluation via joblib (File 1 architecture)
- Rich cycle logic, penalty logging, warm start (File 2 logic)
"""

#%% Imports

from labothappy.machine.examples.ORC.fpi_TC_orc_example import REC_CO2_TC, basic_CO2_TC, Recomp_CO2_TC, Recomp_CO2_TC_1_recup
from labothappy.connector.mass_connector import MassConnector

import numpy as np
from CoolProp.CoolProp import PropsSI
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm
from joblib import Parallel, delayed
import multiprocessing

import warnings
warnings.filterwarnings('ignore')

#%% Top-level parallel evaluation function
# Must be defined at module level for joblib/loky to pickle it.

def system_RC_parallel(x, input_data):
    """
    Evaluate one particle x.
    Non-Recomp: x = [P_high, m_dot, m_dot_HS_fact, m_dot_CS_fact]
    Recomp/Recomp_1_recup: x = [P_high, m_dot, m_dot_HS_fact, m_dot_CS_fact, spliter_frac]
    Returns scalar cost (lower = better); large positive = infeasible.
    """
    warnings.filterwarnings('ignore')

    fluid    = input_data['fluid']
    params   = input_data['params']
    obj      = input_data['obj']
    hs_props = input_data['HSource']
    cs_props = input_data['CSource']
    arch     = input_data['RC_ARCH']

    P_high        = x[0]
    m_dot         = x[1]
    m_dot_HS_fact = x[2]
    m_dot_CS_fact = x[3]          # <-- position [3], commun à toutes les architectures
    m_dot_HS      = m_dot * m_dot_HS_fact
    m_dot_CS      = m_dot * m_dot_CS_fact

    # --- Build connectors ---
    HSource = MassConnector()
    HSource.set_properties(T=hs_props['T'], P=hs_props['P'],
                           fluid=hs_props['fluid'], m_dot=m_dot_HS)

    CSource = MassConnector()
    CSource.set_properties(T=cs_props['T'], P=cs_props['P'],
                           fluid=cs_props['fluid'], m_dot=m_dot_CS)

    # --- Low pressure initial guess ---
    P_sat_CS    = PropsSI('P', 'T', cs_props['T'], 'Q', 0.5, fluid)
    P_crit      = PropsSI('PCRIT', fluid)
    P_low_guess = min(1.3 * P_sat_CS, 0.8 * P_crit)

    # --- Build cycle ---
    try:
        if arch == 'REC':
            RC = REC_CO2_TC(
                HSource, CSource,
                params['PP_gh'], params['PP_rec'],
                params['eta_pp'], params['eta_exp'],
                params['eta_gh'], params['eta_rec'],
                params['PP_cd'], params['SC_cd'],
                P_low_guess, P_high, m_dot,
                DP_h_rec  = params.get('DP_h_rec',  0e5),
                DP_c_rec  = params.get('DP_c_rec',  0e5),
                DP_h_gh   = params.get('DP_h_gh',   0e5),
                DP_c_gh   = params.get('DP_c_gh',   0e5),
                DP_h_cond = params.get('DP_h_cond', 0e5),   
                DP_c_cond = params.get('DP_c_cond', 0e5),  
            )
        elif arch == 'basic':
            RC = basic_CO2_TC(
                HSource, CSource,
                params['PP_gh'], params['eta_pp'],
                params['eta_exp'], params['eta_gh'],
                params['PP_cd'], params['SC_cd'],
                P_low_guess, P_high, m_dot,
                DP_h_gh   = params.get('DP_h_gh',   0e5),
                DP_c_gh   = params.get('DP_c_gh',   0e5),
                DP_h_cond = params.get('DP_h_cond', 0e5),   
                DP_c_cond = params.get('DP_c_cond', 0e5), 
                mute_print_flag=1,
            )

        elif arch == "Recomp":
            spliter_frac = x[4]      # <-- décalé de [3] à [4]

            RC = Recomp_CO2_TC(
                HSource, CSource, 
                params['PP_gh'], params['PP_rec'], 
                params['eta_pp'], params['eta_exp'], params['eta_cp'], 
                params['eta_rec'], params['eta_rec_HT'], params['eta_gh'],
                params['PP_cd'], params['SC_cd'],
                P_low_guess, P_high, m_dot, spliter_frac,
                DP_h_rec  = params.get('DP_h_rec',  0e5),
                DP_c_rec  = params.get('DP_c_rec',  0e5),
                DP_h_gh   = params.get('DP_h_gh',   0e5),
                DP_c_gh   = params.get('DP_c_gh',   0e5),
                DP_h_cond = params.get('DP_h_cond', 0e5),   
                DP_c_cond = params.get('DP_c_cond', 0e5),  
                mute_print_flag=1)
            
        elif arch == "Recomp_1_recup":
            spliter_frac = x[4]      # <-- décalé de [3] à [4]
    
            RC = Recomp_CO2_TC_1_recup(
                HSource, CSource, 
                params['PP_gh'], params['PP_rec'], 
                params['eta_pp'], params['eta_exp'], params['eta_cp'], 
                params['eta_rec'], params['eta_gh'],
                params['PP_cd'], params['SC_cd'],
                P_low_guess, P_high, m_dot, spliter_frac,
                DP_h_rec  = params.get('DP_h_rec',  0e5),
                DP_c_rec  = params.get('DP_c_rec',  0e5),
                DP_h_gh   = params.get('DP_h_gh',   0e5),
                DP_c_gh   = params.get('DP_c_gh',   0e5),
                DP_h_cond = params.get('DP_h_cond', 0e5),   
                DP_c_cond = params.get('DP_c_cond', 0e5),  
                mute_print_flag=1)
        else:
            return 1000.0
    except Exception:
        return 1000.0

    # --- Solve ---
    try:
        RC.solve()
    except Exception as e:
        return 100.0

    if not getattr(RC, 'converged', True):
        return 10.0

    # --- Superheat check at expander inlet ---
    try:
        T_exp_ex  = RC.components['Expander'].model.ex.T
        P_exp_ex  = RC.components['Expander'].model.ex.p
        T_sat_exp = PropsSI('T', 'P', P_exp_ex, 'Q', 1, 'CO2')
        SH_exp    = T_exp_ex - T_sat_exp
    except Exception:
        SH_exp = 50.0  # assume ok if property call fails

    if SH_exp < 0:
        return 100.0

    # --- Power and efficiency ---
    try:
        if arch == "Recomp" or arch == "Recomp_1_recup":
            W_cp = RC.components['Compressor'].model.W.W_dot
        else:
            W_cp = 0
            
        W_exp  = RC.components['Expander'].model.W.W_dot
        W_pump = RC.components['Pump'].model.W.W_dot
        Q_gh   = RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = RC.components['GasHeater'].model.su_H.D
        m_HS_act   = RC.components['GasHeater'].model.su_H.m_dot
        
        rho_CS     = RC.components['Condenser'].model.su_C.D
        m_CS_act   = RC.components['Condenser'].model.su_C.m_dot
        
        W_pump_aux_HS = params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * params.get('eta_pp_aux', 0.8))

        W_pump_aux_CS = params.get('DP_c_cond', 0.5e5) * m_CS_act / \
                     (rho_CS * params.get('eta_pp_aux', 0.8))

        W_dot_net = W_exp - W_pump - W_pump_aux_HS - W_pump_aux_CS - W_cp
        eta       = W_dot_net / Q_gh if Q_gh > 0 else 0.0
        
    except Exception:
        return 1000.0

    # --- Power target penalty ---
    W_obj   = obj.get('W_dot', 1e6)
    rel_err = abs((W_dot_net - W_obj) / W_obj)
    penalty = 500.0 * (rel_err ** 2) if rel_err > 1e-2 else 0.0

    return -eta + penalty


#%% Optimizer Class

class CO2RC_eff_optimizer:

    def __init__(self, fluid):
        self.fluid  = fluid
        self.RC     = None

        self.inputs  = {}
        self.params  = {}
        self.it_var  = {}
        self.obj     = {}

        self._HSource_props = {}
        self._CSource_props = {}

        self.eta         = None
        self.W_dot_net   = None
        self.penalty_log = {}

    # ------------------------------------------------------------------ setters

    def set_inputs(self, **parameters):
        self.inputs.update(parameters)

    def set_parameters(self, **parameters):
        self.params.update(parameters)

    def set_it_var(self, **parameters):
        self.it_var.update(parameters)

    def set_obj(self, **parameters):
        self.obj.update(parameters)

    def set_HSource(self, T, P, fluid, m_dot=1.0):
        self._HSource_props = dict(T=T, P=P, fluid=fluid, m_dot=m_dot)

    def set_CSource(self, T, P, fluid, m_dot=1000.0):
        self._CSource_props = dict(T=T, P=P, fluid=fluid, m_dot=m_dot)

    # ------------------------------------------------------------------ RC build

    def set_RC(self):
        """Builds self.RC from current it_var and source props."""
        HSource = MassConnector()
        HSource.set_properties(
            T     = self._HSource_props['T'],
            P     = self._HSource_props['P'],
            fluid = self._HSource_props['fluid'],
            m_dot = self.it_var['mdot_HS'],
        )

        CSource = MassConnector()
        CSource.set_properties(**self._CSource_props)

        P_sat_CS    = PropsSI('P', 'T', self._CSource_props['T'], 'Q', 0.5, self.fluid)
        P_crit      = PropsSI('PCRIT', self.fluid)
        P_low_guess = min(1.3 * P_sat_CS, 0.8 * P_crit)

        arch = self.params.get('RC_ARCH', 'REC')

        if arch == 'REC':
            self.RC = REC_CO2_TC(
                HSource, CSource,
                self.params['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'],
                self.params['eta_gh'], self.params['eta_rec'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'],
                DP_h_rec  = self.params.get('DP_h_rec',  0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond', 0e5),
                DP_c_cond = self.params.get('DP_c_cond', 0e5),   
                mute_print_flag=1,
            )
            
        elif arch == 'basic':
            self.RC = basic_CO2_TC(
                HSource, CSource,
                self.params['PP_gh'], self.params['eta_pp'], 
                self.params['eta_exp'], self.params['eta_gh'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'],
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond', 0e5),
                DP_c_cond = self.params.get('DP_c_cond', 0e5),  
                mute_print_flag=1,
            )
        elif arch == 'Recomp':
            self.RC = Recomp_CO2_TC(
                HSource, CSource, 
                self.params['PP_gh'], self.params['PP_rec'], 
                self.params['eta_pp'], self.params['eta_exp'], self.params['eta_cp'], 
                self.params['eta_rec'], self.params['eta_rec_HT'], self.params['eta_gh'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'], self.it_var['spliter_frac'],
                DP_h_rec  = self.params.get('DP_h_rec',  0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond', 0e5),
                DP_c_cond = self.params.get('DP_c_cond', 0e5), 
                mute_print_flag=1)
        elif arch == 'Recomp_1_recup':
            self.RC = Recomp_CO2_TC_1_recup(
                HSource, CSource, 
                self.params['PP_gh'], self.params['PP_rec'], 
                self.params['eta_pp'], self.params['eta_exp'], self.params['eta_cp'], 
                self.params['eta_rec'], self.params['eta_gh'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'], self.it_var['spliter_frac'],
                DP_h_rec  = self.params.get('DP_h_rec',  0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0e5),
                DP_c_gh   = self.params.get('DP_c_gh',   0e5),
                DP_h_cond = self.params.get('DP_h_cond', 0e5),
                DP_c_cond = self.params.get('DP_c_cond', 0e5),   
                mute_print_flag=1)
        else:
            raise ValueError("'RC_ARCH' parameter shall be either 'basic', 'REC', 'Recomp', 'Recomp_1_recup'")

    def _log_penalty(self, reason):
        self.penalty_log[reason] = self.penalty_log.get(reason, 0) + 1

    # ------------------------------------------------------------------ final eval

    def _evaluate_final(self, best_pos):
        if self.params['RC_ARCH'] == "Recomp" or self.params['RC_ARCH'] == "Recomp_1_recup":
            P_high, m_dot, m_dot_HS_fact, m_dot_CS_fact, spliter_frac = best_pos
            m_dot_HS = m_dot * m_dot_HS_fact
            m_dot_CS = m_dot * m_dot_CS_fact

            self.it_var['spliter_frac'] = spliter_frac
            self.it_var['P_high']  = P_high
            self.it_var['mdot']    = m_dot
            self.it_var['mdot_HS'] = m_dot_HS
            self.it_var['mdot_CS'] = m_dot_CS

        else:
            P_high, m_dot, m_dot_HS_fact, m_dot_CS_fact = best_pos
            m_dot_HS = m_dot * m_dot_HS_fact
            m_dot_CS = m_dot * m_dot_CS_fact
    
            self.it_var['P_high']  = P_high
            self.it_var['mdot']    = m_dot
            self.it_var['mdot_HS'] = m_dot_HS
            self.it_var['mdot_CS'] = m_dot_CS
    
        # Débit source froide optimisé
        self._CSource_props['m_dot'] = m_dot_CS   # doit précéder set_RC()
    
        # Update source props so set_RC picks up the optimised m_dot_HS
        self._HSource_props['m_dot'] = m_dot_HS
    
        self.set_RC()
        RC = self.RC

        try:
            RC.solve()
        except Exception as e:
            self._log_penalty(f"Final solve exception: {e}")
            self.eta = self.W_dot_net = None
            return

        if not getattr(RC, 'converged', True):
            self._log_penalty("Final solve did not converge")
            self.eta = self.W_dot_net = None
            return

        # Superheat check
        try:
            T_exp_ex  = RC.components['Expander'].model.ex.T
            P_exp_ex  = RC.components['Expander'].model.ex.p
            T_sat_exp = PropsSI('T', 'P', P_exp_ex, 'Q', 1, 'CO2')
            SH_exp    = T_exp_ex - T_sat_exp
        except Exception:
            SH_exp = 50.0

        if SH_exp < 0:
            self._log_penalty(f"Drops in expansion (SH = {SH_exp:.1f} K)")
            self.eta = self.W_dot_net = None
            return
        
        if self.params['RC_ARCH'] == "Recomp" or self.params['RC_ARCH'] == "Recomp_1_recup":
            W_cp = RC.components['Compressor'].model.W.W_dot
        else:
            W_cp = 0
            
        W_exp  = RC.components['Expander'].model.W.W_dot
        W_pump = RC.components['Pump'].model.W.W_dot
        Q_gh   = RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = RC.components['GasHeater'].model.su_H.D
        m_HS_act   = RC.components['GasHeater'].model.su_H.m_dot
        
        rho_CS     = RC.components['Condenser'].model.su_C.D
        m_CS_act   = RC.components['Condenser'].model.su_C.m_dot
        
        self.W_pump_aux_HS = W_pump_aux_HS = self.params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * self.params.get('eta_pp_aux', 0.8))

        self.W_pump_aux_CS = W_pump_aux_CS = self.params.get('DP_c_cond', 0.5e5) * m_CS_act / \
                     (rho_CS * self.params.get('eta_pp_aux', 0.8))

        W_dot_net = W_exp - W_pump - W_pump_aux_HS - W_pump_aux_CS - W_cp

        self.W_dot_net = W_dot_net
        self.eta       = self.W_dot_net / Q_gh if Q_gh > 0 else 0.0

        # self.Q_dot_waste = RC.components['GasHeater'].model.ex_H.m_dot * (
        #     RC.components['GasHeater'].model.ex_H.h
        #     - PropsSI('H', 'T', 273.15 + 15, 'P',
        #               RC.components['GasHeater'].model.ex_H.p,
        #               RC.components['GasHeater'].model.ex_H.fluid)
        # )
        
    # ------------------------------------------------------------------ optimise

    def opt_RC(self, n_jobs = 1, n_particles=100, max_iter=30, patience=None,
               init_pos=None, warm_spread=0.05, warm_fraction=0.5):
        """
        PSO optimisation with parallel particle evaluation via joblib.

        Parameters
        ----------
        n_particles   : swarm size
        max_iter      : maximum PSO iterations
        patience      : early-stop after this many stagnant iterations
                        (default: max_iter // 5)
        init_pos      : 1-D seed array [P_high, m_dot, HS_factor] for warm start,
                        or 2-D (n_particles, 3) matrix, or None for random init.
        warm_spread   : relative noise around the seed (±fraction)
        warm_fraction : fraction of particles initialised near the seed
        """
        if patience is None:
            patience = max(1, max_iter // 5)
        
        if self.params['RC_ARCH'] == "Recomp" or self.params['RC_ARCH'] == "Recomp_1_recup":
            lb = np.array([
                self.params['P_high_min'],
                self.params['m_dot_min'],
                self.params['m_dot_HS_fact_min'],
                self.params['m_dot_CS_fact_min'],   # <-- position [3]
                self.params['spliter_frac_min'],    # <-- décalé en [4]
            ])
            ub = np.array([
                self.params['P_high_max'],
                self.params['m_dot_max'],
                self.params['m_dot_HS_fact_max'],
                self.params['m_dot_CS_fact_max'],   # <-- position [3]
                self.params['spliter_frac_max'],    # <-- décalé en [4]
            ])
        else:   
            lb = np.array([
                self.params['P_high_min'],
                self.params['m_dot_min'],
                self.params['m_dot_HS_fact_min'],
                self.params['m_dot_CS_fact_min'],   # <-- position [3]
            ])
            ub = np.array([
                self.params['P_high_max'],
                self.params['m_dot_max'],
                self.params['m_dot_HS_fact_max'],
                self.params['m_dot_CS_fact_max'],   # <-- position [3]
            ])

        # --- warm start ---
        pso_init_pos = None
        if init_pos is not None:
            seed = np.asarray(init_pos, dtype=float)
            if seed.ndim == 1:
                seed    = np.clip(seed, lb, ub)
                n_warm  = max(1, int(round(warm_fraction * n_particles)))
                n_rand  = n_particles - n_warm
                noise   = np.random.uniform(-warm_spread, warm_spread,
                                            size=(n_warm, len(lb)))
                warm    = np.clip(seed[None, :] * (1.0 + noise), lb, ub)
                rand    = np.random.uniform(lb, ub, size=(n_rand, len(lb)))
                pso_init_pos = np.vstack([warm, rand])
                print(f"  → Warm start: {n_warm}/{n_particles} particles around "
                      f"P={seed[0]/1e5:.1f} bar, ṁ={seed[1]:.1f}, f={seed[2]:.3f}")
            else:
                pso_init_pos = np.clip(seed, lb, ub)

        # --- pack input_data for pickling ---
        input_data = {
            'fluid'   : self.fluid,
            'params'  : self.params,
            'obj'     : self.obj,
            'HSource' : {
                'T'     : self._HSource_props['T'],
                'P'     : self._HSource_props['P'],
                'fluid' : self._HSource_props['fluid'],
            },
            'CSource' : {
                'T'     : self._CSource_props['T'],
                'P'     : self._CSource_props['P'],
                'fluid' : self._CSource_props['fluid'],
                'm_dot' : self._CSource_props.get('m_dot', 1000.0),
            },
            'RC_ARCH' : self.params.get('RC_ARCH', 'REC'),
        }

        # --- parallel objective wrapper ---
        def objective_wrapper(X):
            return np.array(
                Parallel(n_jobs=n_jobs, backend='loky')(
                    delayed(system_RC_parallel)(x, input_data) for x in X
                )
            )

        # --- PSO ---
        optimizer = GlobalBestPSO(
            n_particles = n_particles,
            dimensions  = len(ub),
            options     = {'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds      = (lb, ub),
            init_pos    = pso_init_pos,
        )

        best_cost  = np.inf
        no_improve = 0

        pbar = tqdm(range(max_iter), desc="PSO Optimizing", ncols=80)
        for i in pbar:
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current = optimizer.swarm.best_cost

            if current < best_cost - 1e-3:
                best_cost  = current
                no_improve = 0
            else:
                no_improve += 1

            pbar.set_postfix(best_cost=f"{best_cost:.6f}")

            if no_improve >= patience:
                pbar.set_description("Stopped (no improvement)")
                break

        pbar.close()

        # --- Final evaluation with full diagnostics ---
        self._evaluate_final(optimizer.swarm.best_pos)

        bp = optimizer.swarm.best_pos
        print("\n" + "="*55)
        print("  OPTIMAL RESULT")
        print("="*55)
        print(f"  P_high            : {bp[0]/1e5:.2f}  bar")
        print(f"  m_dot (CO2)       : {bp[1]:.4f}  kg/s")
        print(f"  m_dot_HS_factor   : {bp[2]:.4f}  [-]")
        print(f"  m_dot_HS          : {bp[1]*bp[2]:.4f}  kg/s")
        print(f"  m_dot_CS_factor   : {bp[3]:.4f}  [-]")
        print(f"  m_dot_CS          : {self.it_var.get('mdot_CS', float('nan')):.4f}  kg/s")
        
        if self.W_dot_net is not None:
            print(f"  W_net             : {self.W_dot_net/1e3:.3f}  kW")
            print(f"  Thermal η         : {self.eta*100:.3f}  %")
        else:
            print("  ⚠️  Final solve failed — see penalty log.")
        print("="*55)

        # --- Penalty summary ---
        print("\n" + "="*55)
        print("  PENALTY SUMMARY")
        print("="*55)
        total = sum(self.penalty_log.values())
        if total:
            for reason, count in sorted(self.penalty_log.items(),
                                        key=lambda kv: kv[1], reverse=True):
                print(f"  [{count:4d} | {count/total*100:5.1f}%] : {reason}")
        else:
            print("  No evaluations logged.")
        print("="*55)

        # Flag bad power match
        if self.W_dot_net is not None:
            W_obj    = self.obj.get('W_dot', 1e6)
            rel_err  = abs((self.W_dot_net - W_obj) / W_obj)
            if rel_err > 0.05:
                print(f"  ⚠️  WARNING: W_net error = {rel_err*100:.1f}%")
                self.eta = np.nan

        self.penalty_log = {}
        return optimizer

#%% Main

if __name__ == "__main__":
    
    case_study = "Comparison"
    
    if case_study == "test":
    
        n_cores = multiprocessing.cpu_count()
        
        import matplotlib.pyplot as plt
    
        # ---- sweep ----
        T_vec = np.linspace(100, 350, 6) + 273.15
        # T_vec = np.linspace(150, 150, 1) + 273.15
        ARCH = ['basic', 'REC', 'Recomp_1_recup', 'Recomp']
        
        # T_vec = np.array([350]) + 273.15
    
        eta_vec      = []
        P_high_vec   = []
        m_dot_vec    = []
        m_dot_HS_vec = []
        T_h_ex_vec   = []
        Q_dot_waste  = []
        
        n_MW = 1
        W_dot_test = n_MW*1e6  # 1 MW target
    
        Optimizer = CO2RC_eff_optimizer('CO2')
    
        for T in T_vec:
    
            Optimizer.set_parameters(
                RC_ARCH = 'REC',
    
                eta_pp  = 0.8,
                eta_gh  = 0.95,
                eta_rec = 0.9,
                eta_rec_HT = 0.9, # Recompression Case
                eta_exp = 0.9,
                eta_cp = 0.8,
                eta_pp_aux = 0.8, 

                PP_gh   = 5,
                PP_rec  = 0,
                PP_cd   = 5,
                SC_cd   = 0.1,
    
                DP_h_gh  = 50e3,
                DP_c_gh  = 50e3,
                DP_h_rec = 50e3,
                DP_c_rec = 50e3,
                DP_h_cond  = 50e3,
                DP_c_cond  = 50e3,
    
                P_high_min       = 80e5,
                P_high_max       = 200e5,
                m_dot_min        = 10.0*n_MW,
                m_dot_max        = 100.0*n_MW,
                m_dot_HS_fact_min = 0.1,
                m_dot_HS_fact_max = 2,
                m_dot_CS_fact_min = 1,
                m_dot_CS_fact_max = 20,
                spliter_frac_min = 0,
                spliter_frac_max = 1
            )
            
            if Optimizer.params['RC_ARCH'] == "Recomp" or Optimizer.params['RC_ARCH'] == "Recomp_1_recup":
                Optimizer.set_it_var(P_high=100e5, mdot=20.0, mdot_HS=15.0, spliter_frac = 1)
            else:
                Optimizer.set_it_var(P_high=100e5, mdot=20.0, mdot_HS=15.0)
                
            Optimizer.set_obj(W_dot=W_dot_test)
    
            Optimizer.set_CSource(T=15 + 273.15, P=5e5,  fluid='Water', m_dot=1000.0)
            Optimizer.set_HSource(T=T, P=10e5, fluid='INCOMP::TVP1', m_dot=50.0)
    
            Optimizer.set_RC()
            Optimizer.opt_RC(n_jobs = n_cores - 1, n_particles=50, max_iter=50, patience = 10)
    
            eta_vec.append(Optimizer.eta)
            P_high_vec.append(Optimizer.it_var['P_high'])
            m_dot_vec.append(Optimizer.it_var['mdot'])
            m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
            T_h_ex_vec.append(Optimizer.RC.components['GasHeater'].model.ex_H.T)
            Q_dot_waste.append(getattr(Optimizer, 'Q_dot_waste', None))
    
        # ---- plots ----
        T_C = T_vec - 273.15
    
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
        axes[0, 0].plot(T_C, eta_vec, linewidth=2)
        axes[0, 0].set_title("Efficiency vs Temperature")
        axes[0, 0].set_xlabel("Temperature (°C)")
        axes[0, 0].set_ylabel("Efficiency [-]")
        axes[0, 0].grid(True)
    
        axes[0, 1].plot(T_C, [p / 1e5 for p in P_high_vec], linewidth=2)
        axes[0, 1].set_title("High Pressure vs Temperature")
        axes[0, 1].set_xlabel("Temperature (°C)")
        axes[0, 1].set_ylabel("P_high [bar]")
        axes[0, 1].grid(True)
    
        axes[1, 0].plot(T_C, m_dot_vec, linewidth=2)
        axes[1, 0].set_title("CO2 Mass Flow Rate vs Temperature")
        axes[1, 0].set_xlabel("Temperature (°C)")
        axes[1, 0].set_ylabel("ṁ_CO2 [kg/s]")
        axes[1, 0].grid(True)
    
        axes[1, 1].plot(T_C, m_dot_HS_vec, linewidth=2)
        axes[1, 1].set_title("Heat Source Mass Flow Rate vs Temperature")
        axes[1, 1].set_xlabel("Temperature (°C)")
        axes[1, 1].set_ylabel("ṁ_HS [kg/s]")
        axes[1, 1].grid(True)
    
        plt.tight_layout()
        plt.show()
    
        Optimizer.RC.plot_cycle_Ts()

    elif case_study == "Comparison":
        # -*- coding: utf-8 -*-
        """
        Script de balayage : efficacité (η) vs température de source chaude,
        une courbe par architecture de cycle CO2 (basic, REC, Recomp_1_recup, Recomp).
        
        Pour chaque (architecture, température), on relance l'optimisation PSO
        N_RUNS fois (le PSO étant stochastique, plusieurs essais réduisent le
        risque de rester bloqué dans un optimum local).
        
        Pour chaque condition (architecture, température), on conserve les
        BEST_KEEP (par défaut 5) meilleurs runs valides rencontrés :
          - si moins de BEST_KEEP runs valides ont été trouvés, on les ajoute tous ;
          - une fois BEST_KEEP atteints, un nouveau run ne remplace le moins bon
            de la liste que s'il fait mieux (eta plus élevé).
        
        Ces meilleurs runs sont enregistrés (et réécrits à chaque mise à jour,
        pour être robuste à une interruption) dans un fichier CSV.
        """
        
        import csv
        import multiprocessing
        import numpy as np
        import matplotlib.pyplot as plt
        
        # ---------------------------------------------------------------------
        # Bornes d'optimisation par architecture et par température (°C)
        # ---------------------------------------------------------------------
        BOUNDS = {
            'REC': {
                100.0: dict(P_high_min=1.449e+07, P_high_max=1.580e+07, m_dot_min=69.197, m_dot_max=73.322, m_dot_HS_fact_min=1.5721, m_dot_HS_fact_max=3.6346, m_dot_CS_fact_min=10.6236, m_dot_CS_fact_max=15.4395),
                150.0: dict(P_high_min=1.604e+07, P_high_max=1.704e+07, m_dot_min=33.763, m_dot_max=35.160, m_dot_HS_fact_min=1.1774, m_dot_HS_fact_max=2.8640, m_dot_CS_fact_min=7.7265, m_dot_CS_fact_max=21.3089),
                200.0: dict(P_high_min=1.858e+07, P_high_max=2e+07, m_dot_min=22.302, m_dot_max=22.680, m_dot_HS_fact_min=0.8287, m_dot_HS_fact_max=2.9539, m_dot_CS_fact_min=16.8675, m_dot_CS_fact_max=19.6892),
                250.0: dict(P_high_min=1.636e+07, P_high_max=2e+07, m_dot_min=16.796, m_dot_max=19.942, m_dot_HS_fact_min=0.5621, m_dot_HS_fact_max=5.1976, m_dot_CS_fact_min=3.6010, m_dot_CS_fact_max=22.1168),
                300.0: dict(P_high_min=1.541e+07, P_high_max=2e+07, m_dot_min=13.836, m_dot_max=16.337, m_dot_HS_fact_min=0.5281, m_dot_HS_fact_max=3.3241, m_dot_CS_fact_min=5.9445, m_dot_CS_fact_max=21.5712),
                350.0: dict(P_high_min=1.466e+07, P_high_max=2e+07, m_dot_min=11.571, m_dot_max=18.151, m_dot_HS_fact_min=0.4916, m_dot_HS_fact_max=5.1040, m_dot_CS_fact_min=0.4834, m_dot_CS_fact_max=19.1780),
            },
            'Recomp': {
                100.0: dict(P_high_min=1.307e+07, P_high_max=1.670e+07, m_dot_min=68.991, m_dot_max=94.019, m_dot_HS_fact_min=1.4913, m_dot_HS_fact_max=5.1649, m_dot_CS_fact_min=8.5376, m_dot_CS_fact_max=19.6382),
                150.0: dict(P_high_min=1.043e+07, P_high_max=1.701e+07, m_dot_min=27.424, m_dot_max=90.502, m_dot_HS_fact_min=0.6013, m_dot_HS_fact_max=4.4915, m_dot_CS_fact_min=11.3029, m_dot_CS_fact_max=17.9534),
                200.0: dict(P_high_min=1.060e+07, P_high_max=1.697e+07, m_dot_min=26.044, m_dot_max=56.034, m_dot_HS_fact_min=0.4544, m_dot_HS_fact_max=5.4410, m_dot_CS_fact_min=6.5479, m_dot_CS_fact_max=21.5822),
                250.0: dict(P_high_min=1.269e+07, P_high_max=1.846e+07, m_dot_min=24.450, m_dot_max=36.267, m_dot_HS_fact_min=0.4694, m_dot_HS_fact_max=3.5117, m_dot_CS_fact_min=3.4953, m_dot_CS_fact_max=18.9767),
                300.0: dict(P_high_min=1.483e+07, P_high_max=2e+07, m_dot_min=14.451, m_dot_max=25.908, m_dot_HS_fact_min=0.2716, m_dot_HS_fact_max=3.3641, m_dot_CS_fact_min=7.3641, m_dot_CS_fact_max=16.2130),
                350.0: dict(P_high_min=7.751e+06, P_high_max=2e+07, m_dot_min=7.686, m_dot_max=64.905, m_dot_HS_fact_min=0.2388, m_dot_HS_fact_max=5.1119, m_dot_CS_fact_min=0.0000, m_dot_CS_fact_max=19.7339),
            },
            'Recomp_1_recup': {
                100.0: dict(P_high_min=1.347e+07, P_high_max=1.751e+07, m_dot_min=72.328, m_dot_max=103.851, m_dot_HS_fact_min=1.3625, m_dot_HS_fact_max=4.1721, m_dot_CS_fact_min=3.1581, m_dot_CS_fact_max=22.0434),
                150.0: dict(P_high_min=1.202e+07, P_high_max=1.558e+07, m_dot_min=29.343, m_dot_max=83.315, m_dot_HS_fact_min=0.3673, m_dot_HS_fact_max=5.3512, m_dot_CS_fact_min=6.8280, m_dot_CS_fact_max=19.5688),
                200.0: dict(P_high_min=1.242e+07, P_high_max=1.896e+07, m_dot_min=19.407, m_dot_max=74.857, m_dot_HS_fact_min=0.9736, m_dot_HS_fact_max=5.3729, m_dot_CS_fact_min=7.6422, m_dot_CS_fact_max=16.3157),
                250.0: dict(P_high_min=9.989e+06, P_high_max=2e+07, m_dot_min=16.573, m_dot_max=63.474, m_dot_HS_fact_min=0.0000, m_dot_HS_fact_max=5.1693, m_dot_CS_fact_min=3.3169, m_dot_CS_fact_max=13.2376),
                300.0: dict(P_high_min=1.036e+07, P_high_max=2e+07, m_dot_min=19.989, m_dot_max=26.812, m_dot_HS_fact_min=0.6716, m_dot_HS_fact_max=3.1881, m_dot_CS_fact_min=10.6725, m_dot_CS_fact_max=16.5460),
                350.0: dict(P_high_min=1.057e+07, P_high_max=2e+07, m_dot_min=13.028, m_dot_max=24.634, m_dot_HS_fact_min=1.5937, m_dot_HS_fact_max=4.9266, m_dot_CS_fact_min=4.2053, m_dot_CS_fact_max=16.5487),
            },
            'basic': {
                100.0: dict(P_high_min=1.488e+07, P_high_max=1.557e+07, m_dot_min=68.180, m_dot_max=70.718, m_dot_HS_fact_min=1.8095, m_dot_HS_fact_max=3.4618, m_dot_CS_fact_min=10.1894, m_dot_CS_fact_max=18.7389),
                150.0: dict(P_high_min=1.938e+07, P_high_max=1.995e+07, m_dot_min=33.787, m_dot_max=34.969, m_dot_HS_fact_min=1.0745, m_dot_HS_fact_max=5.2738, m_dot_CS_fact_min=8.1782, m_dot_CS_fact_max=20.2015),
                200.0: dict(P_high_min=1.785e+07, P_high_max=2e+07, m_dot_min=19.508, m_dot_max=44.410, m_dot_HS_fact_min=0.2997, m_dot_HS_fact_max=4.4355, m_dot_CS_fact_min=5.0229, m_dot_CS_fact_max=20.3731),
                250.0: dict(P_high_min=7.654e+06, P_high_max=2e+07, m_dot_min=6.271, m_dot_max=101.155, m_dot_HS_fact_min=0.1, m_dot_HS_fact_max=5.0931, m_dot_CS_fact_min=7.8644, m_dot_CS_fact_max=19.2222),
                300.0: dict(P_high_min=1.156e+07, P_high_max=2e+07, m_dot_min=12.747, m_dot_max=25.952, m_dot_HS_fact_min=1.2133, m_dot_HS_fact_max=5.2711, m_dot_CS_fact_min=0.8459, m_dot_CS_fact_max=20.8587),
                350.0: dict(P_high_min=9.394e+06, P_high_max=2e+07, m_dot_min=11.152, m_dot_max=23.315, m_dot_HS_fact_min=0.4942, m_dot_HS_fact_max=5.5206, m_dot_CS_fact_min=3.9687, m_dot_CS_fact_max=17.6152),
            },
        }
        
        _AVAILABLE_T_BY_ARCH = {arch: sorted(d.keys()) for arch, d in BOUNDS.items()}
        
        def get_bounds(arch, T_C):
            """Renvoie le dict de bornes pour (arch, T_C). Si T_C n'est pas une clé
            exacte (arrondi/erreur flottante), utilise la température disponible
            la plus proche pour cette architecture."""
            arch_bounds = BOUNDS[arch]
            if T_C in arch_bounds:
                return arch_bounds[T_C]
        
            available = _AVAILABLE_T_BY_ARCH[arch]
            T_used = min(available, key=lambda t: abs(t - T_C))
            if abs(T_used - T_C) > 1.0:
                print(f"    [!] Pas de bornes exactes pour {arch} @ {T_C}°C -> utilisation de {T_used}°C")
            return arch_bounds[T_used]
        
        # ---------------------------------------------------------------------
        # Gestion des "N meilleurs runs" par condition (architecture, température)
        # ---------------------------------------------------------------------
        BEST_KEEP = 5  # nombre de meilleurs runs conservés par condition
        CSV_PATH = "best_runs.csv"
        
        FIELDNAMES = [
            "architecture", "T_C", "rank",
            "P_high_Pa", "m_dot_CO2_kg_s", "m_dot_HS_factor", "m_dot_HS_kg_s",
            "m_dot_CS_factor", "m_dot_CS_kg_s", "W_net_W", "eta",
        ]
        
        # best_runs[(arch, T_C)] = liste des BEST_KEEP meilleurs runs (dicts), triée
        # par eta décroissant
        best_runs = {}
        
        
        def extract_run_data(Optimizer, eta):
            """Récupère les grandeurs d'intérêt sur l'objet Optimizer après
            résolution. Adapter les noms d'attributs ci-dessous si l'API de
            CO2RC_eff_optimizer diffère (ex. Optimizer.m_dot vs Optimizer.mdot)."""
            return {
                "P_high_Pa": getattr(Optimizer, "P_high", None),
                "m_dot_CO2_kg_s": getattr(Optimizer, "m_dot", getattr(Optimizer, "mdot", None)),
                "m_dot_HS_factor": getattr(Optimizer, "m_dot_HS_factor", None),
                "m_dot_HS_kg_s": getattr(Optimizer, "m_dot_HS", None),
                "m_dot_CS_factor": getattr(Optimizer, "m_dot_CS_factor", None),
                "m_dot_CS_kg_s": getattr(Optimizer, "m_dot_CS", None),
                "W_net_W": getattr(Optimizer, "W_net", getattr(Optimizer, "W_dot", None)),
                "eta": eta,
            }
        
        
        def update_best_runs(key, new_run, max_keep=BEST_KEEP):
            """Ajoute new_run à la liste des meilleurs runs pour `key`. Si la liste
            dépasse max_keep après ajout, le run le plus faible (eta le plus bas)
            est écarté -- ce qui peut être new_run lui-même s'il n'est pas assez bon."""
            lst = best_runs.setdefault(key, [])
            lst.append(new_run)
            lst.sort(key=lambda r: r["eta"], reverse=True)
            del lst[max_keep:]
        
        
        def write_best_runs_csv(path=CSV_PATH):
            rows = []
            for (arch, T_C), lst in best_runs.items():
                for rank, run in enumerate(lst, start=1):
                    row = {"architecture": arch, "T_C": T_C, "rank": rank}
                    row.update(run)
                    rows.append(row)
            rows.sort(key=lambda r: (r["architecture"], r["T_C"], r["rank"]))
        
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
                writer.writeheader()
                for row in rows:
                    writer.writerow(row)
        
        
        n_cores = multiprocessing.cpu_count()
        
        # ---- sweep ----
        T_vec = np.linspace(100, 350, 6) + 273.15
        ARCH = ['basic', 'REC', 'Recomp_1_recup', 'Recomp']
        
        N_RUNS = 5  # nombre d'optimisations par point (arch, T)
        
        W_dot_test = 1e6  # 1 MW target
        
        # Résultats : dict arch -> liste de eta (ou None) alignée sur T_vec
        results = {arch: [] for arch in ARCH}
        
        # Détail de tous les essais, pour inspection : results_all[arch][i_T] = liste des eta (ou None) des N_RUNS essais
        results_all = {arch: [] for arch in ARCH}
        
        for arch in ARCH:
            print(f"\n{'='*60}")
            print(f"  Architecture : {arch}")
            print(f"{'='*60}")
        
            for T in T_vec:
                T_C = round(T - 273.15, 1)
                print(f"\n--- T = {T_C:.1f} °C ---")
        
                b = get_bounds(arch, T_C)
                condition_key = (arch, T_C)
        
                run_etas = []  # eta valides sur les N_RUNS essais de ce point
        
                for run_idx in range(N_RUNS):
                    print(f"  Run {run_idx+1}/{N_RUNS}...")
        
                    eta_run = None
        
                    try:
                        Optimizer = CO2RC_eff_optimizer('CO2')
        
                        Optimizer.set_parameters(
                            RC_ARCH = arch,
                            eta_pp  = 0.8,
                            eta_gh  = 0.95,
                            eta_rec = 0.9,
                            eta_rec_HT = 0.9,  # utile seulement pour Recomp
                            eta_exp = 0.9,
                            eta_cp  = 0.8,
                            PP_gh   = 5,
                            PP_rec  = 0,
                            PP_cd   = 5,
                            SC_cd   = 0.1,
                            eta_pp_aux = 0.8,
        
                            DP_h_gh  = 50e3,
                            DP_c_gh  = 50e3,
                            DP_h_rec = 50e3,
                            DP_c_rec = 50e3,
                            DP_h_cond  = 50e3,
                            DP_c_cond  = 50e3,
        
                            P_high_min        = b['P_high_min'],
                            P_high_max        = b['P_high_max'],
                            m_dot_min         = b['m_dot_min'],
                            m_dot_max         = b['m_dot_max'],
                            m_dot_HS_fact_min = b['m_dot_HS_fact_min'],
                            m_dot_HS_fact_max = b['m_dot_HS_fact_max'],
                            m_dot_CS_fact_min = b['m_dot_CS_fact_min'],
                            m_dot_CS_fact_max = b['m_dot_CS_fact_max'],
                            spliter_frac_min = 0,
                            spliter_frac_max = 1
                        )
        
                        # Point de départ du solveur : centre des bornes de ce point (arch, T)
                        P_high_init = 0.5 * (b['P_high_min'] + b['P_high_max'])
                        mdot_init = 0.5 * (b['m_dot_min'] + b['m_dot_max'])
                        mdot_HS_init = 0.5 * (b['m_dot_HS_fact_min'] + b['m_dot_HS_fact_max'])
        
                        if arch in ("Recomp", "Recomp_1_recup"):
                            Optimizer.set_it_var(P_high=P_high_init, mdot=mdot_init, mdot_HS=mdot_HS_init, spliter_frac=1)
                        else:
                            Optimizer.set_it_var(P_high=P_high_init, mdot=mdot_init, mdot_HS=mdot_HS_init)
        
                        Optimizer.set_obj(W_dot=W_dot_test)
                        Optimizer.set_CSource(T=15 + 273.15, P=5e5,  fluid='Water', m_dot=1000.0)
                        Optimizer.set_HSource(T=T, P=10e5, fluid='INCOMP::TVP1', m_dot=50.0)
                        Optimizer.set_RC()
                        Optimizer.opt_RC(n_jobs=n_cores - 1, n_particles=50, max_iter=100, patience=10)
        
                        eta = getattr(Optimizer, 'eta', None)
        
                        # Validation : eta doit être un nombre fini et physiquement sensé
                        if eta is not None and np.isfinite(eta) and 0 < eta < 1:
                            eta_run = eta
                            print(f"    -> eta = {eta:.4f}")
        
                            # Mise à jour du top BEST_KEEP pour cette condition (arch, T)
                            run_data = extract_run_data(Optimizer, eta)
                            update_best_runs(condition_key, run_data)
                            write_best_runs_csv()  # réécriture à chaque amélioration -> robuste à une interruption
                        else:
                            print(f"    ⚠️ Pas de solution valide (eta={eta})")
        
                    except Exception as e:
                        print(f"    ⚠️ Échec : {e}")
                        eta_run = None
        
                    run_etas.append(eta_run)
        
                results_all[arch].append(run_etas)
        
                # Meilleur résultat valide sur les N_RUNS essais
                valid_etas = [e for e in run_etas if e is not None]
                best_eta = max(valid_etas) if valid_etas else None
        
                if best_eta is not None:
                    print(f"  ✅ Meilleur eta sur {N_RUNS} essais : {best_eta:.4f}")
                else:
                    print(f"  ⚠️ Aucune solution valide sur {N_RUNS} essais")
        
                results[arch].append(best_eta)
        
        print(f"\n{len(best_runs)} conditions (architecture, T) traitées.")
        print(f"Top {BEST_KEEP} runs par condition enregistrés dans {CSV_PATH}")    
        
        # ---- Affichage récapitulatif ----
        print(f"\n{'='*60}")
        print("  RÉSUMÉ (meilleur eta par point)")
        print(f"{'='*60}")
        for arch in ARCH:
            print(f"{arch:20s} : {results[arch]}")
        
        # ---- Plot ----
        T_C = T_vec - 273.15
        
        plt.figure(figsize=(9, 6))
        
        for arch in ARCH:
            eta_arr = np.array([np.nan if v is None else v for v in results[arch]], dtype=float)
            plt.plot(T_C, eta_arr, marker='o', linewidth=2, label=arch)
        
        plt.title(f"Efficacité du cycle vs température (meilleur de {N_RUNS} essais)")
        plt.xlabel("Température source chaude [°C]")
        plt.ylabel("Efficacité η [-]")
        plt.grid(True)
        plt.legend(title="Architecture")
        plt.tight_layout()
        plt.show()
        
        