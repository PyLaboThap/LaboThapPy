# -*- coding: utf-8 -*-
"""
CO2 Transcritical Rankine Cycle Optimizer
- Parallel PSO evaluation via joblib (File 1 architecture)
- Rich cycle logic, penalty logging, warm start (File 2 logic)
"""

#%% Imports

from labothappy.machine.examples.ORC.fpi_TC_orc_example import REC_CO2_TC, basic_CO2_TC, Recomp_CO2_TC
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
    Evaluate one particle x = [P_high, m_dot, m_dot_HS_factor].
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
    m_dot_HS      = m_dot * m_dot_HS_fact

    # --- Build connectors ---
    HSource = MassConnector()
    HSource.set_properties(T=hs_props['T'], P=hs_props['P'],
                           fluid=hs_props['fluid'], m_dot=m_dot_HS)

    CSource = MassConnector()
    CSource.set_properties(T=cs_props['T'], P=cs_props['P'],
                           fluid=cs_props['fluid'], m_dot=cs_props.get('m_dot', 1000.0))

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
                DP_h_rec  = params.get('DP_h_rec',  1.0e5),
                DP_c_rec  = params.get('DP_c_rec',  2.0e5),
                DP_h_gh   = params.get('DP_h_gh',   0.5e5),
                DP_c_gh   = params.get('DP_c_gh',   2.0e5),
                DP_h_cond = params.get('DP_cond',   1.0e5),
                mute_print_flag=1,
            )
        elif arch == 'basic':
            RC = basic_CO2_TC(
                HSource, CSource,
                params['PP_gh'], params['PP_rec'],
                params['eta_pp'], params['eta_exp'],
                params['eta_gh'],
                params['PP_cd'], params['SC_cd'],
                P_low_guess, P_high, m_dot,
                mute_print_flag=1,
            )
        elif arch == "Recomp":
            spliter_frac = x[3]

            RC = Recomp_CO2_TC(
                HSource, CSource, 
                params['PP_gh'], params['PP_rec'], 
                params['eta_pp'], params['eta_exp'], params['eta_cp'], 
                params['eta_gh'], params['eta_rec'],
                params['PP_cd'], params['SC_cd'],
                P_low_guess, P_high, m_dot, spliter_frac,
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
        if arch == "Recomp":
            W_cp = RC.components['Compressor'].model.W.W_dot
        else:
            W_cp = 0
            
        W_exp  = RC.components['Expander'].model.W.W_dot
        W_pump = RC.components['Pump'].model.W.W_dot
        Q_gh   = RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = RC.components['GasHeater'].model.su_H.D
        m_HS_act   = RC.components['GasHeater'].model.su_H.m_dot
        W_pump_aux = params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * params.get('eta_pp', 0.8))

        W_dot_net = W_exp - W_pump - W_pump_aux - W_cp
        eta       = W_dot_net / Q_gh if Q_gh > 0 else 0.0
        
    except Exception:
        return 1000.0

    # --- Power target penalty ---
    W_obj   = obj.get('W_dot', 1e6)
    rel_err = abs((W_dot_net - W_obj) / W_obj)
    penalty = 500.0 * (rel_err ** 2) if rel_err > 1e-2 else 0.0

    return -eta + penalty


#%% Optimizer Class

class CO2RCOptimizer:

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
                DP_h_rec  = self.params.get('DP_h_rec',  1.0e5),
                DP_c_rec  = self.params.get('DP_c_rec',  2.0e5),
                DP_h_gh   = self.params.get('DP_h_gh',   0.5e5),
                DP_c_gh   = self.params.get('DP_c_gh',   2.0e5),
                DP_h_cond = self.params.get('DP_cond',   1.0e5),
                mute_print_flag=1,
            )
        elif arch == 'basic':
            self.RC = basic_CO2_TC(
                HSource, CSource,
                self.params['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'],
                self.params['eta_gh'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'],
                mute_print_flag=1,
            )
        elif arch == 'Recomp':
            self.RC = Recomp_CO2_TC(
                HSource, CSource, 
                self.params['PP_gh'], self.params['PP_rec'], 
                self.params['eta_pp'], self.params['eta_exp'], self.params['eta_cp'], 
                self.params['eta_gh'], self.params['eta_rec'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, self.it_var['P_high'], self.it_var['mdot'], self.it_var['spliter_frac'],
                mute_print_flag=1)
        else:
            raise ValueError("'RC_ARCH' parameter shall be either 'basic', 'REC', 'Recomp")

    def _log_penalty(self, reason):
        self.penalty_log[reason] = self.penalty_log.get(reason, 0) + 1

    # ------------------------------------------------------------------ final eval

    def _evaluate_final(self, best_pos):
        """
        Re-evaluates the best PSO position with full diagnostics and
        populates self.eta, self.W_dot_net, self.RC, self.it_var.
        Mirrors system_RC_parallel exactly so there is no disconnect.
        """
        if self.params['RC_ARCH'] == "Recomp":
            P_high, m_dot, m_dot_HS_fact, spliter_frac = best_pos
            m_dot_HS = m_dot * m_dot_HS_fact
            self.it_var['spliter_frac'] = spliter_frac
            self.it_var['P_high']  = P_high
            self.it_var['mdot']    = m_dot
            self.it_var['mdot_HS'] = m_dot_HS
            
        else:
            P_high, m_dot, m_dot_HS_fact = best_pos
            m_dot_HS = m_dot * m_dot_HS_fact

            self.it_var['P_high']  = P_high
            self.it_var['mdot']    = m_dot
            self.it_var['mdot_HS'] = m_dot_HS

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
        
        if self.params['RC_ARCH'] == "Recomp":
            W_cp = RC.components['Compressor'].model.W.W_dot
        else:
            W_cp = 0
            
        W_exp  = RC.components['Expander'].model.W.W_dot
        W_pump = RC.components['Pump'].model.W.W_dot
        Q_gh   = RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = RC.components['GasHeater'].model.su_H.D
        m_HS_act   = RC.components['GasHeater'].model.su_H.m_dot
        W_pump_aux = self.params.get('DP_h_gh', 0.5e5) * m_HS_act / \
                     (rho_HS * self.params.get('eta_pp', 0.8))

        self.W_dot_net = W_exp - W_pump - W_pump_aux - W_cp
        self.eta       = self.W_dot_net / Q_gh if Q_gh > 0 else 0.0

        self.Q_dot_waste = RC.components['GasHeater'].model.ex_H.m_dot * (
            RC.components['GasHeater'].model.ex_H.h
            - PropsSI('H', 'T', 273.15 + 15, 'P',
                      RC.components['GasHeater'].model.ex_H.p,
                      RC.components['GasHeater'].model.ex_H.fluid)
        )

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
        
        if self.params['RC_ARCH'] == "Recomp":
            # --- bounds ---
            lb = np.array([
                self.params['P_high_min'],
                self.params['m_dot_min'],
                self.params['m_dot_HS_fact_min'],
                self.params['spliter_frac_min'],
            ])
            ub = np.array([
                self.params['P_high_max'],
                self.params['m_dot_max'],
                self.params['m_dot_HS_fact_max'],
                self.params['spliter_frac_max'],
            ])
        else:   
            # --- bounds ---
            lb = np.array([
                self.params['P_high_min'],
                self.params['m_dot_min'],
                self.params['m_dot_HS_fact_min'],
            ])
            ub = np.array([
                self.params['P_high_max'],
                self.params['m_dot_max'],
                self.params['m_dot_HS_fact_max'],
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

    n_cores = multiprocessing.cpu_count()
    
    import matplotlib.pyplot as plt

    # ---- sweep ----
    T_vec = np.linspace(100, 350, 6) + 273.15
    # T_vec = np.array([200]) + 273.15

    eta_vec      = []
    P_high_vec   = []
    m_dot_vec    = []
    m_dot_HS_vec = []
    T_h_ex_vec   = []
    Q_dot_waste  = []

    W_dot_test = 1e6  # 1 MW target

    Optimizer = CO2RCOptimizer('CO2')

    for T in T_vec:

        Optimizer.set_parameters(
            RC_ARCH = 'Recomp',

            eta_pp  = 0.8,
            eta_gh  = 0.95,
            eta_rec = 0.9,
            eta_exp = 0.9,
            eta_cp = 0.8,

            PP_gh   = 5,
            PP_rec  = 0,
            PP_cd   = 5,
            SC_cd   = 0.1,

            DP_h_gh  = 50e3,
            DP_c_gh  = 50e3,
            DP_h_rec = 50e3,
            DP_c_rec = 50e3,
            DP_cond  = 50e3,

            P_high_min       = 80e5,
            P_high_max       = 200e5,
            m_dot_min        = 10.0,
            m_dot_max        = 100.0,
            m_dot_HS_fact_min = 0.1,
            m_dot_HS_fact_max = 2,
            spliter_frac_min = 0,
            spliter_frac_max = 1
        )
        
        if Optimizer.params['RC_ARCH'] == "Recomp":
            Optimizer.set_it_var(P_high=100e5, mdot=20.0, mdot_HS=15.0, spliter_frac = 1)
        else:
            Optimizer.set_it_var(P_high=100e5, mdot=20.0, mdot_HS=15.0)
            
        Optimizer.set_obj(W_dot=W_dot_test)

        Optimizer.set_CSource(T=15 + 273.15, P=5e5,  fluid='Water', m_dot=1000.0)
        Optimizer.set_HSource(T=T,           P=10e5, fluid='Water', m_dot=50.0)

        Optimizer.set_RC()
        Optimizer.opt_RC(n_jobs = n_cores - 1, n_particles=100, max_iter=50, patience = 10)

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
