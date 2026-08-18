#%% Imports

import time

from labothappy.machine.examples.CO2_Transcritical_Circuits.CO2_Transcritical_circuit import REC_CO2_TC, basic_CO2_TC, Recomp_CO2_TC_1_recup, Recomp_CO2_TC
from labothappy.connector.mass_connector import MassConnector

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from pyswarms.single import GlobalBestPSO
from tqdm import tqdm


#%% Optimizer class

class CO2RCOptimizer(object):

    def __init__(self, fluid):
        self.fluid = fluid
        self.RC = None

        self.inputs = {}
        self.params = {}
        self.it_var = {}
        self.obj    = {}

        self._HSource_props = {}
        self._CSource_props = {}

        self.eta       = None
        self.W_dot_net = None
        self.penalty_log = {}

#%% Set Methods

    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value

    def set_parameters(self, **parameters):
        for key, value in parameters.items():
            self.params[key] = value

    def set_it_var(self, **parameters):
        for key, value in parameters.items():
            self.it_var[key] = value

    def set_obj(self, **parameters):
        for key, value in parameters.items():
            self.obj[key] = value

    def set_HSource(self, T, P, fluid, m_dot):
        self._HSource_props = dict(T=T, P=P, fluid=fluid, m_dot=m_dot)

    def set_CSource(self, T, P, fluid, m_dot):
        self._CSource_props = dict(T=T, P=P, fluid=fluid, m_dot=m_dot)

#%% Cycle Set

    def set_RC(self):
        
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
        
        if self.params['RC_ARCH'] == "RECUP":        
            self.RC = REC_CO2_TC(
                HSource       = HSource,
                CSource       = CSource,
                Pinch_min_GH  = self.params['PP_gh'],
                Pinch_min_REC = self.params['PP_rec'],
                eta_pp        = self.params['eta_pp'],
                eta_exp       = self.params['eta_exp'],
                eta_gh        = self.params['eta_gh'],
                eta_rec       = self.params['eta_rec'],
                PP_cd         = self.params['PP_cd'],
                SC_cd         = self.params['SC_cd'],
                P_low         = P_low_guess,
                P_high        = self.it_var['P_high'],
                m_dot         = self.it_var['mdot'],
                DP_h_rec      = self.params.get('DP_h_rec',  1.0e5),
                DP_c_rec      = self.params.get('DP_c_rec',  2.0e5),
                DP_h_gh       = self.params.get('DP_h_gh',   0.5e5),
                DP_c_gh       = self.params.get('DP_c_gh',   2.0e5),
                DP_h_cond     = self.params.get('DP_h_cond', 1.0e5),
                DP_c_cond     = self.params.get('DP_c_cond', 0.5e5),
                mute_print_flag = 1,
            )
            
        elif self.params['RC_ARCH'] == "BASIC": 
            self.RC = basic_CO2_TC(
                HSource       = HSource,
                CSource       = CSource,
                Pinch_min_GH  = self.params['PP_gh'],
                eta_pp        = self.params['eta_pp'],
                eta_exp       = self.params['eta_exp'],
                eta_gh        = self.params['eta_gh'],
                PP_cd         = self.params['PP_cd'],
                SC_cd         = self.params['SC_cd'],
                P_low         = P_low_guess,
                P_high        = self.it_var['P_high'],
                m_dot         = self.it_var['mdot'],
                DP_h_gh       = self.params.get('DP_h_gh',   0.5e5),
                DP_c_gh       = self.params.get('DP_c_gh',   2.0e5),
                DP_h_cond     = self.params.get('DP_h_cond', 1.0e5),
                DP_c_cond     = self.params.get('DP_c_cond', 0.5e5),
                mute_print_flag = 1,
                )

        elif self.params['RC_ARCH'] == "RECOMP_1_REC": 
            self.RC = Recomp_CO2_TC_1_recup(
                HSource       = HSource,
                CSource       = CSource,
                Pinch_min_GH  = self.params['PP_gh'],
                Pinch_min_REC = self.params['PP_rec'],
                eta_pp        = self.params['eta_pp'],
                eta_cp        = self.params['eta_cp'],
                eta_exp       = self.params['eta_exp'],
                eta_gh        = self.params['eta_gh'],
                eta_rec       = self.params['eta_rec'],
                PP_cd         = self.params['PP_cd'],
                SC_cd         = self.params['SC_cd'],
                P_low         = P_low_guess,
                P_high        = self.it_var['P_high'],
                m_dot         = self.it_var['mdot'],
                DP_h_rec      = self.params.get('DP_h_rec',  1.0e5),
                DP_c_rec      = self.params.get('DP_c_rec',  2.0e5),
                DP_h_gh       = self.params.get('DP_h_gh',   0.5e5),
                DP_c_gh       = self.params.get('DP_c_gh',   2.0e5),
                DP_h_cond     = self.params.get('DP_h_cond', 1.0e5),
                DP_c_cond     = self.params.get('DP_c_cond', 0.5e5),
                spliter_frac  = self.it_var['spliter_frac'],
                mute_print_flag = 1,
            )
            
            
        elif self.params['RC_ARCH'] == "RECOMP": 
            self.RC = Recomp_CO2_TC(
                HSource       = HSource,
                CSource       = CSource,
                Pinch_min_GH  = self.params['PP_gh'],
                Pinch_min_REC = self.params['PP_rec'],
                eta_pp        = self.params['eta_pp'],
                eta_cp        = self.params['eta_cp'],
                eta_exp       = self.params['eta_exp'],
                eta_gh        = self.params['eta_gh'],
                eta_rec       = self.params['eta_rec'],
                PP_cd         = self.params['PP_cd'],
                SC_cd         = self.params['SC_cd'],
                P_low         = P_low_guess,
                P_high        = self.it_var['P_high'],
                m_dot         = self.it_var['mdot'],
                DP_h_rec      = self.params.get('DP_h_rec',  1.0e5),
                DP_c_rec      = self.params.get('DP_c_rec',  2.0e5),
                DP_h_gh       = self.params.get('DP_h_gh',   0.5e5),
                DP_c_gh       = self.params.get('DP_c_gh',   2.0e5),
                DP_h_cond     = self.params.get('DP_h_cond', 1.0e5),
                DP_c_cond     = self.params.get('DP_c_cond', 0.5e5),
                spliter_frac  = self.it_var['spliter_frac'],
                mute_print_flag = 1,
            )

    def log_penalty(self, reason):
        if reason not in self.penalty_log:
            self.penalty_log[reason] = 0
        self.penalty_log[reason] += 1

#%% Cycle Optimization

    def system_RC(self, x):
        
        if "RECOMP" in self.params["RC_ARCH"]:
            self.it_var['P_high']  = x[0]
            self.it_var['mdot']    = x[1]
            self.it_var['mdot_HS'] = x[1] * x[2]
            self.it_var['spliter_frac'] = x[3]
            
        else:
            self.it_var['P_high']  = x[0]
            self.it_var['mdot']    = x[1]
            self.it_var['mdot_HS'] = x[1] * x[2]

        self.set_RC()
        CO2_RC = self.RC

        try:
            CO2_RC.solve(method="wegstein")
        except Exception as e:
            error_msg = str(e)
            low = error_msg.lower()
            if "pinch" in low and "cond" in low:
                reason = "Condenser pinch failure"
            elif "coolprop" in low:
                reason = "CoolProp property error"
            else:
                reason = f"Solver exception: {error_msg[:80]}"
            self.log_penalty(reason)
            self.eta = self.W_dot_net = 0.0
            return 100.0

        if not CO2_RC.converged:
            msg = getattr(CO2_RC, 'error_msg', 'Not converged')
            self.log_penalty(f"Non-convergence: {msg}")
            self.eta = self.W_dot_net = 0.0
            return 10.0

        try:
            T_exp_ex  = CO2_RC.components['Expander'].model.ex.T
            P_exp_ex  = CO2_RC.components['Expander'].model.ex.p
            T_sat_exp = PropsSI('T', 'P', P_exp_ex, 'Q', 1, 'CO2')
            SH_exp    = T_exp_ex - T_sat_exp
        except Exception:
            SH_exp = 50.0

        if SH_exp < 0:
            self.log_penalty(f"Drops in the expansion (SH = {SH_exp:.1f} K < 0)")
            self.eta = self.W_dot_net = 0.0
            return 100.0

        W_exp  = CO2_RC.components['Expander'].model.W.W_dot
        W_pump = CO2_RC.components['Pump'].model.W.W_dot
        Q_gh   = CO2_RC.components['GasHeater'].model.Q.Q_dot

        rho_HS     = CO2_RC.components['GasHeater'].model.su_H.D
        m_HS       = CO2_RC.components['GasHeater'].model.su_H.m_dot
        W_pump_aux = self.params.get('DP_h_gh', 0.5e5) * m_HS / (rho_HS * self.params.get('eta_pp', 0.8))
        
        if "RECOMP" in self.params["RC_ARCH"]:
            W_cp = CO2_RC.components['Compressor'].model.W.W_dot
            self.W_dot_net = W_exp - W_pump - W_pump_aux - W_cp
        else:
            self.W_dot_net = W_exp - W_pump - W_pump_aux
            
        self.eta       = self.W_dot_net / Q_gh if Q_gh > 0 else 0.0

        W_obj   = self.obj.get('W_dot', 10e6)
        rel_err = abs((self.W_dot_net - W_obj) / W_obj)

        if rel_err > 1e-2:
            self.log_penalty("Power out of bounds (> 1% deviation)")
            penalty = 50.0 * (rel_err ** 2)
        else:
            penalty = 0.0
            self.log_penalty("Valid cycle evaluated")

        return -self.eta + penalty

    def opt_RC(self, n_particles=50, max_iter=30, patience=None,
               init_pos=None, warm_spread=0.05, warm_fraction=0.5):
        """
        Launches the PSO with optional warm start.

        init_pos : array-like or None
          - 1D length-n_dim → central seed; warm_fraction of particles around it.
          - 2D (n_particles, n_dim) → used verbatim.
          - None → fully random.
        """
        
        if patience is None:
            patience = max_iter/5
        
        if "RECOMP" in self.params["RC_ARCH"]:
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

        pso_init_pos = None
        if init_pos is not None:
            seed = np.asarray(init_pos, dtype=float)

            if seed.ndim == 1:
                seed = np.clip(seed, lb, ub)
                n_warm = max(1, int(round(warm_fraction * n_particles)))
                n_rand = n_particles - n_warm

                noise = np.random.uniform(-warm_spread, warm_spread,
                                          size=(n_warm, len(lb)))
                warm = seed[None, :] * (1.0 + noise)
                warm = np.clip(warm, lb, ub)

                rand = np.random.uniform(lb, ub, size=(n_rand, len(lb)))
                pso_init_pos = np.vstack([warm, rand])
                print(f"  → Warm start: {n_warm}/{n_particles} particles around "
                      f"P={seed[0]/1e5:.1f} bar, ṁ={seed[1]:.1f}, f={seed[2]:.3f}")
            else:
                pso_init_pos = np.clip(seed, lb, ub)

        def objective_wrapper(X):
            return np.array([self.system_RC(xi) for xi in X])

        optimizer = GlobalBestPSO(
            n_particles = n_particles,
            dimensions  = len(ub),
            options     = {'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds      = (lb, ub),
            init_pos    = pso_init_pos,
        )

        tol        = 1e-3
        no_improve = 0
        best_cost  = np.inf

        for i in tqdm(range(max_iter), desc="PSO Optimization", ncols=80):
            optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current = optimizer.swarm.best_cost

            if current < best_cost - tol:
                best_cost  = current
                no_improve = 0
            else:
                no_improve += 1

            print(f"  [{i+1:02d}] Best cost = {best_cost:.6f}")

            if no_improve >= patience:
                print("  → Early stopping due to stagnation.")
                break

        self.system_RC(optimizer.swarm.best_pos)
        bp = optimizer.swarm.best_pos

        print("\n" + "="*55)
        print(f"  OPTIMAL RESULT : {self._HSource_props['T']-273.15} °C")
        print("="*55)
        print(f"  P_high            : {bp[0]/1e5:.2f}  bar")
        print(f"  m_dot (CO2)       : {bp[1]:.4f}  kg/s")
        print(f"  m_dot_HS_factor   : {bp[2]:.4f}  [-]")
        print(f"  m_dot_HS          : {bp[1]*bp[2]:.4f}  kg/s")
        
        if "RECOMP" in self.params["RC_ARCH"]:
            print(f"  spliter_frac      : {bp[3]:.4f}  [-]")

        print(f"  W_net             : {self.W_dot_net/1e3:.3f}  kW")
        print(f"  Thermal η         : {self.eta*100:.3f}  %")
        print("="*55)

        print("\n" + "="*55)
        print("  PENALTY SUMMARY  ")
        print("="*55)
        
        total_evals = sum(self.penalty_log.values())
        if total_evals > 0:
            sorted_log = sorted(self.penalty_log.items(),
                                key=lambda item: item[1], reverse=True)
            for reason, count in sorted_log:
                percentage = (count / total_evals) * 100
                print(f"  [{count:4d} | {percentage:5.1f}%] : {reason}")
        else:
            print("  No evaluations logged.")
            
        print("="*55)
        
        W_obj = self.obj.get('W_dot', 10e6)
        rel_err_final = abs((self.W_dot_net - W_obj) / W_obj)
        if rel_err_final > 0.05:
            print(f"  ⚠️  WARNING: (W_net error = {rel_err_final*100:.1f}%)")
            self.eta = np.nan

        self.penalty_log = {}

        return

#%% Parameter Setting

Optimizer = CO2RCOptimizer('CO2')

T = 100

if T == 95:
    Optimizer.set_parameters(
        eta_pp = 0.95, # [-]
        eta_gh = 0.95, # [-]
        eta_rec = 0.9, # [-]
        eta_exp = 0.9, # [-]
        PP_cd = 5, # [K]
        PP_gh = 5, # [K]
        PP_rec = 0, # [K]
        SC_cd = 0.1, # [K]
        
        P_high_min = 80*1e5, # [Pa]
        P_high_max = 200*1e5, # [Pa]
        
        m_dot_min = 70, # [kg/s]
        m_dot_max = 90, # [kg/s]
        
        m_dot_HS_fact_min = 0.5, # [-]
        m_dot_HS_fact_max = 1, # [-]
        )
    
    Optimizer.set_it_var(
        P_high = 100*1e5, # [Pa]
        mdot_HS = 70, # [kg/s]
        mdot = 70 # [kg/s]
        )
    
    Optimizer.set_obj(
        W_dot = 1e6, # [W]
        )
    
    Optimizer.CSource.set_properties(
        T = 15 + 273.15, # [K]
        P = 5*1e5, # [Pa]
        fluid = 'Water' # []
        )
    
    Optimizer.HSource.set_properties(
        T = 95 + 273.15, # [K]
        P = 5*1e5, # [Pa]
        fluid = 'Water', # []
        )
    
    Optimizer.set_RC()
    Optimizer.opt_RC()
    
else:
    
    # T_vec = np.linspace(100, 150, 6) + 273.15
    # T_vec = np.linspace(100, 300, 6) + 273.15
    T_vec = np.linspace(200, 300, 2) + 273.15
    eta_vec = []
    P_high_vec = []
    m_dot_vec = []
    m_dot_HS_vec = []
    T_h_ex_vec = []
    
    Opt_dict = {}
    
    for T in T_vec:
        Optimizer = CO2RCOptimizer('CO2')
        
        Optimizer.set_parameters(
            RC_ARCH = "RECOMP_1_REC", # "RECUP", "BASIC", "RECOMP_1_REC", "RECOMP"
            
            eta_pp = 0.8, # [-]
            eta_cp = 0.8, # [-] : Used if recompression
            eta_gh = 0.95, # [-]
            eta_rec = 0.9, # [-]
            eta_exp = 0.9, # [-]
            PP_cd = 5, # [K]
            PP_gh = 5, # [K]
            PP_rec = 0, # [K]
            SC_cd = 0.1, # [K]
            
            P_high_min = 80*1e5, # [Pa]
            P_high_max = 200*1e5, # [Pa]
            
            m_dot_HS_fact_min = 0.4, # [-]
            m_dot_HS_fact_max = 1.3, # [-]
            
            m_dot_min = 5, # [kg/s]
            m_dot_max = 100, # [kg/s]
            
            spliter_frac_min = 0.1, # [-] : Used if recompression
            spliter_frac_max = 1, # [-] : Used if recompression
            )
        
        Optimizer.set_it_var(
            P_high = 140*1e5, # [Pa]
            mdot_HS = 20, # [kg/s]
            mdot = 17, # [kg/s]
            spliter_frac = 0.5, # [-] : Used if recompression
            )
        
        Optimizer.set_obj(
            W_dot = 1e6, # [W]
            )
    
        Optimizer.set_CSource(T=10 + 273.15, P=5e5, fluid='Water', m_dot=5000.0)
        Optimizer.set_HSource(T=T, P=20e5, fluid='Water', m_dot=192.0)
    
        Optimizer.set_RC()
        Optimizer.opt_RC(n_particles=50)

        Opt_dict[str(int(T-273.15))] = Optimizer

        eta_vec.append(Optimizer.eta)
        P_high_vec.append(Optimizer.it_var['P_high'])
        m_dot_vec.append(Optimizer.it_var['mdot'])
        m_dot_HS_vec.append(Optimizer.it_var['mdot_HS'])
        T_h_ex_vec.append(Optimizer.RC.components['GasHeater'].model.ex_H.T)
        
    import matplotlib.pyplot as plt
    
    # Plot 1: Efficiency vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, eta_vec, linewidth=2)
    plt.title("Efficiency vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Efficiency")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Plot 2: P_high vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, P_high_vec, linewidth=2, label='P_high')
    plt.title("High Pressure vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("P_high (units?)")  # Replace with actual units
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot 3: Mass Flow Rate (mdot) vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, m_dot_vec, linewidth=2, label='mdot')
    plt.title("Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Mass Flow Rate (kg/s)")  # Adjust unit if needed
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    # Plot 4: Heat Source Mass Flow Rate (mdot_HS) vs Temperature
    plt.figure(figsize=(8, 5))
    plt.plot(T_vec - 273.15, m_dot_HS_vec, linewidth=2, label='mdot_HS')
    plt.title("Heat Source Mass Flow Rate vs Temperature")
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Heat Source Mass Flow Rate (kg/s)")  # Adjust unit if needed
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
