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

from sizing.turbomachinery.turbine.axial.design_1D.mean_line_axial_turbine_loss_model_design_aungier import AxialTurbineMeanLineDesign
from sizing.heat_exchanger.shell_and_tube.shell_and_tube_sizing_PSO_parallel import ShellAndTubeSizingOpt
from sizing.heat_exchanger.PCHE.PCHE_PSO import PCHESizingOpt

#%% Define parallel system evaluation outside the class

def system_RC_parallel(x, input_data):
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
        RC = REC_CO2_TC(HSource, CSource.T, params['PP_gh'], params['PP_rec'], params['eta_pp'],
                        params['eta_exp'], params['eta_gh'], params['eta_rec'],
                        params['PP_cd'], params['SC_cd'], P_low_guess, x[0], x[1], mute_print_flag=1)

    elif input_data['RC_ARCH'] == 'basic':
        RC = basic_CO2_TC(HSource, CSource.T, params['PP_gh'], params['PP_rec'], params['eta_pp'],
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

        W_dot_net = RC.components['Expander'].model.W_exp.W_dot*0.95 - RC.components['Pump'].model.W_pp.W_dot/0.95 - pp_power/0.95
        eta = W_dot_net / RC.components['GasHeater'].model.Q_dot.Q_dot

        Th_out = RC.components['GasHeater'].model.ex_H.T
        Th_out_obj = 15 + 273.15

        penalty = 0

        if abs((W_dot_net - obj['W_dot'])/obj['W_dot']) > 1e-2:
            penalty = abs((W_dot_net - obj['W_dot'])/obj['W_dot']) * 10

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
            n_particles=40,
            dimensions=3,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        best_cost = np.inf
        no_improve_counter = 0
        patience = 5
        tol = 1e-3
        max_iter = 30
    
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
                self.HSource, self.CSource.T,
                self.params['PP_gh'], self.params['PP_rec'],
                self.params['eta_pp'], self.params['eta_exp'],
                self.params['eta_gh'], self.params['eta_rec'],
                self.params['PP_cd'], self.params['SC_cd'],
                P_low_guess, best_P_high, best_m_dot,
                mute_print_flag=1
            )
    
            self.RC.solve()
    
            # Turbine Sizing
            
            Turb_model = self.RC.components['Expander'].model
            Turb_sizing = self.RC.components['Expander'].sizing = AxialTurbineMeanLineDesign(self.fluid)
            
            Turb_sizing.set_inputs(
                mdot = Turb_model.su.m_dot, # kg/s
                W_dot = Turb_model.W_exp.W_dot, # W : 
                p0_su = Turb_model.su.p, # Pa
                T0_su = Turb_model.su.T, # K
                p_ex = Turb_model.ex.p, # Pa
                )
            
            Turb_sizing.set_parameters(
                Zweifel = 0.8, # [-]
                AR_min = 0.8, # [-]
                r_hub_tip_max = 0.95, # [-]
                r_hub_tip_min = 0.6, # [-]
                Re_bounds = [3*1e6,8*1e6], # [-]
                psi_bounds = [1,1.9], # [-]
                phi_bounds = [0.5,0.8], # [-]
                R_bounds = [0.45,0.55], # [-]
                M_1st_bounds = [0.3, 0.5], # [-]
                r_m_bounds = [0.15, 0.6], # [m]
                # Omega_choices = [500,750,1000,1500,3000], # [RPM] : [500,750,1000,1500,3000]
                damping = 0.3, # [-]
                p_rel_tol = 0.05, # [-]
                delta_tip = 0.4*1e-3, # [m] : tip clearance
                N_lw = 0, # [-] : Number of lashing wires
                D_lw = 0, # [m] : Diameter of lashing wires
                e_blade = 0.002*1e-3, # [m] : blade roughness
                t_TE_o = 0.05, # [-] : trailing edge to throat opening ratio
                t_TE_min = 5*1e-4, # [m]
                )
            
            Turb_sizing.design_parallel(n_jobs=-1)
            
            # ---------------------------------------------------------------------------------------------------------------------------------
            # Recuperator Sizing

            REC_model = self.RC.components['Recuperator'].model
            REC_sizing = self.RC.components['Recuperator'].sizing = PCHESizingOpt()
                        
            REC_sizing.set_inputs(
                # First fluid
                fluid_H = REC_model.su_H.fluid,
                T_su_H = REC_model.su_H.T, # K
                P_su_H = REC_model.su_H.p, # Pa
                m_dot_H = REC_model.su_H.m_dot, # kg/s
        
                # Second fluid
                fluid_C = REC_model.su_C.fluid,
                T_su_C = REC_model.su_C.T, # K
                P_su_C = REC_model.su_C.p, # Pa
                m_dot_C = REC_model.su_C.m_dot, # kg/s  # Make sure to include fluid information
                )
        
            REC_sizing.set_parameters(
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
            
            REC_sizing.set_corr(H_Corr, C_Corr, H_DP, C_DP)
        
            L_x_bounds = np.array([0.1, 1.5])*2
            L_y_bounds = np.array([0.1, 2.3])*2
            L_z_bounds = np.array([0.1, 0.6])*2
        
            REC_sizing.set_bounds(
                alpha = [10,40], # [°]
                D_c = [1*1e-3, 3*1e-3], # [m]
                L_x = L_x_bounds, # [m] : 1.5 limit fixed by Heatric (PCHE manufacturer) : Fluid direction
                L_y = L_y_bounds, # [m] : 2.3 limit for shipping requirements : Vertical direction
                L_z = L_z_bounds, # [m] : 0.6 limit fixed by Heatric (PCHE manufacturer) : Width
                )
            
            Q_dot_cstr = REC_model.Q_dot.Q_dot
            DP_c_cstr = 100*1e3
            DP_h_cstr = 200*1e3
            
            REC_sizing.set_constraints(Q_dot = Q_dot_cstr, DP_h = DP_h_cstr, DP_c = DP_c_cstr)
            REC_sizing.design_parallel()
            
            # ---------------------------------------------------------------------------------------------------------------------------------
            
            # Fan and pump calculations
            DP = 50e3
            rho = self.RC.components['GasHeater'].model.su_H.D
            mdot = self.RC.components['GasHeater'].model.su_H.m_dot
            h_ex = self.RC.components['GasHeater'].model.ex_H.h
    
            h_ex_req = PropsSI('H', 'T', 15 + 273.15, 'P', 101325, 'Water')
            Q_dot_req = mdot * (h_ex - h_ex_req)
            W_dot_fan = 0.2 * Q_dot_req
            W_dot_pp = DP * mdot / (rho * 0.8)
    
            W_net = self.RC.components['Expander'].model.W_exp.W_dot*0.95 - \
                    self.RC.components['Pump'].model.W_pp.W_dot/0.95 - W_dot_pp/0.95
    
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
    
    power = 10*1e6
    
    m_dot_min = power*20/1e6
    m_dot_max = power*70/1e6
    
    import matplotlib.pyplot as plt
    
    T_vec = 130+273.15
    
    # Create optimizer instance
    Optimizer = CO2RCOptimizer('CO2')
    
    # Sweep parameters
    m_dot_HS_fact_min = 0.5
    m_dot_HS_fact_max = 1
    P_high_min = 100 * 1e5
    P_high_max = 180 * 1e5
    
    # Set model parameters
    Optimizer.set_parameters(
        RC_ARCH= 'REC', # 'REC'
        eta_pp=0.8,
        eta_gh=0.95,
        eta_rec=0.8,
        eta_exp=0.9,
        PP_cd=5,
        PP_gh=5,
        PP_rec=0,
        SC_cd=0.1,
        
        DP_h_rec = 100*1e3,
        DP_c_rec = 100*1e3,
        
        DP_h_gh = 100*1e3,
        DP_c_gh = 100*1e3,
        
        DP_h_ev = 100*1e3,
        DP_c_ev = 100*1e3,
        
        P_high_bounds=[P_high_min, P_high_max],
        m_dot_HS_fact_bounds=[m_dot_HS_fact_min,m_dot_HS_fact_max],
        m_dot_bounds=[m_dot_min, m_dot_max],
    )

    # Initial guess
    Optimizer.set_it_var(
        P_high=100e5,
        mdot=17,
        mdot_HS=20
    )

    # Objective
    Optimizer.set_obj(W_dot=10*1e6)

    # Source definitions
    Optimizer.CSource.set_properties(
        T=15 + 273.15,
        P=5e5,
        fluid='Water'
    )

    Optimizer.HSource.set_properties(
        T=T_vec,
        P=10e5,
        fluid='Water'
    )

    # Prepare model and optimize
    Optimizer.opt_RC()
