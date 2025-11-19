# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:12:02 2025

@author: Basile
"""

#%% Imports

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
from sizing.turbomachinery.pump.radial.radial_pump_0D_design import RadialPumpODDesign

import warnings
warnings.filterwarnings('ignore')

#%% Define parallel system evaluation outside the class

def TCO2_rec_comp_sizing(RC):
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Recuperator Sizing
    
    try:

        REC_model = RC.components['Recuperator'].model
        REC_sizing = RC.components['Recuperator'].sizing = PCHESizingOpt()
                    
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
        DP_c_cstr = REC_model.DP_c
        DP_h_cstr = REC_model.DP_h
        
        REC_sizing.set_constraints(Q_dot = Q_dot_cstr, DP_h = DP_h_cstr, DP_c = DP_c_cstr)
        REC_sizing.design_parallel()

    except Exception as e:
        print(f"⚠️ Failed to design Recuperator: {e}")
        return RC, 0
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    # GasHeater Sizing
    try:     
        GH_model = RC.components['GasHeater'].model
        GH_sizing = RC.components['GasHeater'].sizing = ShellAndTubeSizingOpt()
                                 
        GH_sizing.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'])
    
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            # 'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                            #     29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Shell_ID_inch' : [25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2], # [1,2,4,6,8,10]
                            'tube_layout' : [0,45,60]}
    
        GH_sizing.set_choice_vectors(choice_vectors)
    
        GH_sizing.set_max_cycle_prop(T_max_cycle = RC.sources['GH_Water'].properties.T, p_max_cycle = RC.components['Pump'].model.ex.p)
        
        GH_sizing.set_inputs(
            # First fluid
            fluid_H = GH_model.su_H.fluid,
            T_su_H = GH_model.su_H.T, # K
            P_su_H = GH_model.su_H.p, # Pa
            m_dot_H = GH_model.su_H.m_dot, # kg/s
    
            # Second fluid
            fluid_C = GH_model.su_C.fluid,
            T_su_C = GH_model.su_C.T, # K
            P_su_C = GH_model.su_C.p, # Pa
            m_dot_C = GH_model.su_C.m_dot, # kg/s  # Make sure to include fluid information
            )

        GH_sizing.set_parameters(
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
                                n_disc = 50
                              )
    
        H_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}
        C_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Flow_boiling"}
        
        H_DP = "Shell_Kern_DP"
        C_DP = "Gnielinski_DP"
        
        GH_sizing.set_corr(H_Corr, C_Corr, H_DP, C_DP)

        bounds = {
                    "L_shell" : [1,15], # 10],
                    "D_o_inch" : [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
                    "Shell_ID_inch" : [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]],
                    "Tube_pass" : [choice_vectors['Tube_pass'][0], choice_vectors['Tube_pass'][-1]],
                    "tube_layout" : [choice_vectors['tube_layout'][0], choice_vectors['tube_layout'][-1]],
                    "Baffle_cut" : [15, 45]
                    }

        GH_sizing.set_bounds(bounds)
        GH_sizing.set_constraints(Q_dot = GH_model.Q_dot.Q_dot, DP_h = max(GH_model.DP_h, 1e3), DP_c = max(GH_model.DP_c, 1e3))

        global_best_position, global_best_score, best_particle = GH_sizing.opt_size()
        
    except Exception as e:
        print(f"⚠️ Failed to design GasHeater: {e}")
        return RC, 0
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Condenser Sizing
    try: 
    
        CD_model = RC.components['Condenser'].model
        CD_sizing = RC.components['Condenser'].sizing = ShellAndTubeSizingOpt()
                                 
        CD_sizing.set_opt_vars(['D_o_inch', 'L_shell', 'Shell_ID_inch', 'Central_spac', 'Tube_pass', 'tube_layout', 'Baffle_cut'])
    
        choice_vectors = {
                            'D_o_inch' : [0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5],
                            # 'Shell_ID_inch' : [8, 10, 12, 13.25, 15.25, 17.25, 19.25, 21.25, 23.25, 25, 27,        
                            #     29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Shell_ID_inch' : [25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48, 54, 60, 66, 72, 78, 84, 90, 96, 108, 120],
                            'Tube_pass' : [2], # [1,2,4,6,8,10]
                            'tube_layout' : [0,45,60]}
    
        CD_sizing.set_choice_vectors(choice_vectors)
    
        CD_sizing.set_max_cycle_prop(T_max_cycle = RC.sources['CD_Water'].properties.T, p_max_cycle = RC.components['Pump'].model.ex.p)
        
        CD_sizing.set_inputs(
            # First fluid
            fluid_H = CD_model.su_H.fluid,
            T_su_H = CD_model.su_H.T, # K
            P_su_H = CD_model.su_H.p, # Pa
            m_dot_H = CD_model.su_H.m_dot, # kg/s
    
            # Second fluid
            fluid_C = CD_model.su_C.fluid,
            T_su_C = CD_model.su_C.T, # K
            P_su_C = CD_model.su_C.p, # Pa
            m_dot_C = CD_model.su_C.m_dot, # kg/s  # Make sure to include fluid information
            )

        CD_sizing.set_parameters(
                                n_series = 1, # [-]
                                # OPTI -> Oui (regarder le papier pour déterminer ça)
    
                                foul_t = 0.000176, # (m^2 * K/W)
                                foul_s =  0.000176, # (m^2 * K/W)
                                tube_cond = 50, # W/(m*K)
                                Overdesign = 0,
                                
                                Shell_Side = 'C',
    
                                Flow_Type = 'Shell&Tube',
                                H_DP_ON = True,
                                C_DP_ON = True,
                                n_disc = 50
                              )
    
        H_Corr = {"SC" : "Gnielinski", "1P" : "Gnielinski", "2P" : "Thome_Condensation"}
        C_Corr = {"SC" : "Shell_Kern_HTC", "1P" : "Shell_Kern_HTC", "2P" : "Shell_Kern_HTC"}

        H_DP = "Choi_DP"
        C_DP = "Shell_Kern_DP"
        
        CD_sizing.set_corr(H_Corr, C_Corr, H_DP, C_DP)

        bounds = {
                    "L_shell" : [1,15], # 10],
                    "D_o_inch" : [choice_vectors['D_o_inch'][0], choice_vectors['D_o_inch'][-1]],
                    "Shell_ID_inch" : [choice_vectors['Shell_ID_inch'][0], choice_vectors['Shell_ID_inch'][-1]],
                    "Tube_pass" : [choice_vectors['Tube_pass'][0], choice_vectors['Tube_pass'][-1]],
                    "tube_layout" : [choice_vectors['tube_layout'][0], choice_vectors['tube_layout'][-1]],
                    "Baffle_cut" : [15, 45]
                    }

        CD_sizing.set_bounds(bounds)
        CD_sizing.set_constraints(Q_dot = CD_model.Q_dot.Q_dot, DP_h = max(CD_model.DP_h, 1e3), DP_c = max(CD_model.DP_c, 1e3))

        global_best_position, global_best_score, best_particle = CD_sizing.opt_size()
        
    except Exception as e:
        print(f"⚠️ Failed to design Condenser: {e}")
        return RC, 0
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Pump Sizing
    try: 
    
        Pump_model = RC.components['Pump'].model
        Pump_sizing = RC.components['Pump'].sizing = RadialPumpODDesign(RC.fluid)
                              
        Pump_sizing.set_inputs(
            P_su = Pump_model.su.p, # Pa
            P_ex = Pump_model.ex.p, # Pa
            T_su = Pump_model.su.T, # K
                    
            H1 = 0, # m
            H2 = 0, # m
            
            v1 = 0, # m/s
            v2 = 0, # m/s
            
            m_dot = Pump_model.su.m_dot, # kg/s
            )

        Pump_sizing.set_parameters(
            Omega_choices = np.array([750, 1000, 1500, 3000]),
            n_parallel_choices = np.array([1, 2, 3, 4, 5, 6, 7, 8])
            )
        
        Pump_sizing.design()

    except Exception as e:
        print(f"⚠️ Failed to design Pump: {e}")
        return RC, 0
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Turbine Sizing

    try:                
        Turb_model = RC.components['Expander'].model
        Turb_sizing = RC.components['Expander'].sizing = AxialTurbineMeanLineDesign(RC.fluid)
        
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
        
        Turb_sizing.design_parallel(n_jobs=-1, n_particles = 30, max_iter=50)
        
    except Exception as e:
        print(f"⚠️ Failed to design the Turbine: {e}")
        return RC, 0
    
    RC.CAPEX = {
        "Pump" : np.round(Pump_sizing.CAPEX['Total']),
        "GasHeater" : np.round(GH_sizing.CAPEX['Total']),
        "Recuperator" : np.round(REC_sizing.CAPEX['Total']),
        "Expander" : np.round(Turb_sizing.CAPEX['Total']),            
        "Condenser" : np.round(CD_sizing.CAPEX['Total']),    
        }
    
    RC.CAPEX["Total"] = RC.CAPEX["Pump"] + RC.CAPEX["GasHeater"] + RC.CAPEX["Recuperator"] + RC.CAPEX["Expander"] + RC.CAPEX["Condenser"]

    return RC, 1

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
        RC = basic_CO2_TC(HSource, CSource.T, params['PP_gh'], params['PP_rec'], params['eta_pp'],
                          params['eta_exp'], params['eta_gh'], params['PP_cd'], params['SC_cd'],
                          P_low_guess, x[0], x[1], mute_print_flag=1)

    try:
        RC.solve()
    
        if not RC.converged:
            # cost, penalty, eta
            return 1000.0, np.inf, np.nan
    
        DP = 50e3
        rho = RC.components['GasHeater'].model.su_H.D
        mdot = RC.components['GasHeater'].model.su_H.m_dot + RC.components['GasHeater'].model.su_C.m_dot
        
        eta_pp = 0.8
        pp_power = DP * mdot / (rho * eta_pp)
    
        W_dot_net = (RC.components['Expander'].model.W_exp.W_dot * 0.95
                     - RC.components['Pump'].model.W_pp.W_dot / 0.95
                     - pp_power / 0.95)
        eta = W_dot_net / RC.components['GasHeater'].model.Q_dot.Q_dot
    
        penalty = 0.0
    
        if abs((W_dot_net - obj['W_dot'])/obj['W_dot']) > 2e-2:
            penalty += abs((W_dot_net - obj['W_dot'])/obj['W_dot']) * 10
    
        if abs((obj['eta'] - eta)/obj['eta']) > 2e-2:
            penalty += abs((obj['eta'] - eta)/obj['eta']) * 10
            
        # objective = RC.components['GasHeater'].model.Q_dot.Q_dot
        RC.eta = eta
        RC.W_dot_net = W_dot_net
        
        Q_cond = RC.components['Condenser'].model.Q_dot.Q_dot
        RC.components['Condenser'].model.equivalent_effectiveness()
        eta_cond = RC.components['Condenser'].model.epsilon
        
        Q_rec = RC.components['Recuperator'].model.Q_dot.Q_dot
        eta_rec = RC.components['Recuperator'].model.epsilon
    
        Q_gh = RC.components['GasHeater'].model.Q_dot.Q_dot
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

    def evaluate_systems(self):
        
        # invalid_RCs = []
        RC_scores = []
        delta_dicts = []

        for RC in self.potential_RC:
            delta_dicts.append({})
            
            eta_exp = RC.components['Expander'].sizing.eta_is
            eta_pp = RC.components['Pump'].sizing.eta_is
    
            DP_h_gh = RC.components['GasHeater'].sizing.best_particle.DP_h
            DP_c_gh = RC.components['GasHeater'].sizing.best_particle.DP_c
            
            DP_h_cond = RC.components['Condenser'].sizing.best_particle.DP_h
            DP_c_cond = RC.components['Condenser'].sizing.best_particle.DP_c
            
            DP_h_rec = RC.components['Recuperator'].sizing.HX.DP_h
            DP_c_rec = RC.components['Recuperator'].sizing.HX.DP_c
            
            delta_dicts[-1]['eta_exp'] = delta_exp = ((eta_exp - self.params['eta_exp'])/self.params['eta_exp'])**2
            delta_dicts[-1]['eta_pp'] = delta_pp = ((eta_pp - self.params['eta_pp'])/self.params['eta_pp'])**2
            
            delta_dicts[-1]['DP_h_gh'] = delta_h_gh = ((np.max([DP_h_gh, self.params['DP_h_gh']]) - self.params['DP_h_gh'])/self.params['DP_h_gh'])**2
            delta_dicts[-1]['DP_c_gh'] = delta_c_gh = ((np.max([DP_c_gh, self.params['DP_c_gh']]) - self.params['DP_c_gh'])/self.params['DP_c_gh'])**2
            
            delta_dicts[-1]['DP_h_cond'] = delta_h_cond = ((np.max([DP_h_cond, self.params['DP_h_cond']]) - self.params['DP_h_cond'])/self.params['DP_h_cond'])**2
            delta_dicts[-1]['DP_c_cond'] = delta_c_cond = ((np.max([DP_c_cond, self.params['DP_c_cond']]) - self.params['DP_c_cond'])/self.params['DP_c_cond'])**2
                        
            delta_dicts[-1]['DP_h_rec'] = delta_h_rec = ((np.max([DP_h_rec, self.params['DP_h_rec']]) - self.params['DP_h_rec'])/self.params['DP_h_rec'])**2
            delta_dicts[-1]['DP_c_rec'] = delta_c_rec = ((np.max([DP_c_rec, self.params['DP_c_rec']]) - self.params['DP_c_rec'])/self.params['DP_c_rec'])**2
            
            score_current = delta_exp + delta_pp + delta_h_gh + delta_c_gh + delta_h_cond + delta_c_cond + delta_h_rec + delta_c_rec
            
            RC_scores.append(score_current)
    
        # # remove all failing RCs
        # for rc in invalid_RCs:
        #     self.potential_RC.remove(rc) 
                        
        index_of_min = RC_scores.index(np.min(RC_scores))
        
        best_RC = self.potential_RC[index_of_min]
        delta_dict = delta_dicts[index_of_min]
        
        new_params_dict = {}
        
        new_params_dict['eta_exp']   = np.round(best_RC.components['Expander'].sizing.eta_is,3)
        new_params_dict['eta_pp']    = np.round(best_RC.components['Pump'].sizing.eta_is,3)

        new_params_dict['DP_h_gh']   = np.round(best_RC.components['GasHeater'].sizing.best_particle.DP_h)
        new_params_dict['DP_c_gh']   = np.round(best_RC.components['GasHeater'].sizing.best_particle.DP_c)
        
        new_params_dict['DP_h_cond'] = np.round(best_RC.components['Condenser'].sizing.best_particle.DP_h)
        new_params_dict['DP_c_cond'] = np.round(best_RC.components['Condenser'].sizing.best_particle.DP_c)
        
        new_params_dict['DP_h_rec']  = np.round(best_RC.components['Recuperator'].sizing.HX.DP_h)
        new_params_dict['DP_c_rec']  = np.round(best_RC.components['Recuperator'].sizing.HX.DP_c)
        
        return new_params_dict, np.min(RC_scores), delta_dict

    def size_components(self):
        
        i = 0
        n_pos = len(self.top_positions)
        
        self.potential_RC = []
        
        for allowable_position in self.top_positions:
            
            print(f"Component Optimization for top position : {i+1}/{n_pos}")
            i = i + 1
            
            self.it_var['P_high'] = allowable_position['x'][0]
            self.it_var['mdot'] = allowable_position['x'][1]
            self.it_var['mdot_HS'] = best_m_dot_HS = allowable_position['x'][1] * allowable_position['x'][2]
            self.it_var['mdot_CS'] = best_m_dot_CS = allowable_position['x'][1] * allowable_position['x'][7]
            self.it_var['eta_gh'] = allowable_position['x'][3]
            self.it_var['PP_gh'] = allowable_position['x'][4]
            self.it_var['eta_rec'] = allowable_position['x'][5]
            self.it_var['PP_cd'] = allowable_position['x'][6]
        
            self.HSource.set_properties(m_dot=best_m_dot_HS)
            self.CSource.set_properties(m_dot=best_m_dot_CS)
        
            # Estimate low pressure for initialization
            P_sat_T_CSource = PropsSI('P', 'T', self.CSource.T, 'Q', 0.5, self.fluid)
            P_crit_CO2 = PropsSI('PCRIT', self.fluid)
            P_low_guess = min(1.3 * P_sat_T_CSource, 0.8 * P_crit_CO2)
        
            try:
                
                if self.params['RC_ARCH'] == 'REC':
                
                    self.current_RC = REC_CO2_TC(
                        self.HSource, self.CSource, self.it_var['PP_gh'], self.params['PP_rec'], self.params['eta_pp'],
                        self.params['eta_exp'], self.it_var['eta_gh'], self.it_var['eta_rec'], self.it_var['PP_cd'], self.params['SC_cd'],
                        P_low_guess, self.it_var['P_high'], self.it_var['mdot'], DP_h_rec = self.params['DP_h_rec'], DP_c_rec = self.params['DP_c_rec'],
                        DP_h_gh = self.params['DP_h_gh'], DP_c_gh = self.params['DP_c_gh'], DP_h_cond = self.params['DP_h_cond'],
                        DP_c_cond = self.params['DP_c_cond'],
                        mute_print_flag=1
                    )
                            
                self.current_RC.solve()
                
            except Exception as e:
                print(f"⚠️ Failed to solve final RC circuit: {e}")
                self.RC = None
                self.eta = None
    
            # Size Components
            self.current_RC, flag = TCO2_rec_comp_sizing(self.current_RC)
            
            if flag == 1:
                self.potential_RC.append(self.current_RC)

        return

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
                'T': self.CSource.T,
                'P': self.CSource.p,
                'fluid': self.CSource.fluid,
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
                    if pen_i == 0 and np.isfinite(cost):
                        x_disc = discretize(x_i)
                        self.allowable_positions.append({
                            'x': x_disc.copy(),
                            'score': float(cost)
                        })
        
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
        
        while self.criterion == 0 and it < n_it_max:
        
            self.opt_RC(n_jobs=n_jobs, n_particles=n_particles, max_iter=max_iter, patience=patience, tol=tol, ntop = ntop)
            
            #%% 2) Size Components
    
            self.size_components()
    
            #%% 3) From best position update your performance guesses
    
            new_params, best_score, delta_dict = self.evaluate_systems()
            
            self.new_params = new_params
            self.delta_dict = delta_dict
            
            #%% 4) Fine Tune the Optimization with new values
            
            print("\n")
            print("----------------------------------------")
            print(f"New Values - Best Score : {best_score}")
            print("----------------------------------------")
    
            print(f"eta_exp :   {new_params['eta_exp']}   - {delta_dict['eta_exp']*100}")
            print(f"eta_pp :    {new_params['eta_pp']}    - {delta_dict['eta_pp']*100}")
            print(f"DP_h_gh :   {new_params['DP_h_gh']}   - {delta_dict['DP_h_gh']*100}")
            print(f"DP_c_gh :   {new_params['DP_c_gh']}   - {delta_dict['DP_c_gh']*100}")
            print(f"DP_h_cond : {new_params['DP_h_cond']} - {delta_dict['DP_h_cond']*100}")
            print(f"DP_c_cond : {new_params['DP_c_cond']} - {delta_dict['DP_c_cond']*100}")
            print(f"DP_h_rec :  {new_params['DP_h_rec']}  - {delta_dict['DP_h_rec']*100}")
            print(f"DP_c_rec :  {new_params['DP_c_rec']}  - {delta_dict['DP_c_rec']*100}")
            
            self.set_parameters(
                eta_exp    = new_params['eta_exp'],
                eta_pp     = new_params['eta_pp'],
                DP_h_gh    = (new_params['DP_h_gh']+self.params['DP_h_gh'])/2,
                DP_c_gh    = (new_params['DP_c_gh']+self.params['DP_c_gh'])/2,
                DP_h_cond  = (new_params['DP_h_cond']+self.params['DP_h_cond'])/2,
                DP_c_cond  = (new_params['DP_c_cond']+self.params['DP_c_cond'])/2,
                DP_h_rec   = (new_params['DP_h_rec']+self.params['DP_h_rec'])/2,
                DP_c_rec   = (new_params['DP_c_rec']+self.params['DP_c_rec'])/2,
                )
            
            self.criterion = 1

            for key in self.delta_dict:
                value = self.delta_dict[key]
                if value > 1e-3:
                    self.criterion = 0
                    break
            
            it = it + 1

            
        # #%% 1) Import RC
        
        # self.opt_RC(n_jobs=n_jobs, n_particles=n_particles, max_iter=max_iter, patience=patience, tol=tol, ntop = ntop)
        
        # #%% 2) Size Components

        # self.size_components()

        # #%% 3) From best position update your performance guesses

        # new_params, best_score, delta_dict = self.evaluate_systems()
        
        # print("\n")
        # print("----------------------------------------")
        # print(f"New Values - Best Score : {best_score}")
        # print("----------------------------------------")

        # print(f"eta_exp :   {new_params['eta_exp']}   - {delta_dict['eta_exp']*100}")
        # print(f"eta_pp :    {new_params['eta_pp']}    - {delta_dict['eta_pp']*100}")
        # print(f"DP_h_gh :   {new_params['DP_h_gh']}   - {delta_dict['DP_h_gh']*100}")
        # print(f"DP_c_gh :   {new_params['DP_c_gh']}   - {delta_dict['DP_c_gh']*100}")
        # print(f"DP_h_cond : {new_params['DP_h_cond']} - {delta_dict['DP_h_cond']*100}")
        # print(f"DP_c_cond : {new_params['DP_c_cond']} - {delta_dict['DP_c_cond']*100}")
        # print(f"DP_h_rec :  {new_params['DP_h_rec']}  - {delta_dict['DP_h_rec']*100}")
        # print(f"DP_c_rec :  {new_params['DP_c_rec']}  - {delta_dict['DP_c_rec']*100}")
        
        # self.new_params = new_params
        # self.delta_dict = delta_dict
        
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
    m_dot_CS_fact_bounds = [10,20]
    P_high_bounds = np.array([120, 160]) * 1e5
    m_dot_bounds = np.array([30,50])*n_MW
    
    # Discrete Variable choices
    eta_gh_disc = np.arange(0.8,0.99,0.01)
    PP_gh_disc = np.arange(1,10,1)
    eta_rec_disc= np.arange(0.7,0.98,0.02)
    PP_cd_disc = np.arange(1,15,1)
    
    # Set model parameters
    Optimizer.set_parameters(
        RC_ARCH= 'REC', # 'REC'
        
        # Pump
        eta_pp=0.8,
        
        # GasHeater
        DP_h_gh = 100*1e3,
        DP_c_gh = 4*1e5,

        # Recuperator
        PP_rec=0,
        DP_h_rec = 4*1e5,
        DP_c_rec = 2*1e5,
        
        # Expander
        eta_exp=0.9,
        
        # Condenser
        SC_cd=0.1,
        DP_h_cond = 2*1e5,
        DP_c_cond = 100*1e3,
        
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

    Optimizer.cycle_design(ntop = 5)

