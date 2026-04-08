# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:57:09 2026

@author: Basile
"""

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
    # discrete_vars = input_data.get('discrete_vars', {})
    
    # for idx, allowed_vals in discrete_vars.items():
    #     allowed_vals = np.array(allowed_vals, dtype=float)
    #     x[idx] = allowed_vals[np.argmin(np.abs(allowed_vals - x[idx]))]
        
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
    P_low_guess = min(1.1*P_sat_T_CSource, 0.7*P_crit_CO2)
    
    print(x)
    
    if input_data['RC_ARCH'] == 'REC':
        RC = REC_CO2_TC(HSource, CSource, x[4], params['PP_rec'], params['eta_pp'],
                        params['eta_exp'], x[3], x[5], x[6], params['SC_cd'], 
                        P_low_guess, x[0], x[1], DP_h_rec=params['DP_h_rec'], DP_c_rec=params['DP_c_rec'], 
                        DP_h_gh=params['DP_h_gh'], DP_c_gh=params['DP_c_gh'], DP_h_cond=params['DP_h_cond'],
                        DP_c_cond=params['DP_c_cond'], mute_print_flag=0)

    elif input_data['RC_ARCH'] == 'basic':
        RC = basic_CO2_TC(HSource, CSource, params['PP_gh'], params['PP_rec'], params['eta_pp'],
                          params['eta_exp'], params['eta_gh'], params['PP_cd'], params['SC_cd'],
                          P_low_guess, x[0], x[1], mute_print_flag=0)

    RC.solve(method='wegstein')
    
    if not RC.converged:
        # cost, penalty, eta
        return 1000, np.inf, np.nan, RC

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
    
    print(f"x : {x}")
    
    if abs((W_dot_net - obj['W_dot'])/obj['W_dot']) > 2e-2:
        penalty += abs((W_dot_net - obj['W_dot'])/obj['W_dot']) * 100
        print(f"W_dot_net: {W_dot_net}")
        print(f"obj['W_dot']: {obj['W_dot']}")

    if abs((obj['eta'] - eta)/obj['eta']) > 2e-2:
        penalty += abs((obj['eta'] - eta)/obj['eta']) * 100
    
        print(f"eta: {eta}")
        print(f"obj['eta']: {obj['eta']}")
    
    print("------------------------------------\n")
    
    # objective = RC.components['GasHeater'].model.Q.Q_dot
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
    
    objective = 0 # (Q_gh*(-np.log(1-eta_gh)) + Q_rec*(-np.log(1-eta_rec)) + Q_cond*(-np.log(1-eta_cond)))/(Q_cond + Q_rec + Q_gh)
    
    cost = objective + penalty

    # Return all three so we can use them in objective_wrapper
    return cost, penalty, eta, RC

        # T_cold_source = 0.1+273.15
        # T_hot_source = 130+273.15

        # m_dot = 3.00123861e+01

        # eta_is_pp = 0.8
        # eta_is_exp = 0.9
        # eta_gh = 0.95
        # eta_rec = 0.8

        # PPTD_cd = 5
        # SC_cd = 0.1

        # Pinch_min_GH = 5
        # Pinch_min_REC = 0

        # P_high = 140*1e5
        # P_sat_T_CSource = PropsSI('P', 'T', T_cold_source,'Q',0.5,'CO2')
        # P_crit_CO2 = PropsSI('PCRIT','CO2')

        # P_low_guess = min(1.3*P_sat_T_CSource,0.8*P_crit_CO2)   
        
        # T_cold_source = 0.1+273.15
        
        # HSource = MassConnector()
        # HSource.set_properties(fluid = 'Water', T = T_hot_source, p = 5e5, m_dot = m_dot*1)
        
        # CSource = MassConnector()
        # CSource.set_properties(fluid = 'Water', T = T_cold_source, p = 5e5, m_dot = 10000)
        
        # DP_h_rec = 1*1e5
        # DP_c_rec = 2*1e5
        
        # DP_h_gh = 50*1e3
        # DP_c_gh = 2*1e5
        
        # DP_h_cond = 1*1e5
        # DP_c_cond = 50*1e3
        


eta_gh_disc = 0.95
PP_gh_disc = 5    
eta_rec_disc = 0.8
PP_cd_disc = 5    

P_high = 140*1e5
mdot = 30
m_dot_HS_fact = 1
eta_gh = 0.95
Pinch_min_GH = 5
eta_rec = 0.8
PP_cd = 5
m_dot_CS_fact = 200

x = [P_high, mdot, m_dot_HS_fact, eta_gh, Pinch_min_GH, eta_rec, PP_cd, m_dot_CS_fact]

input_data = {
    'fluid': 'CO2',
    'params' : {
        'PP_rec' : 0,
        'eta_pp' : 0.8,
        'eta_exp' : 0.9,
        'SC_cd' : 0.1,
        'DP_h_rec' : 0,
        'DP_c_rec' : 0,
        'DP_h_gh' : 0,
        'DP_c_gh' : 0,
        'DP_h_cond' : 0,
        'DP_c_cond' : 0,
        },
    
    'obj': {
        'W_dot' : 1e6,
        'eta' : 0.1
        },
    'HSource': {
        'T': 273.15+130, # K
        'P': 5*1e5, # Pa
        'fluid': 'Water'
    },
    'CSource': {
        'fluid': 'Water',
        'T': 273.15+0.1, # K
        'P': 5*1e5, # Pa
        # 'm_dot': self.CSource.m_dot # !!!
    },
    'RC_ARCH': 'REC',
    
    'discrete_vars': {
        3: eta_gh_disc,   # eta_gh
        4: PP_gh_disc,    # PP_gh
        5: eta_rec_disc,  # eta_rec
        6: PP_cd_disc,    # PP_cd
    },

    }
    
cost, penalty, eta, RC = system_RC_parallel(x, input_data)

