# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:59:44 2023

@author: Basile
"""
import __init__

import numpy as np
from CoolProp.CoolProp import PropsSI
from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector
from component.steady_state.heat_exchanger.heat_pipe_based.modules.Airflow import nozzle
from component.steady_state.heat_exchanger.heat_pipe_based.modules.P_max_steel_pipes import P_max_adm
from component.steady_state.heat_exchanger.heat_pipe_based.modules.HP_tube_model import thermosyphon_model
from component.base_component import BaseComponent

class HP_HTX(BaseComponent):
    def __init__(self):
        super().__init__()
        self.su_C = MassConnector() # Working fluid supply
        self.su_H = MassConnector() # Secondary fluid supply
        self.ex_C = MassConnector()
        self.ex_H = MassConnector()

        self.Q_dot = HeatConnector()
        self.res = None
        
        "Flue gas side"
        self.su_fg = None
        self.ex_fg = None
        
        "Working fluid side"
        self.su_wf = None
        self.ex_wf = None
        self.evap_type = None

#%% 

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_H_fluid', 'su_H_T', 'su_H_m_dot']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su_C.fluid is not None:
            self.inputs['su_C_fluid'] = self.su_C.fluid
        if self.su_C.h is not None:
            self.inputs['su_C_h'] = self.su_C.h
        if self.su_C.T is not None:
            self.inputs['su_C_T'] = self.su_C.T
        if self.su_C.m_dot is not None:
            self.inputs['su_C_m_dot'] = self.su_C.m_dot

        if self.su_H.fluid is not None:
            self.inputs['su_H_fluid'] = self.su_H.fluid
        if self.su_H.h is not None:
            self.inputs['su_H_h'] = self.su_H.h
        if self.su_H.T is not None:
            self.inputs['su_H_T'] = self.su_H.T
        if self.su_H.m_dot is not None:
            self.inputs['su_H_m_dot'] = self.su_H.m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'su_C_fluid' in self.inputs:
            self.su_C.set_fluid(self.inputs['su_C_fluid'])
        if 'su_C_T' in self.inputs:
            self.su_C.set_T(self.inputs['su_C_T'])
        if 'su_C_h' in self.inputs:
            self.su_C.set_h(self.inputs['su_C_h'])
        if 'su_C_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['su_C_m_dot'])
        if 'su_C_p' in self.inputs:
            self.su_C.set_p(self.inputs['su_C_p'])

        if 'su_H_fluid' in self.inputs:
            self.su_H.set_fluid(self.inputs['su_H_fluid'])
        if 'su_H_T' in self.inputs:
            self.su_H.set_T(self.inputs['su_H_T'])
        if 'su_H_h' in self.inputs:
            self.su_H.set_h(self.inputs['su_H_h'])
        if 'su_H_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['su_H_m_dot'])
        if 'su_H_p' in self.inputs:
            self.su_H.set_p(self.inputs['su_H_p'])

        return['su_C_fluid', 'su_C_h', 'su_C_m_dot', 'su_C_p', 'su_H_fluid', 'su_H_h', 'su_H_m_dot', 'su_H_p']

#%%

    def get_required_parameters(self):
        return [
            'p_CO2', 'p_H2O', 'beta', 'D_o', 't', 'F_r', 'k_pipe', 'geo', 'H_core', 'L_core', 'W_core',
            'coef_evap', 'foul', 'arrang', 'pitch_T', 'pitch_L', 'D_chimney', 'Bank_side'
        ]

    def print_setup(self):
        print("=== Heat Exchanger Setup ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_dot.Q_dot}")

        print("\nInputs:")
        for input in self.get_required_inputs():
            if input in self.inputs:
                print(f"  - {input}: {self.inputs[input]}")
            else:
                print(f"  - {input}: Not set")

        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")

        print("======================")


#%%            

    # With compressor speed and mass flow rate known
    def System(self, h_wf_out_guess, T_wf_out_guess):

        "1) Tube matrix parameters"
        
        L_a = self.params['H_core']/20  # 0.05  # [m] : adiabatic zone length
        L_e = self.params['coef_evap']*(self.params['H_core']-L_a) # 0.75 # np.linspace(0.1, 2, 20)  # 1.3 # 0.756  # [m] : evap zone length
        L_c = (1-self.params['coef_evap'])*(self.params['H_core']-L_a) # 0.75 # np.linspace(0.1, 2, 20)  # 0.75 # 0.17  # [m] : cond zone length
            
        if self.params['arrang'] == 'Inline':
            S_T = self.params['pitch_T'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch (écartement)
            S_L = self.params['pitch_L'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch
        else:
            S_T = self.params['pitch_T'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch (écartement)
            S_L = self.params['pitch_L'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch
    
        N_col = int(np.floor(self.params['L_core']/(S_T))-1)
        
        L_M = 1.08*(S_T*S_L-0.785*self.params['D_o']**2)/self.params['D_o']  # mean beam length
    
        D_i = self.params['D_o'] - 2*self.params['t']  # [m] : Inside diameter
    
        N_tubes_row = np.floor((self.params['W_core'] - 0.012)/S_T)
    
        # Evaporator
        W_casing = self.params['W_core']  # [m] : cross-sectional width of thermosyphon evaporator
        # [m] : cross-section height of thermosyphon evaporator
        h_casing = L_e + 0.016
        A_casing = W_casing*h_casing  # m^2 : cross-section area of thermosyphon evaporator
    
        N_tubes = N_tubes_row*np.ones(N_col)
    
        # Vector for containing
        T_fg = np.zeros(N_col+1)
        P_fg = np.zeros(N_col+1)
        
        T_wf = np.zeros(N_col+1)
        h_wf = np.zeros(N_col+1)
        
        Q_dot_row = np.zeros(N_col)
        Q_dot_tube = np.zeros(N_col)
        T_v_e = np.zeros(N_col)
        
        T_wc_o = np.zeros(N_col)
        T_we_o = np.zeros(N_col)
        
        P_v_e = np.zeros(N_col)
        
        m_dot_v = np.zeros(N_col)
        
        R_tot_tube = np.zeros(N_col)  
        R_c_o = np.zeros(N_col)  
        R_e_o = np.zeros(N_col)  
        
        "2) External fluids data"
    
        "2.1) Used fluids inside and outside the thermosyphon"
    
        fluid_thermosiphon = 'Water'
    
        "2.2) Fume gases"
    
        A_chimney = (np.pi/4)*self.params['D_chimney']**2
    
        # Thermo state related :
        cp_fg_in = PropsSI('C', 'T', self.su_H.T,'P', self.su_H.p, self.su_H.fluid)
    
        # Flow rate related
        V_dot_gas_fume = self.su_H.m_dot / PropsSI('D', 'T', self.su_H.T, 'P', self.su_H.p, self.su_H.fluid)
        u_fg_chimney = V_dot_gas_fume/A_chimney
    
        "Nozzle computation - HTX Inlet"
        (u_gas_fume, T_fg[0], P_fg[0]) = nozzle(self.su_H.m_dot, self.su_H.T, cp_fg_in, u_fg_chimney, self.su_H.p, A_chimney, A_casing)
            
        DP_gas_nozzle = self.su_H.p - P_fg[0]
        DP_gas_nozzle_vect = DP_gas_nozzle
        
        "2.3) Evaporating Working Fluid"
    
        T_sat_wf = PropsSI('T','P',self.su_C.p,'Q',0,self.su_C.fluid)
        
        T_wf[0] = PropsSI('T', 'P',self.su_C.p,'H',h_wf_out_guess, self.su_C.fluid) # [K]
            
        H_casing = L_c + 0.02  # [m] height of oil divergent
        A_casing = H_casing*W_casing
    
        rho_wf_in = PropsSI('D', 'P',self.su_C.p,'H',h_wf_out_guess,self.su_C.fluid)  # [kg/m^3] : oil density at 110°C
        V_dot_wf = self.su_C.m_dot/rho_wf_in  # [m^3/s] : Oil volumetric flow rate
    
        u_wf_in = V_dot_wf/A_casing  # [m/s] : Oil velocity in the casing
        
        "3) Thermosyphon Simulation"
            
        vect_init = (7*101325, 7*101325, 1000, 1000, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 0, 513.15, 635.15, 101325,50) # First initial conditions
        # P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
        DP_air = 0
        
        for i in range(N_col):
            R_tot_tube_mean = 0
            R_e_o_mean = 0
            R_c_o_mean = 0
            
            (Q_dot_tube[i], T_wc_o[i], T_we_o[i], P_v_e[i], T_v_e[i], m_dot_v[i], P_fg[i+1], V_max_air, a_gf, R_tot_tube[i], R_e_o[i], R_c_o[i]) = thermosyphon_model(self.params['beta'], D_i, self.params['D_o'], self.params['F_r'], self.su_C.fluid, self.su_H.fluid, fluid_thermosiphon, self.params['geo'], self.params['k_pipe'], L_a, L_c, L_e, L_M, self.su_H.m_dot, self.su_C.m_dot, N_col, self.params['p_CO2'], self.params['p_H2O'], self.su_C.p, P_fg[i], S_L, S_T, T_wf[i], T_fg[i], u_gas_fume, u_wf_in, A_casing, W_casing, self.params['arrang'], h_wf[i], self.params['foul'], vect_init)  # call of the model for one tube
            
            # if i == 0:
            #     P_in_max = P_v_e[i]
            #     T_in_max = T_v_e[i]
            #     P_max_adm_vect = P_max_adm(D_o*1000, T_we_o[i], 0) # !!! attention au "0"
    
            DP_air = DP_air + P_fg[i]-P_fg[i+1]
    
            Q_dot_row[i] = Q_dot_tube[i]*N_tubes[i]
    
            cp_g = PropsSI('C', 'P', P_fg[i], 'T', T_fg[i], self.su_H.fluid)
            T_fg[i+1] = T_fg[i] - Q_dot_row[i]/(cp_g*self.su_H.m_dot)
            
            if i == 0:
                h_wf[i] = PropsSI('H', 'P', self.su_C.p, 'T', T_wf[i], self.su_C.fluid)
                
            h_wf_row = h_wf[i] - Q_dot_row[i]/self.su_C.m_dot
            h_wf[i+1] = h_wf_row
            
            T_wf[i+1] = PropsSI('T', 'P', self.su_C.p, 'H', h_wf_row, self.su_C.fluid)
            
            R_tot_tube_mean = R_tot_tube_mean + R_tot_tube[i]
            R_e_o_mean = R_e_o_mean + R_e_o[i]
            R_c_o_mean = R_c_o_mean + R_c_o[i]
            
            vect_init = (P_v_e[i],P_v_e[i],1000,Q_dot_tube[i],T_v_e[i],T_v_e[i],T_wc_o[i],T_wc_o[i],T_wc_o[i],T_we_o[i],T_we_o[i],T_we_o[i], 0, T_wf[i+1], T_fg[i+1], P_fg[i+1],R_tot_tube[i])
                        # P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
            
        rho_fg_out = PropsSI('D', 'T',T_fg[-1],'P',P_fg[-1],'air')
        
        V_dot_gas_fume_out = self.su_H.m_dot/rho_fg_out
        u_gas_fume_out = V_dot_gas_fume_out/A_casing
        
        "Nozzle computation - HTX Inlet"
        (u_fg_out, T_fg_out, P_fg_out) = nozzle(self.su_H.m_dot, T_fg[-1], cp_g, u_gas_fume_out, P_fg[-1], A_casing, A_chimney)
        
        # print((P_gas_fume_out - P_gas_fume[-1]))
        
        DP_gas_nozzle_vect = DP_gas_nozzle_vect + (P_fg_out - P_fg[-1])
            
        return T_wf[-1], h_wf[-1], h_wf[0], P_fg_out, T_fg_out
           
    #------------------------------------------------------------------------
    def solve(self, n_it_max, res_tol, C_T_out_guess_1, C_T_out_guess_2):
        
        self.check_calculable()
        self.check_parametrized()
        
        if self.calculable != True:
            print("Supply fluid not fully known.")
            return
      
        if self.parametrized != True:
            print("Heat Exchanger Parameters not fully known.")
            return    
      
        if C_T_out_guess_1 > C_T_out_guess_2:
            temp = C_T_out_guess_1
            C_T_out_guess_2 = C_T_out_guess_1
            C_T_out_guess_2 = temp
        
        if self.params['arrang'] == "inline":
            self.set_parameters(arrang = "Inline")
    
        if self.params['arrang'] == "staggered":
            self.set_parameters(arrang = "Staggered")
        
        C_h_in = PropsSI('H','P', self.su_C.p,'T', self.su_C.T,self.su_C.fluid)
        
        C_T_out_guess = (C_T_out_guess_1 + C_T_out_guess_2)/2
        C_h_out_guess = PropsSI('H','P', self.su_C.p,'T', C_T_out_guess, self.su_C.fluid)    

        n_it = 0
        self.res = 10000

        import warnings

        # Ignore all RuntimeWarnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        while abs(self.res) > res_tol:
                    
            if n_it >= n_it_max:
                print("No convergence in the fixed number of iterations")
                return
                
            (C_T_in, C_h_it, C_h_out, H_P_out, H_T_out) = self.System(C_h_out_guess, C_T_out_guess)
            
            "Compute residuals"
                        
            self.res = C_h_in - C_h_it

            if self.res > 0: # h_real > h_est => T_in > T_in,est => T_out > T_out_est => Rise the value of T_out_guess_1
                C_T_out_guess_1 = PropsSI('T', 'P', self.su_C.p,'H',C_h_out,self.su_C.fluid)
    
            else:
                C_T_out_guess_2 = PropsSI('T', 'P', self.su_C.p,'H',C_h_out,self.su_C.fluid)
                
            C_T_out_guess = (C_T_out_guess_1 + C_T_out_guess_2)/2
            C_h_out_guess = PropsSI('H', 'P', self.su_C.p, 'T', C_T_out_guess,self.su_C.fluid)

            n_it = n_it + 1
            
            # if n_it%5 == 0:
            print("-------------------------")
            print("Iteration : ",n_it, " / ", n_it_max)
            print("T_in_guess :", C_T_in)
            print("T_out_guess :", C_T_out_guess)
            print("C_T_out_guesses : [", C_T_out_guess_1, " ; ", C_T_out_guess_2, "]")
                    
        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_m_dot(self.su_H.m_dot)  # Example mass flow rate [kg]
        self.ex_H.set_T(H_T_out) # Example temperature [K]
        self.ex_H.set_p(H_P_out)  # Example Pressure [Pa]

        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)  # Example mass flow rate [kg]
        self.ex_C.set_T(C_T_in) # Example temperature [K]
        self.ex_C.set_p(self.su_C.p)  # Example Pressure [Pa]
        
        if abs(self.res) < res_tol:
                print("-------------------------")
                print("Success !")
                print("-------------------------")                
                print("Iteration : ",n_it, " / ", n_it_max)
                print("T_in_input :", self.su_wf.T)
                print("T_in_final :", C_T_in)
                print("T_out_final :", C_T_out_guess)
                print("-------------------------")                

