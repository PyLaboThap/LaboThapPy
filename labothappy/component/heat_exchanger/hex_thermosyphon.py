# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:59:44 2023

@author: Basile
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve, minimize

from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

from toolbox.nozzle.airflow import nozzle
from correlations.heat_pipe.HP_h_coeffs import fg_radiative_htc_corr, radiative_htc_fg, external_flow_inline_bank, pool_boiling, external_flow_staggered_bank, ext_conv_boil #, h_cond_Th66, external_flow_finnedtubebank
# from component.heat_exchanger.steady_state.heat_pipe_based.modules.HP_tube_model import thermosyphon_model
from component.base_component import BaseComponent

# from labothappy.component.heat_exchanger.steady_state.heat_pipe_based.modules.HP_internal import Delta_P_v, thermal_res_esdu
#from labothappy.component.heat_exchanger.steady_state.heat_pipe_based.modules.HP_h_coeffs import radiative_htc_fg, external_flow_inline_bank, pool_boiling, external_flow_staggered_bank, ext_conv_boil #, h_cond_Th66, external_flow_finnedtubebank


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

    """
    Functions in this library are used for preliminary computations related to the internal behaviours of the heat pipe. 
        - figures_of_merit : Figures of merit of the fluid inside the heat pipe used in thermal_res_esdu
        - thermal_res_esdu : Internal thermal resistances of the heat pipe
        - Delta_P_v : Pressure loss inside the heat pipe
    """
    
    def figures_of_merit(self, fluid, T_sat):
        """
        ---- Inputs : -------- 
            
        fluid (char) : fluid name
        F_r : Filling ratio of the thermosyphon [/]
        T_sat : Saturation temperature of the fluid [K]
        
        ---- Outputs : --------
        
        phi_cond : Figure of merit for condensation
        phi_boil : Figure of merit for boiling
    
        ---- Reference(s) : --------
        
        Investigation of low Global Warming Potential working fluids for a closed two-phase thermosyphon
        Robert W. MacGregor, Peter A. Kew, David A. Reay
        
        """
        
        if T_sat > 273.15 + 373.6:
            T_sat = 273.15 + 373.6
        
        P_atm = 101325 # [Pa] : atmospheric pressure
        
        "1) Get fluid data"
        
        # Liquid properties
        (P_sat, k_l, rho_l, mu_l, cp_l) = PropsSI(('P','L','D','V','C'), 'T', T_sat, 'Q', 0, fluid) # Sat Pressure, thermal conductivity, density, viscosity, spec. heat capacity
        # Gaseous properties
        rho_v = PropsSI('D', 'T', T_sat, 'Q', 1, fluid) # density
        Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg
        
        "2) Compute figures of merit"
        
        # Condensation
        phi_cond = ((Dh_evap*k_l**3*rho_l**2)/mu_l)**(1/4)
    
        # Boiling
        phi_boil = (0.32*(rho_l**0.65)*(k_l**0.3)*(cp_l**0.7))/(rho_v**(0.25)*Dh_evap**(0.4)*mu_l**(0.1))*(P_sat/P_atm)**(0.23)
        
        return (phi_cond, phi_boil)
    
    def thermal_res_esdu(self, fluid, F_r, D_i, L_c, L_e, Q_dot_radial, T_sat):
        """
        ---- Inputs : -------- 
            
        fluid (char) : fluid name
        F_r : Filling ratio of the thermosyphon [/]
        D_i : Themosyphon internal diameter [m]
        L_c : Length of condensation zone [m]
        L_e : Length of evaporation zone [m]
        Q_dot_radial : Radial heat flux [W]
        T_sat : Saturation temperature of the fluid [K]
        
        ---- Outputs : --------
        
        R_cond : Internal Thermal Resistance for condenser [K/W]
        R_evap : Internal Thermal Resistance for evaporator [K/W] (composed of R_pool : pool boiling resistance 
                                                                   and R_film : Liquid film resistance) 
        Re_f : Reynolds number of the liquid flow in the thermosyphon
        
        ---- Reference(s) : --------
        
        Investigation of low Global Warming Potential working fluids for a closed two-phase thermosyphon
        Robert W. MacGregor, Peter A. Kew, David A. Reay
        
        """
      
        g = 9.81 # [m/s^2] : gravity acceleration constant   
        
        "1) Fluid properties"
    
        mu_l = PropsSI('V', 'T', T_sat, 'Q', 1, fluid) # liquid viscosity [Pa*s]
        Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg
    
        # Reynolds number
        Re_f = 4*Q_dot_radial/(Dh_evap*mu_l*np.pi*D_i)
        
        "2) Condenser resistance"
        
        # Figures of merit
        (phi_cond, phi_boil) = self.figures_of_merit(fluid, T_sat)
        
        # Intermediate resistance
        R_cond_EDSU = 0.235*Q_dot_radial**(1/3)/(D_i**(4/3)*g**(1/3)*L_c*phi_cond**(4/3))
        
        if Re_f > 1300:
            R_cond = R_cond_EDSU*191*Re_f**(-0.733)
        elif Re_f < 50: 
            R_cond = R_cond_EDSU		
        else:
            R_cond = R_cond_EDSU
    
        "3) Evaporator resistance"
        
        # Film resistance
        R_film = R_cond*(L_c/L_e)
        
        # Pool resistance
        R_pool = 1/(phi_boil*g**(0.2)*Q_dot_radial**(0.4)*(np.pi*D_i*L_e)**(0.6))
    
        # Evap resistance
        if R_film > R_pool:
            R_evap = R_pool
        else:
            R_evap = R_pool*F_r + R_film*(1-F_r)
        
        return (R_cond, R_evap, R_pool, R_film, Re_f)
    
    def Delta_P_v(self, fluid, D_v, L_eff, Q_dot_radial, T_sat):
        """
        ---- Inputs : -------- 
            
        fluid (char) : fluid name
        D_v : Vapor diameter (assumed to be the thermosyphon internal diameter) [m]
        L_eff : Effective length of the thermosyphon [m]
        Q_dot_radial : Radial heat flux [W]
        T_sat : Saturation temperature of the fluid [K]
        
        ---- Outputs : --------
        
        DP_v : Pressure drop of fluid inside the thermosyphon [Pa]
        
        ---- Reference(s) : --------
        
        /
        
        """
        
        if T_sat > 273.15 + 370:
            T_sat = 273.15 + 370
        
        "1) Get fluid data"
        
        (P_sat, rho_v, mu_v) = PropsSI(('P','D','V'),'T',T_sat,'Q',0,fluid) # Gas : (Pressure : Pa, density : kg/m^3, viscosity : Pa*s) 
        Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg
    
        A_v = (np.pi*D_v**2)/4 # vapour cross section area
        # u_v = Q_dot_radial/(A_v*Dh_evap*rho_v)
        
        "2) Compute Pressure Difference"
        
        # Reynolds number
        Re_a = (Q_dot_radial*D_v)/(mu_v*A_v*Dh_evap)
        
        if Re_a < 2300:
            DeltaP_v = (32*mu_v*Q_dot_radial*L_eff)/(rho_v*A_v*Dh_evap*D_v**2)
        else:
            DeltaP_v = (0.3164/Re_a**(1/4))*(Q_dot_radial**2 * L_eff)/(2*rho_v*A_v**2*Dh_evap**2*D_v)
            
        return DeltaP_v

#%% 

    def operating_limits(self, fluid, F_r, D_i, L_e, beta, geo, T_sat):
        """
        ---- Inputs : -------- 
            
        fluid (char) : fluid name
        F_r : Filling ratio of the thermosyphon [/]
        D_i : Themosyphon internal diameter [m]
        L_e : Length of evaporation zone [m]
        beta : Thermosyphon inclination angle [rad]
        geo (char) : geometry of the thermosyphon (circular or annular)
        T_sat : Saturation temperature of the fluid [K]
        
        ---- Outputs : --------
        
        Q_dot_ent (Entrainment limit) [W] : Even when there is sufficient liquid present in the thermosyphon to prevent dry-out occuring, 
                                            the overall rate of heat transfer is subject to another limit; this occurs when the rate of 
                                            entrainment of liquid by the vapour prevents the downward flow of liquid.
                                            
        Q_dot_boil (Boiling limit) [W] : Boiling limit occurs when a stable film of vapour is formed between the liquid and the heated 
                                         wall of the evaporator.
            
        Q_dot_son (Sonic limit) [W] : At low operating pressures, the vapour velocity may be appreciate compared with the sonic velocity 
                                      in the vapour.
            
        Q_dot_dry (Dryout limit) [W] : As applied to thermosyphons, the term "dryout" implies that the volume of the liquid fill is not 
                                       sufficient to cover all of the pipe above the pool with a film of liquid. Thus with a vertical pipe,
                                       most of the falling film of liquid would have evaporated before reaching the pool, leaving dry 
                                       patches, with a few rivulets of liquid returning to the pool; with an inclined pipe, dry patches 
                                       would appear at the top of the evaporator. the available evidence suggest that dryout in a vertical 
                                       thermosyphon is avoided if the volume of liquid fill meets the contitions called for in 
                                       ESDU 81038 Section 2.3 .
        
        Q_dot_vap (Vapor pressure limit) [W] : When operating a thermosyphon at a pressure substantially below atmospheric, the pressure 
                                               drop of the vapor may be significant compared to the pressure in the evaporator. Vapor 
                                               pressures are low but necessarily exceed zero at temperature close to the bottom of the 
                                               operational range of a heat pipe. The mimimum vapour pressure, which occurs at the closed 
                                               end of the condenser, can be very small. The pressure drop in the vapour duct, Deltap_v, is 
                                               then constrained by this effectively zero pressure and by the low vapour pressure existing 
                                               at the closed end of the evaporator. Because Deltap_v increases with the overall rate of heat 
                                               transfer, the constraint on Deltap_v requires Q_dot (if not otherwise limited) to be 
                                               approximately limited to a value, called the vapour limit in this text.
        
        ---- Reference(s) : --------
        
        ESDU. Heat Pipes - Performance of Two phase Closed Thermosyphons. London, U.K.: Engineering Sciences Data Unit; 1981."
        Enhancement of heat transport in thermosyphon airpreheater at high temperature with binary working fluid : A case study of TEG–water (F_1 with respect to Bo)
        
        """
        
        g = 9.81 # m/s^2 : gravity acceleration constant
        
        "1) Get fluid data"
        
        (P_sat, rho_v, mu_v, sigma) = PropsSI(('P','D','V','I'),'T',T_sat,'Q',0,fluid) # Gas : (Pressure : Pa, density : kg/m^3, viscosity : Pa*s, surface tension : N/m) 
        (rho_l, mu_l) = PropsSI(('D','V'),'T',T_sat,'Q',1,fluid) # Liquid : (Density : kg/m^3, viscosity : Pa*s)
        Dh_evap = PropsSI('H','T',T_sat,'Q',1,fluid) - PropsSI('H','T',T_sat,'Q',0,fluid) # Evaporation specific enthalpy J/kg
            
        "2) Parameter F_1"
        
        # Bond number (or Eötvös number, dimensionless number measuring the importance of gravitational forces compared to surface tension forces for the movement of liquid front)
        if geo == 'circular':
            Bo = D_i*(g*abs(rho_l-rho_v)/sigma)**(1/2)
        elif geo == 'annular':
            Bo = (D_i/2)*(g*abs(rho_l-rho_v)/sigma)**(1/2)
        else:
            print("Error. Geometry shall be defined to either 'annular' or 'circular'. \n")
            return
        
        # F_1 parameter
        
        if Bo > 11:
            F_1 = 8.2
        else:
            F_1 = -0.0331*Bo**2 + 0.8161*Bo + 3.2134
    
        "3) Parameter F_2"
        
        # Dimensionless pressure parameter
        
        K_p = P_sat/(g*abs(rho_l-rho_v)*sigma)**(1/2)
        
        # F_2 parameter
        
        if K_p <= 4e4:
            F_2 = K_p**(-0.17)
        else:
            F_2 = 0.165    
            
        "4) Parameter F_3 : function of inclination of the pipe (beta in rad)"
        
        if beta == np.pi/2:
            F_3 = 1
        else:
            if Bo>=4:
                F_3 = 0.306075777 + 3.88126365*beta - 8.71941343*beta**2 + 13.4517469*beta**3 - 11.9947367*beta**4 + 5.32966132*beta**5 - 0.93010743*beta**6
            elif Bo >=2:
                F_3 = 0.222344657 + 2.71186727*beta - 6.78785522*beta**2 + 13.1807234*beta**3 - 13.6843299*beta**4 + 6.73969856*beta**5 - 1.26256636*beta**6
            else:
                F_3 = 0.213964282 + 0.922627057*beta + 0.285274233*beta**2 - 0.721665592*beta**3 + 0.235369565*beta**4
        
        "5) Compute limits"
        
        # Entrainment limit (F_1*F_2*F_3 group can be called Kutateladze number)
        Q_dot_ent = F_1*F_2*F_3*Dh_evap*(rho_v**(1/2))*(sigma*g*abs(rho_l-rho_v))**(1/4)
        
        # Boiling limit
        Q_dot_boil = 0.12*rho_v*Dh_evap*((sigma*abs(rho_l-rho_v)*g)/rho_v**2)**(1/4)
            
        # Sonic limit
        Q_dot_son = 0.474*Dh_evap*(rho_v*P_sat)**(1/2)
        
        # Dryout limit
        Q_dot_dry = (2*(Dh_evap**2)*sigma*rho_v)**(3/7)*(rho_l*abs(rho_l-rho_v)*g*Dh_evap/(3*mu_l*L_e))**(1/7)*(F_r/(447*(1-F_r)))**(4/7)
    
        # Vapor pressure limit
        Q_dot_vap = (D_i**2)*Dh_evap*P_sat*rho_v/(64*mu_v*L_e)
        
        return (Q_dot_ent, Q_dot_boil, Q_dot_son, Q_dot_dry, Q_dot_vap)
    
    def equ_syst(self, p, *param):
        
        "System to be solved using fsolve in the thermosyphon_model function, look there for variable descriptions"
    
        # Unknowns and parameters
        (P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_su, T_e_o_ex, P_evap_ex, R_tube_tot) = p # Unknowns
        (fluid, fluid_cd, F_r, D_i, L_e, L_c, L_eff, T_e_o_su, T_c_o_ex, R_axial, A_c_o, A_e_o, R_wc, R_we, S_T, S_L, P_cond_ex, u_wf_in, N_col, D_o, M_dot_wf, L_M, p_CO2, p_H2O, P_evap_su, M_dot_g, u_gas_fume, h_wf, foul, arrang, f_0, f_1, f_2, f_3) = param # Parameters
    
        if np.isnan(T_e_o_ex) or np.isnan(T_we_o): # or np.isnan(h_e_o):
            print("NAN_in")
            print(T_e_o_ex)
            print(T_we_o)
        
        # heat coefficients and resistances for the condenser  
    
        try: 
            h_wf = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,fluid_cd)    
            flag_1_phase = 1
        except:
            flag_1_phase = 0
    
        if flag_1_phase == 1:        
            if arrang == 'Inline':
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_inline_bank(fluid_cd, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
                h_wf_ex = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,fluid_cd)
            else:
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_staggered_bank(fluid_cd, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
        else: 
            h_c_o = ext_conv_boil(D_o, fluid_cd, T_c_o_ex, T_wc_o, u_wf_in)        
        
        R_c_o = (1+foul)/(A_c_o*h_c_o) # condenser external thermal resistance : takes into account the fouling factor
        
        # heat coefficients and resistances for the evaporator         
        h_r_e_o = radiative_htc_fg(p_CO2, p_H2O, (T_e_o_su + T_e_o_ex)/2, T_we_o, L_M,f_0,f_1,f_2,f_3)
        
        if arrang == 'Inline':
            (h_c_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_inline_bank('Air', T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[0:5]
        else: 
            (h_c_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_staggered_bank('Air', T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, N_col, D_o, S_T, S_L)[0:5]
    
        h_e_o = h_r_e_o + h_c_e_o
        R_e_o = (1+foul)/(A_e_o*h_e_o) # condenser external thermal resistance : takes into account the fouling factor
        
        # Pressure drop
    
        f1 = PropsSI('P', 'T', T_v_e, 'Q', 1, fluid) - P_v_e # Saturation temperature at evaporator
        f2 = self.Delta_P_v(fluid, D_i, L_eff, Q_dot_radial, T_v_e) - DP_v
        
        f3 = P_v_e - DP_v - P_v_c # Pressure at condenser side => Impacts condenser saturation temperaturee
        f4 = PropsSI('T', 'P', P_v_c, 'Q', 1, fluid) - T_v_c
        
        # Thermosiphon internal sides thermal resistances
        
        (R_cond_i, R_evap_i, R_pool, R_film, Re_f) = self.thermal_res_esdu(fluid, F_r, D_i, L_c, L_e, Q_dot_radial, T_v_e)
        
        # Heat balances
        
        # Evaporator
        f5 = (T_e_o_su - T_we_o)/R_e_o - Q_dot_axial - Q_dot_radial # Heat transfer from outside fluid to evap wall
        f6 = (T_we_o - T_we_i)/R_we - Q_dot_radial # Heat transfer through evap wall
        f7 = (T_we_i - T_v_e)/R_evap_i - Q_dot_radial # Heat transfer from evap wall to thermosyphon fluid
        
        # Condenser
        f8 = (T_v_c - T_wc_i)/R_cond_i - Q_dot_radial # Heat transfer from thermosyphon fluid to cond wall 
        f9 = (T_wc_i - T_wc_o)/R_wc - Q_dot_radial # Heat transfer through cond wall
        f10 = (T_wc_o - T_c_o_su)/R_c_o - Q_dot_axial - Q_dot_radial # Heat transfer from cond wall to outside fluid 
    
        # Heat transfer through thermosyphon wall in the axial direction
        f11 = (T_we_i+T_we_o)/2 - T_we # average wall temperature of evaporator section
        f12 = (T_wc_i+T_wc_o)/2 - T_wc # average wall temperature of condenser section
        
        f13 = (T_we - T_wc)/R_axial - Q_dot_axial
        
        if flag_1_phase == 1:
            cp_wf = PropsSI('C', 'P', P_cond_ex, 'T', T_c_o_su, fluid_cd) # [J/(kg*K)]
            f14 = T_c_o_su + (Q_dot_axial + Q_dot_radial)/(M_dot_wf*cp_wf) - T_c_o_ex
        else:
            f14 = T_c_o_su - T_c_o_ex        
        
        cp_g = PropsSI('C', 'P', P_evap_su,'T',(T_e_o_ex + T_e_o_su)/2,'Air')
        
        f15 = T_e_o_su - (Q_dot_axial + Q_dot_radial)/(M_dot_g*cp_g) - T_e_o_ex
        f16 = P_evap_su - DP_evap - P_evap_ex  
        
        R_DP = abs((T_v_e - T_v_c)/Q_dot_radial)
        A_in = (np.pi/4)*D_i**2
        
        f17 = ((R_we + R_evap_i+ R_DP + R_cond_i + R_wc)**(-1) + R_axial**(-1))**(-1) - R_tube_tot
    
        return (f1, f2, f3, f4, f5[0], f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17)

    def thermosyphon_model(self, vect_init, T_e_o_su, T_c_o_ex, P_evap_su, P_cond_ex, u_gas_fume, M_dot_g, u_wf_in, M_dot_wf, h_wf):
    
        """
        ---- Inputs : -------- 
    
        beta : Thermosyphon inclination angle [rad]
        D_i : Themosyphon internal diameter [m]
        D_o : Themosyphon external diameter [m]        
        F_r : Filling ratio of the thermosyphon [/]
        fluid (char) : fluid name
        geo : tube geometry
        k_pipe : Thermal conductivity of the pipe [W/(m*K)]
        L_a : Length of adiabatic zone [m]
        L_c : Length of condensation zone [m]
        L_e : Length of evaporation zone [m]
        L_M : Mean beam length of a tube in the configuration [m]
        M_dot_g : Mass flowrate of fume gases in the evaporator [kg/s]
        M_dot_wf : Mass flowrate of working fluid in the condenser [kg/s]
        N_col : Number of tube columns [/]
        p_CO2 : Partial pressure of CO2 in the fume gases [/]
        p_H2O : Partial pressure of H2O in the fume gases [/]
        P_cond_su : Pressure of workingg fluid at the supply of the tube in the condenser side [/]
        P_evap_su : Pressure of fume gases at the supply of the tube in the evaporator side [/]
        S_L : Longitudinal pitch [m]
        S_T : Tranversal pitch [m]
        T_c_o : Temperature of the fluid supply outside the condenser [K]
        T_e_o : Temperature of the fluid supply outside the evaporator (fume gases) [K]
        u_gas_fume : Flow speed of fume gases in the evaporator [m/s]
        u_wf_in : Flow speed of oil in the condenser [m/s]
        A_casing : Casing area [m^2]
        W_casing : Casing width [m]
        arrang : Tube arrangement [-]
        h_wf : Inlet working fluid enthalpy [J/kg]
        foul : Fouling factor [-]
        vect_init : Vector with initial conditions of the equation system
        
        ---- Outputs : --------
        
        Q_dot : Heat transfer rate through the thermpsyphon [W]
        T_we_o : Surface temperature at the outside of the evaporator [K]
        P_v_e : Pressure in the evaporation zone [Pa] 
        T_v_e : Temperature in the evaporation zone [K]
        m_dot_v : Gas mass flowrate in the thermosyphon [kg/s]
        P_cond_ex : Pressure of oil after the tube [Pa]
        P_evap_ex : Pressure of fume gases after the tube [Pa]
        
        (Q_dot : Heat transfer rate through the thermpsyphon [W]
        
        Q_dot_ent : Entrainment limit [W]
        Q_dot_boil : Boiling limit [W]
        Q_dot_son : Sonic limit [W]
        Q_dot_vap : Vapor pressure limit [W]
        (Look in operating_limits function description for more info)
            
        T_we_o : Surface temperature at the outside of the evaporator [K]
        T_we_i : Surface temperature at the inside of the evaporator [K]
        T_v_e : Temperature in the evaporation zone [K]
        T_wc_o : Surface temperature at the outside of the condenser [K]
        T_wc_i : Surface temperature at the outside of the condenser [K]
        T_v_c : Temperature in the condensing zone [K]
        Q_dot_axial : Heat transfer rate in the thermosyphon wall in the axial direction [K]
        R_axial : Heat transfer resistance of the thermosyphon wall in the axial direction [W/K]
        P_v_e : Pressure in the evaporation zone [Pa]
        P_v_c : Pressure in the condensing zone [Pa]
        m_dot_v : Gas mass flowrate in the thermosyphon [kg/s]
        V_dot_v : Gas volume flowrate in the thermosyphon [m^3/s]
        V_v : Gas velocity in the thermosyphon [m/s]
        V_dot_l : Liquid volume flowrate in the thermosyphon [m^3/s])
        
        ---- Reference(s) : --------
        
        - Investigation of low Global Warming Potential working fluids for a closed two-phase thermosyphon
        Robert W. MacGregor, Peter A. Kew, David A. Reay
        
        - ESDU. Heat Pipes - Performance of Two phase Closed Thermosyphons. London, U.K.: Engineering Sciences Data Unit; 1981."
        
        """    
        
        "1) Design parameters"
        
        L_eff = self.L_a + (self.L_e + self.L_c)/2 # Thermosyphon effective length
        
        A_e_o = np.pi*self.L_e*self.params['D_o'] # evaporator ext area
        A_e_i = np.pi*self.L_e*self.D_i # evaporator int area
        A_c_o = np.pi*self.L_c*self.params['D_o'] 
            
        A_c_i = np.pi*self.L_c*self.D_i # condenser int area
        A_axial = (np.pi/4)*(self.params['D_o']**2 - self.D_i**2) # tube annuli cross section area
        
        "2) Determine thermal resistances"
        
        R_we = np.log(self.params['D_o']/self.D_i)/(2*np.pi*self.L_e*self.params['k_pipe']) # evaporator wall thermal resistance
        R_wc = np.log(self.params['D_o']/self.D_i)/(2*np.pi*self.L_c*self.params['k_pipe']) # condenser wall thermal resistance
        R_axial = L_eff/(self.params['k_pipe']*A_axial) # tube axial thermal resistance
        
        "3) System of equations : determine temperatures and heat rates"
        
        x_init = vect_init
        syst_param = (self.HP_fluid, self.su_C.fluid, self.params['F_r'], self.D_i, self.L_e, self.L_c, L_eff, T_e_o_su, T_c_o_ex, R_axial, A_c_o, A_e_o, R_wc, R_we, self.S_T, self.S_L, P_cond_ex, u_wf_in, self.N_col, self.params['D_o'], M_dot_wf, self.L_M, self.params['p_CO2'], self.params['p_H2O'], P_evap_su, M_dot_g, u_gas_fume, h_wf, self.params['foul'], self.params['arrang'], self.f_0, self.f_1, self.f_2, self.f_3)
    
        (sol) = fsolve(self.equ_syst, x0 = x_init, args = syst_param)
        
        (P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_su, T_e_o_ex, P_evap_ex, R_tube_tot) = sol
    
        # Compute air velocities
        V_max_air, a_gas_fume = external_flow_inline_bank('Air', T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, self.N_col, self.params['D_o'], self.S_T, self.S_L)[4:6]
        
        "4) Results computation"
        
        # Phase change heats at both thermosyphon ends
        Dh_evap_e = PropsSI('H', 'P', P_v_e, 'Q', 1, self.HP_fluid)-PropsSI('H', 'P', P_v_e, 'Q', 0, self.HP_fluid)
        Dh_evap_c = PropsSI('H', 'P', P_v_c, 'Q', 1, self.HP_fluid)-PropsSI('H', 'P', P_v_c, 'Q', 0, self.HP_fluid)
        
        # heat transfer coefficients
        try: 
            h_wf = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,self.su_C.fluid)    
            flag_1_phase = 1
        except:
            flag_1_phase = 0
        
        if flag_1_phase == 1:        
            if self.params['arrang'] == 'Inline':
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_inline_bank(self.su_C.fluid, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
                h_wf_ex = PropsSI('H', 'P',P_cond_ex,'T',T_c_o_ex,self.su_C.fluid)
            else:
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_staggered_bank(self.su_C.fluid, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
        
        else: 
            h_c_o = ext_conv_boil(self.params['D_o'], self.su_C.fluid, T_c_o_ex, T_wc_o, u_wf_in)        
            
        if self.params['arrang'] == 'Inline':
            (h_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_inline_bank(self.su_H.fluid, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
        else: 
            (h_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_staggered_bank(self.su_H.fluid, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_gas_fume, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
        
        R_c_o = 1/(A_c_o*h_c_o)
        R_e_o = 1/(A_e_o*h_e_o)
    
        # Results
        
        Q_dot = Q_dot_radial + Q_dot_axial
        
        m_dot_v = Q_dot_radial/Dh_evap_e # gas mass flow rate
        V_dot_v = m_dot_v/PropsSI('D', 'P', P_v_e, 'Q', 1, self.HP_fluid) # gas volumic flow rate
        V_v = V_dot_v/((np.pi*self.D_i**2)/4) # gas velocity in thermosyphon
        
        m_dot_l = Q_dot_radial/Dh_evap_c # gas mass flow rate
        V_dot_l = m_dot_l/PropsSI('D', 'P', P_v_c, 'Q', 0, self.HP_fluid) # gas volumic flow rate
        V_l = V_dot_l/((np.pi*self.D_i**2)/4) # gas velocity in thermosyphon
        
        "5) Operating limits"
        (Q_dot_ent, Q_dot_boil, Q_dot_son, Q_dot_dry, Q_dot_vap) = self.operating_limits(self.HP_fluid, self.params['F_r'], self.D_i, self.L_e, self.params['beta'], self.params['geo'], T_v_e)
    
        return (Q_dot, T_wc_o, T_we_o, P_v_e, T_v_e, m_dot_v, P_evap_ex, V_max_air, a_gas_fume, R_tube_tot, R_e_o, R_c_o)  


#%%            

    # With compressor speed and mass flow rate known
    def System(self, h_wf_out_guess, T_wf_out_guess):

        "1) Tube matrix parameters"
        
        self.L_a = self.params['H_core']/20  # 0.05  # [m] : adiabatic zone length
        self.L_e = self.params['coef_evap']*(self.params['H_core']-self.L_a) # 0.75 # np.linspace(0.1, 2, 20)  # 1.3 # 0.756  # [m] : evap zone length
        self.L_c = (1-self.params['coef_evap'])*(self.params['H_core']-self.L_a) # 0.75 # np.linspace(0.1, 2, 20)  # 0.75 # 0.17  # [m] : cond zone length
            
        if self.params['arrang'] == 'Inline':
            self.S_T = self.params['pitch_T'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch (écartement)
            self.S_L = self.params['pitch_L'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch
        else:
            self.S_T = self.params['pitch_T'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch (écartement)
            self.S_L = self.params['pitch_L'] * self.params['D_o']  # + 2*h_fins [m] : Tube pitch
    
        self.N_col = int(np.floor(self.params['L_core']/(self.S_T))-1)
        
        self.L_M = 1.08*(self.S_T*self.S_L-0.785*self.params['D_o']**2)/self.params['D_o']  # mean beam length
    
        self.D_i = self.params['D_o'] - 2*self.params['t']  # [m] : Inside diameter
    
        self.N_tubes_row = np.floor((self.params['W_core'] - 0.012)/self.S_T)
    
        # Evaporator
        self.W_casing = self.params['W_core']  # [m] : cross-sectional width of thermosyphon evaporator
        # [m] : cross-section height of thermosyphon evaporator
        self.h_casing = self.L_e + 0.016
        self.A_casing = self.W_casing*self.h_casing  # m^2 : cross-section area of thermosyphon evaporator
        
        # Vector for containing
        T_fg = np.zeros(self.N_col+1)
        P_fg = np.zeros(self.N_col+1)
        
        T_wf = np.zeros(self.N_col+1)
        h_wf = np.zeros(self.N_col+1)
        
        Q_dot_row = np.zeros(self.N_col)
        Q_dot_tube = np.zeros(self.N_col)
        T_v_e = np.zeros(self.N_col)
        
        T_wc_o = np.zeros(self.N_col)
        T_we_o = np.zeros(self.N_col)
        
        P_v_e = np.zeros(self.N_col)
        
        m_dot_v = np.zeros(self.N_col)
        
        R_tot_tube = np.zeros(self.N_col)  
        R_c_o = np.zeros(self.N_col)  
        R_e_o = np.zeros(self.N_col)  
        
        "2) External fluids data"
    
        self.HP_fluid = "Water"
    
        "2.1) Flue gases"
    
        A_chimney = (np.pi/4)*self.params['D_chimney']**2
    
        # Thermo state related :
        cp_fg_in = PropsSI('C', 'T', self.su_H.T,'P', self.su_H.p, self.su_H.fluid)
    
        # Flow rate related
        V_dot_gas_fume = self.su_H.m_dot / PropsSI('D', 'T', self.su_H.T, 'P', self.su_H.p, self.su_H.fluid)
        u_fg_chimney = V_dot_gas_fume/A_chimney
    
        "Nozzle computation - HTX Inlet"
        (u_gas_fume, T_fg[0], P_fg[0]) = nozzle(self.su_H.m_dot, self.su_H.T, cp_fg_in, u_fg_chimney, self.su_H.p, A_chimney, self.A_casing)
            
        DP_gas_nozzle = self.su_H.p - P_fg[0]
        DP_gas_nozzle_vect = DP_gas_nozzle
        
        "2.3) Evaporating Working Fluid"
    
        T_sat_wf = PropsSI('T','P',self.su_C.p,'Q',0,self.su_C.fluid)
        
        T_wf[0] = PropsSI('T', 'P',self.su_C.p,'H',h_wf_out_guess, self.su_C.fluid) # [K]
            
        H_casing = self.L_c + 0.02  # [m] height of oil divergent
        self.A_casing = H_casing*self.W_casing
    
        rho_wf_in = PropsSI('D', 'P',self.su_C.p,'H',h_wf_out_guess,self.su_C.fluid)  # [kg/m^3] : oil density at 110°C
        V_dot_wf = self.su_C.m_dot/rho_wf_in  # [m^3/s] : Oil volumetric flow rate
    
        u_wf_in = V_dot_wf/self.A_casing  # [m/s] : Oil velocity in the casing
        
        "3) Thermosyphon Simulation"
            
        vect_init = (7*101325, 7*101325, 1000, 1000, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 0, 513.15, 635.15, 101325,50) # First initial conditions
        # P_v_c, P_v_e, Q_dot_axial, Q_dot_radial, T_v_c, T_v_e, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
        DP_air = 0
        
        for i in range(self.N_col):
            R_tot_tube_mean = 0
            R_e_o_mean = 0
            R_c_o_mean = 0
            
            (Q_dot_tube[i], T_wc_o[i], T_we_o[i], P_v_e[i], T_v_e[i], m_dot_v[i], P_fg[i+1], V_max_air, a_gf, R_tot_tube[i], R_e_o[i], R_c_o[i]) = self.thermosyphon_model(vect_init, T_fg[i], T_wf[i], P_fg[i], self.su_C.p, u_gas_fume, self.su_H.m_dot, u_wf_in, self.su_C.m_dot, h_wf)  # call of the model for one tube
                        
            DP_air = DP_air + P_fg[i]-P_fg[i+1]
    
            Q_dot_row[i] = Q_dot_tube[i]*self.N_tubes_row
    
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
        u_gas_fume_out = V_dot_gas_fume_out/self.A_casing
        
        "Nozzle computation - HTX Inlet"
        (u_fg_out, T_fg_out, P_fg_out) = nozzle(self.su_H.m_dot, T_fg[-1], cp_g, u_gas_fume_out, P_fg[-1], self.A_casing, A_chimney)
        
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
        
        [self.f_0, self.f_1, self.f_2, self.f_3] = fg_radiative_htc_corr()
        
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
                print("T_in_input :", self.su_C.T)
                print("T_in_final :", C_T_in)
                print("T_out_final :", C_T_out_guess)
                print("-------------------------")                

