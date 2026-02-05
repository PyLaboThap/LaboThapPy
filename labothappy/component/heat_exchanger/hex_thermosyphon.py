# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 09:59:44 2023

@author: Basile
"""

import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector

from toolbox.nozzle.airflow import nozzle
from correlations.heat_pipe.HP_h_coeffs import fg_radiative_htc_corr, radiative_htc_fg, external_flow_inline_bank, external_flow_staggered_bank, ext_conv_boil #, h_cond_Th66, external_flow_finnedtubebank

from component.base_component import BaseComponent

class HexThermosyphon(BaseComponent):
    """
    **Component**: Heat Pipe Heat Exchanger (HP_HTX)
    
    **Model**: Thermosyphon-Based Heat Pipe Heat Exchanger with Detailed Thermal and Fluid Dynamics
    
    **Description**:
    
        This model simulates a heat pipe-based heat exchanger (HP_HTX) that uses thermosyphon principles for efficient heat transfer between two fluids. It integrates detailed thermal resistances, flow dynamics, and operating limit calculations, including entrainment, boiling, sonic, dryout, and vapor pressure limits. The model accounts for complex thermal interactions, pressure drops, phase changes, and geometrical configurations.
    
    **Assumptions**:
    
        - Steady-state operation.
        - Uniform flow distribution across tube banks.
        - Neglects transient effects and startup behavior.
        - Fluid properties are retrieved from CoolProp.
        - Correlations for external and internal heat transfer coefficients are used (e.g., ESDU, radiative, convective, pool boiling).
        - Operating limits are checked based on established methods (e.g., ESDU guidelines).
    
    **Connectors**:
    
        su_C (MassConnector): Mass connector for the cold-side working fluid (e.g., oil or secondary fluid).
        su_H (MassConnector): Mass connector for the hot-side fluid (e.g., flue gas).
        ex_C (MassConnector): Mass connector for the cold-side exhaust.
        ex_H (MassConnector): Mass connector for the hot-side exhaust.
        Q_dot (HeatConnector): Heat transfer connector for the total exchanged heat.
    
    **Parameters**:
    
        p_CO2 (float): Partial pressure of CO₂ in the flue gas.
        p_H2O (float): Partial pressure of H₂O in the flue gas.
        beta (float): Inclination angle of the thermosyphon [rad].
        D_o (float): Outer diameter of the thermosyphon tube [m].
        t (float): Tube wall thickness [m].
        F_r (float): Filling ratio of the thermosyphon [-].
        k_pipe (float): Thermal conductivity of the pipe material [W/(m·K)].
        geo (str): Geometry of the thermosyphon ('circular' or 'annular').
        H_core, L_core, W_core (float): Geometric dimensions of the core [m].
        coef_evap (float): Fraction of core height allocated to the evaporator.
        foul (float): Fouling factor for external surfaces [-].
        arrang (str): Tube arrangement ('Inline' or 'Staggered').
        pitch_T, pitch_L (float): Transverse and longitudinal tube pitch ratios [-].
        D_chimney (float): Chimney diameter for flue gas outlet [m].
        Bank_side (str): Bank side configuration (optional/advanced).
    
    **Inputs**:

        P_su_H: Hot suction side pressure. [Pa]

        h_su_H: Hot suction side enthalpy. [J/kg]

        fluid_H: Hot suction side fluid. [-]

        m_dot_H: Hot suction side mass flowrate. [kg/s]

        P_su_C: Cold suction side pressure. [Pa]

        h_su_C: Cold suction side enthalpy. [J/kg]

        fluid_C: Cold suction side fluid. [-]

        m_dot_C: Cold suction side mass flowrate. [kg/s]

    **Ouputs**:

        h_ex_H: Hot exhaust side specific enthalpy. [J/kg]

        P_ex_H: Hot exhaust side pressure. [Pa]

        h_ex_C: Cold exhaust side specific enthalpy. [J/kg]

        P_ex_C: Cold exhaust side pressure. [Pa]
    
        Q_dot: Total heat transfer rate across the exchanger [W].
        
        Operating limits: Entrainment limit, boiling limit, sonic limit, dryout limit, vapor pressure limit [W].
        
        Thermal resistances (R_e_o, R_c_o, R_tube_tot) [K/W].
        
        Intermediate states: Temperatures, pressures, mass flow rates, velocities, and heat transfer rates across discretization segments.
    
    **References**:
    
        - ESDU. Heat Pipes - Performance of Two-Phase Closed Thermosyphons. Engineering Sciences Data Unit; 1981.
        - MacGregor, R.W., Kew, P.A., Reay, D.A. (2005). Investigation of Low Global Warming Potential Working Fluids for a Closed Two-Phase Thermosyphon.
        - Custom correlations for radiative and convective heat transfer, based on external flow arrangements and pool boiling.
    
    """


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
        return['fluid_C', 'h_su_C', 'P_su_C', 'm_dot_C', 'h_su_H', 'P_su_H', 'm_dot_H']

    def get_required_parameters(self):
        return [
            'p_CO2', 'p_H2O', 'beta', 'D_o', 't', 'F_r', 'k_pipe', 'geo', 'H_core', 'L_core', 'W_core',
            'coef_evap', 'foul', 'arrang', 'pitch_T', 'pitch_L', 'D_chimney', 'Bank_side', 'HP_fluid'
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
        self.AS_HP.update(CP.QT_INPUTS, 0, T_sat)
        
        P_sat = self.AS_HP.p()
        k_l = self.AS_HP.conductivity()
        rho_l = self.AS_HP.rhomass()
        mu_l = self.AS_HP.viscosity()
        cp_l = self.AS_HP.cpmass()
        h_l = self.AS_HP.hmass()
        
        # Gaseous properties
        self.AS_HP.update(CP.QT_INPUTS, 1, T_sat)
        rho_v = self.AS_HP.rhomass()
        h_v = self.AS_HP.hmass()

        Dh_evap = h_v - h_l # Evaporation specific enthalpy J/kg
        
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
    
        # Liquid properties
        self.AS_HP.update(CP.QT_INPUTS, 0, T_sat)
        mu_l = self.AS_HP.viscosity()
        h_l = self.AS_HP.hmass()
        
        # Gaseous properties
        self.AS_HP.update(CP.QT_INPUTS, 1, T_sat)
        h_v = self.AS_HP.hmass()
 
        Dh_evap = h_v - h_l # Evaporation specific enthalpy J/kg
    
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
    
    def Delta_P_v(self, D_v, L_eff, Q_dot_radial, T_sat):
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
        
        # Liquid properties
        self.AS_HP.update(CP.QT_INPUTS, 0, T_sat)   
        h_l = self.AS_HP.hmass()
        
        # Gaseous properties
        self.AS_HP.update(CP.QT_INPUTS, 1, T_sat)
        rho_v = self.AS_HP.rhomass()
        h_v = self.AS_HP.hmass()
        mu_v = self.AS_HP.viscosity()
        
        Dh_evap = h_v - h_l # Evaporation specific enthalpy J/kg
    
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
        
                
        # Liquid properties
        self.AS_HP.update(CP.QT_INPUTS, 0, T_sat)
        
        P_sat = self.AS_HP.p()
        rho_l = self.AS_HP.rhomass()
        mu_l = self.AS_HP.viscosity()
        h_l = self.AS_HP.hmass()
        
        # Gaseous properties
        self.AS_HP.update(CP.QT_INPUTS, 1, T_sat)
        rho_v = self.AS_HP.rhomass()
        h_v = self.AS_HP.hmass()
        mu_v = self.AS_HP.viscosity()
        
        # 2-Phase Properties
        self.AS_HP.update(CP.QT_INPUTS, 0.5, T_sat)
        sigma = self.AS_HP.surface_tension()
        
        Dh_evap = h_v - h_l # Evaporation specific enthalpy J/kg
            
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
        (P_HP_c, P_HP_h, Q_dot_axial, Q_dot_radial, T_HP_c, T_HP_h, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_su, T_e_o_ex, P_evap_ex, R_tube_tot) = p # Unknowns
        (fluid, su_C_fluid, F_r, D_i, L_e, L_c, L_eff, T_e_o_su, T_c_o_ex, R_axial, A_c_o, A_e_o, R_wc, R_we, S_T, S_L, P_cond_ex, u_wf_in, N_col, D_o, M_dot_wf, L_M, p_CO2, p_H2O, P_evap_su, M_dot_g, u_h, h_wf, foul, arrang, f_0, f_1, f_2, f_3) = param # Parameters
    
        if np.isnan(T_e_o_ex) or np.isnan(T_we_o): # or np.isnan(h_e_o):
            print("NAN_in")
            print(T_e_o_ex)
            print(T_we_o)
        
        # heat coefficients and resistances for the condenser  
    
        try: 
            self.AS_C.update(CP.PT_INPUTS, P_cond_ex, T_c_o_ex)
            flag_1_phase = 1
        except:
            flag_1_phase = 0
    
        if flag_1_phase == 1:        
            if arrang == 'Inline':
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_inline_bank(self.AS_C, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]                
            else:
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_staggered_bank(self.AS_C, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, N_col, D_o, S_T, S_L)[0:5]
        else: 
            h_c_o = ext_conv_boil(D_o, self.AS_C, T_c_o_ex, T_wc_o, u_wf_in)        
        
        R_c_o = (1+foul)/(A_c_o*h_c_o) # condenser external thermal resistance : takes into account the fouling factor
        
        # heat coefficients and resistances for the evaporator         
        h_r_e_o = radiative_htc_fg(p_CO2, p_H2O, (T_e_o_su + T_e_o_ex)/2, T_we_o, L_M,f_0,f_1,f_2,f_3)
        
        if arrang == 'Inline':
            (h_c_e_o, DP_evap, Nu_evap, Re_evap, V_max_air) = external_flow_inline_bank(self.AS_H, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_h, N_col, D_o, S_T, S_L)[0:5]
        else: 
            (h_c_e_o, DP_evap, Nu_evap, Re_evap, V_max_air) = external_flow_staggered_bank(self.AS_H, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_h, N_col, D_o, S_T, S_L)[0:5]
    
        h_e_o = h_r_e_o + h_c_e_o
        R_e_o = (1+foul)/(A_e_o*h_e_o) # condenser external thermal resistance : takes into account the fouling factor
        
        # Pressure drop
        
        self.AS_HP.update(CP.QT_INPUTS, 1, T_HP_h)
        
        f1 = self.AS_HP.p() - P_HP_h # Saturation temperature at evaporator
        f2 = self.Delta_P_v(D_i, L_eff, Q_dot_radial, T_HP_h) - DP_v

        self.AS_HP.update(CP.PQ_INPUTS, P_HP_c, 1)
        
        f3 = P_HP_h - DP_v - P_HP_c # Pressure at condenser side => Impacts condenser saturation temperaturee
        f4 = self.AS_HP.T() - T_HP_c
        
        # Thermosiphon internal sides thermal resistances
        
        (R_cond_i, R_evap_i, R_pool, R_film, Re_f) = self.thermal_res_esdu(fluid, F_r, D_i, L_c, L_e, Q_dot_radial, T_HP_h)
        
        # Heat balances
        
        # Evaporator
        f5 = (T_e_o_su - T_we_o)/R_e_o - Q_dot_axial - Q_dot_radial # Heat transfer from outside fluid to evap wall
        f6 = (T_we_o - T_we_i)/R_we - Q_dot_radial # Heat transfer through evap wall
        f7 = (T_we_i - T_HP_h)/R_evap_i - Q_dot_radial # Heat transfer from evap wall to thermosyphon fluid
        
        # Condenser
        f8 = (T_HP_c - T_wc_i)/R_cond_i - Q_dot_radial # Heat transfer from thermosyphon fluid to cond wall 
        f9 = (T_wc_i - T_wc_o)/R_wc - Q_dot_radial # Heat transfer through cond wall
        f10 = (T_wc_o - T_c_o_su)/R_c_o - Q_dot_axial - Q_dot_radial # Heat transfer from cond wall to outside fluid 
    
        # Heat transfer through thermosyphon wall in the axial direction
        f11 = (T_we_i+T_we_o)/2 - T_we # average wall temperature of evaporator section
        f12 = (T_wc_i+T_wc_o)/2 - T_wc # average wall temperature of condenser section
        
        f13 = (T_we - T_wc)/R_axial - Q_dot_axial
        
        if flag_1_phase == 1:
            self.AS_C.update(CP.PT_INPUTS, P_cond_ex, T_c_o_su)
            
            cp_wf = self.AS_C.cpmass()
            f14 = T_c_o_su + (Q_dot_axial + Q_dot_radial)/(M_dot_wf*cp_wf) - T_c_o_ex
        else:
            f14 = T_c_o_su - T_c_o_ex        
        
        self.AS_C.update(CP.PT_INPUTS, P_evap_su, (T_e_o_ex + T_e_o_su)/2)
        cp_H = self.AS_H.cpmass() # PrI('C', 'P', P_evap_su,'T',(T_e_o_ex + T_e_o_su)/2,'Air')
        
        f15 = T_e_o_su - (Q_dot_axial + Q_dot_radial)/(M_dot_g*cp_H) - T_e_o_ex
        f16 = P_evap_su - DP_evap - P_evap_ex  
        
        R_DP = abs((T_HP_h - T_HP_c)/Q_dot_radial)
        
        f17 = ((R_we + R_evap_i+ R_DP + R_cond_i + R_wc)**(-1) + R_axial**(-1))**(-1) - R_tube_tot
    
        return (f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17)

    def thermosyphon_model(self, vect_init, T_e_o_su, T_c_o_ex, P_evap_su, P_cond_ex, u_h, M_dot_g, u_wf_in, M_dot_wf, h_wf):
    
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
        u_h : Flow speed of hot fluid [m/s]
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
        P_HP_h : Pressure in the evaporation zone [Pa] 
        T_HP_h : Temperature in the evaporation zone [K]
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
        T_HP_h : Temperature in the evaporation zone [K]
        T_wc_o : Surface temperature at the outside of the condenser [K]
        T_wc_i : Surface temperature at the outside of the condenser [K]
        T_HP_c : Temperature in the condensing zone [K]
        Q_dot_axial : Heat transfer rate in the thermosyphon wall in the axial direction [K]
        R_axial : Heat transfer resistance of the thermosyphon wall in the axial direction [W/K]
        P_HP_h : Pressure in the evaporation zone [Pa]
        P_HP_c : Pressure in the condensing zone [Pa]
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
        A_c_o = np.pi*self.L_c*self.params['D_o']        
        A_axial = (np.pi/4)*(self.params['D_o']**2 - self.D_i**2) # tube annuli cross section area
        
        "2) Determine thermal resistances"
        
        R_we = np.log(self.params['D_o']/self.D_i)/(2*np.pi*self.L_e*self.params['k_pipe']) # evaporator wall thermal resistance
        R_wc = np.log(self.params['D_o']/self.D_i)/(2*np.pi*self.L_c*self.params['k_pipe']) # condenser wall thermal resistance
        R_axial = L_eff/(self.params['k_pipe']*A_axial) # tube axial thermal resistance
        
        "3) System of equations : determine temperatures and heat rates"
        
        x_init = vect_init
        syst_param = (self.AS_HP.fluid_names()[0], self.su_C.fluid, self.params['F_r'], self.D_i, self.L_e, self.L_c, L_eff, T_e_o_su, T_c_o_ex, R_axial, A_c_o, A_e_o, R_wc, R_we, self.S_T, self.S_L, P_cond_ex, u_wf_in, self.N_col, self.params['D_o'], M_dot_wf, self.L_M, self.params['p_CO2'], self.params['p_H2O'], P_evap_su, M_dot_g, u_h, h_wf, self.params['foul'], self.params['arrang'], self.f_0, self.f_1, self.f_2, self.f_3)
    
        (sol) = fsolve(self.equ_syst, x0 = x_init, args = syst_param)
        
        (P_HP_c, P_HP_h, Q_dot_axial, Q_dot_radial, T_HP_c, T_HP_h, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_su, T_e_o_ex, P_evap_ex, R_tube_tot) = sol
    
        # Compute air velocities
        V_max_air, a_gas_fume = external_flow_inline_bank(self.AS_H, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_h, self.N_col, self.params['D_o'], self.S_T, self.S_L)[4:6]
        
        "4) Results computation"
                        
        # heat transfer coefficients
        try: 
            self.AS_C.update(CP.PT_INPUTS, P_cond_ex, T_c_o_ex)
            flag_1_phase = 1
        except:
            flag_1_phase = 0
        
        if flag_1_phase == 1:        
            if self.params['arrang'] == 'Inline':
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_inline_bank(self.AS_C, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
            else:
                (h_c_o,DP_cond,Nu_cond,Re_cond,V_max_oil) = external_flow_staggered_bank(self.AS_C, T_c_o_su, T_c_o_ex, T_wc_o, P_cond_ex, u_wf_in, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
        
        else: 
            h_c_o = ext_conv_boil(self.params['D_o'], self.AS_C, T_c_o_ex, T_wc_o, u_wf_in)        
            
        if self.params['arrang'] == 'Inline':
            (h_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_inline_bank(self.AS_H, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_h, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
        else: 
            (h_e_o, DP_evap, Nu_evap, Re_evap,V_max_air) = external_flow_staggered_bank(self.AS_H, T_e_o_su, T_e_o_ex, T_we_o, P_evap_su, u_h, self.N_col, self.params['D_o'], self.S_T, self.S_L)[0:5]
        
        R_c_o = 1/(A_c_o*h_c_o)
        R_e_o = 1/(A_e_o*h_e_o)
    
        # Results
        
        Q_dot = Q_dot_radial + Q_dot_axial
                
        "5) Operating limits"
        (Q_dot_ent, Q_dot_boil, Q_dot_son, Q_dot_dry, Q_dot_vap) = self.operating_limits(self.AS_HP.fluid_names()[0], self.params['F_r'], self.D_i, self.L_e, self.params['beta'], self.params['geo'], T_HP_h)
    
        return (Q_dot, T_wc_o, T_we_o, P_HP_h, T_HP_h, P_evap_ex, V_max_air, a_gas_fume, R_tube_tot, R_e_o, R_c_o)  


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
        self.T_h = np.zeros(self.N_col+1)
        self.P_h = np.zeros(self.N_col+1)
        
        self.T_c = np.zeros(self.N_col+1)
        self.h_c = np.zeros(self.N_col+1)
        
        self.Q_dot_row = np.zeros(self.N_col)
        Q_dot_tube = np.zeros(self.N_col)
        self.T_HP_h = np.zeros(self.N_col)
        
        self.T_wc_o = np.zeros(self.N_col)
        self.T_we_o = np.zeros(self.N_col)
        
        self.P_HP_h = np.zeros(self.N_col)
                
        self.R_tot_tube = np.zeros(self.N_col)  
        self.R_c_o = np.zeros(self.N_col)  
        self.R_e_o = np.zeros(self.N_col)  
        
        "2) External fluids data"
        
        "2.1) Flue gases"
    
        A_chimney = (np.pi/4)*self.params['D_chimney']**2
    
        # Thermo state related :
        self.AS_H.update(CP.HmassP_INPUTS, self.su_H.h, self.su_H.p)
        cp_h_in = self.AS_H.cpmass()
    
        # Flow rate related
        V_dot_h = self.su_H.m_dot / self.su_H.D
        u_h_chimney = V_dot_h/A_chimney
    
        "Nozzle computation - HTX Inlet"
        (u_h, self.T_h[0], self.P_h[0]) = nozzle(self.su_H.m_dot, self.su_H.T, cp_h_in, u_h_chimney, self.su_H.p, A_chimney, self.A_casing)
            
        DP_gas_nozzle = self.su_H.p - self.P_h[0]
        DP_gas_nozzle_vect = DP_gas_nozzle
        
        "2.3) Evaporating Working Fluid"     
        self.AS_C.update(CP.HmassP_INPUTS, h_wf_out_guess, self.su_C.p)
        self.T_c[0] = self.AS_C.T()
            
        H_casing = self.L_c + 0.02  # [m] height of oil divergent
        self.A_casing = H_casing*self.W_casing
    
        rho_wf_in = self.AS_C.rhomass()  # [kg/m^3] : oil density at 110°C
        V_dot_wf = self.su_C.m_dot/rho_wf_in  # [m^3/s] : Oil volumetric flow rate
    
        u_wf_in = V_dot_wf/self.A_casing  # [m/s] : Oil velocity in the casing
        
        "3) Thermosyphon Simulation"
            
        vect_init = (7*101325, 7*101325, 1000, 1000, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 635.15, 0, 513.15, 635.15, 101325,50) # First initial conditions
        # P_HP_c, P_HP_h, Q_dot_axial, Q_dot_radial, T_HP_c, T_HP_h, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
        DP_h = 0
        
        for i in range(self.N_col):
            R_tot_tube_mean = 0
            R_e_o_mean = 0
            R_c_o_mean = 0
            
            (Q_dot_tube[i], self.T_wc_o[i], self.T_we_o[i], self.P_HP_h[i], self.T_HP_h[i], self.P_h[i+1], V_max_air, a_gf, self.R_tot_tube[i], self.R_e_o[i], self.R_c_o[i]) = self.thermosyphon_model(vect_init, self.T_h[i], self.T_c[i], self.P_h[i], self.su_C.p, u_h, self.su_H.m_dot, u_wf_in, self.su_C.m_dot, self.h_c)  # call of the model for one tube
                        
            DP_h = DP_h + self.P_h[i]-self.P_h[i+1]
    
            self.Q_dot_row[i] = Q_dot_tube[i]*self.N_tubes_row
    
            self.AS_H.update(CP.PT_INPUTS, self.P_h[i], self.T_h[i])
    
            cp_h = self.AS_H.cpmass()
            self.T_h[i+1] = self.T_h[i] - self.Q_dot_row[i]/(cp_h*self.su_H.m_dot)
            
            if i == 0:
                self.AS_C.update(CP.PT_INPUTS, self.su_C.p, self.T_c[i])
                self.h_c[i] = self.AS_C.hmass()
                
            self.h_c[i+1] = self.h_c[i] - self.Q_dot_row[i]/self.su_C.m_dot
            
            try:
                self.AS_C.update(CP.HmassP_INPUTS, self.h_c[i+1], self.su_C.p)
                self.T_c[i+1] = self.AS_C.T()
            except:
                T_current = CP.PropsSI('T', 'H', self.h_c[i+1]-1000, 'P', self.su_C.p, self.su_C.fluid)
                # self.AS_C.update(CP.HmassP_INPUTS, self.h_c[i+1]-10, self.su_C.p)
                self.T_c[i+1] = T_current
                
            R_tot_tube_mean = R_tot_tube_mean + self.R_tot_tube[i]
            R_e_o_mean = R_e_o_mean + self.R_e_o[i]
            R_c_o_mean = R_c_o_mean + self.R_c_o[i]
            
            vect_init = (self.P_HP_h[i],self.P_HP_h[i],1000,Q_dot_tube[i],self.T_HP_h[i],self.T_HP_h[i],self.T_wc_o[i],self.T_wc_o[i],self.T_wc_o[i],self.T_we_o[i],self.T_we_o[i],self.T_we_o[i], 0, self.T_c[i+1], self.T_h[i+1], self.P_h[i+1], self.R_tot_tube[i])
                        # P_HP_c, P_HP_h, Q_dot_axial, Q_dot_radial, T_HP_c, T_HP_h, T_wc, T_wc_i, T_wc_o, T_we, T_we_i, T_we_o, DP_v, T_c_o_ex, T_e_o_ex, P_evap_ex, R_tube_tot
                
        self.AS_H.update(CP.PT_INPUTS, self.P_h[-1], self.T_h[-1])
        rho_h_out = self.AS_H.rhomass()
        
        V_dot_h_out = self.su_H.m_dot/rho_h_out
        u_h_out = V_dot_h_out/self.A_casing
        
        "Nozzle computation - HTX Inlet"
        (u_h_out, T_h_out, P_h_out) = nozzle(self.su_H.m_dot, self.T_h[-1], cp_h, u_h_out, self.P_h[-1], self.A_casing, A_chimney)
        
        # print((P_gas_fume_out - P_gas_fume[-1]))
        
        DP_gas_nozzle_vect = DP_gas_nozzle_vect + (P_h_out - self.P_h[-1])
            
        return self.T_c[-1], self.h_c[-1], self.h_c[0], P_h_out, T_h_out
           
#%%

    def solve(self, Tc_out_guess_1, Tc_out_guess_2, res_tol = 1e-3, n_it_max = 20, print_flag = 0):        
        
        # ---- Initialize -------------
        
        self.check_calculable()
        self.check_parametrized()
        
        if self.calculable != True:
            print("Supply fluid not fully known.")
            return
      
        if self.parametrized != True:
            print("Heat Exchanger Parameters not fully known.")
            return    
      
        self.AS_H  = CP.AbstractState('HEOS', self.su_H.fluid)
        self.AS_C  = CP.AbstractState('HEOS', self.su_C.fluid)
        self.AS_HP = CP.AbstractState('HEOS', self.params['HP_fluid']) 
        
        if Tc_out_guess_1 > Tc_out_guess_2:
            temp = Tc_out_guess_1
            Tc_out_guess_2 = Tc_out_guess_1
            Tc_out_guess_2 = temp
        
        if self.params['arrang'] == "inline":
            self.set_parameters(arrang = "Inline")
    
        if self.params['arrang'] == "staggered":
            self.set_parameters(arrang = "Staggered")
        
        [self.f_0, self.f_1, self.f_2, self.f_3] = fg_radiative_htc_corr()
        
        n_it = 0
        self.res = 10000
            
        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        
        # ---- Guess cold outlet state -------------
        
        Tc_out_guess = (Tc_out_guess_1 + Tc_out_guess_2)/2
        
        self.AS_C.update(CP.PT_INPUTS, self.su_C.p, Tc_out_guess)
        hc_out_guess = self.AS_C.hmass() 
    
        # ---- Bisection Iteration to find cold outlet -------------
    
        while abs(self.res) > res_tol:
                    
            if n_it >= n_it_max:
                print("No convergence in the fixed number of iterations")
                return
                
            (Tc_in, hc_it, hc_out, Ph_out, Th_out) = self.System(hc_out_guess, Tc_out_guess)
            
            "Compute residuals"
                        
            self.res = (self.su_C.h - hc_it)/self.su_C.h

            self.AS_C.update(CP.HmassP_INPUTS, hc_out, self.su_C.p)
    
            if self.res > 0: # h_real > h_est => T_in > T_in,est => T_out > T_out_est => Rise the value of T_out_guess_1            
                Tc_out_guess_1 = self.AS_C.T()
            else:
                Tc_out_guess_2 = self.AS_C.T()
                
            Tc_out_guess = (Tc_out_guess_1 + Tc_out_guess_2)/2
            
            self.AS_C.update(CP.PT_INPUTS, self.su_C.p, Tc_out_guess)

            hc_out_guess = self.AS_C.hmass()
    
            n_it = n_it + 1

            if print_flag == 1:
                print("-------------------------")
                print("Iteration : ",n_it, " / ", n_it_max)
                print("T_in_guess :", Tc_in)
                print("T_out_guess :", Tc_out_guess)
                print("Tc_out_guesses : [", Tc_out_guess_1, " ; ", Tc_out_guess_2, "]")
        
        self.update_connectors(Th_out, Ph_out, Tc_in)
        
        if abs(self.res) < res_tol and print_flag == 1:
                print("-------------------------")
                print("Success !")
                print("-------------------------")                
                print("Iteration : ",n_it, " / ", n_it_max)
                print("T_in_input :", self.su_C.T)
                print("T_in_final :", Tc_in)
                print("T_out_final :", Tc_out_guess)
                print("-------------------------")    
    #%%  
    
    def update_connectors(self, Th_out, Ph_out, Tc_in):
    
        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_m_dot(self.su_H.m_dot)  # Example mass flow rate [kg]
        self.ex_H.set_T(Th_out) # Example temperature [K]
        self.ex_H.set_p(Ph_out)  # Example Pressure [Pa]
    
        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_m_dot(self.su_C.m_dot)  # Example mass flow rate [kg]
        self.ex_C.set_T(Tc_in) # Example temperature [K]
        self.ex_C.set_p(self.su_C.p)  # Example Pressure [Pa]
    
    def plot_disc(self):
        
        plt.plot(self.T_c, 'b', label='Cold fluid')
        plt.plot(self.T_h, 'r', label='Hot fluid')
        
        plt.grid()
        plt.xlabel("Row number")
        plt.ylabel("Temperature [K]")
        plt.legend(loc='best')  # automatically picks a good position
        
        plt.show()

        
        