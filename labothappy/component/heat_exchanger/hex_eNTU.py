# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 14:32:39 2024

@author: Basile
"""

from component.base_component import BaseComponent

# from correlations.convection.pipe_htc import gnielinski_pipe_htc 
from correlations.heat_exchanger.e_NTU import e_NTU
from correlations.convection.plate_htc import han_BPHEX_DP, water_plate_HTC, martin_BPHEX_HTC, muley_manglik_BPHEX_HTC, han_boiling_BPHEX_HTC, han_cond_BPHEX_HTC, thonon_plate_HTC, kumar_plate_HTC, martin_holger_plate_HTC, amalfi_plate_HTC, shah_condensation_plate_HTC
from correlations.convection.pipe_htc import gnielinski_pipe_htc, boiling_curve, horizontal_tube_internal_condensation, horizontal_flow_boiling, flow_boiling_gungor_winterton, Liu_sCO2, Cheng_sCO2, thome_condensation, choi_boiling
from correlations.convection.shell_and_tube_htc import shell_bell_delaware_htc, shell_htc_kern
from correlations.convection.tube_bank_htc import ext_tube_film_condens
from correlations.convection.fins_htc import htc_tube_and_fins
from correlations.convection.printed_circuit_htc import PCHE_Lee, PCHE_conv


from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from scipy.interpolate import interp1d
import numpy as np


from toolbox.geometries.heat_exchanger.c_geometry_HXs_Zorlu import Zorlu_HXs

class HexeNTU(BaseComponent):
    
    """
    **Component**: Heat Exchanger

    **Model**: ε-NTU (Effectiveness - Number of Transfer Units) method.

    **Description**:

        This component models a heat exchanger using the ε-NTU method, a widely used approach for estimating heat transfer performance 
        in steady-state conditions when outlet temperatures are not known a priori. It calculates heat transfer based on fluid properties, 
        flow configuration, and geometry using thermal resistances and heat transfer correlations.

        The model is applicable to various geometries (e.g., pipe-type, plate-type) and requires fluid and geometric properties. 
        The thermal effectiveness is computed via an external ε-NTU correlation, which supports multiple flow configurations 
        (e.g., CounterFlow, ParallelFlow, CrossFlow).

    **Assumptions**:

        - Steady-state operation.
        - No heat loss to the environment.
        - No pressure drop considered (isenthalpic mixing assumed).
        - Thermophysical properties are evaluated at average temperatures.
        - No phase change within the exchanger.

    **Connectors**:

        su_H (MassConnector): Hot fluid inlet connector.
        su_C (MassConnector): Cold fluid inlet connector.
        ex_H (MassConnector): Hot fluid outlet connector.
        ex_C (MassConnector): Cold fluid outlet connector.
        Q_hex (HeatConnector): Connector for total heat transfer rate.

    **Parameters**:
            
        Flow_Type : Flow configuration of the fluid ('CounterFlow', 'CrossFlow', 'Shell&Tube', 'ParallelFlow') [-]
        hex_type: Type of heat exchanger ('Shell&Tube', 'Tube&Fins', 'Plate')
        A_htx: Total heat exchange area [m²]
        L_HTX: Length of the heat exchanger [m]
        V_HTX: Volume of the heat exchanger [m³]
        A_canal_H: Cross-sectional area of hot fluid channels [m²]
        A_canal_C: Cross-sectional area of cold fluid channels [m²]
        D_h_C: Hydraulic diameter of the cold fluid channel [m]
        D_h_H: Hydraulic diameter of the hot fluid channel [m]
        k_plate: Thermal conductivity of the separating plate [W/m.K]
        t_plate: Thickness of the separating plate [m]
        n_plates: Number of plates [-]
        co_pitch: Plate corrugation pitch [m]
        chevron_angle: Plate chevron angle [degrees]
        fouling: Fouling resistance [m².K/W]

    **Inputs**:

        T_su_H: Hot fluid inlet temperature [K]
        P_su_H: Hot fluid inlet pressure [Pa]
        h_su_H: Hot fluid inlet enthalpy [J/kg]
        fluid_H: Hot fluid identifier [-]
        m_dot_H: Hot fluid mass flow rate [kg/s]

        T_su_C: Cold fluid inlet temperature [K]
        P_su_C: Cold fluid inlet pressure [Pa]
        h_su_C: Cold fluid inlet enthalpy [J/kg]
        fluid_C: Cold fluid identifier [-]
        m_dot_C: Cold fluid mass flow rate [kg/s]

    **Outputs**:

        h_ex_H: Hot fluid outlet enthalpy [J/kg]
        P_ex_H: Hot fluid outlet pressure [Pa]
        h_ex_C: Cold fluid outlet enthalpy [J/kg]
        p_ex_C: Cold fluid outlet pressure [Pa]
        Q_dot: Heat transfer rate [W]
"""
    
    class H:
        def __init__(self):
            self.Correlation_1phase = None
            self.Correlation_2phase = None
            self.HeatExchange_Correlation = None  
            self.PressureDrop_Correlation = None
            self.h_liq = None
            self.h_vap = None
            self.h_twophase = None
            self.h_vapwet = None
            self.h_tpdryout = None
            self.h_transcrit = None
            self.f_dp = None

    class C:
        def __init__(self):
            self.Correlation_1phase = None
            self.Correlation_2phase = None
            self.HeatExchange_Correlation = None
            self.PressureDrop_Correlation = None
            self.h_liq = None
            self.h_vap = None
            self.h_twophase = None
            self.h_vapwet = None
            self.h_tpdryout = None
            self.h_transcrit = None
            self.f_dp = None
        
        
        
    def __init__(self, hex_type):
        super().__init__()
        self.su_H = MassConnector()
        self.su_C = MassConnector()

        self.ex_H = MassConnector()
        self.ex_C = MassConnector() # Mass_connector

        self.Q_hex = HeatConnector()
        
        if hex_type == 'Plate' or hex_type == 'Shell&Tube' or hex_type == 'Tube&Fins':
            self.hex_type = hex_type
        else:
            raise ValueError("Heat exchanger types implemented for this model are : 'Plate', 'Shell&Tube', 'Tube&Fins'.")
        self.hex_type = hex_type
        
        self.H = self.H()
        self.C = self.C()
        
        # self.H = HexeNTU.H()
        # self.C = HexeNTU.C()


    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su_H', 'T_su_H', 'm_dot_H', 'fluid_H', 'P_su_C', 'T_su_C', 'm_dot_C', 'fluid_C']

    def get_required_parameters(self):
        """ Returns the list of required parameters to describe the geometry and physical configuration """
        # HEX_TYPE = [self.hex_type]
        self.params['hex_type'] = self.hex_type
        
        if self.hex_type == 'Plate':
            geometry_parameters = ['A_htx', 'L_HTX', 'V_HTX', 'Flow_Type',
                               'A_canal_H', 'A_canal_C', 'D_h_C', 'D_h_H',
                               'k_plate', 't_plate', 'n_plates',
                               'co_pitch', 'chevron_angle', 'fouling']
        elif self.hex_type == 'Shell&Tube':
            geometry_parameters = ['A_htx'] 
            # to be continued
        elif self.hex_type == 'Tube&Fins':
            geometry_parameters = ['A_htx'] 
            # to be continued
        
        
        elif self.hex_type == 'Shell&Tube':
                
            if self.H.Correlation_1phase == "Shell_Bell_Delaware_HTC" or self.C.Correlation_1phase == "Shell_Bell_Delaware_HTC":

                geometry_parameters = ['Baffle_cut', 'D_OTL', 'N_strips', 'Shell_ID', 'Tube_L', 'Tube_OD', 'Tube_pass',
                                    'Tube_t', 'Tubesheet_t', 'central_spacing', 'clear_BS', 'clear_TB',
                                    'cross_passes', 'foul_s', 'foul_t', 'inlet_spacing', 'n_series', 'n_parallel', 
                                    'n_tubes', 'outlet_spacing', 'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']

            if self.H.Correlation_1phase == "Shell_Kern_HTC" or self.C.Correlation_1phase == "Shell_Kern_HTC":

                geometry_parameters = ['Baffle_cut', 'Shell_ID', 'Tube_L', 'Tube_OD', 'Tube_pass','Tube_t', 'central_spacing',
                                    'cross_passes', 'foul_s', 'foul_t', 'n_series', 'n_parallel', 'n_tubes', 'pitch_ratio', 
                                    'tube_cond', 'tube_layout', 'Shell_Side']
        
        elif self.hex_type == 'Tube&Fins':
            
            geometry_parameters = ['A_flow', 'Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type',
                                    'Finned_tube_flag', 'Tube_L', 'Tube_OD',
                                    'Tube_cond', 'Tube_t', 'fouling', 'h', 'k_fin',
                                    'Tube_pass', 'n_rows', 'n_series', 'n_parallel', 'n_tubes', 'pitch', 'pitch_ratio', 'tube_arrang',
                                    'w','Fin_Side']
            
        else:
            print("Unrecognised HEX type. hex_type shall be equal to 'Shell&Tube', 'Tube&Fins' or 'Plate'.")
        
        # return HEX_TYPE + geometry_parameters
        return geometry_parameters
    
    
    # changer et mettre une seule corrélation demandée
    
    # Setting heat transfer coefficient, user defined or from a correlation
    def set_htc(self, htc_type = None, Corr_H = None, Corr_C = None, UD_H_HTC = None, UD_C_HTC = None):
        """
        General Parameters : 
            
            - htc_type : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - Corr_H   : Correlations for hot side
            - Corr_C   : Correlations for cold side
            - UD_H_HTC : User-Defined HTC for hot side
            - UD_C_HTC : User-Defined HTC for cold side
            
        """
        
        
        self.params['htc_type'] = htc_type
        
        self.UD_C_HTC = UD_C_HTC
        self.UD_H_HTC = UD_H_HTC
        self.Corr_H = Corr_H
        self.Corr_C = Corr_C
        
        if htc_type == "User-Defined":
            
            self.H.HeatExchange_Correlation = "User-Defined"
            self.C.HeatExchange_Correlation = "User-Defined"
            
        
        elif htc_type == "Correlation":
            # Type 
            self.H.HeatExchange_Correlation = "Correlation"
            self.C.HeatExchange_Correlation = "Correlation"
        
            if self.hex_type == 'Plate' or self.hex_type == 'Shell&Tube' or self.hex_type == 'Tube&Fins':
                
                self.Corr_C = Corr_C
                self.Corr_H = Corr_H

            if self.C.Correlation_2phase == "Boiling_curve": # Compute the fluid boiling curve beforehand
                try:
                    self.AS_C = CP.AbstractState("BICUBIC&HEOS", self.su_C.fluid)   
                    self.AS_C.update(CP.PQ_INPUTS, self.su_C.p, 0)
                    
                    T_sat = self.AS_C.T()
                    (h_boil, DT_vect) = boiling_curve(self.params['Tube_OD'], self.su_C.fluid, T_sat, self.su_C.p)
                    self.C_f_boiling = interp1d(DT_vect,h_boil)
                except:
                    self.C_f_boiling = interp1d([0,10000],[20000,20000])
    
    
    
    
    # def set_DP(self, DP_type = None, Corr_H = None, Corr_C = None, UD_H_DP = None, UD_C_DP = None):
    def set_DP(self, UD_H_DP = None, UD_C_DP = None):
        """
        General Parameters : 
            
            # - htc_type : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            # - Corr_H   : Correlations for hot side
            # - Corr_C   : Correlations for cold side
            - UD_H_HTC : User-Defined HTC for hot side
            - UD_C_HTC : User-Defined HTC for cold side
            
        """

        # self.params['DP_type'] = DP_type
        
        self.params['UD_H_DP'] = UD_H_DP
        self.params['UD_C_DP'] = UD_C_DP
        
        self.C.Correlation_DP = {}
        self.H.Correlation_DP = {}
        
        # if self.params['DP_type'] is None:
        if self.params['UD_C_DP'] is None and self.params['UD_H_DP'] is None:

            self.H.DP_val = None
            self.C.DP_val = None

        # elif self.params['DP_type'] == 'User-Defined':
        elif self.params['UD_C_DP'] is not None and self.params['UD_H_DP'] is not None:

            self.H.DP_val = UD_H_DP
            self.C.DP_val = UD_C_DP            

        
        # elif self.params['DP_type'] == "Correlation":
            
        #     self.H.Correlation_DP["1P"] = Corr_H["1P"]
        #     if "2P" in Corr_H:
        #         self.H.Correlation_DP["2P"] = Corr_H["2P"]
        #     else:
        #         self.H.Correlation_DP["2P"] = None
                
        #     if "SC" in Corr_H:
        #         self.H.Correlation_DP["SC"] = Corr_H["SC"]
        #     else:
        #         self.H.Correlation_DP["SC"] = None
            
        #     self.C.Correlation_DP["1P"] = Corr_C["1P"]
        #     if "2P" in Corr_C:
        #         self.C.Correlation_DP["2P"] = Corr_C["2P"]
        #     else:
        #         self.C.Correlation_DP["2P"] = None
        
        #     if "SC" in Corr_C:
        #         self.C.Correlation_DP["SC"] = Corr_C["SC"]
        #     else:
        #         self.C.Correlation_DP["SC"] = None
        else:
            # raise ValueError("Error in pressure drop setting. DP_type entry shall either not be set or be set to: \n - 'User-Defined' \n - 'Correlation'")
            raise ValueError("Error in pressure drop setting. 'UD_H_DP' and 'UD_C_DP' shall be set.")
            
        return

    def compute_htc_H(self, k, Pr, T, p, T_w, mu, mu_w, G):
        
        htc_type = self.params.get('htc_type')
        
        
        if htc_type == None and self.hex_type == 'Plate':
            gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['H_Dh'], self.params['l'])
            
        elif htc_type == None and (self.hex_type == 'Shell&Tube' or self.hex_type == 'Tube&Fins'):
            raise ValueError("The htc_type was either not set thanks to 'set_htc' or was set but is not equal to 'User-Defined' or 'Correlation'.")
             
            
        else:
            if self.C.HeatExchange_Correlation == "User-Defined" and self.H.HeatExchange_Correlation == "User-Defined" :
                h_conv = self.UD_H_HTC
                 
            elif self.C.HeatExchange_Correlation == "Correlation" and self.H.HeatExchange_Correlation == "Correlation" :
                if self.hex_type == 'Plate' and (self.H_su.fluid == 'water' or self.H_su.fluid == 'Water'):
                    h_conv = water_plate_HTC(mu, Pr, k, G, self.params['C_Dh'])
                
                
                elif htc_type == "Gnielinski":
                    if self.hex_type == 'Plate':
                        h_conv = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['H_Dh'], self.params['l'])[0]
                    elif self.hex_type == 'Shell&Tube':
                        h_conv = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass'])[0] 
                    elif self.hex_type == 'Tube&Fins':
                        h_conv = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass'])[0] 
                        
                        
                elif htc_type == "Shell_Bell_Delaware_HTC":
                    h_conv = shell_bell_delaware_htc(self.su_H.m_dot, T, T_w, p, self.H_su.fluid, self.params)
                elif htc_type == 'Shell_Kern_HTC':
                    h_conv = shell_htc_kern(self.su_H.m_dot, T_w, T, p, self.AS_H, self.params)     
                elif htc_type == 'Tube_And_Fins':
                    h_conv = htc_tube_and_fins(self.H_su.fluid, self.params, p, self.su_H.h , self.su_H.m_dot, self.params['Fin_type'])[0]
                elif htc_type == 'water_plate_HTC':
                    h_conv = water_plate_HTC(mu, Pr, k, G, self.params['H_Dh'])
                elif htc_type  == 'martin_holger_plate_HTC':
                    h_conv = martin_holger_plate_HTC(mu, Pr, k, self.su_H.m_dot, self.params['H_n_canals'], T, p, self.H_su.fluid, self.params['H_Dh'], self.params['l'], self.params['w'], self.params['amplitude'], self.params['chevron_angle'])            
                
                
                
            elif self.C.HeatExchange_Correlation == None or self.H.HeatExchange_Correlation == None :
                raise ValueError("The htc_type was set but not the rest. Please set 'UD_H_DP' and 'UD_C_DP' if htc_type='User-defined' or 'Corr_C' and 'Corr_H' if htc_type='Correlation'.")
        
        return h_conv
            
    

    def compute_htc_C(self, k, Pr, T, p, T_w, mu, mu_w, G):
        
        htc_type = self.params.get('htc_type')
        
        
        if htc_type == None and self.hex_type == 'Plate':
            gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['C_Dh'], self.params['l'])
            
        elif htc_type == None and (self.hex_type == 'Shell&Tube' or self.hex_type == 'Tube&Fins'):
            raise ValueError("The htc_type was either not set thanks to 'set_htc' or was set but is not equal to 'User-Defined' or 'Correlation'.")
             
            
        else:
            if self.C.HeatExchange_Correlation == "User-Defined" and self.H.HeatExchange_Correlation == "User-Defined" :
                h_conv = self.UD_C_HTC
                 
            elif self.C.HeatExchange_Correlation == "Correlation" and self.H.HeatExchange_Correlation == "Correlation" :
                if self.hex_type == 'Plate' and (self.H_su.fluid == 'water' or self.C_su.fluid == 'Water'):
                    h_conv = water_plate_HTC(mu, Pr, k, G, self.params['C_Dh'])
                
                
                elif htc_type == "Gnielinski":
                    if self.hex_type == 'Plate':
                        h_conv = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['C_Dh'], self.params['l'])[0]
                    elif self.hex_type == 'Shell&Tube':
                        h_conv = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass'])[0] 
                    elif self.hex_type == 'Tube&Fins':
                        h_conv = gnielinski_pipe_htc(mu, Pr, mu_w, k, G, self.params['Tube_OD']-2*self.params['Tube_t'], self.params['Tube_L']*self.params['Tube_pass'])[0] 
                        
                        
                elif htc_type == "Shell_Bell_Delaware_HTC":
                    h_conv = shell_bell_delaware_htc(self.su_C.m_dot, T, T_w, p, self.C_su.fluid, self.params)
                elif htc_type == 'Shell_Kern_HTC':
                    h_conv = shell_htc_kern(self.su_C.m_dot, T_w, T, p, self.AS_C, self.params)     
                elif htc_type == 'Tube_And_Fins':
                    h_conv = htc_tube_and_fins(self.C_su.fluid, self.params, p, self.su_C.h , self.su_C.m_dot, self.params['Fin_type'])[0]
                elif htc_type == 'water_plate_HTC':
                    h_conv = water_plate_HTC(mu, Pr, k, G, self.params['C_Dh'])
                elif htc_type  == 'martin_holger_plate_HTC':
                    h_conv = martin_holger_plate_HTC(mu, Pr, k, self.su_C.m_dot, self.params['C_n_canals'], T, p, self.C_su.fluid, self.params['C_Dh'], self.params['l'], self.params['w'], self.params['amplitude'], self.params['chevron_angle'])            
                
                
                
            elif self.C.HeatExchange_Correlation == None or self.C.HeatExchange_Correlation == None :
                raise ValueError("The htc_type was set but not the rest. Please set 'UD_H_DP' and 'UD_C_DP' if htc_type='User-defined' or 'Corr_C' and 'Corr_H' if htc_type='Correlation'.")
        
        return h_conv
    
    
    
    def G_h_c_computation(self):
        if self.hex_type == 'Plate':          
            G_c = (self.su_C.m_dot/self.params['C_n_canals'])/self.params['C_CS']
            G_h = (self.su_H.m_dot/self.params['H_n_canals'])/self.params['H_CS']
        elif self.hex_type == 'Shell&Tube':
            A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            G_h = (self.params["Tube_pass"]/self.params["n_parallel"])*self.su_H.m_dot/(A_in_one_tube*self.params['n_tubes'])
            G_c = (self.params["Tube_pass"]/self.params["n_parallel"])*self.su_C.m_dot/(A_in_one_tube*self.params['n_tubes'])
        elif self.hex_type == 'Tube&Fins':
            A_in_one_tube = np.pi*((self.params['Tube_OD']-2*self.params['Tube_t'])/2)**2
            G_h = (self.params["Tube_pass"]/self.params["n_parallel"])*(self.su_H.m_dot/self.params['n_tubes'])/A_in_one_tube
            G_c = (self.params["Tube_pass"]/self.params["n_parallel"])*(self.su_C.m_dot/self.params['n_tubes'])/A_in_one_tube
        elif self.hex_type == 'PCHE': 
            A_in_one_channel = np.pi*(self.params['D_c']**2)/8
            G_h = self.su_H.m_dot/(A_in_one_channel*self.params['N_c']*self.params['N_p']*(self.params['R_p']/(1+self.params['R_p'])))
            G_c = self.su_C.m_dot/(A_in_one_channel*self.params['N_c']*self.params['N_p']*(1/(1+self.params['R_p'])))
            
        return G_c, G_h
    
    

    

    def solve(self):
        self.check_calculable()
        self.check_parametrized()
        
        self.AS_H = CP.AbstractState('HEOS', self.su_H.fluid)
        self.AS_C = CP.AbstractState('HEOS', self.su_C.fluid)
        
        
        if self.calculable and self.parametrized:
            
            # Detect Phase change
            # self.detect_phase_change()

            self.AS_H.update(CP.HmassP_INPUTS, self.su_H.h, self.su_H.p)
            cp_h = self.AS_H.cpmass()
            
            self.AS_C.update(CP.HmassP_INPUTS, self.su_H.h, self.su_C.p)
            cp_c = self.AS_C.cpmass()

            
            
            C_h = cp_h*self.su_H.m_dot #Heat capacity rate
            C_c = cp_c*self.su_C.m_dot
            
            C_min = min(C_h, C_c)
            C_max = max(C_h, C_c)
            C_r = C_min/C_max # Heat capacity ratio
                        
            # Calcul de NTU
            T_w = (self.su_H.T + self.su_C.T)/2
            
            self.AS_H.update(CP.PT_INPUTS, self.su_H.p, T_w)
            mu_h_w = self.AS_H.viscosity()
            
            
            self.AS_H.update(CP.PT_INPUTS, self.su_C.p, T_w)
            mu_c_w = self.AS_H.viscosity()
            
        
            self.AS_H.update(CP.HmassP_INPUTS, self.su_H.h, self.su_H.p)
            mu_h = self.AS_H.viscosity()
            Pr_h = self.AS_H.Prandtl()
            k_h = self.AS_H.conductivity()
            
            self.AS_C.update(CP.HmassP_INPUTS, self.su_C.h, self.su_C.p)
            mu_c = self.AS_C.viscosity()
            Pr_c = self.AS_C.Prandtl()
            k_c = self.AS_C.conductivity()

            
            G_c, G_h = self.G_h_c_computation()
                        
            
            
            h_h = self.compute_htc_H(k_h, Pr_h, self.su_H.T, self.su_H.p, T_w, mu_h, mu_h_w, G_h)
            h_c = self.compute_htc_C(k_h, Pr_h, self.su_C.T, self.su_C.p, T_w, mu_h, mu_h_w, G_h)

            
            # h_c = gnielinski_pipe_htc(mu_c, Pr_c, Pr_c, k_c, G_c, self.params['D_h_C'], self.params['L_HTX'])[0]

            AU = (1/(self.params['A_htx']*h_h) + 1/(self.params['A_htx']*h_c) + self.params['t_plate']/(self.params['k_plate']*self.params['A_htx']) + self.params['fouling']/self.params['A_htx'])**(-1)         

            NTU = AU/C_min
                        
            # --- Calculate effectiveness from NTU correlation ---
            eps = e_NTU(NTU, C_r, self.params)

                        
            # --- Estimate maximum heat transfer Q(ideal case with infinite area) ---
            # h_c_Th = PropsSI('H','T',self.su_hot.T,'P',self.su_cold.p,self.su_cold.fluid)
            # h_h_Tc = PropsSI('H','T',self.su_cold.T,'P',self.su_hot.p,self.su_hot.fluid)
            
            self.AS_C.update(CP.PT_INPUTS, self.su_C.p, self.su_H.T)
            h_c_Th = self.AS_C.hmass()
            
            self.AS_H.update(CP.PT_INPUTS, self.su_H.p, self.su_C.T)
            h_h_Tc = self.AS_H.hmass()

            
            # Special cases for incompressibles
            
            if "INCOMP" not in self.su_C.fluid:
                
                self.AS_C.update(CP.PQ_INPUTS, self.su_C.p, 0)
                h_l_cold = self.AS_C.hmass()
                
                self.AS_C.update(CP.PQ_INPUTS, self.su_C.p, 1)
                h_v_cold = self.AS_C.hmass()
                
                DH_pc_c = h_v_cold - h_l_cold
            else:
                DH_pc_c = 0

            
            if "INCOMP" not in self.su_H.fluid:
                
                self.AS_H.update(CP.PQ_INPUTS, self.su_H.p, 0)
                h_l_hot = self.AS_H.hmass()
                
                self.AS_H.update(CP.PQ_INPUTS, self.su_H.p, 1)
                h_v_hot = self.AS_H.hmass()
                
                DH_pc_h = h_v_hot - h_l_hot

            else:
                DH_pc_h = 0
            
            
            
            Qmax_c = self.su_C.m_dot*((h_c_Th - self.su_C.h))
            Qmax_h = self.su_H.m_dot*((self.su_H.h - h_h_Tc))
            
            Qmax = min(Qmax_c, Qmax_h)
            
            Q = eps*Qmax  # Actual heat exchanged
            
            

            # Hypothesis: isenthalpic pressure drop 
            """ Setting the DP """
            
            UD_C_DP = self.params.get('UD_C_DP')
            UD_H_DP = self.params.get('UD_H_DP')
            
            if UD_C_DP is None and UD_H_DP is None:
                self.DP_c = 0
                self.ex_C.p = self.su_C.p - self.DP_c  
                
                self.DP_h = 0
                self.ex_H.p = self.su_H.p - self.DP_h 
                
            # elif self.params['DP_type'] == "User-Defined":
            # elif 'UD_H_DP' in self.params and 'UD_C_DP' in self.params:
            elif UD_C_DP is not None and UD_H_DP is not None:
                self.DP_h = self.H.DP_val
                self.ex_H.p = self.su_H.p - self.DP_h
                    
                self.DP_c = self.C.DP_val        
                self.ex_C.p = self.su_C.p - self.DP_c
            # Correlation-based DP are not considered yet



            # --- Set exhaust states (new enthalpies) and link to connectors ---
            self.ex_H.set_properties(H = self.su_H.h - Q/self.su_H.m_dot, fluid = self.su_H.fluid, m_dot = self.su_H.m_dot, P = self.ex_H.p)
            self.ex_C.set_properties(H = self.su_C.h + Q/self.su_C.m_dot, fluid = self.su_C.fluid, m_dot = self.su_C.m_dot, P = self.ex_C.p)

            self.Q_hex.set_Q_dot(Q)

            self.defined = True
        
        else:
            if not self.calculable:
                print("Input of the component not completely known. Required inputs:")
                for input in self.get_required_inputs():
                    if input not in self.inputs:
                        print(f"  - {input}")
                        
            if not self.parametrized:
                print("Parameters of the component not completely known. Required parameters:")
                for param in self.get_required_parameters():
                    if param not in self.params:
                        print(f"  - {param}")
    
    def print_results(self):
        if self.defined:
            print("=== Heat Exchanger Results ===")
            print(f"  - H_ex: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
            print(f"  - C_ex: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
            print(f"  - Q_dot: {self.Q_hex.Q_dot}")

        else:
            print("Heat Exchanger component is not defined. Ensure it is solved first.")
            
            
# if __name__ == "__main__":
    # geom_obj = Zorlu_HXs()
    # geom_obj.set_parameters("ORC_recuperator")
    
    # HX = HexeNTU()
    # HX.set_parameters(**geom_obj.geom)
    
    # A_c, A_h = HX.A_in(hex_type = 'Plate')
    
    
    