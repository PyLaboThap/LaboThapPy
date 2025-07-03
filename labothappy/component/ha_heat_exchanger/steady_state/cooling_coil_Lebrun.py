# -*- coding: utf-8 -*-
"""
Created on Fri May 30 11:49:21 2025

@author: Samuel Gendebien
"""
import __init__

# External Toolbox 
import numpy as np
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import HAPropsSI

# Connectors
from connector.mass_connector import MassConnector
from connector.heat_connector import HeatConnector
from connector.humid_air_connector import HAConnector

# HTC correlations
from correlations.convection.fins_htc import htc_tube_and_fins
from correlations.convection.pipe_htc import gnielinski_pipe_htc
from correlations.convection.pipe_htc import horizontal_tube_internal_condensation
from correlations.convection.pipe_htc import horizontal_flow_boiling

#Pressure drop correlations
from correlations.pressure_drop.fins_DP import DP_tube_and_fins
from correlations.pressure_drop.pipe_DP import gnielinski_pipe_DP, Muller_Steinhagen_Heck_DP

# Component base frame
from component.base_component import BaseComponent

class CoolingCoilLebrun(BaseComponent):
    
    def __init__(self):
        """
        su : Supply - 'H' : hot
                    - 'C' : cold

        ex : Exhaust - 'H' : hot
                     - 'C' : cold
                    
        Q_dot : Heat connection to the ambient
        
        HTX_Type : Type of HTX - Tube and Fins
        """
        
        super().__init__()
        
        self.su_H = HAConnector()
        self.su_C = MassConnector()
        
        self.ex_H = HAConnector()
        self.ex_C = MassConnector() # Mass_connector
        
        self.Q_dot = HeatConnector()
        self.debug = 0

    #%%    
    
    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        """
        Hot side (humid air) required inputs : 
            
            - Hsu_T or Hsu_h : Hot supply temperature or enthalpy
            - Hsu_w or Hsu_RH or Hsu_dp or Hsu_wb : Hot Humidity ratio, Relative humidity, dew point or wet bulb
            - Hsu_p          : Hot supply pressure
            - Hsu_fluid      : Hot supply fluid
            - Hsu_m_dot      : Hot supply flow rate
            
        Cold side (usually brines) required inputs : 
            
            - Csu_T or Hsu_h : Cold supply temperature or enthalpy
            - Csu_p          : Cold supply pressure
            - Csu_fluid      : Cold supply fluid
            - Csu_m_dot      : Cold supply flow rate
        """
        self.sync_inputs()
        # Return a list of required inputs
        return['Hsu_T', 'Hsu_h', 'Hsu_T', 'Hsu_m_dot', 'Hsu_p', 'Hsu_w', 'Hsu_T', 'Hsu_m_dot', 'Csu_p', 'Csu_T', 'Csu_m_dot', 'Csu_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        # Hot Fluid
        if self.su_H.T is not None:
            self.inputs['Hsu_T'] = self.su_H.T
        elif self.su_H.h is not None:
            self.inputs['Hsu_h'] = self.su_H.h
            
        if self.su_H.w is not None:
            self.inputs['Hsu_w'] = self.su_H.w
        elif self.su_H.RH is not None:
            self.inputs['Hsu_RH'] = self.su_H.RH    
        elif self.su_H.Tdp is not None:
            self.inputs['Hsu_Tdp'] = self.su_H.Tdp
        elif self.su_H.Twb is not None:
            self.inputs['Hsu_Twb'] = self.su_H.Twb
            
        if self.su_H.p is not None:
            self.inputs['Hsu_p'] = self.su_H.p
        if self.su_H.m_dot is not None:
            self.inputs['Hsu_m_dot'] = self.su_H.m_dot
        if self.su_H.V_dot is not None:
            self.inputs['Hsu_V_dot'] = self.su_H.V_dot    
            
        # Cold Fluid                
        if self.su_C.T is not None:
            self.inputs['Csu_T'] = self.su_C.T
        elif self.su_C.h is not None:
            self.inputs['Csu_h'] = self.su_C.h
        if self.su_C.p is not None:
            self.inputs['Csu_p'] = self.su_C.p
        if self.su_C.fluid is not None:
            self.inputs['Csu_fluid'] = self.su_C.fluid
        if self.su_C.m_dot is not None:
            self.inputs['Csu_m_dot'] = self.su_C.m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        # Hot Fluid
        if 'Hsu_T' in self.inputs:
            self.su_H.set_properties(T=self.inputs['Hsu_T'])
        elif 'Hsu_h' in self.inputs:
            self.su_H.set_properties(h=self.inputs['Hsu_h'])
        
        if 'Hsu_w' in self.inputs:
            self.su_H.set_properties(w=self.inputs['Hsu_w'])
        elif 'Hsu_Tdp' in self.inputs:
            self.su_H.set_properties(Tdp=self.inputs['Hsu_Tdp'])
        elif 'Hsu_Twb' in self.inputs:
            self.su_H.set_properties(Tdp=self.inputs['Hsu_Twb'])
        
        if 'Hsu_p' in self.inputs:
            self.su_H.set_properties(p=self.inputs['Hsu_p'])
        if 'Hsu_m_dot' in self.inputs:
            self.su_H.set_properties(m_dot=self.inputs['Hsu_m_dot'])
        if 'Hsu_V_dot' in self.inputs:
            self.su_H.set_properties(V_dot=self.inputs['Hsu_V_dot'])   

        # Cold Fluid
        self.su_C.set_fluid(self.inputs['Csu_fluid'])
        if 'Csu_T' in self.inputs:
            self.su_C.set_T(self.inputs['Csu_T'])
        elif 'Csu_h' in self.inputs:
            self.su_C.set_h(self.inputs['Csu_h'])
        if 'Csu_p' in self.inputs:
            self.su_C.set_p(self.inputs['Csu_p'])
        if 'Csu_m_dot' in self.inputs:
            self.su_C.set_m_dot(self.inputs['Csu_m_dot'])

        return['fluid_wf', 'su_wf_h', 'su_wf_m_dot', 'fluid_sf', 'su_sf_T', 'su_sf_cp', 'su_sf_m_dot']

#%%

    def get_required_parameters(self):
        """
        General Parameters : 
            
            - htc_type  : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - H_DP_ON   : Hot side pressure drop considered or not
            - C_DP_ON   : Cold side pressure drop considered or not
                   
        Geometry Parameters depend on specific geometry python files.
            
        """

        general_parameters = ['H_DP_ON', 'C_DP_ON']
            
        geometry_parameters = ['A_finned', 'A_flow', 'A_in_tot', 'A_out_tot', 'A_unfinned',
                                'B_V_tot', 'Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type',
                                'Finned_tube_flag', 'L', 'T_V_tot', 'Tube_L', 'Tube_OD',
                                'Tube_cond', 'Tube_t', 'fouling', 'h', 'k_fin',
                                'Tube_pass', 'n_rows', 'n_tubes', 'pitch', 'pitch_ratio', 'tube_arrang',
                                'w','Fin_Side']
        
        return general_parameters + geometry_parameters

#%%
            
    def print_setup(self):
        print("=== Cooling coil Setup ===")
        print("Connectors:")
        print(f"  - H_su: fluid={self.su_H.fluid}, T={self.su_H.T}, w={self.su.w} p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - C_su: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")

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
    def set_debug(self):
        self.debug = 1
        return

