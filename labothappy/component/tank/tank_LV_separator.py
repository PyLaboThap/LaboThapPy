# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:09:50 2025

@author: Samuel Gendebien
"""

#%% IMPORTS

"EXTERNAL IMPORTS"
 
# External Toolbox 
import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.interpolate import interp1d
import copy


# Ensures that the file can be imported from other library directories
import __init__
 
"INTERNAL IMPORTS"

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector
from component.base_component import BaseComponent

# !!! Put what you need to import here from the library


#%% CLASS DEFINITION : # !!! The name of the class shall be written in the same fashion as hereunder

class TankLVSeparator(BaseComponent):
    """
        **Component**: Tank separator

        **Model**:  Classical way to model a LV separator 
        
        **Reference**: /

        **Descritpion**:

            - Method of modeling: This code defines an LV Separator that takes an input stream and splits it into liquid and vapor components based on quality (x) and temperature.
            - Use of the model: The model can used in thermodynamic cycle
            - Specific outputs: mass flow rate in both liquid and vapor state and exhaust conditions for both lines 
            - Pressure drop not modeled
            - Geometry not required
            - Rating of the models level of complexity, precision, robustness, adaptability : simple
                    
        **Assumptions**:

            - No work, no heat losses
                    
        **Connectors**:
            su (MassConnector): Mass connector for the supply conditions
            ex_v (MassConnector): Mass connector for the vapor line
            ex_l (MassConnector): Mass connector for the liquid line


        **Parameters**:
            
            No parameters are expected

        **Inputs**:    
            
            x_su          : Supply quality [-]
            T_su or h_su  : Supply temperature or enthalpy [K] or [J/kg]
            P_su          : Supply pressure [Pa]
            fluid         : Supply fluid [string]
            m_dot         : Supply flow rate [kg/s]
                
        **Ouputs**:
            
            Vapor line
            h_ex_v       : Exhaust enthalpy at the vapor line [J/kg]
            P_ex_v       : Exhaust pressure at the vapor line [Pa]
            fluid_v      : Exhaust fluid at the vapor line [string]
            m_dot_v      : Exhaust flow rate at the vapor line [kg/s]
            
            Vapor line
            h_ex_l       : Exhaust enthalpy at the liquid line [J/kg]
            P_ex_l       : Exhaust pressure at the liquid line [Pa]
            fluid_l      : Exhaust fluid at the liquid line [string]
            m_dot_l      : Exhaust flow rate at the liquid line [kg/s]
    
    """

    def __init__(self):
        """
        su : Supply 
        ex_v: Exhaust vapor line
        ex_l : Exhaust liquid line
        """
        
        # !!! This nomenclature shall be common to all HX models
        
        super().__init__() # Imports the stuff from BaseComponent (important to know what is inside)
        
        self.su = MassConnector()
        self.ex_l = MassConnector()
        self.ex_v = MassConnector () 
            
    def get_required_inputs(self): # Used in check_calculable (in BaseComponent) to see if all of the required inputs are set
        """
        Supply required inputs : 
            - x_su          : Supply quality [-]
            - T_su          : Supply temperature [K]
            - h_su          : Supply enthalpy [J/kg]
            - P_su          : Supply pressure [Pa]
            - fluid         : Supply fluid [string]
            - m_dot         : Supply flow rate [kg/s]
            
        """
        # The inputs should be defined in order to define the state of the fluid (i.e x_su should be defined at saturation)
        return['x_su', 'T_su', 'h_su', 'P_su', 'm_dot', 'fluid'] 
    
   
    def get_required_parameters(self): # Used in check_parametrized (in BaseComponent) to see if all of the required parameters are set
        """
        No required parameters for the model 
        """
        # print ('No required parameters')
        return []

           

    def update_connectors(self, x_ex_l, T_ex_l, h_ex_l, p_ex_l, m_dot_l, x_ex_v, T_ex_v, h_ex_v, p_ex_v, m_dot_v):
        """Update mass connectors after solving"""
        self.ex_l.set_fluid(self.su.fluid)
        self.ex_l.set_x(x_ex_l)
        self.ex_l.set_T(T_ex_l)
        self.ex_l.set_h(h_ex_l) 
        self.ex_l.set_p(p_ex_l)
        self.ex_l.set_m_dot(m_dot_l)  
    
        self.ex_v.set_fluid(self.su.fluid)
        self.ex_v.set_x(x_ex_v)
        self.ex_v.set_T(T_ex_v)
        self.ex_l.set_h(h_ex_v) 
        self.ex_v.set_p(p_ex_v)
        self.ex_v.set_m_dot(m_dot_v)  

    def solve(self):    
        """
        Solve the LV Separator model.
        """
    
        self.check_calculable()
        self.check_parametrized()
    
        if not self.calculable:
            print("Component not calculable, check input")
            return
    
        if not self.parametrized:
            print("Component not parametrized, check parameters")      
            return
    
        self.AS=CP.AbstractState('HEOS', self.su.fluid) 
        
        # Extract input properties
        x_su = self.su.x  # Quality (0 = liquid, 1 = vapor)
        P_su = self.su.p
        T_su = self.su.T
        m_dot_su = self.su.m_dot
    
        # Get saturation temperature at supply pressure
        self.AS.update(CP.PQ_INPUTS,P_su, 0.5)
        T_sat =self.AS.T ()
    
        if 0 <= x_su <= 1:  # Two-phase mixture
            m_dot_l = (1 - x_su) * m_dot_su
            m_dot_v = x_su * m_dot_su
        else:
            if T_su > T_sat:  # Vapor phase
                m_dot_v = m_dot_su
                m_dot_l = 0
            else:  # Liquid phase
                m_dot_l = m_dot_su
                m_dot_v = 0
                
                # Set exhaust enthalpy
        x_ex_v=1
        T_ex_v=T_su
        self.AS.update(CP.PQ_INPUTS, P_su, x_ex_v)
        h_ex_v = self.AS.hmass ()
        P_ex_v=P_su
        
        x_ex_l=0
        self.AS.update(CP.PQ_INPUTS, P_su, x_ex_l)
        h_ex_l = self.AS.hmass ()
        T_ex_l=T_su
        P_ex_l=P_su
        
                
        # Update mass connectors
        self.update_connectors(x_ex_l, T_ex_l, h_ex_l, P_ex_l, m_dot_l, x_ex_v, T_ex_v, h_ex_v, P_ex_v, m_dot_v)
        self.solved = True  # Mark as solved
        
    # !!! These shall fit the output of your model
    def print_results(self):
        print("=== Liquid vapor separator ===")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex_liq: fluid={self.ex_l.fluid}, T={self.ex_l.T}, p={self.ex_l.p}, m_dot_l={self.ex_l.m_dot}")
        print(f"  - ex_vap: fluid={self.ex_v.fluid}, T={self.ex_v.T}, p={self.ex_v.p}, m_dot_v={self.ex_v.m_dot}")
        print("=========================")

    # !!! Again, what is in the print shall correspond to what is in the __init__() method
    def print_states_connectors(self):
        print("=== Liquid vapor separator ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex_liq: fluid={self.ex_l.fluid}, T={self.ex_l.T}, p={self.ex_l.p}, m_dot={self.ex_l.m_dot}")
        print(f"  - ex_vap: fluid={self.ex_v.fluid}, T={self.ex_v.T}, p={self.ex_v.p}, m_dot={self.ex_v.m_dot}")
        print("======================")

