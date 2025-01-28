# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:38:18 2025

@author: Basile
"""

"""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

READ THIS: 
    
    All tasks have been marked with !!! or # !!! symbols and can be seen in the scrollbar of some IDEs.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

# !!! Put what you need to import here

# Ensures that the file can be imported from other library directories
import __init__
 
"INTERNAL IMPORTS"

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector
from component.base_component import BaseComponent

# !!! Put what you need to import here from the library


#%% CLASS DEFINITION : # !!! The name of the class shall be written in the same fashion as hereunder

class HeatExchangerMB(BaseComponent):
    
#%% DOCSTRING : #!!! IMPORTANT AS IT AUTOMATIZES PART OF THE DOCUMENTATION

    """
        Component: Turbomachine

        Model:  !!! Write where the model comes from

        **Descritpion**:

                !!! Write a succint description of the model 
                    - Method of modelling
                    - Use of the model
                    - Specific outputs
                    - Whether pressure drops are or can be taken into account
                    - Application to geometries
                    - Rating of the models level of complexity, precision, robustness, adaptability
                    
        **Assumptions**:

                !!! Write the assumptions linked to the model

        **Connectors**:

            su (MassConnector): Mass connector for the suction side.

            ex (MassConnector): Mass connector for the exhaust side.

            W_dot (WorkConnector): Work connector for the power transfered to or taken from the fluid

            (Q_dot (HeatConnector): Heat connector for the heat transfer between the fluid and the ambient)

                !!! These connectors shall be present in the class __init__() method with the same name
                
                !!! Explain if there is anything specific to the model concerning connectors (if another connector is required for example)

        **Parameters**:
            
                !!! Detail what parameters are necessary for the model operation
                
                !!! The geometry could be considered as a parameter to not detail everything here

                !!! Explain if there is anything specific to the model concerning parameters (if some are for example correlation dependent)

                Parameters required by the model are listed by the class 'get_required_parameters' method
                The parameter dictionnary can be accessed inside the model using self.params['Name of the parameter']

        **Inputs**:

                !!! Write all the required inputs for the model to be able to compute. Here are some nomenclature examples:        

                su_p: Hot suction side pressure. [Pa]
    
                su_h: Hot suction side enthalpy. [J/kg]
    
                su_fluid: Hot suction side fluid. [-]
    
                su_m_dot: Hot suction side mass flowrate. [kg/s]
                
                N: Rotating speed. [RPM]

        **Ouputs**:

                !!! Write all the expected outputs of the model once solved. Here are some nomenclature examples:                  

                ex_h: Hot exhaust side specific enthalpy. [J/kg]
    
                ex_p: Hot exhaust side pressure. [Pa]
    
                ex_h: Cold exhaust side specific enthalpy. [J/kg]
    
                ex_p: Cold exhaust side pressure. [Pa]
    
                W_dot: Power transfered to/taken from the fluid [W]
                
                W_dot_el: Electrical power produced/consumed [W]
    """

#%% SUBCLASSES 

    """
    !!! If subclasses are needed, prefered for the model, define them here and 
    """
    
#%% INIT METHODD: # !!! WRITE A SHORT DOCSTRING 

    def __init__(self):
        """
        su : Supply 

        ex : Exhaust 
                    
        W_dot : Work connection to the fluid
        """
        
        # !!! This nomenclature shall be common to all HX models
        
        super().__init__() # Imports the stuff from BaseComponent (important to know what is inside)
        
        self.su = MassConnector()
        
        self.ex = MassConnector()
                
        self.W_dot = WorkConnector()

        # !!! WRITE WHAT IS MISSING FOR YOUR MODEL HERE 

#%% INPUTS AND PARAMETERS HANDLING METHODS
    
    def get_required_inputs(self): # Used in check_calculable (in BaseComponent) to see if all of the required inputs are set
        """
        Hot side required inputs : 
            
            - su_T or su_h : Supply temperature or enthalpy [K] or [J/kg]
            - su_p          : Supply pressure [Pa]
            - su_fluid      : Supply fluid [string]
            - su_m_dot      : Supply flow rate [kg/s]
            
            - N : Rotating speed [RPM]
        """
        self.sync_inputs()
        
        # !!! Here under is an example for heat exchangers, write in the same fashion the required inputs of your model
        
        # Return a list of required inputs
        return['su_p', 'su_T', 'su_m_dot', 'su_fluid', 'N'] # !!! This nomenclature shall be common to all HX models
    
    def sync_inputs(self): # Makes the link between inputs and MassConnector() in object supply and exhaust attributes
        """Synchronize the inputs dictionary with the connector states."""

        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        elif self.su.h is not None:
            self.inputs['su_h'] = self.su.h
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.m_dot is not None:
            self.inputs['su_m_dot'] = self.su.m_dot
        if self.W_dot.N is not None:
            self.inputs['N'] = self.W_dot.N

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        # Hot Fluid
        self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        elif 'su_h' in self.inputs:
            self.su.set_h(self.inputs['su_h'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'su_m_dot' in self.inputs:
            self.su.set_m_dot(self.inputs['su_m_dot'])
        if 'N' in self.inputs:
            self.W_dot.N = self.inputs['N']

    def get_required_parameters(self): # Used in check_parametrized (in BaseComponent) to see if all of the required parameters are set
        """
        Required parameters could be chosen depending on the geometry, type of the model. 
        
        Here are some examples, (these might not all be useful, a filtering might be necessary)
            
        """
        return [
            'D_inlet', 'N_turb_rated', 'turb_voltage', 'turb_phases', 'eta_max_motor', 
            'W_dot_el_rated', 'eta_m', 'eta_is_coefs', 'eta_is_coefs_red', 'A_th'
        ]
    
    def print_setup(self): 
        print("=== Heat Exchanger Setup ===")
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
           
#%% INTERNAL METHODS # !!! DON'T FORGET DOCSTRINGS

    """
    !!! IF YOUR MODEL REQUIRES OTHER METHODS FOR ITS SOLVING, WRITE THEM HERE 
    
    (REALLY SPECIFIC TO THIS MODEL. IF IT MIGHT BE COMMON TO OTHER MODELS, THEN IT IS EITHER IN CORRELATIONS OR IN TOOLBOX)
    
    """
    def update_connectors(self): # !!! What is inside the parentheses shall fit what is computed in your model 
        """Update the connectors with the calculated values."""
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_h(self.h_ex)
        self.ex.set_p(self.P_ex)
        self.ex.set_m_dot(self.su.m_dot)

        self.W_dot.set_W_dot(self.W_dot_fluid)

#%% SOLVE METHOD # !!! DON'T FORGET DOCSTRINGS
    
    def solve(self):    
        """
        # !!! DOCSTRING: DESCRIBE THE METHOD HERE IF ANYTHING SPECIAL
        """
        
        self.check_calculable()
        self.check_parametrized()
        
        if not self.calculable:
            print("Component not calculable, check input")
            return
            
        if not self.parametrized:
            print("Component not parametrized, check parameters")      
            return
            
        # !!! WRITE YOUR SOLVING CODE HERE
        
        # !!! DO NOT FORGET TO UPDATE EXHAUST MASS CONNECTOR VALUES AS WELL AS HEAT CONNECTOR (OR WORK CONNECTOR IF THERE IS)
        
        self.update_connectors()
        self.solved = True # Shall be done only if the solving was successful, otherwise, solved shall stay at False value
        return

#%% PRINT THE RESULTS 

    # !!! These shall fit the output of your model
    def print_results(self):
        print("=== Solar Collector Results ===")
        print(f"  - ex_h: {self.ex.h} [J/kg]")
        print(f"  - ex_T: {self.ex.T} [K]")
        print(f"  - W_dot: {self.W_dot.W_dot} [W]")
        print("=========================")

    # !!! Again, what is in the print shall correspond to what is in the __init__() method
    def print_states_connectors(self):
        print("=== Solar Collector States ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")
        print(f"  - W_dot: {self.W_dot.W_dot}")
        print("======================")

