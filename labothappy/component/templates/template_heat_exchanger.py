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
from connector.heat_connector import HeatConnector
from component.base_component import BaseComponent

# !!! Put what you need to import here from the library


#%% CLASS DEFINITION : # !!! The name of the class shall be written in the same fashion as hereunder

class HeatExchangerMB(BaseComponent):
    
#%% DOCSTRING : #!!! IMPORTANT AS IT AUTOMATIZES PART OF THE DOCUMENTATION

    """
        Component: Heat Exchanger

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

            su_H (MassConnector): Mass connector for the hot suction side.
            su_C (MassConnector): Mass connector for the cold suction side.

            ex_H (MassConnector): Mass connector for the hot exhaust side.
            ex_C (MassConnector): Mass connector for the cold exhaust side.

            Q_HX (HeatConnector): Heat connector for the heat transfer between the fluids

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

                su_H_p: Hot suction side pressure. [Pa]
    
                su_H_h: Hot suction side enthalpy. [J/kg]
    
                su_H_fluid: Hot suction side fluid. [-]
    
                su_H_m_dot: Hot suction side mass flowrate. [kg/s]
    
                su_C_p: Cold suction side pressure. [Pa]
    
                su_C_h: Cold suction side enthalpy. [J/kg]
    
                su_C_fluid: Cold suction side fluid. [-]
    
                su_C_m_dot: Cold suction side mass flowrate. [kg/s]

        **Ouputs**:

                !!! Write all the expected outputs of the model once solved. Here are some nomenclature examples:                  

                ex_H_h: Hot exhaust side specific enthalpy. [J/kg]
    
                ex_H_p: Hot exhaust side pressure. [Pa]
    
                ex_C_h: Cold exhaust side specific enthalpy. [J/kg]
    
                ex_C_p: Cold exhaust side pressure. [Pa]
    
                Q_dot: Heat transfer rate [W]
    
                M_H : Hot fluid charge [kg]
    
                M_C : Cold fluid charge [kg]
    """

#%% SUBCLASSES 

    """
    !!! If subclasses are needed, prefered for the model, define them here and 
    """
    
#%% INIT METHODD: # !!! WRITE A SHORT DOCSTRING 

    def __init__(self):
        """
        su : Supply - 'H' : hot
                    - 'C' : cold

        ex : Exhaust - 'H' : hot
                     - 'C' : cold
                    
        Q_HX : Heat connection between fluids
        
        HTX_Type : Type of HTX - Plate 
                               - Shell and Tube
                               - Tube and Fins
        """
        
        # !!! This nomenclature shall be common to all HX models
        
        super().__init__() # Imports the stuff from BaseComponent (important to know what is inside)
        
        self.su_H = MassConnector()
        self.su_C = MassConnector()
        
        self.ex_H = MassConnector()
        self.ex_C = MassConnector() 
                
        self.Q_HX = HeatConnector()
        
        # !!! WRITE WHAT IS MISSING FOR YOUR MODEL HERE 

#%% INPUTS AND PARAMETERS HANDLING METHODS
    
    def get_required_inputs(self): # Used in check_calculable (in BaseComponent) to see if all of the required inputs are set
        """
        Hot side required inputs : 
            
            - Hsu_T or Hsu_h : Hot supply temperature or enthalpy
            - Hsu_p          : Hot supply pressure
            - Hsu_fluid      : Hot supply fluid
            - Hsu_m_dot      : Hot supply flow rate
            
        Cold side required inputs : 
            
            - Csu_T or Hsu_h : Cold supply temperature or enthalpy
            - Csu_p          : Cold supply pressure
            - Csu_fluid      : Cold supply fluid
            - Csu_m_dot      : Cold supply flow rate
        """
        self.sync_inputs()
        
        # !!! Here under is an example for heat exchangers, write in the same fashion the required inputs of your model
        
        # Return a list of required inputs
        return['Hsu_p', 'Hsu_T', 'Hsu_m_dot', 'Hsu_fluid', 'Csu_p', 'Csu_T', 'Csu_m_dot', 'Csu_fluid'] # !!! This nomenclature shall be common to all HX models
    
    def sync_inputs(self): # Makes the link between inputs and MassConnector() in object supply and exhaust attributes
        """Synchronize the inputs dictionary with the connector states."""
        # Hot Fluid
        if self.su_H.T is not None:
            self.inputs['Hsu_T'] = self.su_H.T
        elif self.su_H.h is not None:
            self.inputs['Hsu_h'] = self.su_H.h
        if self.su_H.p is not None:
            self.inputs['Hsu_p'] = self.su_H.p
        if self.su_H.fluid is not None:
            self.inputs['Hsu_fluid'] = self.su_H.fluid
        if self.su_H.m_dot is not None:
            self.inputs['Hsu_m_dot'] = self.su_H.m_dot
            
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
        self.su_H.set_fluid(self.inputs['Hsu_fluid'])
        if 'Hsu_T' in self.inputs:
            self.su_H.set_T(self.inputs['Hsu_T'])
        elif 'Hsu_h' in self.inputs:
            self.su_H.set_h(self.inputs['Hsu_h'])
        if 'Hsu_p' in self.inputs:
            self.su_H.set_p(self.inputs['Hsu_p'])
        if 'Hsu_m_dot' in self.inputs:
            self.su_H.set_m_dot(self.inputs['Hsu_m_dot'])

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


    def get_required_parameters(self): # Used in check_parametrized (in BaseComponent) to see if all of the required parameters are set
        """
        General Parameters : 
            
            - Flow_Type : Flow configuration of the fluid ('CounterFlow', 'CrossFlow', 'Shell&Tube', 'ParallelFlow')
            - htc_type  : Heat Transfer coefficient type ('User-Defined' or 'Correlation')
            - H_DP_ON   : Hot side pressure drop considered or not
            - C_DP_ON   : Cold side pressure drop considered or not
            - n_disc    : number of discretizations
        
        Geometry Parameters depend on specific geometry python files.
            
        """
        general_parameters = ['Flow_Type','htc_type', 'H_DP_ON', 'C_DP_ON','n_disc']
        
        
        # !!! Here is a complex example for a geometry dependent model that can adapt to many HX technologies
        
        if self.HTX_Type == 'Plate':     
            geometry_parameters = ['A_c', 'A_h', 'h', 'l', 'l_v',
                                   'C_CS', 'C_Dh', 'C_V_tot', 'C_canal_t', 'C_n_canals',
                                   'H_CS', 'H_Dh', 'H_V_tot', 'H_canal_t', 'H_n_canals',
                                   'casing_t', 'chevron_angle', 'fooling', 
                                   'n_plates', 'plate_cond', 'plate_pitch_co', 't_plates', 'w']

        elif self.HTX_Type == 'Shell&Tube':
                
            geometry_parameters = ['A_eff', 'Baffle_cut', 'D_OTL', 'N_strips', 'S_V_tot',
                                   'Shell_ID', 'T_V_tot', 'Tube_L', 'Tube_OD', 'Tube_pass',
                                   'Tube_t', 'Tubesheet_t', 'central_spacing', 'clear_BS', 'clear_TB',
                                   'cross_passes', 'foul_s', 'foul_t', 'inlet_spacing', 'n_series',
                                   'n_tubes', 'outlet_spacing', 'pitch_ratio', 'tube_cond', 'tube_layout', 'Shell_Side']
        
        elif self.HTX_Type == 'Tube&Fins':
            
            geometry_parameters = ['A_finned', 'A_flow', 'A_in_tot', 'A_out_tot', 'A_unfinned',
                                   'B_V_tot', 'Fin_OD', 'Fin_per_m', 'Fin_t', 'Fin_type',
                                   'Finned_tube_flag', 'L', 'T_V_tot', 'Tube_L', 'Tube_OD',
                                   'Tube_cond', 'Tube_t', 'fouling', 'h', 'k_fin',
                                   'n_passes', 'n_rows', 'n_tubes', 'pitch', 'pitch_ratio', 'tube_arrang',
                                   'w','Fin_Side']
        
        return general_parameters + geometry_parameters # !!! In the end, what is important is that all your parameters are here
    
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
        self.ex_H.set_fluid(self.su_H.fluid)
        self.ex_H.set_h(self.h_ex)
        self.ex_H.set_p(self.P_ex)
        self.ex_H.set_m_dot(self.su_H.m_dot)

        self.ex_C.set_fluid(self.su_C.fluid)
        self.ex_C.set_h(self.h_ex)
        self.ex_C.set_p(self.P_ex)
        self.ex_C.set_m_dot(self.su_C.m_dot)

        self.Q_HX.set_Q_dot(self.Q_dot)

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
        print("=== Expander Results ===")
        print(f"  - Hex_h: {self.ex_H.h} [J/kg]")
        print(f"  - Hex_T: {self.ex_H.T} [K]")
        print(f"  - Cex_h: {self.ex_C.h} [J/kg]")
        print(f"  - Cex_T: {self.ex_C.T} [K]")
        print(f"  - Q_dot: {self.Q_HX.Q_dot} [W]")
        print("=========================")

    # !!! Again, what is in the print shall correspond to what is in the __init__() method
    def print_states_connectors(self):
        print("=== Heat Exchanger States ===")
        print("Connectors:")
        print(f"  - su_C: fluid={self.su_C.fluid}, T={self.su_C.T}, p={self.su_C.p}, m_dot={self.su_C.m_dot}")
        print(f"  - su_H: fluid={self.su_H.fluid}, T={self.su_H.T}, p={self.su_H.p}, m_dot={self.su_H.m_dot}")
        print(f"  - ex_C: fluid={self.ex_C.fluid}, T={self.ex_C.T}, p={self.ex_C.p}, m_dot={self.ex_C.m_dot}")
        print(f"  - ex_H: fluid={self.ex_H.fluid}, T={self.ex_H.T}, p={self.ex_H.p}, m_dot={self.ex_H.m_dot}")
        print(f"  - Q_dot: {self.Q_HX.Q_dot}")
        print("======================")

