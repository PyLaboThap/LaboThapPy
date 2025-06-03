# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""
import sys
import os

# Get the absolute path of the directory that contains the script (simulation_model.py)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Determine the project root directory (which contains both 'connector' and 'component')
project_root = os.path.abspath(os.path.join(current_dir, '..', '..', '..'))

# Add the project root to sys.path if it's not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)
###########################################

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from component.base_component import BaseComponent
import numpy as np

class Spliter(BaseComponent):
    
    class geom():
            pass 

    def __init__(self, geom=None, outlet_repartition=None):
        
        super().__init__()
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.solved = False
        
        "Geometry sub-object"
        self.geom = geom # parameters 
        
        if self.geom is not None:
            self.outlet_repartition = self.geom.flow_coef
        else: 
            if outlet_repartition is not None:
                self.outlet_repartition = outlet_repartition
            else:
                raise ValueError("'Spliter' model requires to set an array for its 'outlet_repartition' input")
        
        "Input"
        self.su = MassConnector()

        "Outputs"
        for i in range(len(self.outlet_repartition)):
            outlet_num = i + 1
            setattr(self, f"ex_{outlet_num}", MassConnector())
                
#%%    

    def get_required_inputs(self):
            self.sync_inputs()
            # Return a list of required inputs
            return ['su_p', 'su_T', 'su_m_dot', 'su_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        if self.su.h is not None:
            self.inputs['su_h'] = self.su.h
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.su.m_dot is not None:
            self.inputs['su_m_dot'] = self.su.m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        if 'su_h' in self.inputs:
            self.su.set_T(self.inputs['su_h'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'su_m_dot' in self.inputs:
            self.su.set_m_dot(self.inputs['su_m_dot'])

    def get_required_parameters(self):
        return [
        ]
    
    def print_setup(self):
        print("=== Pump Setup ===")
       
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

    def solve(self):
        
        "1) Compute output"
        
        if self.check_parametrized():
            if self.check_calculable():
                tolerance = 1e-3
                if abs(sum(self.outlet_repartition) - 1) < tolerance:
                    for i in range(len(self.outlet_repartition)):
                        outlet_num = i + 1
                        connector = getattr(self, f"ex_{outlet_num}")
                        connector.set_fluid(self.su.fluid)
                        connector.set_p(self.su.p)
                        connector.set_h(self.su.h)
                        connector.set_m_dot(self.su.m_dot*self.outlet_repartition[i])
            
                    self.solved = True
                else:
                    raise ValueError(f"'Spliter' total outlet repartition is too far from 1 (tolerance : {tolerance}).")
            else:
                raise ValueError("'Spliter' is not fully calculable.")
        else:
            raise ValueError("'Spliter' is not fully parametrized.")
        
        return
        
        
            