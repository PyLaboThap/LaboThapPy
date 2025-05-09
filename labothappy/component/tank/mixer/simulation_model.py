# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI
from component.base_component import BaseComponent
import numpy as np

class Mixer(BaseComponent):
    
    class geom():
            pass 

    def __init__(self, n_inlets=None):
        
        super().__init__()
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.solved = False
        
        "Supply"
        
        if n_inlets is not None:
            self.n_inlets = n_inlets
            
            for i in range(n_inlets):
                inlet_num = i + 1
                setattr(self, f"su_{inlet_num}", MassConnector())
        else:
            raise ValueError("'Mixer' model requires to set a value for its 'n_inlet' input")
        
        "Exhaust"
        self.ex = MassConnector()

                
#%%    

    def get_required_inputs(self):
            self.sync_inputs()
            # Return a list of required inputs
            
            required_inputs = []
            
            for i in range(self.n_inlets):
                inlet_num = i + 1 
                required_inputs.append(f"su_{inlet_num}_p")
                required_inputs.append(f"su_{inlet_num}_T")
                required_inputs.append(f"su_{inlet_num}_m_dot")
                required_inputs.append(f"su_{inlet_num}_fluid")
            
            return required_inputs
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        
        for i in range(self.n_inlets):
            inlet_num = i + 1 
            if getattr(self, f"su_{inlet_num}").fluid is not None:
                self.inputs[f"su_{inlet_num}_fluid"] = getattr(self, f"su_{inlet_num}").fluid
            if getattr(self, f"su_{inlet_num}").T is not None:
                self.inputs[f"su_{inlet_num}_T"] = getattr(self, f"su_{inlet_num}").T 
            if getattr(self, f"su_{inlet_num}").h is not None:
                self.inputs[f"su_{inlet_num}_h"] = getattr(self, f"su_{inlet_num}").h
            if getattr(self, f"su_{inlet_num}").p is not None:
                self.inputs[f"su_{inlet_num}_p"] = getattr(self, f"su_{inlet_num}").p
            if getattr(self, f"su_{inlet_num}").m_dot is not None:
                self.inputs[f"su_{inlet_num}_m_dot"] = getattr(self, f"su_{inlet_num}").m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        for i in range(self.n_inlets):
            inlet_num = i + 1 
            connector = getattr(self, f"su_{inlet_num}")
            # Update the connectors based on the new inputs
            if f"su_{inlet_num}_fluid" in self.inputs:   
                connector.set_fluid(self.inputs[f"su_{inlet_num}_fluid"])
            if f"su_{inlet_num}_T" in self.inputs:
                connector.set_T(self.inputs[f"su_{inlet_num}_T"])
            if f"su_{inlet_num}_h" in self.inputs:
                connector.set_T(self.inputs[f"su_{inlet_num}_h"])
            if f"su_{inlet_num}_p" in self.inputs:
                connector.set_p(self.inputs[f"su_{inlet_num}_p"])
            if f"su_{inlet_num}_m_dot" in self.inputs:
                connector.set_m_dot(self.inputs[f"su_{inlet_num}_m_dot"])

    def get_required_parameters(self):
        return []
    
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

    def are_pressures_close(self,pressures, tolerance):
        """Check if all values are close to each other within the specified tolerance."""
        reference_value = np.mean(pressures)
        for pressure in pressures:
            if abs(pressure - reference_value) > tolerance:
                return False
        return True

    def solve(self):
        
        
        
        "1) Compute output"
                
        pressures = np.zeros(self.n_inlets)
        m_dot = np.zeros(self.n_inlets) 
        fluids = [] 
        mean_h = 0
        
        for i in range(self.n_inlets):
            inlet_num = i + 1
            connector = getattr(self, f"su_{inlet_num}")
            pressures[i] = connector.p 
            m_dot[i] = connector.m_dot
            fluids.append(connector.fluid)
            try:
                mean_h += connector.h*connector.m_dot
            except:
                mean_h += PropsSI('H', 'T', connector.T, 'P', connector.p, connector.fluid)*connector.m_dot
        
        mean_h = mean_h/(sum(m_dot))
        mean_p = np.mean(pressures)
        
        if len(set(fluids)) == 1: # All fluids are the same
            tolerance = 100
            if self.are_pressures_close(pressures, tolerance):
                self.ex.set_fluid(self.su_1.fluid)
                self.ex.set_p(mean_p)
                self.ex.set_m_dot(sum(m_dot))
                self.ex.set_h(mean_h)
                
                self.solved = True
            else:
                self.solved = False
                return
                # raise ValueError(f"Mixing different pressure flows (difference higher than tolerance = {tolerance} Pa) in 'Mixer'")
        else:
            raise ValueError("Mixing different fluids in 'Mixer'")

        
        
        return