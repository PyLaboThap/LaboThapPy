# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile - Marie 
"""

from connector.mass_connector import MassConnector
from CoolProp.CoolProp import PropsSI

from component.base_component import BaseComponent
import numpy as np         
            
            
class Mixer(BaseComponent):
    """
    Component: Tank

    Model: Mixer

    Reference: /

    **Description**:

        This model simulates a perfectly mixed tank where multiple inlet streams of the same fluid are combined into a single outlet stream. 
        It computes the outlet specific enthalpy and mass flow rate based on an energy balance with no heat losses and no accumulation.


    **Assumptions**:

        - Steady-state operation.
        - Perfect mixing (single outlet enthalpy and pressure).
        - All inlet streams must contain the same working fluid.
        - All inlet pressures must be within a fixed tolerance (default: 100 Pa).
        - No heat exchange with the environment (adiabatic mixing).
        - No phase change during mixing (inherent in the use of enthalpy balance).


    **Connectors**:

        su (MassConnector): Mass connector for the suction side ((su_1, su_2, ..., su_n).

        ex (MassConnector): Mass connector for the exhaust side.
        
        
    **Parameters**: 
        
        n_inlets : Number of inlet streams to be mixed [-]
        

    **Inputs**:

        For each inlet stream `i` from 1 to n_inlets:
            
        p_su_i: Suction side pressure. [Pa]

        T_su_i: Suction side temperature. [K]
        
        m_dot_su_i: Suction side mass flow flow rate. [kg/s]

        fluid: Working fluid [-]
    

    **Ouputs**:

        p_ex: Outlet pressure [Pa] (mean of inlet pressures)

        h_ex: Exhaust side specific enthalpy. [J/kg]
            
    """


    def __init__(self, n_inlets=None):
        
        super().__init__()
        
        "Status variables"
        self.calculable = None
        self.parametrized = None
        self.solved = False
        
        "Supply"
        if n_inlets is not None: 
            self.n_inlets = n_inlets
            
            self.params['n_inlets'] = n_inlets
            
            for i in range(n_inlets): # Create inlet connectors dynamically based on number of inlets
                inlet_num = i + 1
                setattr(self, f"su_{inlet_num}", MassConnector())
        else:
            raise ValueError("'Mixer' model requires to set a value for its 'n_inlets' input")
        
        "Exhaust"
        self.ex = MassConnector() # Create single exhaust connector

                
#%%    

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
    
        """
        Required inputs for the Mixer:
        
            For each inlet i:
                - T_su_i or h_su_i : Inlet temperature or specific enthalpy [K or J/kg]
                - p_su_i          : Inlet pressure [Pa]
                - fluid_su_i       : Inlet fluid [-]
                - m_dot_su_i       : Inlet mass flow rate [kg/s]
        """

        self.sync_inputs()
        # Return a list of required inputs
        
        required_inputs = []
        
        for i in range(self.n_inlets):
            inlet_num = i + 1 
            required_inputs.append(f"p_su_{inlet_num}")
            required_inputs.append(f"T_su_{inlet_num}")
            required_inputs.append(f"m_dot_su_{inlet_num}")
            required_inputs.append(f"fluid_su_{inlet_num}")
        
        return required_inputs
    
    

    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        
        for i in range(self.n_inlets):
            inlet_num = i + 1 
            
            if getattr(self, f"su_{inlet_num}").T is not None:
                self.inputs[f"T_su_{inlet_num}"] = getattr(self, f"su_{inlet_num}").T 
            elif getattr(self, f"su_{inlet_num}").h is not None:
                self.inputs[f"h_su_{inlet_num}"] = getattr(self, f"su_{inlet_num}").h
            if getattr(self, f"su_{inlet_num}").p is not None:
                self.inputs[f"p_su_{inlet_num}"] = getattr(self, f"su_{inlet_num}").p
            if getattr(self, f"su_{inlet_num}").fluid is not None:
                self.inputs[f"fluid_su_{inlet_num}"] = getattr(self, f"su_{inlet_num}").fluid
            if getattr(self, f"su_{inlet_num}").m_dot is not None:
                self.inputs[f"m_dot_su_{inlet_num}"] = getattr(self, f"su_{inlet_num}").m_dot

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.
        
        # Update the connectors based on the new inputs
        for i in range(self.n_inlets):
            inlet_num = i + 1 
            connector = getattr(self, f"su_{inlet_num}")
            # Update the connectors based on the new inputs
            if f"fluid_su_{inlet_num}" in self.inputs:   
                connector.set_fluid(self.inputs[f"fluid_su_{inlet_num}"])
            if f"T_su_{inlet_num}" in self.inputs:
                connector.set_T(self.inputs[f"T_su_{inlet_num}"])
            elif f"h_su_{inlet_num}" in self.inputs:
                connector.set_h(self.inputs[f"h_su_{inlet_num}"])
            if f"p_su_{inlet_num}" in self.inputs:
                connector.set_p(self.inputs[f"p_su_{inlet_num}"])
            if f"m_dot_su_{inlet_num}" in self.inputs:
                connector.set_m_dot(self.inputs[f"m_dot_su_{inlet_num}"])

    def get_required_parameters(self):
        
        """        
        General Parameters : 
            
            - n_inlets : Number of inlets to the tank [-]
        
        Geometry Parameters depend on specific geometry python files.
        
        """
        return ['n_inlets'] 

    
    
    def print_setup(self):
        print("=== Tank Mixer Setup ===")
        print("Connectors:")
        for i in range(1, self.geom.n_inlet + 1):
            connector = getattr(self, f"inlet_{i}")
            print(f"  - inlet_{i}: fluid={connector.fluid}, T={connector.T}, p={connector.p}, m_dot={connector.m_dot}")


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
        """
        Solves the mixer by performing an enthalpy balance across inlets.
        
        Checks for:
            - Same fluid in all inlets
            - Pressure consistency across inlets
        
        Sets:
            - Outlet pressure as mean of inlet pressures
            - Outlet enthalpy as mass-weighted average
            - Outlet mass flow rate as sum of all inlet flow rates
        """
                
        self.check_calculable()
        self.check_parametrized()
                        
        
        if not self.calculable:
            print("Component not calculable, check input")
            
        if not self.parametrized:
            print("Component not parametrized, check parameters") 
        
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
            
            if connector.h is not None:
                mean_h += connector.h * connector.m_dot
            elif connector.T is not None and connector.p is not None:
                h = PropsSI('H', 'T', connector.T, 'P', connector.p, connector.fluid)
                mean_h += h * connector.m_dot
            else:
                raise ValueError(f"Missing enthalpy or temperature/pressure for inlet {inlet_num}")
                
        # Final enthalpy and pressure computation
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
                print("Pressure difference between inlets exceeds tolerance.")
                return
        else:
            raise ValueError("Mixing different fluids in 'Mixer'")

        return
    
    
    
    
    
    