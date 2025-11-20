# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:24 2024

@author: Basile
"""


from connector.mass_connector import MassConnector
from component.base_component import BaseComponent
import numpy as np

class TankSpliter(BaseComponent):
    
    """
    **Component**: Spliter
    
    **Model**: Mass Flow Splitting Model
    
    **Description**:
    
        The `Spliter` component divides an incoming mass flow into multiple outlet streams based on specified mass flow fractions. It assumes no change in fluid properties (pressure, temperature, enthalpy) across the splitter. The model is useful for distributing mass flows in a system model where the thermodynamic state remains unchanged.
    
    **Assumptions**:
    
        - Steady-state operation.
        - No heat or work interactions (adiabatic, no shaft work).
        - Ideal splitting with no pressure or temperature drops.
        - The sum of the outlet flow fractions must equal 1 (within a small numerical tolerance).
    
    **Connectors**:
    
        su (MassConnector): Mass connector for the inlet supply stream.
        ex_1, ex_2, ..., ex_n (MassConnector): Mass connectors for each outlet stream, created based on the number of elements in `outlet_repartition`.
    
    **Parameters**:
    
        outlet_repartition (list of float): Fractions of the total inlet mass flow allocated to each outlet. Must sum to 1 within a small tolerance.
    
    **Inputs**:
    
        fluid (str): Inlet fluid.
        T_su (float): Inlet temperature [K].
        h_su (float): Inlet specific enthalpy [J/kg].
        P_su (float): Inlet pressure [Pa].
        m_dot (float): Inlet mass flow rate [kg/s].
    
    **Outputs**:
    
        For each outlet connector (ex_1, ex_2, ..., ex_n):
            - fluid: Same as inlet fluid.
            - P: Same as inlet pressure [Pa].
            - h: Same as inlet enthalpy [J/kg].
            - m_dot: Mass flow rate for the outlet, computed as `su_m_dot * outlet_repartition[i]` [kg/s].
    
    """


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
                
  

    def get_required_inputs(self):
        # Return a list of required inputs
        return ['P_su', 'T_su', 'm_dot', 'fluid']


    def get_required_parameters(self):
        return []

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
        
        
            