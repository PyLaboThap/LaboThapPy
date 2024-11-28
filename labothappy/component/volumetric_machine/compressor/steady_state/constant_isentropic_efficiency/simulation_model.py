from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

class CompressorCstEff(BaseComponent):
    """
    Component: Compressor

    Model: Constant isentropic efficiency

    **Descritpion**:

        This model determines the exhqust specific enthalpy and the exhaust temperature of a compressor. This model can be used for on-design models of systems.

    **Assumptions**:

        - Steady-state operation.
        - Isentropic efficiency stays constant for all the conditions.

    **Connectors**:

        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

    **Parameters**:

        eta_is: Isentropic efficiency. [-]

    **Inputs**:

        su_p: Suction side pressure. [Pa]

        su_T: Suction side temperature. [K]

        ex_p: Exhaust side pressure. [Pa]

        su_fluid: Suction side fluid. [-]

    **Ouputs**:

        ex_h: Exhaust side specific enthalpy. [J/kg] 

        ex_T: Exhaust side temperature. [K]
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector() # Mass_connector for the suction side
        self.ex = MassConnector() # Mass_connector for the exhaust side

    def get_required_inputs(self):
            self.sync_inputs()
            # Return a list of required inputs
            return ['su_p', 'su_T', 'ex_p', 'su_fluid']
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])

    def get_required_parameters(self):
        return [
            'eta_is',
        ]
    
    def print_setup(self):
        print("=== Compressor Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")


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

    def solve(self):
        if not (self.calculable and self.parametrized):
            self.solved = False
            print("CompressorCstEff could not be solved. It is not calculable and/or not parametrized")
            return
        try:
            h_ex_is = PropsSI('H', 'P', self.ex.p, 'S', self.su.s, self.su.fluid)
            h_ex = self.su.h + (h_ex_is - self.su.h) / self.params['eta_is']

            self.update_connectors(h_ex)

            self.solved = True
        except Exception as e:
            print(f"Error: {e}")
            self.solved = False
            return
    
    def update_connectors(self, h_ex):
        self.ex.set_h(h_ex)

    def print_results(self):
        print("=== Expander Results ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}, p={self.su.p}, m_dot={self.su.m_dot}")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T}, p={self.ex.p}, m_dot={self.ex.m_dot}")

        print("\nResults:")
        print(f"  - h_ex: {self.ex.h}")
        print(f"  - T_ex: {self.ex.T}")
        print("=========================")

    def print_states_connectors(self):
        print("=== Expander Results ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")

