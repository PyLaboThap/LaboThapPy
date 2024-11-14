# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: elise neven
"""

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from component.steadystate.heat_exchanger.pinch_cst.simulation_model import HXPinchCst
from component.steadystate.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.steadystate.pump.constant_efficiency.simulation_model import PumpCstEff

from CoolProp.CoolProp import PropsSI

class Circuit:
    class Component:
        def __init__(self, name, model, fluid=None):
            self.name = name
            self.model = model
            self.previous = {}  # Dictionary to store preceding components by connector
            self.next = {}  # Dictionary to store succeeding components by connector
            self.fluid = fluid

        def add_previous(self, connector, component):
            self.previous[connector] = component 

        def add_next(self, connector, component):
            self.next[connector] = component

        def link(self, output_connector, target_component, input_connector):
            # Determine the type of connector based on the output connector
            connector_type = output_connector.split('-')[0]

            if connector_type == "m":  # Mass connector
                connector = MassConnector(fluid=self.fluid)
            elif connector_type == "q":  # Heat connector
                connector = HeatConnector()
            elif connector_type == "w":  # Work connector
                connector = WorkConnector()
            else:
                raise ValueError(f"Unknown connector type: {connector_type}")

            # Set the connector in both the source and target components
            setattr(self.model, output_connector.split('-')[1], connector)
            setattr(target_component.model, input_connector.split('-')[1], connector)

            # Add the connection to the tracking dictionaries
            self.add_next(output_connector, target_component)
            target_component.add_previous(input_connector, self)

            print(f"Linked {self.name}.{output_connector} to {target_component.name}.{input_connector}")

        def set_properties(self, connector_name, **kwargs):
            # Set properties for a specific connector
            connector = getattr(self.model, connector_name)
            connector.set_properties(**kwargs)

        def solve(self):
            # Solve the component if it is calculable
            self.model.check_calculable()
            if self.model.calculable:
                self.model.solve()
            else:
                print(f"{self.name} is not calculable")
                self.model.print_states_connectors()

    # class Source():
    #     def __init__(self, name, target_component, input_port):
    #         self.name = name
    #         self.properties = MassConnector()
    #         self.next = {}
    #         self.previous = {}
    #         self.link(target_component, input_port)
    #     def set_properties(self, port_name, **kwargs):
    #         port = getattr(self.model, port_name)
    #         port.set_properties(**kwargs)
    #     def link(self, target_component, input_port):
    #         connector_type = input_port.split('-')[0]
 
    #         if connector_type != "m":  # Mass connector
    #             print("Source shall be connected by a mass connector")
    #             return
    #         else:
    #             setattr(target_component.model, input_port.split('-')[1], self.properties) # Voir si ça fait juste référence ou si ça crée un nouvel objet    
    #             self.next[target_component.name] = target_component.model
    #             target_component.add_previous(input_port, self)
 
    # class Sink():
    #     def __init__(self, name, target_component, output_port):
    #         self.name = name
    #         self.properties = MassConnector()
    #         self.previous = {}
    #         self.next = {}
    #         self.link(target_component, output_port)
 
    #     def set_properties(self, port_name, **kwargs):
    #         port = getattr(self.model, port_name)
    #         port.set_properties(**kwargs)
 
    #     def link(self, target_component, output_port):
    #         connector_type = output_port.split('-')[0]
    #         if connector_type != "m":  # Mass connector
    #             print("Source shall be connected by a mass connector")
    #             return
    #         else:                
    #             self.previous[target_component.name] = target_component.model
    #             target_component.add_next(output_port, self)


    def __init__(self, fluid=None):
        self.components = {}  # Store components using a dictionary for easy access
        self.fixed_properties = {}
        self.guesses = {}
        self.parameters = {}
        self.fluid = fluid

    def add_component(self, model, name):
        # Add a component to the cycle
        component = Circuit.Component(name, model, self.fluid)
        self.components[name] = component

    def get_component(self, name):
        # Retrieve a component by name
        if name in self.components:
            return self.components[name]
        raise ValueError(f"Component '{name}' not found")

    def link_components(self, component1_name, output_connector, component2_name, input_connector):
        # Link two components through specified connectors
        component1 = self.get_component(component1_name)
        component2 = self.get_component(component2_name)
        component1.link(output_connector, component2, input_connector)

    def set_cycle_properties(self, **kwargs):
        # Set properties for a specific component and connector
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        component.set_properties(connector_name, **kwargs)

        # Update the fixed properties for the cycle
        if target not in self.fixed_properties:
            self.fixed_properties[target] = {}
        self.fixed_properties[target].update(kwargs)

    def set_cycle_guess(self, **kwargs):
        # Set initial guesses for a specific component and connector
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        
        component.set_properties(connector_name, **kwargs)

        # Update the guesses for the cycle
        if target not in self.guesses:
            self.guesses[target] = {}
        self.guesses[target].update(kwargs)

    def set_cycle_parameters(self, **kwargs):
        # Set parameters for the cycle
        self.parameters.update(kwargs)
