# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: elise neven
"""
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import numpy as np

from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
from labothappy.connector.heat_connector import HeatConnector

from labothappy.machine.circuit import Circuit
from machine.boundary_conditions.mass_source import MassSource
from machine.boundary_conditions.mass_sink import MassSink

from labothappy.component.heat_exchanger.steady_state.epsilon_NTU.simulation_model import HXeNTU
from labothappy.component.volumetric_machine.expander.steady_state.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from labothappy.component.pump.steady_state.constant_efficiency.simulation_model import PumpCstEff

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

            # Set the connector in both the source and target component
            if hasattr(self.model, output_connector.split('-')[1]):
                setattr(self.model, output_connector.split('-')[1], connector)
            else:
                raise ValueError(f"Component '{self.name}' does not have a '{output_connector.split('-')[1]}' port")
            
            if hasattr(target_component.model, input_connector.split('-')[1]):
                setattr(target_component.model, input_connector.split('-')[1], connector)
            else:
                raise ValueError(f"Target component '{target_component.name}' does not have a '{input_connector.split('-')[1]}' port")
            
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

    class Source():
        def __init__(self, name, properties, target_component, next_comp_input_port):
            self.name = name
            self.properties = properties
            self.next = {}
            self.link(target_component, next_comp_input_port)

        def set_properties(self, **kwargs):
            self.properties.set_properties(**kwargs)

        def link(self, target_component, next_comp_input_port):
            connector_type = next_comp_input_port.split('-')[0]
 
            if connector_type != "m":  # Mass connector
                print("Source shall be connected by a mass connector")
                return
            else:
                setattr(target_component.model, next_comp_input_port.split('-')[1], self.properties) # Voir si ça fait juste référence ou si ça crée un nouvel objet    
                self.next[target_component.name] = target_component.model
                target_component.add_previous(next_comp_input_port, self)
                print(f"Linked source {self.name} to {target_component.name}.{next_comp_input_port}")


    class Sink():
        def __init__(self, name, target_component, prev_comp_output_port):
            self.name = name
            self.properties = MassConnector()
            self.previous = {}
            self.link(target_component, prev_comp_output_port)
 
        def set_properties(self, port_name, **kwargs):
            port = getattr(self.model, port_name)
            port.set_properties(**kwargs)
 
        def link(self, target_component, output_port):
            connector_type = output_port.split('-')[0]
            if connector_type != "m":  # Mass connector
                print("Source shall be connected by a mass connector")
                return
            else:                
                self.previous[target_component.name] = target_component.model
                target_component.add_next(output_port, self)

    def __init__(self, fluid=None):
        
        # Building Blocks
        self.components = {}  # Store components using a dictionary for easy access
        self.sources = {}
        self.sinks = {}

        # Properties and guesses
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

    def add_source(self, name, connector, next_comp, next_comp_input_port):
        # Add a source to the cycle
        source = Circuit.Source(name, connector, next_comp, next_comp_input_port)
        self.sources[name] = source

    def get_source(self, name):
        # Retrieve a component by name
        if name in self.sources:
            return self.sources[name]
        raise ValueError(f"Source '{name}' not found")

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

    def set_source_properties(self, **kwargs):
        # Set properties for a specific source
        target = kwargs.pop('target')
        source = self.get_source(target)

        source.set_properties(**kwargs)

        # Update the fixed properties for the cycle
        if target not in self.fixed_properties:
            self.fixed_properties[target] = {}
        self.fixed_properties[target].update(kwargs)

        return 

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

if __name__ == "__main__":

    ORC = Circuit()
    
    # Create components
    Pump = PumpCstEff()
    Evaporator = HXeNTU()
    Condenser = HXeNTU()
    Expander = ExpanderCstEff()

    # Heat transfer area  
    A_htx = 0.752 # m^2
    
    # HTX dimensions
    L = 0.393 # [m] : length
    w = 0.243 # [m] : width
    h = 0.0446 # [m] : height
    l_v = 0.324 # [m] : length between ports
    casing_t = 0.005 # [m] : casing thickness # !!! arbitrary
    
    # Number and thickness of plates
    n_plates = 10
    t_plates = 0.0008 # [m] # !!! arbitrary
    
    # Fooling factor
    fooling = 98.51/1000 # (m^2*K)/W

    # Number of canals
    C_n_canals = 5
    H_n_canals = 4
    
    # Plate values 
    plate_cond = 45 # [W/(m*K)] : plate conduction
    plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
    chevron_angle = 20*np.pi/180 # !!! arbitrary

    # Total volume of each part
    V_tot = 0.9*1e-3 # [m^3]

    # Canal thickness
    C_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*C_n_canals)
    H_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*H_n_canals)

    # Total Canal Surface
    C_CS = C_canal_t*(w-2*casing_t)*C_n_canals
    H_CS = H_canal_t*(w-2*casing_t)*H_n_canals

    # Dh : hydraulic diameter
    C_Dh = (4*C_canal_t*w)/(2*C_canal_t+2*w)
    H_Dh = (4*H_canal_t*w)/(2*H_canal_t+2*w)

    # Set parameters for components
    Pump.set_parameters(eta_is=0.8)

    Evaporator.set_parameters(
        A_htx=0.752, L_HTX=0.393, V_HTX=0.9*1e-3, Flow_Type = 'CounterFlow',
        A_canal_h=H_CS, A_canal_c=C_CS, D_h=H_Dh, 
        k_plate=45, t_plate=0.0008, n_plates = 10,
        co_pitch=0.005, chevron_angle=20*np.pi/180, fouling=98.51/1000
    )
    
    Condenser.set_parameters(NTU=1.5, Cmin=1000, hot_fluid='R1233ZD', cold_fluid='R1233ZD')
    Expander.set_parameters(eta_is=0.8)

    # Add components to the cycle
    ORC.add_component(Pump, "Pump") # :!\ Est-ce que les noms des composants sont importants?
    ORC.add_component(Evaporator, "Evaporator")
    ORC.add_component(Condenser, "Condenser")
    ORC.add_component(Expander, "Expander")
    
    # Link components
    ORC.link_components("Pump", "m-ex", "Evaporator", "m-su_cold")
    ORC.link_components("Evaporator", "m-ex_cold", "Expander", "m-su")
    ORC.link_components("Expander", "m-ex", "Condenser", "m-su_hot")
    ORC.link_components("Condenser", "m-ex_hot", "Pump", "m-su")

    Hot_oil_source = MassConnector()
    ORC.add_source("Hot_Oil", Hot_oil_source, ORC.components["Evaporator"], "m-su_hot")
    ORC.set_source_properties(T=100 + 273.15, fluid='INCOMP::T66', m_dot=0.4, target='Hot_Oil', P = 4e5)


