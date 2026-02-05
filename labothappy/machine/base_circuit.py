# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: basile chaudoir
"""
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import numpy as np

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

class BaseCircuit:
    #%%
    class Component:
        def __init__(self, name, model, fluid=None):
            self.name = name
            self.model = model
            self.previous = {} # Dictionary to store preceding components by connector
            self.next = {} # Dictionary to store succeeding components by connector
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

            # print(f"Linked {self.name}.{output_connector} to {target_component.name}.{input_connector}")

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

    #%%

    class Source():
        
        # Heat, Work, Mass ?? 
        
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
                # print(f"Linked source {self.name} to {target_component.name}.{next_comp_input_port}")

    #%%

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

#%%

    def __init__(self, fluid=None):
        
        # Building Blocks
        self.components = {}  # Store components using a dictionary for easy access
        self.sources = {}
        self.sinks = {}

        # Properties and guesses
        self.print_flag = 1
        self.fluid = fluid
        self.parameters = {}
        self.converged = False

#%% Component related methods

    def add_component(self, model, name):
        
        # Check if print shall be muted
        if not self.print_flag:
            model.mute_print()

        # Add a component to the cycle        
        component = BaseCircuit.Component(name, model, self.fluid)
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
        
#%% Source related methods

    def add_source(self, name, connector, next_comp, next_comp_input_port):
        # Add a source to the cycle
        source = BaseCircuit.Source(name, connector, next_comp, next_comp_input_port)
        self.sources[name] = source

    def get_source(self, name):
        # Retrieve a component by name
        if name in self.sources:
            return self.sources[name]
        raise ValueError(f"Source '{name}' not found")

#%%
    def set_cycle_parameters(self, **kwargs):
        # Set parameters for the cycle
        self.parameters.update(kwargs)

#%% Print related methods

    def print_states(self):
        # Print Source states
        print("\n")
        
        if self.sources != {}:
        
            print("---------------------------")
            print("---------------------------")
    
            for source in self.sources:
                source_item = self.sources[source]
                print(f"Source: {source}")
                print("---------------------------")
                source_item.properties.print_resume()
                print("---------------------------")

        if self.sinks != {}:

            print("\n")
            print("---------------------------")
            print("---------------------------")
    
            for sink in self.sinks:
                sink_item = self.sinks[sink]
                print(f"Sink: {sink}")
                print("---------------------------")
                sink_item.properties.print_resume()
                print("---------------------------")


        if self.components != {}:

            print("\n")
            print("---------------------------")
            print("---------------------------")
    
            for component in self.components:
                component_item = self.components[component]
                
                print(f"Component: {component} inlet")
                print("---------------------------")
                
                for inlet in component_item.previous:
                    print(f"{inlet}:")
                    print("\n")
                    input_name = inlet.split('-')[1]
                    # Dynamically access the attribute using getattr
                    connector_in = getattr(self.components[component].model, input_name, None)
                    
                    if connector_in:  # Check if the connector exists (not None)
                        connector_in.print_resume()  # Assuming `print_resume` is a method of the connector
                        
                    print("---------------------------")

                print(f"Component: {component} outlet")
                print("---------------------------")
                
                for outlet in component_item.next:
                    print(f"{outlet}:")
                    print("\n")

                    output_name = outlet.split('-')[1]
                    
                    # Dynamically access the attribute using getattr
                    connector_out = getattr(self.components[component].model, output_name, None)
                    
                    if connector_out:  # Check if the connector exists (not None)
                        connector_out.print_resume()  # Assuming `print_resume` is a method of the connector
                        
                    print("---------------------------")
                    
        return

    def mute_print(self):
        self.print_flag = 0
        
        for component_name in self.components:
            self.components[component_name].model.mute_print()
        
        return
    
#%% Ph-Plot related methods

    def get_p_h(self, component_name, h, p, start_flag):
        
        if start_flag == 0:
            if component_name == self.solve_start_components[0]:
                return h, p
        
        component = self.components[component_name]
        component_model = component.model
        
        for input_port in component.previous:
            
            connector_type, input_port_name = input_port.split("-")
            
            if connector_type == 'm':
                connector = getattr(self.components[component_name].model, input_port_name, None)
            
            if connector.fluid == self.fluid:
                h.append(connector.h)
                p.append(connector.p)                    
        
        for output_port in component.next:
            
            connector_type, output_port_name = output_port.split("-")
            
            if connector_type == 'm':
                connector = getattr(self.components[component_name].model, output_port_name, None)
            
                if connector.fluid == self.fluid:
                    h.append(connector.h)
                    p.append(connector.p)        
                    
                    next_comp_name = self.components[component_name].next[output_port].name
                    
                    self.get_p_h(next_comp_name, h, p, 0)
                
        return h, p

    def ph_plot(self):
        
        h = []
        p = []
        
        h, p = self.get_p_h(self.solve_start_components[0], h, p, 1)
    
        plt.figure()
    
        for i in range(len(h)):
            plt.plot(h,p)
    
        return

#%% Ts-Plot related methods

    def plot_cycle_Ts(self, saturation_curve = True):
        
        fig = plt.figure()
        
        if saturation_curve:
            def generate_saturation_curve(fluid, n_points=100):
                """
                Generates saturation curve arrays (T, s_liq, s_vap) for the fluid+suffix.
                """
            
                fluid_name = fluid
            
                # Get saturation temperature range
                T_crit = PropsSI('TCRIT', fluid_name)
                T_triple = PropsSI('Ttriple', fluid_name)
            
                # Avoid extremely low T
                T_min = max(T_triple, 0.1 * T_crit)
                T_max = 1 * T_crit  # avoid critical point
                T_sat = np.linspace(T_min, T_max, n_points)
            
                s_liq = np.zeros_like(T_sat)
                s_vap = np.zeros_like(T_sat)
            
                for i, T in enumerate(T_sat):
                    try:
                        s_liq[i] = PropsSI('S', 'T', T, 'Q', 0, fluid_name)  # saturated liquid entropy
                        s_vap[i] = PropsSI('S', 'T', T, 'Q', 1, fluid_name)  # saturated vapor entropy
                    except:
                        s_liq[i] = np.nan
                        s_vap[i] = np.nan
            
                return T_sat, s_liq, s_vap
            
            T_sat, s_liq, s_vap = generate_saturation_curve(self.fluid)
            
            plt.plot(s_liq, T_sat, 'k--')  # saturated liquid
            plt.plot(s_vap, T_sat, 'k--')  # saturated vapor
        
        for comp in self.components:
            model = self.components[comp].model
            
            su_C_flag = 0
            su_H_flag = 0
            
            if hasattr(model, "su_C"):
                su_C_flag = 1
                
            if hasattr(model, "su_H"):
                su_H_flag = 1
            
            if su_C_flag + su_H_flag > 0: # semi or total HX
                if su_C_flag + su_H_flag == 1:
                    
                    if su_C_flag: # Semi HX Case
                        fig = model.plot_Ts(fig = fig, choose_HX_side = 'C')
                    else:
                        fig = model.plot_Ts(fig = fig, choose_HX_side = 'H')

                else:
                    
                    fluid_H = getattr(model, "su_H").fluid
                    fluid_C = getattr(model, "su_C").fluid
                    
                    if fluid_H == fluid_C:
                        dominant_side = None
                    elif fluid_C == self.fluid:
                        dominant_side = 'C'  
                    else:
                        dominant_side = 'H'
                    
                    fig = model.plot_Ts(fig = fig, choose_HX_side = dominant_side)
                    
            else:
                fig = model.plot_Ts(fig = fig)
        
        return
