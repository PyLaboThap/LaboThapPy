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

class Circuit:
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

    class Residual_variable():
        def __init__(self, target, variable, tolerance):

            self.name = target + "-" + variable
            self.target = target 
            self.variable = variable
            self.value = None
            self.tolerance = tolerance
            self.converged = False
            
        def check_convergence(self, new_value):
            
            if (abs(new_value - self.value)/self.value) < self.tolerance:
                self.converged = True
                return True
            else: 
                self.set_value(new_value)
                self.converged = False
                return False
            return
            
        def set_value(self, value):
            self.value = value
            return

    class Guess():
        def __init__(self, target, variable, value):
    
            self.name = target + "-" + variable
            self.target = target 
            self.variable = variable
            self.value = value
            
        def get_value(self, value):
            return self.value
            
    class Fixed_Property():
        def __init__(self, target, variable, value):
    
            self.name = target + "-" + variable
            self.target = target 
            self.variable = variable
            self.value = value
            
        def get_value(self,value):
            return self.value

    class Iteration_variable():
        def __init__(self, target, variable, objective, objective_value, tol, rel, damping_factor):
            self.target = target 
            self.variable = variable
            self.objective = objective
            self.objective_value = objective_value
            self.tol = tol
            self.converged = False
            self.damping_factor = damping_factor
            
            if rel == 1 or rel == -1:
                self.rel = rel
            else:
                raise ValueError("'rel' value for 'Iteration_variable' class shall be either 1 or -1.")
            
        def check(self, target_value, objective_connector, var):
            if abs((target_value - self.objective_value)/self.objective_value) < self.tol:
                self.converged = True
                return None
            else:
                self.converged = False 
                
                if var == 'SC':
                    T_sat_desired = objective_connector.T + self.objective_value
                    target_variable_desired_value = PropsSI(self.variable.upper(), 'T', T_sat_desired, 'Q', 0.5, objective_connector.fluid)
                    
                    return target_variable_desired_value          
                    
        
#%%

    def __init__(self, fluid=None):
        
        # Building Blocks
        self.components = {}  # Store components using a dictionary for easy access
        self.sources = {}
        self.sinks = {}

        # Properties and guesses
        self.fluid = fluid
        self.fixed_properties = {}
        self.parameters = {}

        # Solve related variables
        self.solve_start_components = []
        self.guesses = {}
        self.guess_update = False
        self.res_vars = {}
        self.it_vars = []
        self.converged = False

#%% Component related methods

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

#%% Source related methods

    def add_source(self, name, connector, next_comp, next_comp_input_port):
        # Add a source to the cycle
        source = Circuit.Source(name, connector, next_comp, next_comp_input_port)
        self.sources[name] = source

    def get_source(self, name):
        # Retrieve a component by name
        if name in self.sources:
            return self.sources[name]
        raise ValueError(f"Source '{name}' not found")

#%% Variable set related methods

    def set_fixed_properties(self, **kwargs):

        # Set properties for a specific component and connector
        target = kwargs.pop('target')
                
        try:
            SC = kwargs.pop('SC')
            fix_prop = Circuit.Fixed_Property(target, 'SC', SC)
            
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop

            component_name, connector_name = target.split(':')
            connector = getattr(self.components[component_name].model, connector_name, None)
            
            if connector.p is not None:
                T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
                component = self.get_component(component_name)
                component.set_properties(connector_name, T = T_sat - SC)
                
        except: 
            pass
        
        try:
            SH = kwargs.pop('SH')
            fix_prop = Circuit.Fixed_Property(target, 'SH', SH)
            
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop

            component_name, connector_name = target.split(':')
            connector = getattr(self.components[component_name].model, connector_name, None)
            
            if connector.p is not None:
                T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
                component = self.get_component(component_name)
                component.set_properties(connector_name, T = T_sat + SH)
                
        except:
            pass

        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)

        for arg_name in kwargs:
            value = kwargs[arg_name]
            fix_prop = Circuit.Fixed_Property(target, arg_name, value)
            
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop
                        
            if arg_name == 'P' or arg_name == 'p':
                if target + '-SC' in self.guesses:
                    T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat - self.guesses[target + '-SC'].value)

                if target + '-SC' in self.fixed_properties:
                    T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat - self.fixed_properties[target + '-SC'].value)
                        
                if target + '-SH' in self.guesses:
                    T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat + self.guesses[target + '-SH'].value)

                if target + '-SH' in self.fixed_properties:
                    T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat + self.fixed_properties[target + '-SH'].value)

            
            if arg_name == 'T':
                if target + '-SC' in self.guesses:
                    T_sat = fix_prop.value + self.guesses[target + '-SC']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)

                if target + '-SC' in self.fixed_properties:
                    T_sat = fix_prop.value + self.fixed_properties[target + '-SC']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)
            
                if target + '-SH' in self.guesses:
                    T_sat = fix_prop.value - self.guesses[target + '-SH']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)

                if target + '-SH' in self.fixed_properties:
                    T_sat = fix_prop.value - self.fixed_properties[target + '-SH']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)
            
        component.set_properties(connector_name, **kwargs)

        component.model.check_calculable()

        return

    def set_source_properties(self, **kwargs):
        # Set properties for a specific source
        target = kwargs.pop('target')
        source = self.get_source(target)

        source.set_properties(**kwargs)

        for arg_name in kwargs:
            value = kwargs[arg_name]
            fix_prop = Circuit.Fixed_Property(target, arg_name, value)
            
            # Update the fixed properties for the cycle
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop

        return 

#%%

    def set_cycle_guess(self, external_set = False, **kwargs):
        # Set properties for a specific component and connector
        target = kwargs.pop('target')
                
        try:
            SC = kwargs.pop('SC')
            guess = Circuit.Guess(target, 'SC', SC)
            
            if self.guess_update == True or guess.name not in self.guesses or external_set == True:
                self.guesses[guess.name] = guess

            component_name, connector_name = target.split(':')
            connector = getattr(self.components[component_name].model, connector_name, None)
            
            if connector.p is not None:
                T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
                component = self.get_component(component_name)
                component.set_properties(connector_name, T = T_sat - SC)
                
        except: 
            pass
        
        try:
            SH = kwargs.pop('SH')
            guess = Circuit.Guess(target, 'SH', SH)
            
            if self.guess_update == True or guess.name not in self.guesses or external_set == True:
                self.guesses[guess.name] = guess

            component_name, connector_name = target.split(':')
            connector = getattr(self.components[component_name].model, connector_name, None)
            
            if connector.p is not None:
                T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
                component = self.get_component(component_name)
                component.set_properties(connector_name, T = T_sat + SH)
                
        except:
            pass

        for arg_name in kwargs:
            value = kwargs[arg_name]
            guess = Circuit.Guess(target, arg_name, value)

            component_name, connector_name = target.split(':')
            connector = getattr(self.components[component_name].model, connector_name, None)
            
            if self.guess_update == True or guess.name not in self.guesses or external_set == True:
                self.guesses[guess.name] = guess
                        
            if arg_name.upper() == 'P':
                if target + '-SC' in self.guesses:
                    T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat - self.guesses[target + '-SC'].value)

                if target + '-SC' in self.fixed_properties:
                    T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat - self.fixed_properties[target + '-SC'].value)
                        
                if target + '-SH' in self.guesses:
                    T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat + self.guesses[target + '-SH'].value)

                if target + '-SH' in self.fixed_properties:
                    T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat + self.fixed_properties[target + '-SH'].value)

            if arg_name == 'T':
                if target + '-SC' in self.guesses:
                    T_sat = guess.value + self.guesses[target + '-SC']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)

                if target + '-SC' in self.fixed_properties:
                    T_sat = guess.value + self.fixed_properties[target + '-SC']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)
            
                if target + '-SH' in self.guesses:
                    T_sat = guess.value - self.guesses[target + '-SH']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)

                if target + '-SH' in self.fixed_properties:
                    T_sat = guess.value - self.fixed_properties[target + '-SH']
                    p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, p = p_sat)

        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        component.set_properties(connector_name, **kwargs)

        component.model.check_calculable()

        return

    def set_iteration_variable(self, tol = 1e-2, target = None, variable = None, objective = None, rel = 1, damping_factor = 0.5):
        """
        Parameters
        ----------
        tol : double, optional
            Tolerance for checking iteration objective value is at its targetted value. The default is 1e-2.
        target : string
            Iteration variable component and port, shall already in guesses.
        variable : string
            Iteration physical property.
        objective : string
            Objective component, port and property, shall already be present in the fixed_variables.
        rel : 1 or -1
            Relationship bewteen iteration variable and objective variable.
            1 if they are correlated
            -1 if they have inverse correlation
        """
        
        if objective in self.fixed_properties:
            it_var = Circuit.Iteration_variable(target, variable, objective, self.fixed_properties[objective].value, tol, rel, damping_factor)
            self.it_vars.append(it_var) 
        else:
            raise ValueError(f"{objective} not part of the fixed_properties. Cannot be defined as an objective")
        
        return

#%%

    def print_guesses(self):
        for guess in self.guesses:
            print(f"{guess}: {self.guesses[guess].value}")
        
        return

    def print_fixed_properties(self):
        for fix_prop in self.fixed_properties:
            print(f"{fix_prop}: {self.fixed_properties[fix_prop].value}")
        
        return

#%%

    def set_cycle_parameters(self, **kwargs):
        # Set parameters for the cycle
        self.parameters.update(kwargs)
        
    def set_residual_variable(self, target, variable, tolerance):
        
        res_var = Circuit.Residual_variable(target, variable, tolerance)
        self.res_vars[res_var.name] = res_var
        return

    def print_res_vars(self):
        for res_var in self.res_vars:
            print(f"{res_var}: {self.res_vars[res_var].value}")
        
        return

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

    def get_T_s_p(self, component_name, T, s, p, start_flag):
        
        if start_flag == 0:
            if component_name == self.solve_start_components[0]:
                return T, s, p
        
        component = self.components[component_name]
        component_model = component.model
        
        for input_port in component.previous:
            
            connector_type, input_port_name = input_port.split("-")
            
            if connector_type == 'm':
                connector = getattr(self.components[component_name].model, input_port_name, None)
            
            if connector.fluid == self.fluid:
                T.append(connector.T)
                s.append(connector.s)                    
                p.append(connector.p)                    
        
        for output_port in component.next:
            
            connector_type, output_port_name = output_port.split("-")
            
            if connector_type == 'm':
                connector = getattr(self.components[component_name].model, output_port_name, None)
            
                if connector.fluid == self.fluid:
                    T.append(connector.T)
                    s.append(connector.s)        
                    p.append(connector.p)   
                    
                    next_comp_name = self.components[component_name].next[output_port].name
                    
                    self.get_T_s_p(next_comp_name, T, s, p, 0)
                
        return T, s, p

    def Ts_plot(self):
        
        T = []
        s = []
        p = []
        
        T, s, p = self.get_T_s_p(self.solve_start_components[0], T, s, p, 1)
    
        plt.figure()
            
        for i in range(len(T)):

            print(i)
            
            T_sat = PropsSI('T', 'P', p[i], 'Q', 0.5, self.fluid)

            s_sat_0 = PropsSI('S', 'P', p[i], 'Q', 0, self.fluid)  # Liquid
            s_sat_1 = PropsSI('S', 'P', p[i], 'Q', 1, self.fluid)

            if i < len(T) - 1:
                
                if (T[i] >= T_sat and T[i+1] <= T_sat):
                    
                    plt.plot([s[i], s_sat_1], [T[i], T_sat])
                    plt.plot([s_sat_1, s_sat_0], [T_sat, T_sat])
                    plt.plot([s_sat_0, s[i+1]], [T_sat, T[i+1]])
                                        
                elif (T[i] <= T_sat and T[i+1] >= T_sat):
                    
                    plt.plot([s[i], s_sat_0], [T[i], T_sat])
                    plt.plot([s_sat_0, s_sat_1], [T_sat, T_sat])
                    plt.plot([s_sat_1, s[i+1]], [T_sat, T[i+1]])

                else:
                    plt.plot([s[i], s[i+1]], [T[i],T[i+1]])
    
            else:
                plt.plot([s[i], s[i+1]], [T[i],T[i+1]])
                
        plt.axis([min(s)-100, max(s)+100, min(T)-10, max(T)+10])
        plt.show()

        return

#%% Solve related methods

    def check_all_component_solved(self):
        for component in self.components:
            if self.components[component].model.solved == False:
                return False
        
        return True

    def reset_solved_marker(self):
        for component in self.components:
            self.components[component].model.solved = False
            
        return

    def recursive_solve(self, component_name):
        # print(f"Recursive Solve : {component_name}")
                
        component = self.get_component(component_name)
        component_model = component.model

        # if component_name == "Evaporator":
        #     print(component_model.su_C.T)
        #     print(component_model.su_C.h)
        #     print(component_model.su_C.p)
        #     print(component_model.su_C.m_dot)

        # if component_name == self.solve_start_component:
        if component_model.solved:
            # print("Back to a solved component.")
            return

        component_model.check_calculable()
        component_model.check_parametrized()
        
        if component_model.parametrized:
            if component_model.calculable:        
                component_model.solve()
                # print(f"Component '{component_name}' solved.")
            else:
                # print(f"Component '{component_name}' not calculable.")
                return
        else:
            raise ValueError(f"Component '{component_name}' not parametrized.")
        
        if component_name == 'Mixer':
            print(f"su_1 m_dot : {component_model.su_1.m_dot}")
            print(f"su_2 m_dot : {component_model.su_2.m_dot}")

            print(f"----------------------")
            print(f"ex m_dot : {component_model.ex.m_dot}")

            print(f"----------------------------------------------")

        
        for output_port in component.next:
            next_comp_name = component.next[output_port].name
            self.recursive_solve(next_comp_name)
        
        return

    def solve(self):
        
        # print("Solve Start")
        
        # print("Find first calculable component")
        
        for component in self.components:     
            if self.components[component].model.check_calculable():
                self.solve_start_components.append(component)
        
        # print(f"First calculable components: {self.solve_start_components}")
        
        # print("First iteration")
        
        for start_component in self.solve_start_components:
            self.recursive_solve(start_component)
            if self.check_all_component_solved():
                break
        
        if self.check_all_component_solved():
            print("All components solved in first iteration.")
            pass
        else:
            print("Not all components were solved in first iteration.")
            return
        
        # print("Set residual variables")
        
        for res_var_name in self.res_vars:
            component_res_var, rest = res_var_name.split(":")
            port_res_var, var_res_var = rest.split("-")

            connector = getattr(self.components[component_res_var].model, port_res_var)
            value = float(getattr(connector, var_res_var))
                                    
            self.res_vars[res_var_name].set_value(value)
        
        # print("Print guesses")
        # print("-------------")

        # self.print_guesses()
        
        # print("-------------")

        self.guess_update = True

        i=0
        n_it_max = 100

        # m_dot_su_pp = []
        # m_dot_ex_pp = []
        # m_dot_ex_cd = []
        # m_dot_su_spli = []
        # m_dot_ex_mix = []
        # P_cd = []        

        h_in_cp = []
        h_in_gc = []
        h_in_vv = []
        h_in_ev = []
        res_ev = []
        
        while i < n_it_max:
            
            # h_in_cp.append(self.components['Compressor'].model.su.p)
            # h_in_gc.append(self.components['GasCooler'].model.su_H.p)
            # h_in_vv.append(self.components['Valve'].model.su.p)
            # h_in_ev.append(self.components['Evaporator'].model.su_C.p)
            # res_ev.append(self.components['Evaporator'].model.res)
            
            # P_cd.append(self.components['Pump'].model.su.p)
            
            for it_var in self.it_vars:
                obj_comp, rest = it_var.objective.split(":")
                port, var = rest.split("-")
                
                if var == "SC":
                    connector = getattr(self.components[obj_comp].model, port)
                    T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                    
                    if connector.T <= T_sat:
                        SC = T_sat - connector.T               
                        desired_value = it_var.check(SC, connector, var)
                        
                        if desired_value == None:
                            delta = 0
                        else:
                            target_comp, port = it_var.target[0].split(":")
                            target_connector = getattr(self.components[target_comp].model, port)
                            target_current_value = getattr(target_connector, it_var.variable)
                            
                            new_target_value = (desired_value - target_current_value)*it_var.damping_factor + target_current_value
                            
                            for target in it_var.target:
                                kwargs_dict = {'target': target,
                                                it_var.variable : new_target_value}
                        
                                self.set_cycle_guess(**kwargs_dict)
                            

            for guess in self.guesses:
                component_name, rest = guess.split(":")
                port, var = rest.split("-")

                if var != "SC" and var != "SH":    
                    connector = getattr(self.components[component_name].model, port)
                    value = float(getattr(connector, var))
                    
                    kwargs_dict = {'target': self.guesses[guess].target,
                                    self.guesses[guess].variable : value}
            
                    self.set_cycle_guess(**kwargs_dict)
                
                else:
                    if var == "SC":
                        connector = getattr(self.components[component_name].model, port)    
                        T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                        
                        SC_val = max(T_sat - connector.T,0)
                                                
                        kwargs_dict = {'target': self.guesses[guess].target,
                                        self.guesses[guess].variable : SC_val}
                
                        self.set_cycle_guess(**kwargs_dict) 

                    if var == "SH":
                        connector = getattr(self.components[component_name].model, port)    
                        T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                        
                        SH_val = max(connector.T-T_sat,0)
                                                
                        kwargs_dict = {'target': self.guesses[guess].target,
                                        self.guesses[guess].variable : SH_val}
                
                        self.set_cycle_guess(**kwargs_dict) 
                    
                    
                # self.set_fixed_properties()
                        
            for fix_prop in self.fixed_properties:
                if ":" in fix_prop:
                    component_name, rest = guess.split(":")
                    port, var = rest.split("-")
                    
                    if var == "SC":
                        pass
                        
                    else:                        
                        kwargs_dict = {'target': self.fixed_properties[fix_prop].target,
                                        self.fixed_properties[fix_prop].variable : self.fixed_properties[fix_prop].value}
                        
                        # print(self.fixed_properties)
                        self.set_fixed_properties(**kwargs_dict)
                        # print(self.fixed_properties)

                else:
                    pass

            # print("-------------")    
            # print("Print guesses")
            # print("-------------")

            # self.print_guesses()
    
            # print("-------------")

            # print("Print residual vars")
            # print("-------------")

            # self.print_res_vars()

            # print("-------------")

            self.reset_solved_marker()
    
            for start_component in self.solve_start_components:
                self.recursive_solve(start_component)
                if self.check_all_component_solved():
                    break       
    
            # print("Check res_vars convergence and set new values if needed")

            self.converged = True

            for res_var in self.res_vars:    
                component_name, rest = res_var.split(":")
                port, var = rest.split("-")
    
                connector = getattr(self.components[component_name].model, port)
                value = float(getattr(connector, var))
                
                self.res_vars[res_var].check_convergence(value)
                
                if not self.res_vars[res_var].converged:
                    self.converged = False
                # print(self.res_vars[res_var].converged)
            
            for it_var in self.it_vars:    
                if not it_var.converged:
                    self.converged = False            
            
            # m_dot_su_pp.append(self.components['Pump'].model.su.m_dot)
            # m_dot_ex_pp.append(self.components['Pump'].model.ex.m_dot)
            # m_dot_su_spli.append(self.components['Spliter'].model.su.m_dot)
            # m_dot_ex_mix.append(self.components['Mixer'].model.ex.m_dot)
            # m_dot_ex_cd.append(self.components['Condenser'].model.ex_H.m_dot)
        
            if self.converged:
            #     plt.figure()
                
            #     plt.plot(m_dot_su_pp, 'r')
            #     plt.plot(m_dot_ex_pp, 'g',  marker='o')
            #     plt.plot(m_dot_su_spli, 'b')
            #     plt.plot(m_dot_ex_mix, 'k',  marker='o')
            #     plt.plot(m_dot_ex_cd, 'orange',  marker='o')
                
            #     plt.legend(['su_pp', 'ex_pp', 'su_spli', 'ex_mix', 'ex_cd'])
                
            #     plt.figure()
            #     plt.plot(P_cd, 'r')
                
            #     plt.legend(['P_cd'])

                # plt.figure()

                # plt.plot(h_in_cp, 'r')
                # plt.plot(h_in_gc, 'g',  marker='o')
                # plt.plot(h_in_vv, 'b')
                # plt.plot(h_in_ev, 'k',  marker='o')                            

                # plt.legend(['su_cp', 'su_gc', 'su_vv', 'su_ev'])

                # plt.show()

                print(f"Solver successfully converged in {i+1} iteration(s) !")
                return
            
            i = i + 1
        
        plt.figure()

        plt.plot(h_in_cp, 'r')
        plt.plot(h_in_gc, 'g',  marker='o')
        plt.plot(h_in_vv, 'b')
        plt.plot(h_in_ev, 'k',  marker='o')                            

        plt.legend(['su_cp', 'su_gc', 'su_vv', 'su_ev'])

        plt.show()

        # plt.figure()
        # plt.plot(res_ev, 'r',  marker='o')                            
        # plt.show()

        print(f"Solver failed to converge in {n_it_max} iterations.")
        
        # self.print_res_vars()
        
        # plt.figure()
        
        # plt.plot(m_dot_su_pp, 'r')
        # plt.plot(m_dot_ex_pp, 'g',  marker='o')
        # plt.plot(m_dot_su_spli, 'b')
        # plt.plot(m_dot_ex_mix, 'k',  marker='o')
        # plt.plot(m_dot_ex_cd, 'orange',  marker='o')
        
        # plt.legend(['su_pp', 'ex_pp', 'su_spli', 'ex_mix', 'ex_cd'])
        
        # plt.figure()
        # plt.plot(P_cd, 'r')
        
        # plt.legend(['P_cd'])
        
        return
    
# # -*- coding: utf-8 -*-
# """
# Created on Wed Jul 07 11:47:52 2024
    
# @author: basile chaudoir
# """
# from CoolProp.CoolProp import PropsSI

# import matplotlib.pyplot as plt
# import numpy as np

# from connector.mass_connector import MassConnector
# from connector.work_connector import WorkConnector
# from connector.heat_connector import HeatConnector

# class Circuit:
#     #%%
    
#     class Component:
#         def __init__(self, name, model, fluid=None):
#             self.name = name
#             self.model = model
#             self.previous = {} # Dictionary to store preceding components by connector
#             self.next = {} # Dictionary to store succeeding components by connector
#             self.fluid = fluid

#         def add_previous(self, connector, component):
#             self.previous[connector] = component 

#         def add_next(self, connector, component):
#             self.next[connector] = component

#         def link(self, output_connector, target_component, input_connector):
#             # Determine the type of connector based on the output connector
#             connector_type = output_connector.split('-')[0]

#             if connector_type == "m":  # Mass connector
#                 connector = MassConnector(fluid=self.fluid)
#             elif connector_type == "q":  # Heat connector
#                 connector = HeatConnector()
#             elif connector_type == "w":  # Work connector
#                 connector = WorkConnector()
#             else:
#                 raise ValueError(f"Unknown connector type: {connector_type}")

#             # Set the connector in both the source and target component
#             if hasattr(self.model, output_connector.split('-')[1]):
#                 setattr(self.model, output_connector.split('-')[1], connector)
#             else:
#                 raise ValueError(f"Component '{self.name}' does not have a '{output_connector.split('-')[1]}' port")
            
#             if hasattr(target_component.model, input_connector.split('-')[1]):
#                 setattr(target_component.model, input_connector.split('-')[1], connector)
#             else:
#                 raise ValueError(f"Target component '{target_component.name}' does not have a '{input_connector.split('-')[1]}' port")
            
#             # Add the connection to the tracking dictionaries
#             self.add_next(output_connector, target_component)
#             target_component.add_previous(input_connector, self)

#             print(f"Linked {self.name}.{output_connector} to {target_component.name}.{input_connector}")

#         def set_properties(self, connector_name, **kwargs):
#             # Set properties for a specific connector
#             connector = getattr(self.model, connector_name)
#             connector.set_properties(**kwargs)

#         def solve(self):
#             # Solve the component if it is calculable
#             self.model.check_calculable()
#             if self.model.calculable:
#                 self.model.solve()
#             else:
#                 print(f"{self.name} is not calculable")
#                 self.model.print_states_connectors()

#     #%%

#     class Source():
#         def __init__(self, name, properties, target_component, next_comp_input_port):
#             self.name = name
#             self.properties = properties
#             self.next = {}
#             self.link(target_component, next_comp_input_port)

#         def set_properties(self, **kwargs):
#             self.properties.set_properties(**kwargs)

#         def link(self, target_component, next_comp_input_port):
#             connector_type = next_comp_input_port.split('-')[0]
 
#             if connector_type != "m":  # Mass connector
#                 print("Source shall be connected by a mass connector")
#                 return
#             else:
#                 setattr(target_component.model, next_comp_input_port.split('-')[1], self.properties) # Voir si ça fait juste référence ou si ça crée un nouvel objet    
#                 self.next[target_component.name] = target_component.model
#                 target_component.add_previous(next_comp_input_port, self)
#                 print(f"Linked source {self.name} to {target_component.name}.{next_comp_input_port}")

#     #%%

#     class Sink():
#         def __init__(self, name, target_component, prev_comp_output_port):
#             self.name = name
#             self.properties = MassConnector()
#             self.previous = {}
#             self.link(target_component, prev_comp_output_port)
 
#         def set_properties(self, port_name, **kwargs):
#             port = getattr(self.model, port_name)
#             port.set_properties(**kwargs)
 
#         def link(self, target_component, output_port):
#             connector_type = output_port.split('-')[0]
#             if connector_type != "m":  # Mass connector
#                 print("Source shall be connected by a mass connector")
#                 return
#             else:                
#                 self.previous[target_component.name] = target_component.model
#                 target_component.add_next(output_port, self)

#     #%%

#     class Residual_variable():
#         def __init__(self, target, variable, tolerance):

#             self.name = target + "-" + variable
#             self.target = target 
#             self.variable = variable
#             self.value = None
#             self.tolerance = tolerance
#             self.converged = False
            
#         def check_convergence(self, new_value):
            
#             if (abs(new_value - self.value)/self.value) < self.tolerance:
#                 self.converged = True
#                 return True
#             else: 
#                 self.set_value(new_value)
#                 self.converged = False
#                 return False
#             return
            
#         def set_value(self, value):
#             self.value = value
#             return

#     class Guess():
#         def __init__(self, target, variable, value):
    
#             self.name = target + "-" + variable
#             self.target = target 
#             self.variable = variable
#             self.value = value
            
#         def get_value(self, value):
#             return self.value
            
#     class Fixed_Property():
#         def __init__(self, target, variable, value):
    
#             self.name = target + "-" + variable
#             self.target = target 
#             self.variable = variable
#             self.value = value
            
#         def get_value(self,value):
#             return self.value

#     class Iteration_variable():
#         def __init__(self, target, variable, objective, objective_value, tol, rel, damping_factor):
#             self.target = target 
#             self.variable = variable
#             self.objective = objective
#             self.objective_value = objective_value
#             self.tol = tol
#             self.converged = False
#             self.damping_factor = damping_factor
            
#             if rel == 1 or rel == -1:
#                 self.rel = rel
#             else:
#                 raise ValueError("'rel' value for 'Iteration_variable' class shall be either 1 or -1.")
            
#         def check(self, target_value, objective_connector, var):
#             if abs((target_value - self.objective_value)/self.objective_value) < self.tol:
#                 self.converged = True
#                 return None
#             else:
#                 self.converged = False 
                
#                 if var == 'SC':
#                     T_sat_desired = objective_connector.T + self.objective_value
#                     target_variable_desired_value = PropsSI(self.variable.upper(), 'T', T_sat_desired, 'Q', 0.5, objective_connector.fluid)

#                     return target_variable_desired_value          
                    
        
# #%%

#     def __init__(self, fluid=None):
        
#         # Building Blocks
#         self.components = {}  # Store components using a dictionary for easy access
#         self.sources = {}
#         self.sinks = {}

#         # Properties and guesses
#         self.fluid = fluid
#         self.fixed_properties = {}
#         self.parameters = {}

#         # Solve related variables
#         self.solve_start_components = []
#         self.guesses = {}
#         self.guess_update = False
#         self.res_vars = {}
#         self.it_vars = []
#         self.converged = False

# #%% Component related methods

#     def add_component(self, model, name):
#         # Add a component to the cycle
#         component = Circuit.Component(name, model, self.fluid)
#         self.components[name] = component

#     def get_component(self, name):
#         # Retrieve a component by name
#         if name in self.components:
#             return self.components[name]
#         raise ValueError(f"Component '{name}' not found")

#     def link_components(self, component1_name, output_connector, component2_name, input_connector):
#         # Link two components through specified connectors
#         component1 = self.get_component(component1_name)
#         component2 = self.get_component(component2_name)
#         component1.link(output_connector, component2, input_connector)

# #%% Source related methods

#     def add_source(self, name, connector, next_comp, next_comp_input_port):
#         # Add a source to the cycle
#         source = Circuit.Source(name, connector, next_comp, next_comp_input_port)
#         self.sources[name] = source

#     def get_source(self, name):
#         # Retrieve a component by name
#         if name in self.sources:
#             return self.sources[name]
#         raise ValueError(f"Source '{name}' not found")

# #%% Variable set related methods

#     def set_fixed_properties(self, **kwargs):

#         # Set properties for a specific component and connector
#         target = kwargs.pop('target')
                
#         try:
#             SC = kwargs.pop('SC')
#             fix_prop = Circuit.Fixed_Property(target, 'SC', SC)
            
#             if fix_prop.name not in self.fixed_properties:
#                 self.fixed_properties[fix_prop.name] = fix_prop

#             component_name, connector_name = target.split(':')
#             connector = getattr(self.components[component_name].model, connector_name, None)
            
#             if connector.p is not None:
#                 T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
#                 component = self.get_component(component_name)
#                 component.set_properties(connector_name, T = T_sat - SC)
                
#         except: 
#             pass
        
#         try:
#             SH = kwargs.pop('SH')
#             fix_prop = Circuit.Fixed_Property(target, 'SH', SH)
            
#             if fix_prop.name not in self.fixed_properties:
#                 self.fixed_properties[fix_prop.name] = fix_prop

#             component_name, connector_name = target.split(':')
#             connector = getattr(self.components[component_name].model, connector_name, None)
            
#             if connector.p is not None:
#                 T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
#                 component = self.get_component(component_name)
#                 component.set_properties(connector_name, T = T_sat + SH)
                
#         except:
#             pass

#         component_name, connector_name = target.split(':')
#         component = self.get_component(component_name)

#         for arg_name in kwargs:
#             value = kwargs[arg_name]
#             fix_prop = Circuit.Fixed_Property(target, arg_name, value)
            
#             if fix_prop.name not in self.fixed_properties:
#                 self.fixed_properties[fix_prop.name] = fix_prop
                        
#             if arg_name == 'P' or arg_name == 'p':
#                 if target + '-SC' in self.guesses:
#                     T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat - self.guesses[target + '-SC'].value)

#                 if target + '-SC' in self.fixed_properties:
#                     T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat - self.fixed_properties[target + '-SC'].value)
                        
#                 if target + '-SH' in self.guesses:
#                     T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat + self.guesses[target + '-SH'].value)

#                 if target + '-SH' in self.fixed_properties:
#                     T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat + self.fixed_properties[target + '-SH'].value)

            
#             if arg_name == 'T':
#                 if target + '-SC' in self.guesses:
#                     T_sat = fix_prop.value + self.guesses[target + '-SC']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)

#                 if target + '-SC' in self.fixed_properties:
#                     T_sat = fix_prop.value + self.fixed_properties[target + '-SC']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)
            
#                 if target + '-SH' in self.guesses:
#                     T_sat = fix_prop.value - self.guesses[target + '-SH']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)

#                 if target + '-SH' in self.fixed_properties:
#                     T_sat = fix_prop.value - self.fixed_properties[target + '-SH']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, self.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)
            
#         component.set_properties(connector_name, **kwargs)

#         component.model.check_calculable()

#         return

#     def set_source_properties(self, **kwargs):
#         # Set properties for a specific source
#         target = kwargs.pop('target')
#         source = self.get_source(target)

#         source.set_properties(**kwargs)

#         for arg_name in kwargs:
#             value = kwargs[arg_name]
#             fix_prop = Circuit.Fixed_Property(target, arg_name, value)
            
#             # Update the fixed properties for the cycle
#             if fix_prop.name not in self.fixed_properties:
#                 self.fixed_properties[fix_prop.name] = fix_prop

#         return 

# #%%

#     def set_cycle_guess(self, **kwargs):
#         # Set properties for a specific component and connector
#         target = kwargs.pop('target')
#         try:
#             SC = kwargs.pop('SC')
#             guess = Circuit.Guess(target, 'SC', SC)
            
#             if self.guess_update == True or guess.name not in self.guesses:
#                 self.guesses[guess.name] = guess

#             component_name, connector_name = target.split(':')
#             connector = getattr(self.components[component_name].model, connector_name, None)
            
#             if connector.p is not None:
#                 T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
#                 component = self.get_component(component_name)
#                 component.set_properties(connector_name, T = T_sat - SC)
                
#         except: 
#             pass
        
#         try:
#             SH = kwargs.pop('SH')
#             guess = Circuit.Guess(target, 'SH', SH)
            
#             if self.guess_update == True or guess.name not in self.guesses:
#                 self.guesses[guess.name] = guess

#             component_name, connector_name = target.split(':')
#             connector = getattr(self.components[component_name].model, connector_name, None)
            
#             if connector.p is not None:
#                 T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                
#                 component = self.get_component(component_name)
#                 component.set_properties(connector_name, T = T_sat + SH)
                
#         except:
#             pass

#         for arg_name in kwargs:
#             value = kwargs[arg_name]
#             guess = Circuit.Guess(target, arg_name, value)
            
#             if self.guess_update == True or guess.name not in self.guesses:
#                 self.guesses[guess.name] = guess
                        
#             if arg_name == 'P':
#                 if target + '-SC' in self.guesses:
#                     T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat - self.guesses[target + '-SC'].value)

#                 if target + '-SC' in self.fixed_properties:
#                     T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat - self.fixed_properties[target + '-SC'].value)
                        
#                 if target + '-SH' in self.guesses:
#                     T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat + self.guesses[target + '-SH'].value)

#                 if target + '-SH' in self.fixed_properties:
#                     T_sat = PropsSI('T', 'P', guess.value, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, T = T_sat + self.fixed_properties[target + '-SH'].value)

#             if arg_name == 'T':
#                 if target + '-SC' in self.guesses:
#                     T_sat = guess.value + self.guesses[target + '-SC']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)

#                 if target + '-SC' in self.fixed_properties:
#                     T_sat = guess.value + self.fixed_properties[target + '-SC']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)
            
#                 if target + '-SH' in self.guesses:
#                     T_sat = guess.value - self.guesses[target + '-SH']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)

#                 if target + '-SH' in self.fixed_properties:
#                     T_sat = guess.value - self.fixed_properties[target + '-SH']
#                     p_sat = PropsSI('P', 'T', T_sat, 'Q', 0.5, connector.fluid)
                    
#                     component = self.get_component(component_name)
#                     component.set_properties(connector_name, p = p_sat)

#         component_name, connector_name = target.split(':')
#         component = self.get_component(component_name)
#         component.set_properties(connector_name, **kwargs)

#         component.model.check_calculable()
#         return

#     def set_iteration_variable(self, tol = 1e-2, target = None, variable = None, objective = None, rel = 1, damping_factor = 0.5):
#         """
#         Parameters
#         ----------
#         tol : double, optional
#             Tolerance for checking iteration objective value is at its targetted value. The default is 1e-2.
#         target : string
#             Iteration variable component and port, shall already in guesses.
#         variable : string
#             Iteration physical property.
#         objective : string
#             Objective component, port and property, shall already be present in the fixed_variables.
#         rel : 1 or -1
#             Relationship bewteen iteration variable and objective variable.
#             1 if they are correlated
#             -1 if they have inverse correlation
#         """
        
#         if objective in self.fixed_properties:
#             it_var = Circuit.Iteration_variable(target, variable, objective, self.fixed_properties[objective].value, tol, rel, damping_factor)
#             self.it_vars.append(it_var) 
#         else:
#             raise ValueError(f"{objective} not part of the fixed_properties. Cannot be defined as an objective")
        
#         return

# #%%

#     def print_guesses(self):
#         for guess in self.guesses:
#             print(f"{guess}: {self.guesses[guess].value}")
        
#         return

#     def print_fixed_properties(self):
#         for fix_prop in self.fixed_properties:
#             print(f"{fix_prop}: {self.fixed_properties[fix_prop].value}")
        
#         return

# #%%

#     def set_cycle_parameters(self, **kwargs):
#         # Set parameters for the cycle
#         self.parameters.update(kwargs)
        
#     def set_residual_variable(self, target, variable, tolerance):
        
#         res_var = Circuit.Residual_variable(target, variable, tolerance)
#         self.res_vars[res_var.name] = res_var
#         return

#     def print_res_vars(self):
#         for res_var in self.res_vars:
#             print(f"{res_var}: {self.res_vars[res_var].value}")
        
#         return

# #%% Print related methods

#     def print_states(self):
#         # Print Source states
#         print("\n")
        
#         if self.sources != {}:
        
#             print("---------------------------")
#             print("---------------------------")
    
#             for source in self.sources:
#                 source_item = self.sources[source]
#                 print(f"Source: {source}")
#                 print("---------------------------")
#                 source_item.properties.print_resume()
#                 print("---------------------------")

#         if self.sinks != {}:

#             print("\n")
#             print("---------------------------")
#             print("---------------------------")
    
#             for sink in self.sinks:
#                 sink_item = self.sinks[sink]
#                 print(f"Sink: {sink}")
#                 print("---------------------------")
#                 sink_item.properties.print_resume()
#                 print("---------------------------")


#         if self.components != {}:

#             print("\n")
#             print("---------------------------")
#             print("---------------------------")
    
#             for component in self.components:
#                 component_item = self.components[component]
                
#                 print(f"Component: {component} inlet")
#                 print("---------------------------")
                
#                 for inlet in component_item.previous:
#                     print(f"{inlet}:")
#                     print("\n")
#                     input_name = inlet.split('-')[1]
#                     # Dynamically access the attribute using getattr
#                     connector_in = getattr(self.components[component].model, input_name, None)
                    
#                     if connector_in:  # Check if the connector exists (not None)
#                         connector_in.print_resume()  # Assuming `print_resume` is a method of the connector
                        
#                     print("---------------------------")

#                 print(f"Component: {component} outlet")
#                 print("---------------------------")
                
#                 for outlet in component_item.next:
#                     print(f"{outlet}:")
#                     print("\n")

#                     output_name = outlet.split('-')[1]
                    
#                     # Dynamically access the attribute using getattr
#                     connector_out = getattr(self.components[component].model, output_name, None)
                    
#                     if connector_out:  # Check if the connector exists (not None)
#                         connector_out.print_resume()  # Assuming `print_resume` is a method of the connector
                        
#                     print("---------------------------")
                    
#         return

# #%% Ph-Plot related methods

#     def get_p_h(self, component_name, h, p, start_flag):
        
#         if start_flag == 0:
#             if component_name == self.solve_start_component:
#                 return h, p
        
#         component = self.components[component_name]
#         component_model = component.model
        
#         for input_port in component.previous:
            
#             connector_type, input_port_name = input_port.split("-")
            
#             if connector_type == 'm':
#                 connector = getattr(self.components[component_name].model, input_port_name, None)
            
#             if connector.fluid == self.fluid:
#                 h.append(connector.h)
#                 p.append(connector.p)                    
        
#         for output_port in component.next:
            
#             connector_type, output_port_name = output_port.split("-")
            
#             if connector_type == 'm':
#                 connector = getattr(self.components[component_name].model, output_port_name, None)
            
#                 if connector.fluid == self.fluid:
#                     h.append(connector.h)
#                     p.append(connector.p)        
                    
#                     next_comp_name = self.components[component_name].next[output_port].name
                    
#                     self.get_p_h(next_comp_name, h, p, 0)
                
#         return h, p

#     def ph_plot(self):
        
#         h = []
#         p = []
        
#         h, p = self.get_p_h(self.solve_start_component, h, p, 1)
    
#         plt.figure()
    
#         for i in range(len(h)):
#             plt.plot(h,p)
    
#         return

# #%% Ts-Plot related methods

#     def get_T_s_p(self, component_name, T, s, p, start_flag):
        
#         if start_flag == 0:
#             if component_name == self.solve_start_components[0]:
#                 return T, s, p
        
#         component = self.components[component_name]
#         component_model = component.model
        
#         for input_port in component.previous:
            
#             connector_type, input_port_name = input_port.split("-")
            
#             if connector_type == 'm':
#                 connector = getattr(self.components[component_name].model, input_port_name, None)
            
#             if connector.fluid == self.fluid:
#                 T.append(connector.T)
#                 s.append(connector.s)                    
#                 p.append(connector.p)                    
        
#         for output_port in component.next:
            
#             connector_type, output_port_name = output_port.split("-")
            
#             if connector_type == 'm':
#                 connector = getattr(self.components[component_name].model, output_port_name, None)
            
#                 if connector.fluid == self.fluid:
#                     T.append(connector.T)
#                     s.append(connector.s)        
#                     p.append(connector.p)   
                    
#                     next_comp_name = self.components[component_name].next[output_port].name
                    
#                     self.get_T_s_p(next_comp_name, T, s, p, 0)
                
#         return T, s, p

#     def Ts_plot(self):
        
#         T = []
#         s = []
#         p = []
        
#         T, s, p = self.get_T_s_p(self.solve_start_components[0], T, s, p, 1)
    
#         plt.figure()
            
#         for i in range(len(T)):

#             print(i)
            
#             T_sat = PropsSI('T', 'P', p[i], 'Q', 0.5, self.fluid)

#             s_sat_0 = PropsSI('S', 'P', p[i], 'Q', 0, self.fluid)  # Liquid
#             s_sat_1 = PropsSI('S', 'P', p[i], 'Q', 1, self.fluid)

#             if i < len(T) - 1:
                
#                 if (T[i] >= T_sat and T[i+1] <= T_sat):
                    
#                     plt.plot([s[i], s_sat_1], [T[i], T_sat])
#                     plt.plot([s_sat_1, s_sat_0], [T_sat, T_sat])
#                     plt.plot([s_sat_0, s[i+1]], [T_sat, T[i+1]])
                                        
#                 elif (T[i] <= T_sat and T[i+1] >= T_sat):
                    
#                     plt.plot([s[i], s_sat_0], [T[i], T_sat])
#                     plt.plot([s_sat_0, s_sat_1], [T_sat, T_sat])
#                     plt.plot([s_sat_1, s[i+1]], [T_sat, T[i+1]])

#                 else:
#                     plt.plot([s[i], s[i+1]], [T[i],T[i+1]])
    
#             else:
#                 plt.plot([s[i], s[i+1]], [T[i],T[i+1]])
                
#         plt.axis([min(s)-100, max(s)+100, min(T)-10, max(T)+10])
#         plt.show()

#         return

# #%% Solve related methods

#     def check_all_component_solved(self):
#         for component in self.components:
#             if self.components[component].model.solved == False:
#                 return False
        
#         return True

#     def reset_solved_marker(self):
#         for component in self.components:
#             self.components[component].model.solved = False
            
#         return

#     def recursive_solve(self, component_name):
#         # print(f"Recursive Solve : {component_name}")
                
#         component = self.get_component(component_name)
#         component_model = component.model

#         # if component_name == self.solve_start_component:
#         if component_model.solved:
#             # print("Back to a solved component.")
#             return

#         component_model.check_calculable()
#         component_model.check_parametrized()
        
#         if component_model.parametrized:
#             if component_model.calculable:        
#                 component_model.solve()
#                 # print(f"Component '{component_name}' solved.")
#             else:
#                 # print(f"Component '{component_name}' not calculable.")
#                 return
#         else:
#             raise ValueError(f"Component '{component_name}' not parametrized.")
        
#         for output_port in component.next:
#             next_comp_name = component.next[output_port].name
#             self.recursive_solve(next_comp_name)
        
#         return

#     def solve(self):
        
#         # print("Solve Start")
        
#         # print("Find first calculable component")
        
#         for component in self.components:     
#             if self.components[component].model.check_calculable():
#                 self.solve_start_components.append(component)
        
#         # print(f"First calculable components: {self.solve_start_components}")
        
#         # print("First iteration")
#         for start_component in self.solve_start_components:
#             self.recursive_solve(start_component)
#             if self.check_all_component_solved():
#                 break
        
#         if self.check_all_component_solved():
#             # print("All components solved in first iteration.")
#             pass
#         else:
#             raise ValueError("Not all components were solved in first iteration.")
        
#         # print("Set residual variables")
        
#         for res_var_name in self.res_vars:
#             component_res_var, rest = res_var_name.split(":")
#             port_res_var, var_res_var = rest.split("-")

#             connector = getattr(self.components[component_res_var].model, port_res_var)
#             value = float(getattr(connector, var_res_var))
                                    
#             self.res_vars[res_var_name].set_value(value)
        
#         # print("Print guesses")
#         # print("-------------")

#         # self.print_guesses()
        
#         # print("-------------")

#         self.guess_update = True

#         i=0
#         n_it_max = 50

#         m_dot_su_pp = []
#         m_dot_ex_pp = []
#         m_dot_ex_cd = []
#         m_dot_su_spli = []
#         m_dot_ex_mix = []
#         P_cd = []

#         while i < n_it_max:
            
# #            P_cd.append(self.components['Pump'].model.su.p)
            
#             for guess in self.guesses:
#                 component_name, rest = guess.split(":")
#                 port, var = rest.split("-")

#                 if var != "SC" and var != "SH":    
#                     connector = getattr(self.components[component_name].model, port)
#                     value = float(getattr(connector, var))
                    
#                     kwargs_dict = {'target': self.guesses[guess].target,
#                                     self.guesses[guess].variable : value}
            
#                     self.set_cycle_guess(**kwargs_dict)
                    
#                 # self.set_fixed_properties()

#             for it_var in self.it_vars:
#                 obj_comp, rest = it_var.objective.split(":")
#                 port, var = rest.split("-")
                
#                 if var == "SC":
#                     connector = getattr(self.components[obj_comp].model, port)
#                     T_sat = PropsSI('T', 'P', connector.p, 'Q', 0.5, connector.fluid)
                    
#                     if connector.T <= T_sat:
#                         SC = T_sat - connector.T               
#                         desired_value = it_var.check(SC, connector, var)
                        
#                         if desired_value == None:
#                             delta = 0
#                         else:
#                             target_comp, port = it_var.target[0].split(":")
#                             target_connector = getattr(self.components[target_comp].model, port)
#                             target_current_value = getattr(target_connector, it_var.variable)
                            
#                             new_target_value = (desired_value - target_current_value)*it_var.damping_factor + target_current_value
                            
#                             for target in it_var.target:
#                                 kwargs_dict = {'target': target,
#                                                 it_var.variable : new_target_value}
                        
#                                 self.set_cycle_guess(**kwargs_dict)
                            
                        
#             for fix_prop in self.fixed_properties:
#                 if ":" in fix_prop:
#                     component_name, rest = guess.split(":")
#                     port, var = rest.split("-")
                    
#                     if var == "SC":
#                         pass
                        
#                     else:                        
#                         kwargs_dict = {'target': self.fixed_properties[fix_prop].target,
#                                         self.fixed_properties[fix_prop].variable : self.fixed_properties[fix_prop].value}
                        
#                         # print(self.fixed_properties)
#                         self.set_fixed_properties(**kwargs_dict)
#                         # print(self.fixed_properties)

#                 else:
#                     pass

#             # print("-------------")    
#             # print("Print guesses")
#             # print("-------------")

#             # self.print_guesses()
    
#             # print("-------------")

#             # print("Print residual vars")
#             # print("-------------")

#             # self.print_res_vars()

#             # print("-------------")

#             self.reset_solved_marker()
    
#             for start_component in self.solve_start_components:
#                 self.recursive_solve(start_component)
#                 if self.check_all_component_solved():
#                     break       
    
#             # print("Check res_vars convergence and set new values if needed")

#             self.converged = True

#             for res_var in self.res_vars:    
#                 component_name, rest = res_var.split(":")
#                 port, var = rest.split("-")
    
#                 connector = getattr(self.components[component_name].model, port)
#                 value = float(getattr(connector, var))
                
#                 self.res_vars[res_var].check_convergence(value)
                
#                 if not self.res_vars[res_var].converged:
#                     self.converged = False
#                 # print(self.res_vars[res_var].converged)
            
#             for it_var in self.it_vars:    
#                 if not it_var.converged:
#                     self.converged = False            
            
#             # m_dot_su_pp.append(self.components['Pump'].model.su.m_dot)
#             # m_dot_ex_pp.append(self.components['Pump'].model.ex.m_dot)
#             # m_dot_su_spli.append(self.components['Spliter'].model.su.m_dot)
#             # m_dot_ex_mix.append(self.components['Mixer'].model.ex.m_dot)
#             # m_dot_ex_cd.append(self.components['Condenser'].model.ex_H.m_dot)
        
#             if self.converged:
#                 plt.figure()
                
#                 plt.plot(m_dot_su_pp, 'r')
#                 plt.plot(m_dot_ex_pp, 'g',  marker='o')
#                 plt.plot(m_dot_su_spli, 'b')
#                 plt.plot(m_dot_ex_mix, 'k',  marker='o')
#                 plt.plot(m_dot_ex_cd, 'orange',  marker='o')
                
#                 plt.legend(['su_pp', 'ex_pp', 'su_spli', 'ex_mix', 'ex_cd'])
                
#                 plt.figure()
#                 plt.plot(P_cd, 'r')
                
#                 plt.legend(['P_cd'])
                
#                 print(f"Solver successfully converged in {i} iterations !")
#                 return
            
#             i = i + 1
        
#         print(f"Solver failed to converge in {n_it_max} iterations.")
        
#         # self.print_res_vars()
        
#         plt.figure()
        
#         plt.plot(m_dot_su_pp, 'r')
#         plt.plot(m_dot_ex_pp, 'g',  marker='o')
#         plt.plot(m_dot_su_spli, 'b')
#         plt.plot(m_dot_ex_mix, 'k',  marker='o')
#         plt.plot(m_dot_ex_cd, 'orange',  marker='o')
        
#         plt.legend(['su_pp', 'ex_pp', 'su_spli', 'ex_mix', 'ex_cd'])
        
#         plt.figure()
#         plt.plot(P_cd, 'r')
        
#         plt.legend(['P_cd'])
        
#         return
    
    