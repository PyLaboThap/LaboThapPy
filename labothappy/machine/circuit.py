# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: elise neven
"""
from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import numpy as np

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from machine.boundary_conditions.mass_source import MassSource
from machine.boundary_conditions.mass_sink import MassSink

from component.heat_exchanger.steady_state.moving_boundary.charge_sensitive.simulation_model_elise import HeatExchangerMB
from component.heat_exchanger.steady_state.moving_boundary.charge_sensitive.modules.geometry_plate_hx_swep import PlateGeomSWEP
from component.volumetric_machine.expander.steady_state.semi_empirical.simulation_model import ExpanderSE
from component.volumetric_machine.pump.steady_state.polynomial_efficiency.simulation_model_Pout import PumpPolyEff
from component.tank.Spliter.simulation_model import Spliter
from component.tank.Mixer.simulation_model import Mixer

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
                print(f"Linked source {self.name} to {target_component.name}.{next_comp_input_port}")

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

    # class Guess():
    #     def __init__(self, target, variable, value):
    
    #         self.name = target + "-" + variable
    #         self.target = target 
    #         self.variable = variable
    #         self.value = None
            
    #     def get_value(self, value):
    #         return self.value
            

    class Fixed_Property():
        def __init__(self, target, variable, value):
    
            self.name = target + "-" + variable
            self.target = target 
            self.variable = variable
            self.value = value
            
        def get_value(self,value):
            return self.value

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
        self.solve_start_component = None
        self.guesses = {}
        self.res_vars = {}

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
        
        for arg_name in kwargs:
            value = kwargs[arg_name]
            fix_prop = Circuit.Fixed_Property(target, arg_name, value)
            
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop
                        
            if arg_name == 'P':
                if target + '-SC' in self.fixed_properties:
                    T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, connector.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat - self.fixed_properties[target + '-SC'].value)
            
                if target + '-SH' in self.fixed_properties:
                    T_sat = PropsSI('T', 'P', fix_prop.value, 'Q', 0.5, connector.fluid)
                    
                    component = self.get_component(component_name)
                    component.set_properties(connector_name, T = T_sat + self.fixed_properties[target + '-SH'].value)
            
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
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
        
    def set_residual_variable(self, target, variable, tolerance):
        
        res_var = Circuit.Residual_variable(target, variable, tolerance)
        self.res_vars[res_var.name] = res_var
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
            if component_name == self.solve_start_component:
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
        
        h, p = self.get_p_h(self.solve_start_component, h, p, 1)
    
        plt.figure()
    
        for i in range(len(h)):
            plt.plot(h,p)
    
        print(f"h: {h}")
        print(f"p: {p}")
      
        return

#%% Ts-Plot related methods

    def get_T_s_p(self, component_name, T, s, p, start_flag):
        
        if start_flag == 0:
            if component_name == self.solve_start_component:
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
        
        T, s, p = self.get_T_s_p(self.solve_start_component, T, s, p, 1)
    
        plt.figure()
        
        print(len(T))
    
        for i in range(len(T)):

            print(i)
            
            T_sat = PropsSI('T', 'P', p[i], 'Q', 0.5, self.fluid)

            s_sat_0 = PropsSI('S', 'P', p[i], 'Q', 0, self.fluid)  # Liquid
            s_sat_1 = PropsSI('S', 'P', p[i], 'Q', 1, self.fluid)

            print(s_sat_0)
            print(s_sat_1)

            if i < len(T) - 1:

                print("<")
                
                print(f"T[i]: {T[i]}")
                print(f"T_sat: {T_sat}")
                print(f"T[i+1]: {T[i+1]}")
                
                if (T[i] >= T_sat and T[i+1] <= T_sat):
                    
                    plt.plot([s[i], s_sat_1], [T[i], T_sat])
                    plt.plot([s_sat_1, s_sat_0], [T_sat, T_sat])
                    plt.plot([s_sat_0, s[i+1]], [T_sat, T[i+1]])
                    
                    print("Plot done 1")
                    
                elif (T[i] <= T_sat and T[i+1] >= T_sat):
                    
                    plt.plot([s[i], s_sat_0], [T[i], T_sat])
                    plt.plot([s_sat_0, s_sat_1], [T_sat, T_sat])
                    plt.plot([s_sat_1, s[i+1]], [T_sat, T[i+1]])

                    print("Plot done 2")

                else:
                    plt.plot([s[i], s[i+1]], [T[i],T[i+1]])

                    print("Plot done 3")
    
            else:
                plt.plot([s[i], s[i+1]], [T[i],T[i+1]])

                print("Plot done 4")
        
        print(f"T: {T}")
        print(f"s: {s}")
        
        print(max(s))
        
        plt.axis([min(s)-100, max(s)+100, min(T)-10, max(T)+10])
        plt.show()

        return

#%% Solve related methods

    def recursive_solve(self, component_name):
        print(f"Recursive Solve : {component_name}")
                
        component = self.get_component(component_name)
        component_model = component.model

        if component_name == self.solve_start_component:
            if component_model.solved:
                print("Back to first component, exit from iteration")
                return

        component_model.check_calculable()
        component_model.check_parametrized()
        
        if component_model.parametrized:
            if component_model.calculable:        
                component_model.solve()
                print(f"Component '{component_name}' solved.")
            else:
                print(f"Component '{component_name}' not calculable.")
                return
        else:
            raise ValueError(f"Component '{component_name}' not parametrized.")
        
        for output_port in component.next:
            next_comp_name = component.next[output_port].name
            self.recursive_solve(next_comp_name)
        
        return

    def solve(self):
        
        print("Solve Start")
        
        print("Find first calculable component")
        
        for component in self.components:     
            if self.components[component].model.calculable:
                self.solve_start_component = component
                break
        
        print(f"First calculable component: {self.solve_start_component}")
        
        
        print("First iteration")
        self.recursive_solve(self.solve_start_component)
        
        print("Set residual variables")
        
        for res_var_name in self.res_vars:
            component_res_var, rest = res_var_name.split(":")
            port_res_var, var_res_var = rest.split("-")

            connector = getattr(self.components[component_res_var].model, port_res_var)
            value = float(getattr(connector, var_res_var))
                                    
            self.res_vars[res_var_name].set_value(value)
        
        return

#%% DEFINE MODELS

if __name__ == "__main__":

    ORC = Circuit('R1233zd(E)')
    
    # Create components
    Expander_1 = ExpanderSE()
    Expander_2 = ExpanderSE()
    Expander_3 = ExpanderSE()
    
    Spliter_expander = Spliter(outlet_repartition = [1/3,1/3,1/3])
    Mixer_expander = Mixer(n_inlets = 3)
    
    Condenser = HeatExchangerMB('Plate')
    Pump = PumpPolyEff()
    Evaporator = HeatExchangerMB('Plate')

#%% CONDENSER PARAMETERS

    # Condenser
    condenser_geom = PlateGeomSWEP()
    condenser_geom.set_parameters("P200THx140/1P_Condenser")

    Condenser.set_parameters(
        # Set the geometry of the condenser
        A_c=condenser_geom.A_c, A_h=condenser_geom.A_h, h=condenser_geom.h, l=condenser_geom.l, l_v=condenser_geom.l_v, w_v=condenser_geom.w_v,
        C_CS=condenser_geom.C_CS, C_Dh=condenser_geom.C_Dh, C_V_tot=condenser_geom.C_V_tot, C_canal_t=condenser_geom.C_canal_t, C_n_canals=condenser_geom.C_n_canals,
        H_CS=condenser_geom.H_CS, H_Dh=condenser_geom.H_Dh, H_V_tot=condenser_geom.H_V_tot, H_canal_t=condenser_geom.H_canal_t, H_n_canals=condenser_geom.H_n_canals,
        casing_t=condenser_geom.casing_t, chevron_angle=condenser_geom.chevron_angle, fooling=condenser_geom.fooling,
        n_plates=condenser_geom.n_plates, plate_cond=condenser_geom.plate_cond, plate_pitch_co=condenser_geom.plate_pitch_co, t_plates=condenser_geom.t_plates, w=condenser_geom.w,
        amplitude=condenser_geom.amplitude, phi=condenser_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)

    Corr_H = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_cond_BPHEX"}
    Corr_C = {"1P" : "water_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}

    # Set the pressure drop correlations of the condenser
    Condenser.set_DP()
    
    # Set the heat transfer coefficients correlations of the condenser           
    Condenser.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)
    
#%% EVAPORATOR PARAMETERS

    # Evaporator
    evaporator_geom = PlateGeomSWEP()
    evaporator_geom.set_parameters("P200THx140/1P_Evaporator")

    Evaporator.set_parameters(
        # Set the geometry of the evaporator
        A_c=evaporator_geom.A_c, A_h=evaporator_geom.A_h, h=evaporator_geom.h, l=evaporator_geom.l, l_v=evaporator_geom.l_v, w_v=evaporator_geom.w_v,
        C_CS=evaporator_geom.C_CS, C_Dh=evaporator_geom.C_Dh, C_V_tot=evaporator_geom.C_V_tot, C_canal_t=evaporator_geom.C_canal_t, C_n_canals=evaporator_geom.C_n_canals, 
        H_CS=evaporator_geom.H_CS, H_Dh=evaporator_geom.H_Dh, H_V_tot=evaporator_geom.H_V_tot, H_canal_t=evaporator_geom.H_canal_t, H_n_canals=evaporator_geom.H_n_canals,
        casing_t=evaporator_geom.casing_t, chevron_angle=evaporator_geom.chevron_angle, fooling=evaporator_geom.fooling, 
        n_plates=evaporator_geom.n_plates, plate_cond=evaporator_geom.plate_cond, plate_pitch_co=evaporator_geom.plate_pitch_co, t_plates=evaporator_geom.t_plates, w=evaporator_geom.w, 
        amplitude=evaporator_geom.amplitude, phi=evaporator_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)
    
    Corr_H = {"1P" : "water_plate_HTC", "2P" : "Han_cond_BPHEX"}
    Corr_C = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}

    # Set the pressure drop correlations of the condenser
    Evaporator.set_DP()
    
    # Set the heat transfer coefficients correlations of the condenser           
    Evaporator.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)

#%% EXPANDER PARAMETERS AND INPUTS

    N_exp = 6000
    T_amb = 293

    Expander_1.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)

    Expander_1.set_inputs(N_rot=N_exp, T_amb=T_amb)

    Expander_2.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)

    Expander_2.set_inputs(N_rot=N_exp, T_amb=T_amb)
    
    Expander_3.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)

    Expander_3.set_inputs(N_rot=N_exp, T_amb=T_amb)

#%% PUMP INPUTS

    Pump.set_inputs(N_pp = 2500)

#%% ADD AND LINK COMPONENTS
    ORC.add_component(Expander_1, "Expander_1")
    ORC.add_component(Expander_2, "Expander_2")
    ORC.add_component(Expander_3, "Expander_3")

    ORC.add_component(Mixer_expander, "Mixer")
    
    ORC.add_component(Condenser, "Condenser")
    ORC.add_component(Pump, "Pump") # :!\ Est-ce que les noms des composants sont importants?
    ORC.add_component(Evaporator, "Evaporator")

    ORC.add_component(Spliter_expander, "Spliter")

    # Link components+
    ORC.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
    ORC.link_components("Evaporator", "m-ex_C", "Spliter", "m-su")

    ORC.link_components("Spliter", "m-ex_1", "Expander_1", "m-su")
    ORC.link_components("Spliter", "m-ex_2", "Expander_2", "m-su")
    ORC.link_components("Spliter", "m-ex_3", "Expander_3", "m-su")

    ORC.link_components("Expander_1", "m-ex", "Mixer", "m-su_1")
    ORC.link_components("Expander_2", "m-ex", "Mixer", "m-su_2")
    ORC.link_components("Expander_3", "m-ex", "Mixer", "m-su_3")

    ORC.link_components("Mixer", "m-ex", "Condenser", "m-su_H")
    ORC.link_components("Condenser", "m-ex_H", "Pump", "m-su")

#%% SOURCES AND SINKS

    Hot_oil_source = MassConnector()
    ORC.add_source("EV_Water", Hot_oil_source, ORC.components["Evaporator"], "m-su_H")
    ORC.set_source_properties(T=90 + 273.15, fluid='Water', m_dot=3, target='EV_Water', P = 2e5)
    
    CD_water_source = MassConnector()
    ORC.add_source("CD_Water", CD_water_source, ORC.components["Condenser"], "m-su_C")
    ORC.set_source_properties(T=15 + 273.15, fluid='Water', m_dot=2.6, target='CD_Water', P = 2e5)

#%% CYCLE FIXED VARIABLES

    ORC.set_fixed_properties(target='Pump:su', SC = 5)

    ORC.set_fixed_properties(target='Spliter:su', P = 4*1e5, T = 340, m_dot = 0.5)
    ORC.set_fixed_properties(target='Expander_1:ex', P = 1.1*1e5)

    ORC.set_fixed_properties(target='Expander_2:ex', P = 1.1*1e5)
    
    ORC.set_fixed_properties(target='Expander_3:ex', P = 1.1*1e5)

#%% CYCLE RESIDUAL VARIABLES
    ORC.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
    ORC.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)

    ORC.solve()

    # ORC.print_states()
