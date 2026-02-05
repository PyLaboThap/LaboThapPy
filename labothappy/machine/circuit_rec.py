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

from machine.base_circuit import BaseCircuit

class RecursiveCircuit(BaseCircuit):

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
            
            self.res = (abs(new_value - self.value)/self.value)
            
            if self.res < self.tolerance:
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
        def __init__(self, target, variable, objective, objective_value, tol, rel, damping_factor, cycle):
            self.target = target 
            self.variable = variable
            self.objective = objective
            self.objective_value = objective_value
            self.tol = tol
            self.converged = False
            self.damping_factor = damping_factor
            self.cycle = cycle
            
            if rel == 1 or rel == -1:
                self.rel = rel
            else:
                raise ValueError("'rel' value for 'Iteration_variable' class shall be either 1 or -1.")
            
        def find_backward_path(self):
            return
            
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
        
        super().__init__()
        
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

    def set_source_properties(self, **kwargs):
        # Set properties for a specific source
        target = kwargs.pop('target')
        source = self.get_source(target)

        source.set_properties(**kwargs)

        for arg_name in kwargs:
            value = kwargs[arg_name]
            fix_prop = RecursiveCircuit.Fixed_Property(target, arg_name, value)
            
            # Update the fixed properties for the cycle
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop

        return 


#%% Variable set related methods

    def set_fixed_properties(self, **kwargs):

        # Set properties for a specific component and connector
        target = kwargs.pop('target')
                
        try:
            SC = kwargs.pop('SC')
            fix_prop = RecursiveCircuit.Fixed_Property(target, 'SC', SC)
            
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
            fix_prop = RecursiveCircuit.Fixed_Property(target, 'SH', SH)
            
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
            fix_prop = RecursiveCircuit.Fixed_Property(target, arg_name, value)
            
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

#%%

    def set_cycle_guess(self, external_set = False, **kwargs):
        # Set properties for a specific component and connector
        target = kwargs.pop('target')
                
        try:
            SC = kwargs.pop('SC')
            guess = RecursiveCircuit.Guess(target, 'SC', SC)
            
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
            guess = RecursiveCircuit.Guess(target, 'SH', SH)
            
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
            guess = RecursiveCircuit.Guess(target, arg_name, value)

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

    def set_iteration_variable(self, tol = 1e-2, target = None, variable = None, objective = None, rel = 1, damping_factor = 0.5, cycle = None):
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
            it_var = RecursiveCircuit.Iteration_variable(target, variable, objective, self.fixed_properties[objective].value, tol, rel, damping_factor, self)
            self.it_vars.append(it_var) 
        else:
            if "Link" in objective:
                it_var = RecursiveCircuit.Iteration_variable(target, variable, objective, None, tol, rel, damping_factor, self)
                self.it_vars.append(it_var) 
            else:
                raise ValueError(f"{objective} not part of the fixed_properties. Cannot be defined as an non-linking objective.")
        
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
        
    def set_residual_variable(self, target, variable, tolerance= 1e-3):
        res_var = RecursiveCircuit.Residual_variable(target, variable, tolerance)
        self.res_vars[res_var.name] = res_var
        return

    def print_res_vars(self):
        for res_var in self.res_vars:
            print(f"{res_var}: {self.res_vars[res_var].value}")
        
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
        if self.print_flag:
            print(f"----------------------------------")
            print(f"Component : {component_name}")
                
        component = self.get_component(component_name)
        component_model = component.model

        # !!!
        if component_name == 'Expander_1':
            print(f"P_su_exp : {component_model.su.p}")

        if component_model.solved:
            return

        component_model.check_calculable()
        component_model.check_parametrized()
        
        if component_model.parametrized:
            if component_model.calculable:        
                component_model.solve()
            else:
                return
        else:
            raise ValueError(f"Component '{component_name}' not parametrized.")
        
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
        
        for start_component in self.solve_start_components:
            self.recursive_solve(start_component)
            if self.check_all_component_solved():
                break
        
        if self.check_all_component_solved():
            if self.print_flag:
                print("All components solved in first iteration.")
            pass
        else:
            if self.print_flag:
                print("Not all components were solved in first iteration.")
            return
        
        # print("Set residual variables")
        
        for res_var_name in self.res_vars:
            component_res_var, rest = res_var_name.split(":")
            port_res_var, var_res_var = rest.split("-")

            connector = getattr(self.components[component_res_var].model, port_res_var)
            value = float(getattr(connector, var_res_var))
                                    
            self.res_vars[res_var_name].set_value(value)

        self.guess_update = True

        i=0
        n_it_max = 30   
        
        while i < n_it_max:
 
            for it_var in self.it_vars:
                
                if "Link" in it_var.objective:                    
                    _, obj_comp_name, port_var = it_var.objective.split(":")
                    port_name, variable = port_var.split("-")
                    
                    connector_obj = getattr(self.components[obj_comp_name].model, port_name)
                    value_obj = getattr(connector_obj, variable)
                    
                    for target in it_var.target:
                        kwargs_dict = {'target': target,
                                        it_var.variable : value_obj}
                        
                        self.set_cycle_guess(**kwargs_dict)
                        
                    it_var.converged = True

                else:
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
                    
                    
                        
            for fix_prop in self.fixed_properties:
                if ":" in fix_prop:
                    component_name, rest = guess.split(":")
                    port, var = rest.split("-")
                    
                    if var == "SC":
                        pass
                        
                    else:                        
                        kwargs_dict = {'target': self.fixed_properties[fix_prop].target,
                                        self.fixed_properties[fix_prop].variable : self.fixed_properties[fix_prop].value}
                        
                        self.set_fixed_properties(**kwargs_dict)

                else:
                    pass

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
            
            for it_var in self.it_vars:    
                if not it_var.converged:
                    self.converged = False            
        
            if not self.check_all_component_solved():
                self.converged = False
        
            if self.converged:
                if self.print_flag:
                    print(f"Solver successfully converged in {i+2} iteration(s) !")
                return
            
            i = i + 1
        
        if self.print_flag:
            print(f"Solver failed to converge in {n_it_max} iterations.")
        
        return
    