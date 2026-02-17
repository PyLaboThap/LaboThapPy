# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: basile chaudoir
"""

from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import numpy as np

from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
from labothappy.connector.heat_connector import HeatConnector

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
            
            if self.value is not None and (abs(new_value - self.value)/self.value) < self.tolerance:
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
            
        def find_link_DP(self, components):
                        
            start = self.target[0].split(":")[0]
            end = self.objective.split(":")[1]
            
            start_comp = components[start]
            
            from collections import deque

            def get_suffix(port_name):
                if port_name is None:
                    return ""
                if port_name.endswith("_C"):
                    return "_C"
                if port_name.endswith("_H"):
                    return "_H"
                return ""
            
            def shortest_path(start_component, end_name):
                """
                Returns path as list of tuples:
                (component_name, entrance_port, exit_port_of_previous_component)
                The exit_port is the port on the previous component used to reach this component.
                """
                queue = deque()
                # Start: no entrance, no exit
                queue.append((start_component, [[start_component.name, None, None]]))
                
                visited = set([start_component])
            
                while queue:
                    current, path = queue.popleft()
            
                    if current.name == end_name:
                        return path
            
                    # Suffix filtering
                    _, entrance_port_current, _ = path[-1]
                    entrance_suffix = get_suffix(entrance_port_current)
            
                    for prev_exit_port, next_comp in current.next.items():
                        if next_comp in visited:
                            continue
            
                        # Entrance port of next component pointing back to current
                        entrance_port_next = next(
                            (port for port, prev in next_comp.previous.items() if prev is current),
                            None
                        )
            
                        # Only consider next components whose exit port matches suffix if needed
                        if entrance_suffix in ("_C", "_H") and not prev_exit_port.endswith(entrance_suffix):
                            continue
            
                        visited.add(next_comp)
            
                        # Store:
                        # - component name (next component)
                        # - entrance port on next component
                        # - exit port on previous component (current)
                        queue.append(
                            (next_comp, path + [[next_comp.name, entrance_port_next, prev_exit_port]])
                        )
            
                return None

            path = shortest_path(start_comp, end)
            
            for i in range(1,len(path)):
                path_elem = path[i]
                prev_path_elem = path[i-1]
                                
                prev_path_elem[-1] = path_elem[-1]
                
                if i == len(path)-1:
                    path_elem[-1] = None
            
            def compute_delta_P(path, components):
                """
                path: list of [component_name, entrance_port, exit_port_of_prev]
                components: dict of component_name -> component object
                
                Returns: dict of component_name -> delta_P
                         Only for components with both entrance and exit ports not None
                         Uses DP_c if port has '_C', DP_h if port has '_H'
                """
                delta_Ps = {}
            
                for comp_name, entrance_port, exit_port in path:
                    # Skip components with None ports
                    if entrance_port is None or exit_port is None:
                        continue
            
                    comp = components[comp_name].model  # use model to access DP_c / DP_h
            
                    # Determine suffix
                    suffix = ""
                    if entrance_port.endswith("_C") or exit_port.endswith("_C"):
                        suffix = "_C"
                    elif entrance_port.endswith("_H") or exit_port.endswith("_H"):
                        suffix = "_H"
            
                    # Assign delta_P based on suffix
                    if suffix == "_C":
                        delta_Ps[comp_name] = getattr(comp, "DP_c", 0)
                    elif suffix == "_H":
                        delta_Ps[comp_name] = getattr(comp, "DP_h", 0)
                    else:
                        # fallback if no suffix
                        delta_Ps[comp_name] = 0
            
                return delta_Ps
            
            self.delta_Ps = compute_delta_P(path, components)
            
            self.DP = sum(self.delta_Ps.values())
            
            return self.DP
            
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

            # self.set_residual_variable(target = target, variable=arg_name)

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

#%% Plot 

    import matplotlib.pyplot as plt
    from typing import List
    
    def plot_thermo_cycle(self, axes: str = 'Ts'):
        """
        Plot thermo cycle for each component in self.components.
    
        axes:
            'Ts' → Temperature–entropy diagram
            'Ph' → Pressure–enthalpy diagram
        """
        plt.figure()
    
        # --- Axis configuration ---
        if axes == 'Ts':
            x_attr, y_attr = 's', 'T'
            xlabel, ylabel = 's', 'T'
            title = 'Thermodynamic cycle (T–s)'
        elif axes == 'ph':
            x_attr, y_attr = 'h', 'p'
            xlabel, ylabel = 'h', 'p'
            title = 'Thermodynamic cycle (P–h)'
        else:
            raise ValueError(f"Unsupported axes mode: {axes}")
            
        self.groups_comp = {}
        
        for component in self.components:
            model = self.components[component].model
    
            prefix_supply = "su"
            prefix_exhaust = "ex"
    
            supply_connectors = []
            exhaust_connectors = []
    
            supply_suffixes = []
            exhaust_suffixes = []
    
            # 1) collect connector attribute names and suffixes
            for attr in vars(model):
                if attr.startswith(prefix_supply):
                    supply_connectors.append(attr)
                    supply_suffixes.append(attr.removeprefix(prefix_supply))
                elif attr.startswith(prefix_exhaust):
                    exhaust_connectors.append(attr)
                    exhaust_suffixes.append(attr.removeprefix(prefix_exhaust))
    
            # only keep suffixes that exist in both sets
            common_suffixes = sorted(set(supply_suffixes) & set(exhaust_suffixes))
    
            # Build groups: one group per common suffix
            groups = []
            for suffix in common_suffixes:
                groups.append({
                    "supply": [getattr(model, prefix_supply + suffix)],
                    "exhaust": [getattr(model, prefix_exhaust + suffix)],
                })
    
            # Collect unmatched connectors
            unmatched_supply = [
                s for s in supply_connectors
                if s.removeprefix(prefix_supply) not in common_suffixes
            ]
            unmatched_exhaust = [
                e for e in exhaust_connectors
                if e.removeprefix(prefix_exhaust) not in common_suffixes
            ]
    
            if unmatched_supply or unmatched_exhaust:
                groups.append({
                    "supply": [getattr(model, c) for c in unmatched_supply],
                    "exhaust": [getattr(model, c) for c in unmatched_exhaust]
                })
    
            for group in groups:                
                new_group = []
                
                for connector in group['supply']:
                    dic = {
                        x_attr : getattr(connector, x_attr),
                        y_attr : getattr(connector, y_attr)
                        }
                    
                    new_group.append(dic)

                # print(groups[group])
                group['supply'] = new_group
                
                new_group = []
                for connector in group['exhaust']:
                    dic = {
                        x_attr : getattr(connector, x_attr),
                        y_attr : getattr(connector, y_attr)
                        }
                    
                    new_group.append(dic)

                # print(groups[group])
                group['exhaust'] = new_group
    
            self.groups_comp[component] = groups
    
        for component, groups in self.groups_comp.items():
            for group in groups:
                supply_points = group['supply']
                exhaust_points = group['exhaust']
        
                # Determine how many segments we have
                n = max(len(supply_points), len(exhaust_points))
        
                # Helper: pad shorter list with last point
                def pad_points(points, target_len):
                    if not points:
                        return [{'x': None, 'y': None}] * target_len
                    if len(points) >= target_len:
                        return points[:target_len]
                    return points + [points[-1]] * (target_len - len(points))
        
                supply_points = pad_points(supply_points, n)
                exhaust_points = pad_points(exhaust_points, n)
        
                # Plot a line between each supply-exhaust pair
                for su, ex in zip(supply_points, exhaust_points):
                    xs = [su[x_attr], ex[x_attr]]
                    ys = [su[y_attr], ex[y_attr]]
        
                    if None not in xs and None not in ys:
                        plt.plot(xs, ys, 'b-o')  # you can change color/style
    
        # --- Plot saturation curve using CoolProp ---
        fluid = self.fluid
        
        if fluid:
            if axes == 'Ts':
                # T–s saturation curve
                P_min = PropsSI('Pcrit', fluid) * 0.1
                P_max = PropsSI('Pcrit', fluid)
                P_vals = np.logspace(np.log10(P_min), np.log10(P_max), 200)
                s_liq, s_vap, T_liq, T_vap = [], [], [], []
        
                for P in P_vals:
                    try:
                        T_liq.append(PropsSI('T', 'P', P, 'Q', 0, fluid))
                        T_vap.append(PropsSI('T', 'P', P, 'Q', 1, fluid))
                        s_liq.append(PropsSI('S', 'P', P, 'Q', 0, fluid))
                        s_vap.append(PropsSI('S', 'P', P, 'Q', 1, fluid))
                    except:
                        continue
        
                plt.plot(s_liq, T_liq, 'k--', label='Saturation liquid')
                plt.plot(s_vap, T_vap, 'k--', label='Saturation vapor')
        
            elif axes == 'Ph':
                # P–h saturation curve
                T_min = PropsSI('Ttriple', fluid)
                T_max = PropsSI('Tcrit', fluid)
                T_vals = np.linspace(T_min, T_max, 200)
                P_liq, P_vap, h_liq, h_vap = [], [], [], []
        
                for T in T_vals:
                    try:
                        P_liq.append(PropsSI('P', 'T', T, 'Q', 0, fluid))
                        P_vap.append(PropsSI('P', 'T', T, 'Q', 1, fluid))
                        h_liq.append(PropsSI('H', 'T', T, 'Q', 0, fluid))
                        h_vap.append(PropsSI('H', 'T', T, 'Q', 1, fluid))
                    except:
                        continue
        
                plt.plot(h_liq, P_liq, 'k--', label='Saturation liquid')
                plt.plot(h_vap, P_vap, 'k--', label='Saturation vapor')

        plt.show()

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

        if component_model.solved:
            return

        component_model.check_calculable()
        component_model.check_parametrized()
        
        if component_model.parametrized:
            if component_model.calculable:    
                save = {}
                for next_connector in component.next:
                    type_connector, connector_name = next_connector.split("-")
                    if type_connector == "m":
                        connector = getattr(component_model, connector_name)
                        save[next_connector] = {'p' : connector.p, 'h' : connector.h}
                    
                component_model.solve()
                
                save_new = {}
                for next_connector in component.next:
                    type_connector, connector_name = next_connector.split("-")
                    if type_connector == "m":
                        connector = getattr(component_model, connector_name)
                        save_new[next_connector] = {'p' : connector.p, 'h' : connector.h}
                
                tol = 1e-2
                
                for connector in save:
                    for prop in save[connector]:
                        value_prev = save[connector][prop]
                        value_new = save_new[connector][prop]
                        
                        if value_prev is not None:
                            delta = abs((value_new - value_prev)/value_prev)
                            # print(f"delta : {delta} for {connector} - {prop}")
                            if delta > tol:
                                component.next[connector].model.solved = False
                            
                        
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
                print("All components solved in setup iteration.")
            pass
        else:
            # for start_component in self.solve_start_components:
            #     self.recursive_solve(start_component)
            #     if self.check_all_component_solved():
            #         break
            
            if self.print_flag:
                self.blocking_comp = []
                for component in self.components:     
                    if not self.components[component].model.check_calculable():
                        self.blocking_comp.append(component)
                
                print("Not all components were solved in setup iteration.")
                print(f"{self.blocking_comp} is/are not calculable.")
            return
        
        # print("Set residual variables")
        
        for res_var_name in self.res_vars:
            component_res_var, rest = res_var_name.split(":")
            port_res_var, var_res_var = rest.split("-")

            connector = getattr(self.components[component_res_var].model, port_res_var)
            value = float(getattr(connector, var_res_var))
                                    
            self.res_vars[res_var_name].set_value(value)

        self.guess_update = True
        
        self.convergence_frames.append(self.plot_cycle_Ts(plot_auto = False))

        i=0
        n_it_max = 30   
        
        while i < n_it_max:
            if self.print_flag:
                print(f"Iteration {i+1}")
            
            for it_var in self.it_vars:
                
                if "Link" in it_var.objective:        
                    
                    _, obj_comp_name, port_var = it_var.objective.split(":")
                    port_name, variable = port_var.split("-")

                    if variable == 'p':
                        gain = it_var.find_link_DP(self.components)
                    else:
                        gain = 0
                        
                    connector_obj = getattr(self.components[obj_comp_name].model, port_name)
                    value_obj = getattr(connector_obj, variable)
                    
                    for target in it_var.target:
                        kwargs_dict = {'target': target,
                                        it_var.variable : value_obj+gain}
                                                
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
                    print(f"Solver successfully converged in {i+1} iteration(s) !")
                    
                return
            
            self.convergence_frames.append(self.plot_cycle_Ts(plot_auto=False))
            
            i = i + 1

        # plt.figure()
        # plt.plot(res_ev, 'r',  marker='o')                            
        # plt.show()
        
        if self.print_flag:
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
    
    
