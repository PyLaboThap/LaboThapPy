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
            self.prev_value = None
            self.tolerance = tolerance
            self.converged = False

        def check_convergence(self, new_value):
            
            if self.value is not None and (abs(new_value - self.value)/self.value) < self.tolerance:
                self.set_prev_value(self.value)
                self.set_value(new_value)
                self.converged = True
                return True
            else: 
                self.set_prev_value(self.value)
                self.set_value(new_value)
                self.converged = False
                return False
            
            return
            
        def set_value(self, value):
            self.value = value
            return

        def set_prev_value(self, value):
            self.prev_value = value
            return
        
    class Guess():
        def __init__(self, target, variable, value):
    
            self.name = target + "-" + variable
            self.target = target 
            self.variable = variable
            self.value = value

            # History tracking: list of (input_value, output_value) pairs
            # input_history[i]  = x^i  (the guess fed into iteration i)
            # output_history[i] = f(x^i) (the connector value that came out)
            self.input_history  = []   # x^0, x^1, x^2, ...
            self.output_history = []   # f(x^0), f(x^1), f(x^2), ...

        # ------------------------------------------------------------------
        # History helpers
        # ------------------------------------------------------------------

        def record_input(self, value):
            """Call this BEFORE solving, to record the guess that was fed in."""
            self.input_history.append(float(value))

        def record_output(self, value):
            """Call this AFTER solving, to record the connector value that came out."""
            self.output_history.append(float(value))

        @property
        def prev_input(self):
            """Last recorded input (x^{n-1}), or None if fewer than 2 entries."""
            if len(self.input_history) >= 2:
                return self.input_history[-2]
            return None

        @property
        def prev_output(self):
            """Last recorded output (f(x^{n-1})), or None if fewer than 2 entries."""
            if len(self.output_history) >= 2:
                return self.output_history[-2]
            return None

        @property
        def last_input(self):
            """Most recent recorded input (x^n), or None."""
            if self.input_history:
                return self.input_history[-1]
            return None

        @property
        def last_output(self):
            """Most recent recorded output (f(x^n)), or None."""
            if self.output_history:
                return self.output_history[-1]
            return None

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
            
                    comp = components[comp_name].model # use model to access DP_c / DP_h
            
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
            if abs((target_value - self.objective_value) / self.objective_value) < self.tol:
                self.converged = True
                return None
        
            self.converged = False
        
            # Build a temporary connector at the desired state to read off the iteration variable
            temp = MassConnector(fluid=objective_connector.fluid)
            temp.set_properties(**{var: self.objective_value, 'T': objective_connector.T})
            return getattr(temp, self.variable)  
        
#%%

    def __init__(self, fluid=None):
        
        super().__init__()
        
        # Building Blocks
        self.components = {}
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
        self.solving_order = []
        
        # Convergence tolerance for automatic connector state tracking
        self.convergence_tolerance = 1e-4
        
        # Internal storage for previous connector states {comp_name: {port_name: {p, h}}}
        self._prev_connector_states = {}

        # Wegstein acceleration settings (only used when solve(method='Wegstein') is called)
        # SC and SH guesses always use successive substitution regardless of method.
        self.wegstein_q_min = -5.0   # clipping bounds for the acceleration factor q
        self.wegstein_q_max =  0.0

    # ------------------------------------------------------------------
    # NEW: snapshot and check helpers for global p/h convergence
    # ------------------------------------------------------------------

    def _snapshot_connector_states(self):
        """
        Record current p and h for every MassConnector output port of every
        component.  Only ports whose values are not None are stored.
        Returns a nested dict {comp_name: {port_name: {'p': ..., 'h': ...}}}.
        """
        snapshot = {}
        for comp_name, comp in self.components.items():
            snapshot[comp_name] = {}
            for port_key in comp.next:
                type_connector, port_name = port_key.split("-", 1)
                if type_connector != "m":
                    continue
                connector = getattr(comp.model, port_name, None)
                if connector is None:
                    continue
                p_val = getattr(connector, 'p', None)
                h_val = getattr(connector, 'h', None)
                if p_val is not None and h_val is not None:
                    snapshot[comp_name][port_name] = {'p': float(p_val), 'h': float(h_val)}
        return snapshot

    def _check_connector_convergence(self, current_snapshot):
        """
        Compare current connector states against the previous snapshot.
        Returns (all_converged: bool, messages: list[str]).
        A relative tolerance self.convergence_tolerance is applied to both p and h.
        Guard against near-zero denominators with a small epsilon.
        """
        epsilon = 1e-10
        all_converged = True
        messages = []

        for comp_name, ports in current_snapshot.items():
            prev_ports = self._prev_connector_states.get(comp_name, {})
            for port_name, vals in ports.items():
                prev_vals = prev_ports.get(port_name)
                if prev_vals is None:
                    # No previous value yet — cannot declare convergence on first iteration
                    all_converged = False
                    continue
                for prop in ('p', 'h'):
                    v_new = vals[prop]
                    v_old = prev_vals[prop]
                    rel_err = abs(v_new - v_old) / max(abs(v_old), epsilon)
                    if rel_err > self.convergence_tolerance:
                        all_converged = False
                        messages.append(
                            f"Connector not converged: {comp_name}:{port_name}-{prop} "
                            f"| rel_err={rel_err:.2e} (tol={self.convergence_tolerance:.2e})"
                        )

        return all_converged, messages

    # ------------------------------------------------------------------

    @staticmethod
    def _wegstein_update(x_prev, x_curr, f_prev, f_curr, q_min=-5.0, q_max=0.0):
        """
        Compute the next iterate using Wegstein's acceleration method.

        Parameters
        ----------
        x_prev : float  — guess fed in at iteration n-1  (x^{n-1})
        x_curr : float  — guess fed in at iteration n    (x^n)
        f_prev : float  — circuit output at iteration n-1 (f(x^{n-1}))
        f_curr : float  — circuit output at iteration n   (f(x^n))
        q_min  : float  — lower clip bound for acceleration factor q (default -5)
        q_max  : float  — upper clip bound for acceleration factor q (default  0)

        Returns
        -------
        float — accelerated next guess x^{n+1}

        Notes
        -----
        When history is insufficient (None inputs) the method falls back to
        plain successive substitution, i.e. returns f_curr unchanged.
        The acceleration factor q is clipped to [q_min, q_max] for stability.
        q = 0  →  pure successive substitution.
        q < 0  →  extrapolation (acceleration).
        """
        # Fall back to successive substitution when history is not yet available
        if x_prev is None or f_prev is None:
            return f_curr

        denominator = x_curr - x_prev
        if abs(denominator) < 1e-12:
            return f_curr  # avoid division by zero

        q = (f_curr - f_prev) / denominator
        q = max(q_min, min(q_max, q))  # clip for stability

        if abs(q - 1.0) < 1e-12:
            return f_curr  # avoid division by zero in the Wegstein formula

        return (q * x_curr - f_curr) / (q - 1.0)

    # ------------------------------------------------------------------

    def set_source_properties(self, **kwargs):
        target = kwargs.pop('target')
        source = self.get_source(target)

        source.set_properties(**kwargs)

        for arg_name in kwargs:
            value = kwargs[arg_name]
            fix_prop = RecursiveCircuit.Fixed_Property(target, arg_name, value)
            
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop

        return 

#%% Variable set related methods

    def set_fixed_properties(self, **kwargs):
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
    
        for arg_name, value in kwargs.items():
            fix_prop = RecursiveCircuit.Fixed_Property(target, arg_name, value)
            if fix_prop.name not in self.fixed_properties:
                self.fixed_properties[fix_prop.name] = fix_prop
    
        component.set_properties(connector_name, **kwargs)
        component.model.check_calculable()

        return

#%%

    def set_cycle_guess(self, external_set=False, **kwargs):
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
    
        for arg_name, value in kwargs.items():
            guess = RecursiveCircuit.Guess(target, arg_name, value)
    
            if self.guess_update or guess.name not in self.guesses or external_set:
                if guess.name in self.guesses:
                    # Preserve history when updating an existing guess
                    guess.input_history  = self.guesses[guess.name].input_history
                    guess.output_history = self.guesses[guess.name].output_history
                self.guesses[guess.name] = guess
    
        component.set_properties(connector_name, **kwargs)
        component.model.check_calculable()

        return

    def set_iteration_variable(self, tol = 1e-2, target = None, variable = None, objective = None, rel = 1, damping_factor = 0.5, cycle = None):
        
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

    def print_guess_history(self):
        """Print the full input/output history for every guess variable."""
        for name, guess in self.guesses.items():
            print(f"\n--- {name} ---")
            print(f"  {'Iter':<6} {'Input (x)':<20} {'Output f(x)':<20}")
            n = max(len(guess.input_history), len(guess.output_history))
            for i in range(n):
                x   = guess.input_history[i]  if i < len(guess.input_history)  else "-"
                fx  = guess.output_history[i] if i < len(guess.output_history) else "-"
                x_s  = f"{x:.6g}"  if isinstance(x,  float) else x
                fx_s = f"{fx:.6g}" if isinstance(fx, float) else fx
                print(f"  {i:<6} {x_s:<20} {fx_s:<20}")

    def print_fixed_properties(self):
        for fix_prop in self.fixed_properties:
            print(f"{fix_prop}: {self.fixed_properties[fix_prop].value}")
        
        return

    def reset_input_values(self):
        
        for comp in self.components:
            comp_model = self.components[comp].model
            comp_model.reset_inputs()

        if self.print_flag:
            print(f"----------------------------------")
            print(f"Reset component inputs")

#%%

    def set_cycle_parameters(self, **kwargs):
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

        if component_model.solved:
            if self.print_flag:
                print(f"Component '{component_name}' already solved.")
            return

        component_model.check_calculable()
        component_model.check_parametrized()
        
        if component_model.parametrized:
            if component_model.calculable:    
                component_model.solve()
                
                if component_name not in self.solving_order:
                    self.solving_order.append(component_name)
                        
            else:
                if self.print_flag:
                    print(f"Component '{component_name}' not calculable.")
                return
        else:
            raise ValueError(f"Component '{component_name}' not parametrized.")
        
        for output_port in component.next:
            next_comp_name = component.next[output_port].name
            self.recursive_solve(next_comp_name)
        
        return

    def solve(self, max_iter=30, method='successive_substitution', tol = 1e-6):
        """
        Solve the circuit.

        Parameters
        ----------
        max_iter : int, optional
            Maximum number of iterations. Default is 30.
        method : str, optional
            Convergence method for tear stream guess updates. Options:
            - 'successive_substitution' (default) : plain fixed-point update
            - 'Wegstein'                           : Wegstein acceleration

            Note: SC and SH guesses always use successive substitution
            regardless of the chosen method, as they involve nonlinear
            pressure transformations.
        """
        _VALID_METHODS = ('successive_substitution', 'Wegstein')
        if method not in _VALID_METHODS:
            raise ValueError(
                f"Unknown method '{method}'. "
                f"Valid options are: {_VALID_METHODS}"
            )

        use_wegstein = (method == 'Wegstein')

        if self.print_flag:
            print(f"Solve method: {method}")

        # Find first calculable components and determine solving order
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
        else:
            if self.print_flag:
                self.blocking_comp = []
                for component in self.components:     
                    if not self.components[component].model.check_calculable():
                        self.blocking_comp.append(component)
                print("Not all components were solved in setup iteration.")
                print(f"{self.blocking_comp} is/are not calculable.")
            return
        
        # Legacy res_vars initialisation (kept for backward compatibility)
        for res_var_name in self.res_vars:
            component_res_var, rest = res_var_name.split(":")
            port_res_var, var_res_var = rest.split("-")
            connector = getattr(self.components[component_res_var].model, port_res_var)
            value = float(getattr(connector, var_res_var))
            self.res_vars[res_var_name].set_value(value)

        # Take initial snapshot of all connector states before iteration starts
        self._prev_connector_states = self._snapshot_connector_states()

        self.guess_update = True
        
        i = 0
        
        while i < max_iter:
            self.messages = []
            
            if self.print_flag:
                print("\n")
                print(f"###########################")
                print(f"Iteration {i+1}")
                print(f"###########################")

            # --- Iteration variables update ---
            for it_var in self.it_vars:
                
                if "Link" in it_var.objective:        
                    _, obj_comp_name, port_var = it_var.objective.split(":")
                    port_name, variable = port_var.split("-")
            
                    gain = it_var.find_link_DP(self.components) if variable == 'p' else 0
                        
                    connector_obj = getattr(self.components[obj_comp_name].model, port_name)
                    value_obj = getattr(connector_obj, variable)
                    
                    for target in it_var.target:
                        kwargs_dict = {'target': target, it_var.variable: value_obj + gain}
                        self.set_cycle_guess(**kwargs_dict)
                        
                    it_var.converged = True
            
                else:
                    obj_comp, rest = it_var.objective.split(":")
                    port, var = rest.split("-")
            
                    connector = getattr(self.components[obj_comp].model, port)
                    current_value = getattr(connector, var)  # works for SC, SH, p, T, h, ...
            
                    if current_value is not None:
                        desired_value = it_var.check(current_value, connector, var)
            
                        if desired_value is not None:
                            target_comp, target_port = it_var.target[0].split(":")
                            target_connector = getattr(self.components[target_comp].model, target_port)
                            target_current_value = getattr(target_connector, it_var.variable)
            
                            new_target_value = (desired_value - target_current_value) * it_var.damping_factor + target_current_value
            
                            for target in it_var.target:
                                kwargs_dict = {'target': target, it_var.variable: new_target_value}
                                self.set_cycle_guess(**kwargs_dict)

            # --- Fixed properties re-enforcement ---
            for fix_prop in self.fixed_properties:
                if ":" in fix_prop:
                    kwargs_dict = {
                        'target': self.fixed_properties[fix_prop].target,
                        self.fixed_properties[fix_prop].variable: self.fixed_properties[fix_prop].value
                    }
                    self.set_fixed_properties(**kwargs_dict)
            
            # --- Guess update from current connector states ---
            for guess_name, guess in self.guesses.items():
                component_name, rest = guess_name.split(":")
                port, var = rest.split("-")
            
                connector = getattr(self.components[component_name].model, port)
                f_curr = float(getattr(connector, var))
                x_curr = guess.value
            
                # Record aligned (input, output) pair
                guess.record_input(x_curr)
                guess.record_output(f_curr)
            
                if use_wegstein and var not in ("SC", "SH"):
                    next_value = RecursiveCircuit._wegstein_update(
                        x_prev = guess.prev_input,
                        x_curr = x_curr,
                        f_prev = guess.prev_output,
                        f_curr = f_curr,
                        q_min  = self.wegstein_q_min,
                        q_max  = self.wegstein_q_max,
                    )
                    if self.print_flag:
                        print(f"  Wegstein [{guess_name}]: x_curr={x_curr:.4g}, "
                              f"f_curr={f_curr:.4g} → next={next_value:.4g}")
                else:
                    next_value = f_curr
            
                kwargs_dict = {'target': guess.target, guess.variable: next_value}
                self.set_cycle_guess(**kwargs_dict)
            
            # --- Component solve sweep ---
            self.reset_solved_marker()
        
            for component_name in self.solving_order:
                if self.print_flag:
                    print(f"----------------------------------")
                    print(f"Component : {component_name}")
                    
                comp_model = self.components[component_name]
                comp_model.solve()
        
            # --- Convergence check ---
            self.converged = True

            # 1. Global p/h convergence across all connector output states
            current_snapshot = self._snapshot_connector_states()
            connector_converged, connector_messages = self._check_connector_convergence(current_snapshot)
            
            if not connector_converged:
                self.converged = False
                self.messages.extend(connector_messages)
            
            # Update snapshot for next iteration
            self._prev_connector_states = current_snapshot

            # 2. Iteration variable convergence (unchanged)
            for it_var in self.it_vars:    
                if not it_var.converged:
                    self.converged = False        
                    self.messages.append(f"Iteration variable tolerance not satisfied : {it_var.target}-{it_var.variable}.")

            # 3. All components solved
            if not self.check_all_component_solved():
                self.converged = False
                self.messages.append("Not all components solved.")
                    
            if self.converged:
                if self.print_flag:
                    print(f"Solver successfully converged in {i+1} iteration(s) !")
                return
            
            if self.plot_flag:
                self.convergence_frames.append(self.plot_cycle_Ts(plot_auto=False))
                        
            i += 1
        
        if self.print_flag:
            print(f"Solver failed to converge in {max_iter} iterations.")
        
        return