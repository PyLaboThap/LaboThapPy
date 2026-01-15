# -*- coding: utf-8 -*-
"""
Created on Thu Jan 08 2026
    
@author: elise neven
@email: elise.neven@uliege.be
"""


from machine.base_circuit import BaseCircuit

import numpy as np
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI

class IterativeCircuit(BaseCircuit):
    def __init__(self, fluid=None):
        
        super().__init__()

        # Properties and inputs
        self.fluid = fluid
        self.inputs = {}

        # Variables vectors/dictionary
        self.res_vars = [] # List of residual variables
        self.it_vars = {} # Dictionary of iteration variables

        # Solve related variables
        self.solve_start_components = [] # Ordered list of components to solve
        self.converged = False # Convergence flag
        self.state_cache = { 
            "pre": {}, # State cache before solving
            "post": {} # State cache after solving
        }

    def set_source_properties(self, **kwargs):
        # Set properties for a specific source
        target = kwargs.pop('target')
        source = self.get_source(target)

        source.set_properties(**kwargs)

    def set_cycle_input(self, **kwargs):
        target = kwargs.pop("target")

        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)

        # 1. Store inputs
        for var, value in kwargs.items():
            input_name = f"{target}-{var.upper()}"
            self.inputs[input_name] = {
                "target": target,
                "component": component_name,
                "connector": connector_name,
                "variable": var,
                "value": value
            }

        # 2. Apply values to component model
        component.set_properties(connector_name, **kwargs)

    def set_iteration_variable(self, target, variable, guess=None, tolerance=1e-6):

        if isinstance(target, str):
            targets = [target]
        else:
            targets = list(target)

        var = variable
        name = f"{'+'.join(targets)}-{var}"

        entries = []

        for tgt in targets:
            component_name, connector_name = tgt.split(":")
            component = self.get_component(component_name)
            connector = getattr(component.model, connector_name)

            if not hasattr(connector, var):
                raise AttributeError(f"{tgt} has no attribute '{var}'")

            setter_name = f"set_{var}"
            if not hasattr(connector, setter_name):
                raise AttributeError(f"{tgt} has no setter '{setter_name}()'")

            entries.append({
                "target": tgt,
                "component": component_name,
                "connector": connector_name,
                "variable": var,
                "connector_obj": connector,
                "setter": getattr(connector, setter_name)
            })

        if guess is None:
            raise ValueError(
                f"Initial guess for iteration variable '{name}' is None"
            )

        self.it_vars[name] = {
            "name": name,
            "variable": var,
            "entries": entries,
            "x0": float(guess),
            "tolerance": tolerance
        }

    def set_residual_variable(self, pre_target, post_target, variable, tolerance, scale=None):
        self.res_vars.append({
            "pre_target": pre_target,
            "post_target": post_target,
            "variable": variable,
            "tol": tolerance,
            "scale": scale
        })

    def _default_residual_scale(self, variable):
        default_scales = {
            "p": 1e5,     # Pa
            "T": 1.0,     # K
            "h": 1e5,     # J/kg
            "s": 1e3,     # J/kg/K
            "m": 1e-2,    # kg/s
            "q": 1e3,     # W
            "w": 1e3      # W
        }

        return default_scales.get(variable, 1.0)


    def _get_iteration_vector(self):
        """
        Returns initial solver vector x0
        """
        return np.array([it["x0"] for it in self.it_vars.values()])
    
    def _solve_circuit(self, x):

        # 1. Apply iteration variables
        self._apply_iteration_vector(x)

        # 2. Snapshot PRE-solve values
        self.state_cache["pre"].clear()
        for rv in self.res_vars:
            key = (rv["pre_target"], rv["variable"])
            self.state_cache["pre"][key] = self._read_variable(
                rv["pre_target"], rv["variable"]
            )

        # 3. Solve components (one full circuit pass)
        if not self.solve_start_components:
            self._build_solve_order()

        for name in self.solve_start_components:
            self.components[name].solve()

        # 4. Snapshot POST-solve values
        self.state_cache["post"].clear()
        for rv in self.res_vars:
            key = (rv["post_target"], rv["variable"])
            self.state_cache["post"][key] = self._read_variable(
                rv["post_target"], rv["variable"]
            )

        # 5. Compute residuals
        residuals = []

        for rv in self.res_vars:
            pre_key = (rv["pre_target"], rv["variable"])
            post_key = (rv["post_target"], rv["variable"])

            raw = (
                self.state_cache["post"][post_key]
                - self.state_cache["pre"][pre_key]
            )

            scale = (
                rv["scale"]
                if rv["scale"] is not None
                else self._default_residual_scale(rv["variable"])
            )

            residuals.append(raw / scale)
        return np.array(residuals, dtype=float)


    
    def _apply_iteration_vector(self, x):
        """
        Apply solver vector x to all iteration variables
        """
        for value, it in zip(x, self.it_vars.values()):
            for entry in it["entries"]:
                entry["setter"](value)

    def _read_variable(self, target, variable):
        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)
        connector = getattr(component.model, connector_name)
        return getattr(connector, variable)
    
    def _build_solve_order(self):
        names = list(self.components.keys())

        start_name = self._find_start_component()
        start_index = names.index(start_name)

        # Circular reordering
        ordered_names = names[start_index:] + names[:start_index]

        self.solve_start_components = ordered_names

    def _find_start_component(self):
        for name, comp in self.components.items():
            if comp.model.check_calculable():
                return name
        raise RuntimeError("No calculable component found to start the circuit solve.")

    def solve(self):

        # Initial guess vector
        x0 = self._get_iteration_vector()

        # Solve
        sol = fsolve(self._solve_circuit, x0)

        # Apply final solution
        self._apply_iteration_vector(sol)

        # Convergence check
        residuals = self._solve_circuit(sol)
        self.converged = np.all(np.abs(residuals) < 1.0)


        return sol

        