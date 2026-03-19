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
        self.source_inputs = {}

        # Variables vectors/dictionary
        self.res_vars = [] # List of residual variables
        self.it_vars = {} # Dictionary of iteration variables

        # Solve related variables
        self.solve_start_components = [] # Ordered list of components to solve
        self.solving_order = []          # DFS-determined solving order
        self.converged = False # Convergence flag

    def set_source_properties(self, **kwargs):

        target = kwargs.pop("target")

        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)

        connector = getattr(component.model, connector_name)

        # Store inputs
        for var, value in kwargs.items():

            name = f"{target}-{var.upper()}"

            self.source_inputs[name] = {
                "target": target,
                "component": component_name,
                "connector": connector_name,
                "variable": var,
                "value": value
            }

        # Apply immediately
        connector.set_properties(**kwargs)

    def _apply_source_inputs(self):

        for inp in self.source_inputs.values():

            component = self.get_component(inp["component"])

            connector = getattr(component.model, inp["connector"])

            variable = inp["variable"]

            value = inp["value"]

            connector.set_properties(**{variable: value})

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

    def set_residual_variable(self, target, tolerance, target_value=None, scale=None):
        """
        Define a residual variable.

        Two modes depending on whether target_value is provided:

        Closure mode (target_value=None):
            F = value_after_sweep - value_before_sweep
            The residual closes the cycle: the calculated value must match
            the value that was imposed as a guess before the sweep.

        Target mode (target_value provided):
            F = calculated_value - target_value
            The residual drives a connector variable toward a fixed external value.
            If variable is "SC" or "SH", the residual is computed on enthalpy
            internally via CoolProp.

        Parameters
        ----------
        target : str
            "ComponentName:connector_name-variable"  e.g. "Condenser:ex_H-p"
        tolerance : float
            Convergence tolerance (kept for reference).
        target_value : float, optional
            External target value. If None, closure mode is used.
        scale : float, optional
            Normalisation factor. Defaults to _default_residual_scale(variable).
        """
        connector_part, variable = target.rsplit("-", 1)
        component_name, connector_name = connector_part.split(":")

        self.res_vars.append({
            "target": target,
            "component": component_name,
            "connector": connector_name,
            "variable": variable,
            "target_value": float(target_value) if target_value is not None else None,
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
    
    def _clear_all_connectors(self):
        """
        Reset all connectors in the circuit before a new iteration.
        Keeps the fluid but removes thermodynamic states.
        """

        for comp in self.components.values():
            model = comp.model

            for attr in dir(model):
                obj = getattr(model, attr)

                if hasattr(obj, "clear_state"):  # MassConnector
                    obj.clear_state()


    def _get_iteration_vector(self):
        """
        Returns initial solver vector x0
        """
        return np.array([it["x0"] for it in self.it_vars.values()])

    # -----------------------------------------------------------------------
    # DFS solve order  (ported from circuit_mix)
    # -----------------------------------------------------------------------

    def check_all_component_solved(self) -> bool:
        return all(c.model.solved for c in self.components.values())

    def reset_solved_marker(self):
        for c in self.components.values():
            c.model.solved = False

    def recursive_solve(self, component_name: str):
        component       = self.get_component(component_name)
        component_model = component.model
        if component_model.solved:
            return
        component_model.check_calculable()
        component_model.check_parametrized()
        if not component_model.parametrized:
            raise ValueError(f"Component '{component_name}' not parametrized.")
        if component_model.calculable:
            component_model.solve()
            if component_name not in self.solving_order:
                self.solving_order.append(component_name)
        else:
            return
        for output_port in component.next:
            self.recursive_solve(component.next[output_port].name)

    def _build_solve_order(self):
        """
        Determine solving order via DFS traversal (ported from circuit_mix).
        Populates self.solving_order and self.solve_start_components.
        """
        # Apply initial guesses so components can evaluate check_calculable()
        x0 = self._get_iteration_vector()
        self._apply_iteration_vector(x0)

        for component_name in self.components:
            if self.components[component_name].model.check_calculable():
                self.solve_start_components.append(component_name)

        for start in self.solve_start_components:
            self.recursive_solve(start)
            if self.check_all_component_solved():
                break

        if not self.check_all_component_solved():
            blocking = [c for c in self.components
                        if not self.components[c].model.check_calculable()]
            raise RuntimeError(
                f"Setup failed. Blocking components (not calculable): {blocking}"
            )

        # Reset solved markers — the DFS solve above was only for ordering
        self.reset_solved_marker()

    def _solve_circuit(self, x):

        self.n_it += 1

        # 1️⃣ Reset all connectors
        self._clear_all_connectors()

        # 2️⃣ Apply iteration variables + inputs
        self._apply_iteration_vector(x)

        # 3️⃣ Snapshot PRE-solve values (used by closure residuals)
        pre_snapshot = {}
        for rv in self.res_vars:
            if rv["target_value"] is None:
                key = (rv["component"], rv["connector"], rv["variable"])
                pre_snapshot[key] = self._read_variable(
                    f"{rv['component']}:{rv['connector']}", rv["variable"]
                )
        self._pre_snapshot = pre_snapshot

        # 4️⃣ Solve components in DFS order
        if not self.solving_order:
            self._build_solve_order()

        self.reset_solved_marker()
        for name in self.solving_order:
            self.components[name].solve()

        # 5️⃣ Compute residuals
        residuals = []

        for rv in self.res_vars:
            var      = rv["variable"]
            comp_key = f"{rv['component']}:{rv['connector']}"

            # --- Determine target ---
            if rv["target_value"] is None:
                # Closure mode: target is the value imposed before the sweep
                target_val = self._read_variable_pre(rv)
            elif var in ("SC", "SH"):
                # SC/SH target: convert to enthalpy via CoolProp
                p = self._read_variable(comp_key, "p")
                T_sat = PropsSI('T', 'P', p, 'Q', 0.5, self.fluid)
                sign  = -1 if var == "SC" else 1
                T_target = T_sat + sign * rv["target_value"]
                target_val = PropsSI('H', 'P', p, 'T', T_target, self.fluid)
                var = "h"   # residual is now on h
            else:
                target_val = rv["target_value"]

            calculated = self._read_variable(comp_key, var)

            scale = (
                rv["scale"]
                if rv["scale"] is not None
                else self._default_residual_scale(var)
            )

            residuals.append((calculated - target_val) / scale)

            
        return np.array(residuals, dtype=float)

    def _apply_cycle_inputs(self):
        for inp in self.inputs.values():

            component = self.get_component(inp["component"])

            connector = inp["connector"]

            variable = inp["variable"]

            value = inp["value"]

            component.set_properties(connector, **{variable: value})
        
    def _apply_iteration_vector(self, x):

        # 1️⃣ reapply source inputs
        self._apply_source_inputs()

        # 2️⃣ reapply cycle inputs
        self._apply_cycle_inputs()

        # 3️⃣ apply iteration variables
        for value, it in zip(x, self.it_vars.values()):
            for entry in it["entries"]:
                entry["setter"](value)

    def _read_variable_pre(self, rv):
        """Return the pre-sweep snapshot value for a closure residual."""
        key = (rv["component"], rv["connector"], rv["variable"])
        return self._pre_snapshot.get(key, 0.0)

    def _read_variable(self, target, variable):
        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)
        connector = getattr(component.model, connector_name)
        return getattr(connector, variable)

    def solve(self):
        self.n_it = 0

        # Build DFS solving order before handing off to fsolve
        if not self.solving_order:
            self._build_solve_order()

        # Reset solved markers so _solve_circuit can re-solve cleanly
        self.reset_solved_marker()

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