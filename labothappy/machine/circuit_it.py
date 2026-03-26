# -*- coding: utf-8 -*-
"""
Created on Thu Jan 08 2026
    
@author: elise neven
@email: elise.neven@uliege.be
"""


from machine.base_circuit import BaseCircuit

from scipy.optimize import fsolve, root
import numpy as np
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
        self.cycle_guesses = {} # Initial guesses — updated after each sweep

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

    def set_cycle_guess(self, **kwargs):
        """
        Set an initial value on a connector for DFS bootstrapping only.

        Unlike set_cycle_input, this value is NOT re-enforced at each
        iteration — it is applied once to help the DFS find a starting
        component, then left free for the sweep to overwrite.

        Unlike set_iteration_variable, it is NOT part of the solver's
        unknown vector — fsolve does not control it.

        Use this for connector properties that need an initial value to
        make the circuit calculable, but whose converged value will be
        determined naturally by the cycle propagation.

        Parameters
        ----------
        target : str
            "ComponentName:connector_name"
        **kwargs : float
            One or more variable=value pairs (e.g. SH=34, p=1e5).
        """
        target = kwargs.pop("target")
        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)

        # Store for reapplication after each _clear_all_connectors
        for variable, value in kwargs.items():
            name = f"{target}-{variable}"
            self.cycle_guesses[name] = {
                "target":    target,
                "component": component_name,
                "connector": connector_name,
                "variable":  variable,
                "value":     value
            }

        # Apply immediately for DFS bootstrapping
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

    def set_residual_variable(self, target, tolerance, closes=None, target_value=None, scale=None):
        """
        Define a residual variable.

        Three modes:

        Self-closure (closes=None, target_value=None):
            F = value_after_sweep(target) - value_before_sweep(target)
            Works when the target connector is recalculated by the sweep
            (e.g. an enthalpy that propagates through the cycle).

        Cross-closure (closes provided, target_value=None):
            F = value_after_sweep(target) - value_before_sweep(closes)
            Use when the two ends of the loop live on different connectors
            (e.g. pressure must match between Condenser:ex_H and ExpansionValve:su).

        Target mode (target_value provided):
            F = calculated_value(target) - target_value
            Drives a connector variable toward a fixed external value.
            If variable is "SC" or "SH", residual is computed on enthalpy
            internally via CoolProp.

        Parameters
        ----------
        target : str
            "ComponentName:connector_name-variable" — POST-sweep read location.
        tolerance : float
            Convergence tolerance (kept for reference).
        closes : str, optional
            "ComponentName:connector_name-variable" — PRE-sweep read location.
            Only used in cross-closure mode.
        target_value : float, optional
            External target value. If None, closure mode is used.
        scale : float, optional
            Normalisation factor. Defaults to _default_residual_scale(variable).
        """
        connector_part, variable = target.rsplit("-", 1)
        component_name, connector_name = connector_part.split(":")

        # Parse closes string if provided
        if closes is not None:
            closes_connector_part, closes_variable = closes.rsplit("-", 1)
            closes_component, closes_connector = closes_connector_part.split(":")
        else:
            closes_component  = component_name
            closes_connector  = connector_name
            closes_variable   = variable

        self.res_vars.append({
            "target":           target,
            "component":        component_name,
            "connector":        connector_name,
            "variable":         variable,
            "closes_component": closes_component,
            "closes_connector": closes_connector,
            "closes_variable":  closes_variable,
            "target_value":     float(target_value) if target_value is not None else None,
            "tol":              tolerance,
            "scale":            scale
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

                if hasattr(obj, "reset"):  # MassConnector
                    obj.reset()


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
        if self.print_flag:
            print(f"Solving {component_name}")
        
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
        if self.print_flag:
            print("Starting solving order definition")
        
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

        else:
            if self.print_flag:
                print("Setup completed !")
        
        # Reset solved markers — the DFS solve above was only for ordering
        self.reset_solved_marker()

    def _solve_circuit(self, x):
                        
        self.n_it += 1

        if self.print_flag:            
            print(f"="*30)
            print(f"Iteration : {self.n_it}")
            print(f"="*30)

        # 1️⃣ Reset all connectors
        self._clear_all_connectors()

        # 2️⃣ Apply iteration variables + inputs
        self._apply_iteration_vector(x)

        # 3️⃣ Snapshot PRE-solve values (used by closure residuals)
        pre_snapshot = {}
        for rv in self.res_vars:
            if rv["target_value"] is None:
                key = (rv["closes_component"], rv["closes_connector"], rv["closes_variable"])
                pre_snapshot[key] = self._read_variable(
                    f"{rv['closes_component']}:{rv['closes_connector']}", rv["closes_variable"]
                )
        self._pre_snapshot = pre_snapshot

        # 4️⃣ Solve components in DFS order
        if not self.solving_order:
            self._build_solve_order()

        self.reset_solved_marker()
        for name in self.solving_order: 

            if self.print_flag:
                print(f"Solving {name}")
                
            self.components[name].solve()

        # 5️⃣ Update cycle guesses from calculated connector values
        self._update_cycle_guesses()

        # 6️⃣ Compute residuals
        residuals = []

        print(self.components['Evaporator'].model.su_C.T)

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
            
            res = (calculated - target_val) / scale

            residuals.append(res)
            
        if self.print_flag:            
            print(f"residuals : {residuals}")
            
            
        return np.array(residuals, dtype=float)

    def _apply_cycle_guesses(self):
        """Re-apply stored guess values after _clear_all_connectors."""
        for g in self.cycle_guesses.values():
            component = self.get_component(g["component"])
            component.set_properties(g["connector"], **{g["variable"]: g["value"]})

    def _update_cycle_guesses(self):
        """
        After a sweep, read back the calculated values from the connectors
        and update the stored guesses — so the next iteration starts from
        the latest computed state rather than the original initial value.
        """
        for g in self.cycle_guesses.values():
            component = self.get_component(g["component"])
            connector = getattr(component.model, g["connector"])
            calculated = getattr(connector, g["variable"], None)
            if calculated is not None:
                g["value"] = calculated

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

        # 2️⃣ reapply cycle guesses (updated after each sweep)
        self._apply_cycle_guesses()

        # 3️⃣ reapply cycle inputs (fixed — overrides guesses if overlap)
        self._apply_cycle_inputs()

        # 4️⃣ apply iteration variables
        for value, it in zip(x, self.it_vars.values()):
            for entry in it["entries"]:
                entry["setter"](value)

    def _read_variable_pre(self, rv):
        """Return the pre-sweep snapshot value for a closure residual."""
        key = (rv["closes_component"], rv["closes_connector"], rv["closes_variable"])
        return self._pre_snapshot.get(key, 0.0)

    def _read_variable(self, target, variable):
        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)
        connector = getattr(component.model, connector_name)
        return getattr(connector, variable)
    
    # ── private solver backends ─────────────────────────────────────────────────
    
    def _solve_fsolve(self, x0):
        """Original fsolve backend — unchanged behaviour."""
        return fsolve(self._solve_circuit, x0)
    
    
    def _solve_root(self, x0, root_method='hybr'):
        """
        scipy.optimize.root backend.
    
        Recommended methods:
          'hybr'    — Powell's hybrid (same engine as fsolve, more diagnostics)
          'lm'      — Levenberg-Marquardt (robust to near-singular Jacobians)
          'broyden1'— quasi-Newton, good when the Jacobian is expensive
          'krylov'  — matrix-free, efficient for large systems
        """
        result = root(self._solve_circuit, x0, method=root_method,
                      options={'maxfev': 2000, 'xtol': 1e-8})
        if not result.success:
            print(f"[root/{root_method}] did not converge: {result.message}")
        return result.x
    
    
    def _solve_newton(self, x0, tol=1e-8, max_iter=100, eps=1e-6):
        """
        Baseline Newton-Raphson with finite-difference Jacobian.
    
        Useful as a reference / debugging baseline — the iteration is
        fully transparent and easy to instrument.
        """
        x = np.array(x0, dtype=float)
        n = len(x)
    
        for it in range(max_iter):
            self.n_it += 1
            F = np.array(self._solve_circuit(x), dtype=float)
    
            norm = np.linalg.norm(F)
            if norm < tol:
                print(f"[newton] converged in {it+1} iterations (|F|={norm:.2e})")
                return x
    
            # Finite-difference Jacobian  J[i,j] = dF_i/dx_j
            J = np.zeros((n, n))
            for j in range(n):
                x_fwd = x.copy()
                x_fwd[j] += eps
                F_fwd = np.array(self._solve_circuit(x_fwd), dtype=float)
                J[:, j] = (F_fwd - F) / eps
    
            # Newton step  J·dx = -F
            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                print("[newton] singular Jacobian — aborting")
                return x
    
            x = x + dx
    
        print(f"[newton] reached max_iter={max_iter} without converging "
              f"(|F|={np.linalg.norm(np.array(self._solve_circuit(x))):.2e})")
        return x

    def solve(self, method='fsolve'):
        """
        Solve the circuit.
        
        Parameters
        ----------
        method : str
            'fsolve'   — original scipy fsolve (default, backward-compatible)
            'root'     — scipy.optimize.root (more methods/diagnostics available)
            'newton'   — baseline Newton-Raphson with finite-difference Jacobian
        """
        root_methods = {'hybr', 'lm', 'broyden1', 'broyden2', 'anderson', 'krylov', 'df-sane'}
        
        self.n_it = 0
        if not self.solving_order:
            self._build_solve_order()
        self.reset_solved_marker()
        x0 = self._get_iteration_vector()
    
        # try:
        if method == 'fsolve':
            sol = self._solve_fsolve(x0)

        elif method == 'newton':
            sol = self._solve_newton(x0)

        elif method in root_methods:
            sol = self._solve_root(x0, root_method=method)

        else:
            raise ValueError(f"Unknown method '{method}'. "
                             f"Choose from: 'fsolve', 'newton', a root method: {root_methods}.")

        # except Exception as e:
        #     print(f"Error during solving ({method}): {e}")
        #     self.converged = False
        #     return
    
        self._apply_iteration_vector(sol)
        residuals = self._solve_circuit(sol)
        self.converged = np.all(np.abs(residuals) < 1.0)
        
        return sol
