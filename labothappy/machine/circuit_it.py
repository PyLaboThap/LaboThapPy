# # -*- coding: utf-8 -*-
# """
# Created on Thu Jan 08 2026
    
# @author: elise neven
# @email: elise.neven@uliege.be
# """


# from machine.base_circuit import BaseCircuit

# from scipy.optimize import fsolve, root
# import numpy as np
# from CoolProp.CoolProp import PropsSI

# class IterativeCircuit(BaseCircuit):
#     def __init__(self, fluid=None):
        
#         super().__init__()

#         # Properties and inputs
#         self.fluid = fluid
#         self.inputs = {}
#         self.source_inputs = {}

#         # Variables vectors/dictionary
#         self.res_vars = [] # List of residual variables
#         self.it_vars = {} # Dictionary of iteration variables

#         # Solve related variables
#         self.solve_start_components = [] # Ordered list of components to solve
#         self.solving_order = []          # DFS-determined solving order
#         self.converged = False # Convergence flag
#         self.cycle_guesses = {} # Initial guesses — updated after each sweep

#     def set_source_properties(self, **kwargs):

#         target = kwargs.pop("target")

#         component_name, connector_name = target.split(":")
#         component = self.get_component(component_name)

#         connector = getattr(component.model, connector_name)

#         # Store inputs
#         for var, value in kwargs.items():

#             name = f"{target}-{var.upper()}"

#             self.source_inputs[name] = {
#                 "target": target,
#                 "component": component_name,
#                 "connector": connector_name,
#                 "variable": var,
#                 "value": value
#             }

#     #     # Apply immediately
#     #     connector.set_properties(**kwargs)

#     def _apply_source_inputs(self):

#         for inp in self.source_inputs.values():

#             component = self.get_component(inp["component"])

#             connector = getattr(component.model, inp["connector"])

#             variable = inp["variable"]

#             value = inp["value"]

#             connector.set_properties(**{variable: value})

#     def set_cycle_input(self, **kwargs):
#         target = kwargs.pop("target")

#         component_name, connector_name = target.split(":")
#         component = self.get_component(component_name)

#         # 1. Store inputs
#         for var, value in kwargs.items():
#             input_name = f"{target}-{var.upper()}"
#             self.inputs[input_name] = {
#                 "target": target,
#                 "component": component_name,
#                 "connector": connector_name,
#                 "variable": var,
#                 "value": value
#             }

#         # 2. Apply values to component model
#         component.set_properties(connector_name, **kwargs)

#     def set_cycle_guess(self, **kwargs):
#         """
#         Set an initial value on a connector for DFS bootstrapping only.

#         Unlike set_cycle_input, this value is NOT re-enforced at each
#         iteration — it is applied once to help the DFS find a starting
#         component, then left free for the sweep to overwrite.

#         Unlike set_iteration_variable, it is NOT part of the solver's
#         unknown vector — fsolve does not control it.

#         Use this for connector properties that need an initial value to
#         make the circuit calculable, but whose converged value will be
#         determined naturally by the cycle propagation.

#         Parameters
#         ----------
#         target : str
#             "ComponentName:connector_name"
#         **kwargs : float
#             One or more variable=value pairs (e.g. SH=34, p=1e5).
#         """
#         target = kwargs.pop("target")
#         component_name, connector_name = target.split(":")
#         component = self.get_component(component_name)

#         # Store for reapplication after each _clear_all_connectors
#         for variable, value in kwargs.items():
#             name = f"{target}-{variable}"
#             self.cycle_guesses[name] = {
#                 "target":    target,
#                 "component": component_name,
#                 "connector": connector_name,
#                 "variable":  variable,
#                 "value":     value
#             }

#         # Apply immediately for DFS bootstrapping
#         component.set_properties(connector_name, **kwargs)

#     def set_iteration_variable(self, target, variable, guess=None, tolerance=1e-6):

#         if isinstance(target, str):
#             targets = [target]
#         else:
#             targets = list(target)

#         var = variable
#         name = f"{'+'.join(targets)}-{var}"

#         entries = []

#         for tgt in targets:
#             component_name, connector_name = tgt.split(":")
#             component = self.get_component(component_name)
#             connector = getattr(component.model, connector_name)

#             if not hasattr(connector, var):
#                 raise AttributeError(f"{tgt} has no attribute '{var}'")

#             setter_name = f"set_{var}"
#             if not hasattr(connector, setter_name):
#                 raise AttributeError(f"{tgt} has no setter '{setter_name}()'")

#             entries.append({
#                 "target": tgt,
#                 "component": component_name,
#                 "connector": connector_name,
#                 "variable": var,
#                 "connector_obj": connector,
#                 "setter": getattr(connector, setter_name)
#             })

#         if guess is None:
#             raise ValueError(
#                 f"Initial guess for iteration variable '{name}' is None"
#             )

#         self.it_vars[name] = {
#             "name": name,
#             "variable": var,
#             "entries": entries,
#             "x0": float(guess),
#             "tolerance": tolerance
#         }

#     def set_residual_variable(self, target, tolerance, closes=None, target_value=None, scale=None):
#         """
#         Define a residual variable.

#         Three modes:

#         Self-closure (closes=None, target_value=None):
#             F = value_after_sweep(target) - value_before_sweep(target)
#             Works when the target connector is recalculated by the sweep
#             (e.g. an enthalpy that propagates through the cycle).

#         Cross-closure (closes provided, target_value=None):
#             F = value_after_sweep(target) - value_before_sweep(closes)
#             Use when the two ends of the loop live on different connectors
#             (e.g. pressure must match between Condenser:ex_H and ExpansionValve:su).

#         Target mode (target_value provided):
#             F = calculated_value(target) - target_value
#             Drives a connector variable toward a fixed external value.
#             If variable is "SC" or "SH", residual is computed on enthalpy
#             internally via CoolProp.

#         Parameters
#         ----------
#         target : str
#             "ComponentName:connector_name-variable" — POST-sweep read location.
#         tolerance : float
#             Convergence tolerance (kept for reference).
#         closes : str, optional
#             "ComponentName:connector_name-variable" — PRE-sweep read location.
#             Only used in cross-closure mode.
#         target_value : float, optional
#             External target value. If None, closure mode is used.
#         scale : float, optional
#             Normalisation factor. Defaults to _default_residual_scale(variable).
#         """
#         connector_part, variable = target.rsplit("-", 1)
#         component_name, connector_name = connector_part.split(":")

#         # Parse closes string if provided
#         if closes is not None:
#             closes_connector_part, closes_variable = closes.rsplit("-", 1)
#             closes_component, closes_connector = closes_connector_part.split(":")
#         else:
#             closes_component  = component_name
#             closes_connector  = connector_name
#             closes_variable   = variable

#         self.res_vars.append({
#             "target":           target,
#             "component":        component_name,
#             "connector":        connector_name,
#             "variable":         variable,
#             "closes_component": closes_component,
#             "closes_connector": closes_connector,
#             "closes_variable":  closes_variable,
#             "target_value":     float(target_value) if target_value is not None else None,
#             "tol":              tolerance,
#             "scale":            scale
#         })

#     def _default_residual_scale(self, variable):
#         default_scales = {
#             "p": 1e5,     # Pa
#             "T": 1.0,     # K
#             "h": 1e5,     # J/kg
#             "s": 1e3,     # J/kg/K
#             "m": 1e-2,    # kg/s
#             "q": 1e3,     # W
#             "w": 1e3      # W
#         }

#         return default_scales.get(variable, 1.0)
    
#     def _clear_all_connectors(self):
#         """
#         Reset all connectors in the circuit before a new iteration.
#         Keeps the fluid but removes thermodynamic states.
#         """

#         for comp in self.components.values():
#             model = comp.model

#             for attr in dir(model):
#                 obj = getattr(model, attr)

#                 if hasattr(obj, "reset"):  # MassConnector
#                     obj.reset()


#     def _get_iteration_vector(self):
#         """
#         Returns initial solver vector x0
#         """
#         return np.array([it["x0"] for it in self.it_vars.values()])

#     # -----------------------------------------------------------------------
#     # DFS solve order  (ported from circuit_mix)
#     # -----------------------------------------------------------------------

#     def check_all_component_solved(self) -> bool:
#         return all(c.model.solved for c in self.components.values())

#     def reset_solved_marker(self):
#         for c in self.components.values():
#             c.model.solved = False

#     def recursive_solve(self, component_name: str):
#         if self.print_flag:
#             print(f"Solving {component_name}")
        
#         component       = self.get_component(component_name)
#         component_model = component.model
#         if component_model.solved:
#             return
#         component_model.check_calculable()
#         component_model.check_parametrized()
#         if not component_model.parametrized:
#             raise ValueError(f"Component '{component_name}' not parametrized.")
#         if component_model.calculable:
#             component_model.solve()
#             if component_name not in self.solving_order:
#                 self.solving_order.append(component_name)
#         else:
#             return
#         for output_port in component.next:
#             self.recursive_solve(component.next[output_port].name)

#     def _build_solve_order(self):
#         """
#         Determine solving order via DFS traversal (ported from circuit_mix).
#         Populates self.solving_order and self.solve_start_components.
#         """
#         if self.print_flag:
#             print("Starting solving order definition")
        
#         # Apply initial guesses so components can evaluate check_calculable()
#         x0 = self._get_iteration_vector()
#         self._apply_iteration_vector(x0)

#         for component_name in self.components:
#             if self.components[component_name].model.check_calculable():
#                 self.solve_start_components.append(component_name)

#         for start in self.solve_start_components:
#             self.recursive_solve(start)
#             if self.check_all_component_solved():
#                 break

#         if not self.check_all_component_solved():
#             blocking = [c for c in self.components
#                         if not self.components[c].model.check_calculable()]
#             raise RuntimeError(
#                 f"Setup failed. Blocking components (not calculable): {blocking}"
#             )

#         else:
#             if self.print_flag:
#                 print("Setup completed !")
        
#         # Reset solved markers — the DFS solve above was only for ordering
#         self.reset_solved_marker()

#     def _solve_circuit(self, x):
                        
#         self.n_it += 1

#         self.x_log.append(x)

#         if self.print_flag:            
#             print(f"="*30)
#             print(f"Iteration : {self.n_it}")
#             print(f"="*30)

#         # 1️⃣ Reset all connectors
#         self._clear_all_connectors()

#         # 2️⃣ Apply iteration variables + inputs
#         self._apply_iteration_vector(x)

#         # 3️⃣ Snapshot PRE-solve values (used by closure residuals)
#         pre_snapshot = {}
#         for rv in self.res_vars:
#             if rv["target_value"] is None:
#                 key = (rv["closes_component"], rv["closes_connector"], rv["closes_variable"])
#                 pre_snapshot[key] = self._read_variable(
#                     f"{rv['closes_component']}:{rv['closes_connector']}", rv["closes_variable"]
#                 )
#         self._pre_snapshot = pre_snapshot

#         # 4️⃣ Solve components in DFS order
#         if not self.solving_order:
#             self._build_solve_order()

#         self.reset_solved_marker()
#         for name in self.solving_order: 

#             if self.print_flag:
#                 print(f"Solving {name}")
                
#             self.components[name].solve()

#         # 5️⃣ Update cycle guesses from calculated connector values
#         self._update_cycle_guesses()

#         # 6️⃣ Compute residuals
#         residuals = []

#         for rv in self.res_vars:
#             var      = rv["variable"]
#             comp_key = f"{rv['component']}:{rv['connector']}"

#             # --- Determine target ---
#             if rv["target_value"] is None:
#                 # Closure mode: target is the value imposed before the sweep
#                 target_val = self._read_variable_pre(rv)
                                
#             elif var in ("SC", "SH"):
#                 # SC/SH target: convert to enthalpy via CoolProp
#                 p = self._read_variable(comp_key, "p")
#                 T_sat = PropsSI('T', 'P', p, 'Q', 0.5, self.fluid)
#                 sign  = -1 if var == "SC" else 1
#                 T_target = T_sat + sign * rv["target_value"]
#                 target_val = PropsSI('H', 'P', p, 'T', T_target, self.fluid)
#                 var = "h"   # residual is now on h
#             else:
#                 target_val = rv["target_value"]

#             calculated = self._read_variable(comp_key, var)

#             scale = (
#                 rv["scale"]
#                 if rv["scale"] is not None
#                 else self._default_residual_scale(var)
#             )
            
#             res = (calculated - target_val) / scale

#             residuals.append(float(res))
            
#         if self.print_flag:            
#             print(f"residuals : {residuals}")
            
#         # self.print_states()
        
#         return np.array(residuals, dtype=float)

#     def _apply_cycle_guesses(self):
#         """Re-apply stored guess values after _clear_all_connectors."""
#         for g in self.cycle_guesses.values():
#             component = self.get_component(g["component"])
#             component.set_properties(g["connector"], **{g["variable"]: g["value"]})

#     def _update_cycle_guesses(self):
#         """
#         After a sweep, read back the calculated values from the connectors
#         and update the stored guesses — so the next iteration starts from
#         the latest computed state rather than the original initial value.
#         """
#         for g in self.cycle_guesses.values():
#             component = self.get_component(g["component"])
#             connector = getattr(component.model, g["connector"])
#             calculated = getattr(connector, g["variable"], None)
#             if calculated is not None:
#                 g["value"] = calculated

#     def _apply_cycle_inputs(self):
#         for inp in self.inputs.values():

#             component = self.get_component(inp["component"])

#             connector = inp["connector"]

#             variable = inp["variable"]

#             value = inp["value"]

#             component.set_properties(connector, **{variable: value})
        
#     def _apply_iteration_vector(self, x):

#         # 1️⃣ reapply source inputs
#         self._apply_source_inputs()

#         # 2️⃣ reapply cycle guesses (updated after each sweep)
#         self._apply_cycle_guesses()

#         # 3️⃣ reapply cycle inputs (fixed — overrides guesses if overlap)
#         self._apply_cycle_inputs()

#         # 4️⃣ apply iteration variables
#         for value, it in zip(x, self.it_vars.values()):
#             for entry in it["entries"]:
#                 entry["setter"](value)

#     def _read_variable_pre(self, rv):
#         """Return the pre-sweep snapshot value for a closure residual."""
#         key = (rv["closes_component"], rv["closes_connector"], rv["closes_variable"])
#         return self._pre_snapshot.get(key, 0.0)

#     def _read_variable(self, target, variable):
#         component_name, connector_name = target.split(":")
#         component = self.get_component(component_name)
#         connector = getattr(component.model, connector_name)
#         return getattr(connector, variable)
    
#     # ── private solver backends ─────────────────────────────────────────────────
    
#     def _solve_fsolve(self, x0):
#         """Original fsolve backend — unchanged behaviour."""
#         return fsolve(self._solve_circuit, x0)
    
    
#     def _solve_root(self, x0, root_method='hybr'):
#         """
#         scipy.optimize.root backend.
    
#         Recommended methods:
#           'hybr'    — Powell's hybrid (same engine as fsolve, more diagnostics)
#           'lm'      — Levenberg-Marquardt (robust to near-singular Jacobians)
#           'broyden1'— quasi-Newton, good when the Jacobian is expensive
#           'krylov'  — matrix-free, efficient for large systems
#         """
#         result = root(self._solve_circuit, x0, method=root_method,
#                       options={'maxfev': 2000, 'xtol': 1e-5})
#         if not result.success:
#             print(f"[root/{root_method}] did not converge: {result.message}")
#         return result.x
    
    
#     def _solve_newton(self, x0, tol=1e-5, max_iter=100, eps=1e-6):
#         """
#         Baseline Newton-Raphson with finite-difference Jacobian.
    
#         Useful as a reference / debugging baseline — the iteration is
#         fully transparent and easy to instrument.
#         """
#         x = np.array(x0, dtype=float)
#         n = len(x)
        
        
#         for it in range(max_iter):
#             self.n_it += 1
#             F = np.array(self._solve_circuit(x), dtype=float)
#             # print(f"Evap P_sat : {self.components['Evaporator'].model.ex_C.p}")

#             norm = np.linalg.norm(F)
#             if norm < tol:
#                 if self.print_flag:
#                     print(f"[newton] converged in {it+1} iterations (|F|={norm:.2e})")
#                 return x
    
#             # Finite-difference Jacobian  J[i,j] = dF_i/dx_j
#             J = np.zeros((n, n))
#             for j in range(n):
#                 x_fwd = x.copy()
#                 x_fwd[j] += eps
#                 F_fwd = np.array(self._solve_circuit(x_fwd), dtype=float)
#                 J[:, j] = (F_fwd - F) / eps
#                 # print(f"Evap P_sat : {self.components['Evaporator'].model.su_C.p}")
                
#             # Newton step  J·dx = -F
#             try:
#                 dx = np.linalg.solve(J, -F)
#                 # print(f"Evap P_sat : {self.components['Evaporator'].model.su_C.p}")
#             except np.linalg.LinAlgError:
#                 print("[newton] singular Jacobian — aborting")
#                 self.converged = False
#                 return x
    
#             x = x + dx
    
#             # print(f"x : {x}")
        
        
#         print(f"[newton] reached max_iter={max_iter} without converging "
#               f"(|F|={np.linalg.norm(np.array(self._solve_circuit(x))):.2e})")
#         return x

#     def solve(self, method='fsolve'):
#         """
#         Solve the circuit.
        
#         Parameters
#         ----------
#         method : str
#             'fsolve'   — original scipy fsolve (default, backward-compatible)
#             'root'     — scipy.optimize.root (more methods/diagnostics available)
#             'newton'   — baseline Newton-Raphson with finite-difference Jacobian
#         """
#         root_methods = {'hybr', 'lm', 'broyden1', 'broyden2', 'anderson', 'krylov', 'df-sane'}
        
#         self.n_it = 0
#         if not self.solving_order:
#             self._build_solve_order()
#         self.reset_solved_marker()
#         x0 = self._get_iteration_vector()
        
#         self.x_log = [x0]
        
#         try:
#             if method == 'fsolve':
#                 sol = self._solve_fsolve(x0)
    
#             elif method == 'newton':
#                 sol = self._solve_newton(x0)
    
#             elif method in root_methods:
#                 sol = self._solve_root(x0, root_method=method)
    
#             else:
#                 raise ValueError(f"Unknown method '{method}'. "
#                                  f"Choose from: 'fsolve', 'newton', a root method: {root_methods}.")

#         except Exception as e:
#             print(f"Error during solving ({method}): {e}")
#             self.converged = False
#             return
    
#         self._apply_iteration_vector(sol)
#         residuals = self._solve_circuit(sol)
#         self.converged = np.all(np.abs(residuals) < 1.0)
        
#         return sol

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 08 2026
Revised 2026-11-12:
    - residual function is now deterministic (cycle guesses frozen during a solver run)
    - real convergence check (per-residual tolerances + solver exit flag)
    - controllable finite-difference step (fd_step) for noisy components
    - residuals returned/printed as plain floats
    - debugging tools (states, energy balance, residual tables, determinism check, scans)

Changes are tagged [CHANGED] or [NEW] in the code.

@author: elise neven
@email: elise.neven@uliege.be
"""

import copy
from contextlib import contextmanager

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
        self.res_vars = []  # List of residual variables
        self.it_vars = {}   # Dictionary of iteration variables

        # Solve related variables
        self.solve_start_components = []  # Ordered list of components to solve
        self.solving_order = []           # DFS-determined solving order
        self.converged = False            # Convergence flag
        self.cycle_guesses = {}           # Initial guesses (frozen during a solver run)

        # [NEW] diagnostics
        self.history = []          # one entry per residual evaluation: {"it", "x", "res"}
        self.last_residuals = []   # detailed info on the last residual evaluation
        self.guess_drift = {}      # applied vs recalculated value of each cycle guess
        self.solver_message = ""   # exit message of the last solver run
        self.fd_step = None        # relative finite-difference step (None = scipy default)
        self.n_it = 0

    # =======================================================================
    # Inputs, guesses, iteration and residual variables
    # =======================================================================

    def set_source_properties(self, **kwargs):

        target = kwargs.pop("target")

        component_name, connector_name = target.split(":")

        for var, value in kwargs.items():
            name = f"{target}-{var.upper()}"
            self.source_inputs[name] = {
                "target": target,
                "component": component_name,
                "connector": connector_name,
                "variable": var,
                "value": value
            }

    def _apply_source_inputs(self):
        for inp in self.source_inputs.values():
            component = self.get_component(inp["component"])
            connector = getattr(component.model, inp["connector"])
            connector.set_properties(**{inp["variable"]: inp["value"]})

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

        [CHANGED] The guess is now held FIXED while the solver runs, so that the
        residual function F(x) depends only on x. If the guess sits on a
        connector that is read before the sweep recalculates it, it acts as a
        hidden input: print_cycle_guesses() will flag it. Turn such a guess
        into an iteration variable + closure residual, or use
        solve(update_guesses=True) to update it in an outer loop.

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

        for variable, value in kwargs.items():
            name = f"{target}-{variable}"
            self.cycle_guesses[name] = {
                "target":    target,
                "component": component_name,
                "connector": connector_name,
                "variable":  variable,
                "value":     value,
                "initial":   value,   # [NEW] kept to allow reset_cycle_guesses()
            }

        # Apply immediately for DFS bootstrapping
        component.set_properties(connector_name, **kwargs)

    def reset_cycle_guesses(self):
        """[NEW] Restore all cycle guesses to the values given by the user."""
        for g in self.cycle_guesses.values():
            g["value"] = g["initial"]

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

        Cross-closure (closes provided, target_value=None):
            F = value_after_sweep(target) - value_before_sweep(closes)

        Target mode (target_value provided):
            F = calculated_value(target) - target_value
            If variable is "SC" or "SH", residual is computed on enthalpy
            internally via CoolProp.

        Parameters
        ----------
        target : str
            "ComponentName:connector_name-variable" — POST-sweep read location.
        tolerance : float
            [CHANGED] Now USED for the convergence check. It is expressed in the
            physical unit of the compared variable (Pa for p, J/kg for h, K for T,
            kg/s for m_dot...). For SC/SH targets the comparison is made on
            enthalpy, so the tolerance is in J/kg.
        closes : str, optional
            "ComponentName:connector_name-variable" — PRE-sweep read location.
        target_value : float, optional
            External target value. If None, closure mode is used.
        scale : float, optional
            Normalisation factor. Defaults to _default_residual_scale(variable).
        """
        connector_part, variable = target.rsplit("-", 1)
        component_name, connector_name = connector_part.split(":")

        if closes is not None:
            closes_connector_part, closes_variable = closes.rsplit("-", 1)
            closes_component, closes_connector = closes_connector_part.split(":")
        else:
            closes_component = component_name
            closes_connector = connector_name
            closes_variable = variable

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
            "m_dot": 1e-2,
            "q": 1e3,     # W
            "w": 1e3      # W
        }
        return default_scales.get(variable, 1.0)

    # =======================================================================
    # Connector handling
    # =======================================================================

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
        """Returns initial solver vector x0"""
        return np.array([it["x0"] for it in self.it_vars.values()], dtype=float)

    def _apply_cycle_guesses(self):
        """Re-apply stored guess values after _clear_all_connectors."""
        for g in self.cycle_guesses.values():
            component = self.get_component(g["component"])
            component.set_properties(g["connector"], **{g["variable"]: g["value"]})

    def _update_cycle_guesses(self):
        """
        [CHANGED] Only called BETWEEN solver runs (outer loop), never inside
        _solve_circuit. Reads back the calculated values and stores them as
        the new guesses. Returns the largest relative change.
        """
        max_change = 0.0
        for g in self.cycle_guesses.values():
            component = self.get_component(g["component"])
            connector = getattr(component.model, g["connector"])
            calculated = getattr(connector, g["variable"], None)
            if calculated is not None:
                old = g["value"]
                max_change = max(max_change, abs(calculated - old) / max(abs(old), 1e-12))
                g["value"] = float(calculated)
        return max_change

    def _record_guess_drift(self):
        """[NEW] Store applied vs recalculated value of each cycle guess."""
        self.guess_drift = {}
        for name, g in self.cycle_guesses.items():
            connector = getattr(self.get_component(g["component"]).model, g["connector"])
            self.guess_drift[name] = (g["value"], getattr(connector, g["variable"], None))

    def _apply_cycle_inputs(self):
        for inp in self.inputs.values():
            component = self.get_component(inp["component"])
            component.set_properties(inp["connector"], **{inp["variable"]: inp["value"]})

    def _apply_iteration_vector(self, x):

        # 1. reapply source inputs
        self._apply_source_inputs()

        # 2. reapply cycle guesses (frozen during a solver run)
        self._apply_cycle_guesses()

        # 3. reapply cycle inputs (fixed — overrides guesses if overlap)
        self._apply_cycle_inputs()

        # 4. apply iteration variables
        for value, it in zip(x, self.it_vars.values()):
            for entry in it["entries"]:
                entry["setter"](value)

    def _read_variable_pre(self, rv):
        """Return the pre-sweep snapshot value for a closure residual."""
        key = (rv["closes_component"], rv["closes_connector"], rv["closes_variable"])
        # [CHANGED] no silent fallback to 0.0
        if key not in self._pre_snapshot:
            raise KeyError(f"No pre-sweep snapshot for {key}")
        return self._pre_snapshot[key]

    def _read_variable(self, target, variable):
        component_name, connector_name = target.split(":")
        component = self.get_component(component_name)
        connector = getattr(component.model, connector_name)
        value = getattr(connector, variable)
        # [CHANGED] explicit error instead of a cryptic TypeError later
        if value is None:
            raise ValueError(f"{target}-{variable} is not set (None)")
        return value

    # =======================================================================
    # DFS solve order
    # =======================================================================

    def check_all_component_solved(self) -> bool:
        return all(c.model.solved for c in self.components.values())

    def reset_solved_marker(self):
        for c in self.components.values():
            c.model.solved = False

    def recursive_solve(self, component_name: str):
        if self.print_flag:
            print(f"Solving {component_name}")

        component = self.get_component(component_name)
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
        Determine solving order via DFS traversal.
        Populates self.solving_order and self.solve_start_components.
        """
        if self.print_flag:
            print("Starting solving order definition")

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
        elif self.print_flag:
            print("Setup completed !")

        self.reset_solved_marker()

    # =======================================================================
    # Residual function
    # =======================================================================

    def _solve_circuit(self, x):
        """
        One sweep through the cycle. Returns the scaled residual vector.

        [CHANGED] Pure function of x: cycle guesses are NOT updated here any
        more, so two calls with the same x give the same residuals (as long as
        the components themselves have no memory).
        """
        x = np.asarray(x, dtype=float)
        self.n_it += 1

        if self.print_flag:
            print("=" * 30)
            print(f"Iteration : {self.n_it}")
            print("=" * 30)

        # 1. Reset all connectors
        self._clear_all_connectors()

        # 2. Apply inputs, guesses and iteration variables
        self._apply_iteration_vector(x)

        # 3. Snapshot PRE-solve values (used by closure residuals)
        pre_snapshot = {}
        for rv in self.res_vars:
            if rv["target_value"] is None:
                key = (rv["closes_component"], rv["closes_connector"], rv["closes_variable"])
                pre_snapshot[key] = self._read_variable(
                    f"{rv['closes_component']}:{rv['closes_connector']}", rv["closes_variable"]
                )
        self._pre_snapshot = pre_snapshot

        # 4. Solve components in DFS order
        if not self.solving_order:
            self._build_solve_order()

        self.reset_solved_marker()
        for name in self.solving_order:
            if self.print_flag:
                print(f"Solving {name}")
            self.components[name].solve()

        # 5. [CHANGED] record (but do not apply) how far each guess drifted
        self._record_guess_drift()

        # 6. Compute residuals
        residuals, details = self._compute_residuals()
        self.last_residuals = details
        self.history.append({"it": self.n_it, "x": x.copy(), "res": residuals.copy()})

        if self.print_flag:
            # [CHANGED] .tolist() -> plain floats, no "np.float64(...)"
            print(f"residuals : {residuals.tolist()}")

        return residuals

    def _compute_residuals(self):
        """[NEW] Returns (scaled residual array, list of detail dicts)."""
        residuals = []
        details = []

        for rv in self.res_vars:
            var = rv["variable"]
            comp_key = f"{rv['component']}:{rv['connector']}"

            if rv["target_value"] is None:
                target_val = self._read_variable_pre(rv)
            elif var in ("SC", "SH"):
                p = self._read_variable(comp_key, "p")
                T_sat = PropsSI('T', 'P', p, 'Q', 0.5, self.fluid)
                sign = -1 if var == "SC" else 1
                T_target = T_sat + sign * rv["target_value"]
                target_val = PropsSI('H', 'P', p, 'T', T_target, self.fluid)
                var = "h"
            else:
                target_val = rv["target_value"]

            calculated = self._read_variable(comp_key, var)

            scale = rv["scale"] if rv["scale"] is not None else self._default_residual_scale(var)

            diff = float(calculated - target_val)
            res = diff / scale
            residuals.append(res)

            details.append({
                "target": rv["target"],
                "compared": var,
                "calculated": float(calculated),
                "reference": float(target_val),
                "diff": diff,
                "scaled": res,
                "tol": rv["tol"],
                "ok": abs(diff) <= rv["tol"],
            })

        return np.array(residuals, dtype=float), details

    # =======================================================================
    # Solver backends — each returns (x, success, message)
    # =======================================================================

    def _epsfcn(self):
        """fsolve/hybr/lm use step ≈ sqrt(epsfcn)·|x|, so epsfcn = fd_step²."""
        return None if self.fd_step is None else self.fd_step ** 2

    def _solve_fsolve(self, x0):
        # [CHANGED] full_output to get the real exit flag
        sol, info, ier, msg = fsolve(self._solve_circuit, x0,
                                     full_output=True, epsfcn=self._epsfcn())
        return sol, ier == 1, msg

    def _solve_root(self, x0, root_method='hybr'):
        """
        scipy.optimize.root backend.
          'hybr'    — Powell's hybrid (same engine as fsolve)
          'lm'      — Levenberg-Marquardt (robust to near-singular Jacobians)
          'broyden1'— quasi-Newton
          'krylov'  — matrix-free
        """
        # [CHANGED] options adapted to each method (avoids "unknown option" warnings)
        if root_method == 'hybr':
            options = {'maxfev': 2000, 'xtol': 1e-5}
        elif root_method == 'lm':
            options = {'maxiter': 2000, 'xtol': 1e-5}
        else:
            options = {'maxiter': 200}
        if root_method in ('hybr', 'lm') and self.fd_step is not None:
            options['eps'] = self._epsfcn()

        result = root(self._solve_circuit, x0, method=root_method, options=options)
        return result.x, bool(result.success), str(result.message)

    def _solve_newton(self, x0, tol=1e-6, max_iter=50):
        """
        Newton-Raphson with a finite-difference Jacobian.
        [CHANGED] relative FD step (iteration variables can be ~1e5 Pa) and no
        double counting of n_it.
        """
        step = self.fd_step if self.fd_step is not None else 1e-6
        x = np.array(x0, dtype=float)
        n = len(x)

        for it in range(max_iter):
            F = self._solve_circuit(x)
            norm = np.linalg.norm(F)
            if norm < tol:
                return x, True, f"converged in {it + 1} Newton steps (|F|={norm:.2e})"

            J = np.zeros((n, n))
            for j in range(n):
                h = step * max(abs(x[j]), 1.0)
                x_fwd = x.copy()
                x_fwd[j] += h
                J[:, j] = (self._solve_circuit(x_fwd) - F) / h

            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                return x, False, "singular Jacobian"

            x = x + dx

        return x, False, f"max_iter={max_iter} reached (|F|={norm:.2e})"

    def _run_backend(self, method, x0):
        root_methods = {'hybr', 'lm', 'broyden1', 'broyden2', 'anderson', 'krylov', 'df-sane'}
        if method == 'fsolve':
            return self._solve_fsolve(x0)
        if method == 'newton':
            return self._solve_newton(x0)
        if method in root_methods:
            return self._solve_root(x0, root_method=method)
        raise ValueError(f"Unknown method '{method}'. "
                         f"Choose from: 'fsolve', 'newton', a root method: {root_methods}.")

    def solve(self, method='fsolve', fd_step=None, update_guesses=False,
              max_outer=10, guess_rtol=1e-4):
        """
        Solve the circuit.

        Parameters
        ----------
        method : str
            'fsolve', 'newton', or a scipy.optimize.root method ('hybr', 'lm', ...).
        fd_step : float, optional
            [NEW] Relative finite-difference step for the Jacobian. Use 1e-5..1e-4
            when components (e.g. moving-boundary HX) return slightly noisy results.
        update_guesses : bool
            [NEW] If True, run an outer loop: solve with frozen cycle guesses,
            update the guesses from the result, solve again... until the guesses
            change by less than guess_rtol (relative) or max_outer is reached.
        """
        if fd_step is not None:
            self.fd_step = fd_step

        self.n_it = 0
        self.history = []
        self.converged = False

        if not self.solving_order:
            self._build_solve_order()
        self.reset_solved_marker()

        x = self._get_iteration_vector()
        self.x_log = [x]

        n_outer = max_outer if update_guesses else 1
        success = False

        for k in range(n_outer):
            try:
                x, success, msg = self._run_backend(method, x)
            except Exception as e:
                self.solver_message = f"Error during solving ({method}): {e}"
                print(self.solver_message)
                return None

            self.solver_message = msg
            if not update_guesses:
                break

            change = self._update_cycle_guesses()
            if self.print_flag:
                print(f"[outer {k + 1}] max relative change of cycle guesses: {change:.2e}")
            if change < guess_rtol:
                break

        # Final evaluation at the solution, then the REAL convergence test
        self._solve_circuit(x)
        tol_ok = all(d["ok"] for d in self.last_residuals)
        self.converged = bool(success and tol_ok)

        status = "CONVERGED" if self.converged else "NOT CONVERGED"
        print(f"[{method}] {status} after {self.n_it} evaluations — {self.solver_message.strip()}")
        if not self.converged:
            self.print_residuals()

        return x

    # =======================================================================
    # [NEW] Debugging tools
    # =======================================================================

    @contextmanager
    def _quiet(self):
        """Silence prints and keep history/n_it untouched during debug evaluations."""
        old_flag = getattr(self, "print_flag", False)
        old_hist, old_n = self.history, self.n_it
        self.print_flag = False
        self.history = []
        try:
            yield
        finally:
            self.print_flag = old_flag
            self.history, self.n_it = old_hist, old_n

    @staticmethod
    def _is_mass_connector(obj):
        return all(hasattr(obj, a) for a in ("p", "h", "m_dot", "reset"))

    def _mass_connectors(self, component_name):
        """List of (name, connector) of a component, supply ports first."""
        model = self.get_component(component_name).model
        conns = [(n, o) for n, o in vars(model).items() if self._is_mass_connector(o)]

        def order(item):
            n = item[0]
            return (0 if n.startswith("su") else 1 if n.startswith("ex") else 2, n)

        return sorted(conns, key=order)

    @staticmethod
    def _fmt(v, factor=1.0, offset=0.0, spec="10.3f"):
        if v is None:
            return f"{'-':>10}"
        try:
            return format(v * factor + offset, spec)
        except (TypeError, ValueError):
            return f"{str(v):>10}"

    @staticmethod
    def _phase_info(p, h, fluid):
        """Phase label with subcooling / superheat / quality."""
        if p is None or h is None or fluid is None:
            return "-"
        try:
            h_l = PropsSI("H", "P", p, "Q", 0, fluid)
            h_v = PropsSI("H", "P", p, "Q", 1, fluid)
            T_l = PropsSI("T", "P", p, "Q", 0, fluid)
            T_v = PropsSI("T", "P", p, "Q", 1, fluid)
            T = PropsSI("T", "P", p, "H", h, fluid)
        except Exception:
            return "supercritical / n.a."
        if h < h_l:
            return f"liquid, SC = {T_l - T:.2f} K"
        if h > h_v:
            return f"vapour, SH = {T - T_v:.2f} K"
        return f"TWO-PHASE, x = {(h - h_l) / (h_v - h_l):.3f}"

    def print_component_states(self, component_name):
        """Print supply/exhaust states of one component."""
        print(f"\n--- {component_name} ---")
        print(f"{'port':<8}{'fluid':<12}{'p [bar]':>10}{'T [°C]':>10}"
              f"{'h [kJ/kg]':>10}{'m [kg/s]':>10}   phase")
        for name, c in self._mass_connectors(component_name):
            p, T, h = getattr(c, "p", None), getattr(c, "T", None), getattr(c, "h", None)
            m, fluid = getattr(c, "m_dot", None), getattr(c, "fluid", None)
            print(f"{name:<8}{str(fluid):<12}"
                  f"{self._fmt(p, 1e-5)}{self._fmt(T, 1, -273.15, '10.2f')}"
                  f"{self._fmt(h, 1e-3, 0, '10.2f')}{self._fmt(m, 1, 0, '10.4f')}"
                  f"   {self._phase_info(p, h, fluid)}")

    def print_all_states(self):
        """Print supply/exhaust states of every component, in solving order."""
        order = self.solving_order or list(self.components)
        for name in order:
            self.print_component_states(name)

    def _side_balances(self, component_name):
        """Enthalpy flow change m·(h_ex - h_su) per side ('', '_H', '_C', ...)."""
        sides = {}
        for cname, c in self._mass_connectors(component_name):
            if cname.startswith("su"):
                side, key = cname[2:], "su"
            elif cname.startswith("ex"):
                side, key = cname[2:], "ex"
            else:
                continue
            s = sides.setdefault(side, {"fluid": getattr(c, "fluid", None)})
            m, h = getattr(c, "m_dot", None), getattr(c, "h", None)
            s["m_" + key] = m
            s["H_" + key] = m * h if (m is not None and h is not None) else None
        for s in sides.values():
            H_su, H_ex = s.get("H_su"), s.get("H_ex")
            s["dH"] = H_ex - H_su if (H_su is not None and H_ex is not None) else None
        return sides

    def print_energy_balance(self):
        """
        Per-component enthalpy flow change, and closure of the working-fluid loop.
        For the working fluid, the sum of m·(h_ex - h_su) over all components is
        exactly the mismatch at the links: it must be ~0 for a closed cycle.
        """
        print(f"\n{'component':<14}{'side':<6}{'fluid':<12}{'m_su':>9}{'m_ex':>9}{'ΔH [W]':>14}")
        total_wf, missing = 0.0, False
        order = self.solving_order or list(self.components)
        for name in order:
            for side, s in self._side_balances(name).items():
                print(f"{name:<14}{side or '-':<6}{str(s['fluid']):<12}"
                      f"{self._fmt(s.get('m_su'), 1, 0, '9.4f')[-9:]}"
                      f"{self._fmt(s.get('m_ex'), 1, 0, '9.4f')[-9:]}"
                      f"{self._fmt(s['dH'], 1, 0, '14.1f')}")
                if str(s["fluid"]).lower() == str(self.fluid).lower():
                    if s["dH"] is None:
                        missing = True
                    else:
                        total_wf += s["dH"]
        note = " (some sides missing)" if missing else ""
        print(f"\nWorking-fluid loop closure error: {total_wf:.1f} W  (should be ≈ 0){note}")
        return total_wf

    def print_iteration_variables(self, x=None):
        """Initial guess and current value of each iteration variable."""
        if x is None:
            x = self.history[-1]["x"] if self.history else self._get_iteration_vector()
        print(f"\n{'iteration variable':<45}{'x0':>14}{'current':>14}")
        for val, it in zip(x, self.it_vars.values()):
            print(f"{it['name']:<45}{it['x0']:>14.6g}{val:>14.6g}")

    def print_residuals(self):
        """Detailed table of the last residual evaluation."""
        print(f"\n{'residual':<38}{'var':<5}{'calculated':>14}{'reference':>14}"
              f"{'diff':>13}{'scaled':>11}{'tol':>10}  ok")
        for d in self.last_residuals:
            print(f"{d['target']:<38}{d['compared']:<5}{d['calculated']:>14.6g}"
                  f"{d['reference']:>14.6g}{d['diff']:>13.4g}{d['scaled']:>11.4g}"
                  f"{d['tol']:>10.2g}  {'✓' if d['ok'] else '✗'}")

    def print_cycle_guesses(self, rtol=1e-6):
        """
        Applied value vs value recalculated by the sweep for each cycle guess.
        A guess that differs from its recalculated value was used as an input
        that the cycle does not agree with: it is a hidden loop. Make it an
        iteration variable with a closure residual (or use update_guesses=True).
        """
        if not self.guess_drift:
            print("No cycle guesses (or no evaluation yet).")
            return
        print(f"\n{'cycle guess':<40}{'applied':>14}{'recalculated':>14}  status")
        for name, (applied, calc) in self.guess_drift.items():
            if calc is None:
                status = "not recalculated"
            elif abs(calc - applied) <= rtol * max(abs(applied), 1e-12):
                status = "ok"
            else:
                status = "⚠ HIDDEN INPUT — cycle disagrees"
            calc_s = f"{calc:>14.6g}" if calc is not None else f"{'-':>14}"
            print(f"{name:<40}{applied:>14.6g}{calc_s}  {status}")

    def check_determinism(self, x=None, n_repeat=3, tol=1e-9):
        """
        Evaluate the residuals n_repeat times at the same x. If they differ, the
        solver cannot work (finite-difference Jacobian is meaningless).
        """
        x = np.asarray(x if x is not None else self._get_iteration_vector(), dtype=float)
        with self._quiet():
            R = np.array([self._solve_circuit(x) for _ in range(n_repeat)])
        spread = np.max(np.abs(R - R[0]), axis=0)
        print(f"\nDeterminism check ({n_repeat} evaluations at the same x):")
        for rv, r, s in zip(self.res_vars, R.T, spread):
            print(f"  {rv['target']:<38} values {np.array2string(r, precision=6)}  spread {s:.2e}")
        ok = bool(np.all(spread <= tol))
        print("  -> deterministic ✓" if ok else
              "  -> NOT deterministic ✗ (something keeps state between calls)")
        return ok

    def scan_iteration_variable(self, which, values, x=None, plot=True):
        """
        Evaluate the residuals while varying one iteration variable (others fixed).
        Shows discontinuities, flat regions and where the residuals cross zero.

        which : name (key of it_vars) or index of the iteration variable
        values : iterable of values to test
        """
        names = list(self.it_vars)
        j = names.index(which) if isinstance(which, str) else int(which)
        base = np.asarray(x if x is not None else self._get_iteration_vector(), dtype=float)
        values = np.asarray(list(values), dtype=float)

        R = []
        with self._quiet():
            for v in values:
                xx = base.copy()
                xx[j] = v
                try:
                    R.append(self._solve_circuit(xx))
                except Exception as e:
                    print(f"  {names[j]} = {v:.6g}: failed ({e})")
                    R.append(np.full(len(self.res_vars), np.nan))
        R = np.array(R)

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            for k, rv in enumerate(self.res_vars):
                ax.plot(values, R[:, k], "o-", ms=3, label=rv["target"])
            ax.axhline(0, color="k", lw=0.8)
            ax.set_xlabel(names[j])
            ax.set_ylabel("scaled residual [-]")
            ax.legend()
            ax.grid(alpha=0.3)
            plt.show()

        return values, R

    def plot_history(self):
        """Residuals and iteration variables vs evaluation number (last solve)."""
        import matplotlib.pyplot as plt
        if not self.history:
            print("No history (run solve() first).")
            return
        its = [h["it"] for h in self.history]
        R = np.array([h["res"] for h in self.history])
        X = np.array([h["x"] for h in self.history])

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
        for k, rv in enumerate(self.res_vars):
            ax1.semilogy(its, np.abs(R[:, k]), "o-", ms=3, label=rv["target"])
        ax1.set_ylabel("|scaled residual|")
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3)
        for k, name in enumerate(self.it_vars):
            ax2.plot(its, X[:, k] / X[0, k], "o-", ms=3, label=name)
        ax2.set_ylabel("x / x0")
        ax2.set_xlabel("evaluation")
        ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3)
        plt.tight_layout()
        plt.show()

    def debug_report(self):
        """Everything at once, for the last evaluated state."""
        print("\n" + "=" * 70 + "\nDEBUG REPORT\n" + "=" * 70)
        print(f"Converged: {self.converged} — {self.solver_message.strip()}")
        print(f"Residual evaluations: {self.n_it}")
        self.print_iteration_variables()
        self.print_residuals()
        self.print_cycle_guesses()
        self.print_all_states()
        self.print_energy_balance()
