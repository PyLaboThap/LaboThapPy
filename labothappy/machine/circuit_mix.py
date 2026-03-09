# -*- coding: utf-8 -*-
"""
circuit_mix.py
==============
Unified thermodynamic circuit solver combining the best of circuit_rec and
circuit_it.

Architecture
------------
- Component ordering  : DFS recursive traversal (from circuit_rec)
- Problem definition  : set_cycle_guess / set_fixed_properties /
                        set_iteration_variable (from circuit_rec)
- Convergence check   : global p/h snapshot on all MassConnector output ports
                        (from circuit_rec)
- Numerical methods   : pluggable via solve(method=...)

    Substitution-based (no residual vector needed):
        'successive_substitution'   plain fixed-point
        'wegstein'                  Wegstein acceleration (1 variable at a time)

    Newton-based (residual vector built automatically from guesses):
        'fsolve'                    scipy fsolve  — Newton-Powell hybrid
        'hybr'                      scipy root / hybr — same algorithm
        'lm'                        Levenberg-Marquardt
        'broyden1'                  Broyden quasi-Newton
        'anderson'                  Anderson mixing
        'krylov'                    Newton-Krylov (large systems)
        'df-sane'                   derivative-free

Notes
-----
For Newton-based methods the guess variables are automatically used as the
iteration vector x.  The residual is  F(x) = f(x) - x  where f(x) is the
value read back from the connectors after one full sequential sweep.
The user does NOT need to define residual variables manually.

Authors
-------
Based on circuit_rec.py  (Basile Chaudoir, 2024)
     and circuit_it.py   (Elise Neven, 2026)
Unified by: <your name>
"""

from __future__ import annotations

import warnings
import numpy as np
from scipy.optimize import fsolve, root

from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
from labothappy.connector.heat_connector import HeatConnector

from machine.base_circuit import BaseCircuit


# ---------------------------------------------------------------------------
# Substitution-based methods supported natively
# ---------------------------------------------------------------------------
_SUBST_METHODS = frozenset(['successive_substitution', 'wegstein'])

# Newton-based methods delegated to scipy
_NEWTON_METHODS = frozenset(['fsolve', 'hybr', 'lm', 'broyden1',
                              'broyden2', 'anderson', 'krylov',
                              'diagbroyden', 'linearmixing', 'df-sane'])

_ALL_METHODS = _SUBST_METHODS | _NEWTON_METHODS


class MixedCircuit(BaseCircuit):
    """
    Unified sequential-modular circuit solver with pluggable numerical methods.

    Usage
    -----
    >>> circuit = MixedCircuit(fluid='R1233ZDE')
    >>> # ... add components, connect ports ...
    >>> circuit.set_fixed_properties(target='Pump:su', p=1e5, fluid='R1233ZDE')
    >>> circuit.set_cycle_guess(target='Pump:su', m_dot=0.1)
    >>> circuit.solve(method='wegstein')          # or 'anderson', 'lm', etc.
    """

    # -----------------------------------------------------------------------
    # Inner classes (identical to circuit_rec)
    # -----------------------------------------------------------------------

    class Guess:
        def __init__(self, target: str, variable: str, value: float):
            self.name    = f"{target}-{variable}"
            self.target  = target
            self.variable = variable
            self.value   = value
            self.input_history:  list[float] = []
            self.output_history: list[float] = []

        def record_input(self, v: float):  self.input_history.append(float(v))
        def record_output(self, v: float): self.output_history.append(float(v))

        @property
        def prev_input(self):
            return self.input_history[-2]  if len(self.input_history)  >= 2 else None
        @property
        def prev_output(self):
            return self.output_history[-2] if len(self.output_history) >= 2 else None

    class Fixed_Property:
        def __init__(self, target: str, variable: str, value: float):
            self.name     = f"{target}-{variable}"
            self.target   = target
            self.variable = variable
            self.value    = value

    class Iteration_variable:
        def __init__(self, target, variable, objective,
                     objective_value, tol, rel, damping_factor, cycle):
            self.target          = target
            self.variable        = variable
            self.objective       = objective
            self.objective_value = objective_value
            self.tol             = tol
            self.converged       = False
            self.damping_factor  = damping_factor
            self.cycle           = cycle
            if rel not in (1, -1):
                raise ValueError("'rel' must be 1 or -1.")
            self.rel = rel

        def find_link_DP(self, components):
            """BFS path finder — identical to circuit_rec implementation."""
            from collections import deque

            start = self.target[0].split(":")[0]
            end   = self.objective.split(":")[1]
            start_comp = components[start]

            def get_suffix(port_name):
                if port_name is None: return ""
                if port_name.endswith("_C"): return "_C"
                if port_name.endswith("_H"): return "_H"
                return ""

            def shortest_path(start_component, end_name):
                queue   = deque([(start_component, [[start_component.name, None, None]])])
                visited = {start_component}
                while queue:
                    current, path = queue.popleft()
                    if current.name == end_name:
                        return path
                    _, entrance_port_current, _ = path[-1]
                    entrance_suffix = get_suffix(entrance_port_current)
                    for prev_exit_port, next_comp in current.next.items():
                        if next_comp in visited:
                            continue
                        entrance_port_next = next(
                            (p for p, prev in next_comp.previous.items() if prev is current), None)
                        if entrance_suffix in ("_C", "_H") and \
                                not prev_exit_port.endswith(entrance_suffix):
                            continue
                        visited.add(next_comp)
                        queue.append(
                            (next_comp, path + [[next_comp.name, entrance_port_next, prev_exit_port]]))
                return None

            path = shortest_path(start_comp, end)
            for i in range(1, len(path)):
                path[i - 1][-1] = path[i][-1]
                if i == len(path) - 1:
                    path[i][-1] = None

            delta_Ps = {}
            for comp_name, entrance_port, exit_port in path:
                if entrance_port is None or exit_port is None:
                    continue
                comp   = components[comp_name].model
                suffix = ""
                if entrance_port.endswith("_C") or exit_port.endswith("_C"):
                    suffix = "_C"
                elif entrance_port.endswith("_H") or exit_port.endswith("_H"):
                    suffix = "_H"
                if suffix == "_C":   delta_Ps[comp_name] = getattr(comp, "DP_c", 0)
                elif suffix == "_H": delta_Ps[comp_name] = getattr(comp, "DP_h", 0)
                else:                delta_Ps[comp_name] = 0

            self.delta_Ps = delta_Ps
            self.DP       = sum(delta_Ps.values())
            return self.DP

        def check(self, target_value, objective_connector, var):
            if abs((target_value - self.objective_value) / self.objective_value) < self.tol:
                self.converged = True
                return None
            self.converged = False
            temp = MassConnector(fluid=objective_connector.fluid)
            temp.set_properties(**{var: self.objective_value, 'T': objective_connector.T})
            return getattr(temp, self.variable)

    # -----------------------------------------------------------------------
    # Constructor
    # -----------------------------------------------------------------------

    def __init__(self, fluid: str | None = None):
        super().__init__()

        # Topology
        self.components: dict = {}
        self.sources:    dict = {}
        self.sinks:      dict = {}

        # Problem definition
        self.fluid            = fluid
        self.fixed_properties: dict = {}
        self.parameters:       dict = {}
        self.guesses:          dict = {}
        self.it_vars:          list = []
        self.res_vars:         dict = {}   # kept for backward compat

        # Solve state
        self.guess_update          = False
        self.solve_start_components: list = []
        self.solving_order:          list = []
        self.converged               = False
        self.messages:               list = []
        self._prev_connector_states: dict = {}

        # Convergence tolerance (p/h snapshot)
        self.convergence_tolerance = 1e-4

        # Wegstein clipping bounds
        self.wegstein_q_min = -5.0
        self.wegstein_q_max =  0.0

    # -----------------------------------------------------------------------
    # Connector snapshot helpers  (identical to circuit_rec)
    # -----------------------------------------------------------------------

    def _snapshot_connector_states(self) -> dict:
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
                    snapshot[comp_name][port_name] = {
                        'p': float(p_val), 'h': float(h_val)}
        return snapshot

    def _check_connector_convergence(self, current_snapshot) -> tuple[bool, list]:
        epsilon      = 1e-10
        all_converged = True
        messages     = []
        for comp_name, ports in current_snapshot.items():
            prev_ports = self._prev_connector_states.get(comp_name, {})
            for port_name, vals in ports.items():
                prev_vals = prev_ports.get(port_name)
                if prev_vals is None:
                    all_converged = False
                    continue
                for prop in ('p', 'h'):
                    v_new   = vals[prop]
                    v_old   = prev_vals[prop]
                    rel_err = abs(v_new - v_old) / max(abs(v_old), epsilon)
                    if rel_err > self.convergence_tolerance:
                        all_converged = False
                        messages.append(
                            f"Not converged: {comp_name}:{port_name}-{prop} "
                            f"rel_err={rel_err:.2e}")
        return all_converged, messages

    # -----------------------------------------------------------------------
    # Static numerical helpers
    # -----------------------------------------------------------------------

    @staticmethod
    def _wegstein_update(x_prev, x_curr, f_prev, f_curr,
                         q_min=-5.0, q_max=0.0) -> float:
        """Wegstein acceleration — falls back to SS when history is insufficient."""
        if x_prev is None or f_prev is None:
            return f_curr
        denom = x_curr - x_prev
        if abs(denom) < 1e-12:
            return f_curr
        q = (f_curr - f_prev) / denom
        q = max(q_min, min(q_max, q))
        if abs(q - 1.0) < 1e-12:
            return f_curr
        return (q * x_curr - f_curr) / (q - 1.0)

    # -----------------------------------------------------------------------
    # Problem-definition API  (mirrors circuit_rec)
    # -----------------------------------------------------------------------

    def set_source_properties(self, **kwargs):
        target = kwargs.pop('target')
        source = self.get_source(target)
        source.set_properties(**kwargs)
        for arg_name, value in kwargs.items():
            fp = MixedCircuit.Fixed_Property(target, arg_name, value)
            if fp.name not in self.fixed_properties:
                self.fixed_properties[fp.name] = fp

    def set_fixed_properties(self, **kwargs):
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        for arg_name, value in kwargs.items():
            fp = MixedCircuit.Fixed_Property(target, arg_name, value)
            if fp.name not in self.fixed_properties:
                self.fixed_properties[fp.name] = fp
        component.set_properties(connector_name, **kwargs)
        component.model.check_calculable()

    def set_cycle_guess(self, external_set: bool = False, **kwargs):
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        for arg_name, value in kwargs.items():
            guess = MixedCircuit.Guess(target, arg_name, value)
            if self.guess_update or guess.name not in self.guesses or external_set:
                if guess.name in self.guesses:
                    guess.input_history  = self.guesses[guess.name].input_history
                    guess.output_history = self.guesses[guess.name].output_history
                self.guesses[guess.name] = guess
        component.set_properties(connector_name, **kwargs)
        component.model.check_calculable()

    def set_iteration_variable(self, tol=1e-2, target=None, variable=None,
                               objective=None, rel=1, damping_factor=0.5,
                               cycle=None):
        if objective in self.fixed_properties:
            it_var = MixedCircuit.Iteration_variable(
                target, variable, objective,
                self.fixed_properties[objective].value,
                tol, rel, damping_factor, self)
        elif "Link" in objective:
            it_var = MixedCircuit.Iteration_variable(
                target, variable, objective, None,
                tol, rel, damping_factor, self)
        else:
            raise ValueError(
                f"'{objective}' not in fixed_properties and not a Link objective.")
        self.it_vars.append(it_var)

    def set_cycle_parameters(self, **kwargs):
        self.parameters.update(kwargs)

    def set_residual_variable(self, target, variable, tolerance=1e-3):
        """Kept for backward compatibility — not used for convergence checking."""
        pass

    # -----------------------------------------------------------------------
    # DFS solver  (identical to circuit_rec)
    # -----------------------------------------------------------------------

    def check_all_component_solved(self) -> bool:
        return all(c.model.solved for c in self.components.values())

    def reset_solved_marker(self):
        for c in self.components.values():
            c.model.solved = False

    def recursive_solve(self, component_name: str):
        if self.print_flag:
            print(f"--- Component : {component_name}")
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
            if self.print_flag:
                print(f"  '{component_name}' not calculable.")
            return
        for output_port in component.next:
            self.recursive_solve(component.next[output_port].name)

    # -----------------------------------------------------------------------
    # Core sequential sweep  (shared by ALL methods)
    # -----------------------------------------------------------------------

    def _sequential_sweep(self):
        """One full pass through all components in DFS order."""
        self.reset_solved_marker()
        for component_name in self.solving_order:
            if self.print_flag:
                print(f"--- Component : {component_name}")
            self.components[component_name].solve()

    # -----------------------------------------------------------------------
    # Iteration variable and fixed property enforcement
    # -----------------------------------------------------------------------

    def _enforce_iteration_variables(self):
        """Update guesses based on iteration variable objectives."""
        for it_var in self.it_vars:
            if "Link" in it_var.objective:
                _, obj_comp_name, port_var = it_var.objective.split(":")
                port_name, variable = port_var.split("-")
                gain         = it_var.find_link_DP(self.components) if variable == 'p' else 0
                connector_obj = getattr(self.components[obj_comp_name].model, port_name)
                value_obj    = getattr(connector_obj, variable)
                for target in it_var.target:
                    self.set_cycle_guess(target=target, **{it_var.variable: value_obj + gain})
                it_var.converged = True
            else:
                obj_comp, rest = it_var.objective.split(":")
                port, var      = rest.split("-")
                connector      = getattr(self.components[obj_comp].model, port)
                current_value  = getattr(connector, var)
                if current_value is not None:
                    desired_value = it_var.check(current_value, connector, var)
                    if desired_value is not None:
                        target_comp, target_port = it_var.target[0].split(":")
                        target_connector         = getattr(
                            self.components[target_comp].model, target_port)
                        target_current_value = getattr(target_connector, it_var.variable)
                        new_target_value = (
                            (desired_value - target_current_value)
                            * it_var.damping_factor + target_current_value)
                        for target in it_var.target:
                            self.set_cycle_guess(
                                target=target, **{it_var.variable: new_target_value})

    def _enforce_fixed_properties(self):
        """Re-enforce all fixed properties on their connectors."""
        for fix_prop_name, fix_prop in self.fixed_properties.items():
            if ":" in fix_prop_name:
                self.set_fixed_properties(
                    target=fix_prop.target,
                    **{fix_prop.variable: fix_prop.value})

    # -----------------------------------------------------------------------
    # Substitution-based iteration step
    # -----------------------------------------------------------------------

    def _substitution_step(self, use_wegstein: bool):
        """
        Read connector outputs, compute next guess via SS or Wegstein,
        write back to connectors.
        """
        for guess_name, guess in self.guesses.items():
            component_name, rest = guess_name.split(":")
            port, var            = rest.split("-")
            connector            = getattr(self.components[component_name].model, port)
            f_curr               = float(getattr(connector, var))
            x_curr               = guess.value

            guess.record_input(x_curr)
            guess.record_output(f_curr)

            if use_wegstein and var not in ("SC", "SH"):
                next_value = MixedCircuit._wegstein_update(
                    x_prev=guess.prev_input,
                    x_curr=x_curr,
                    f_prev=guess.prev_output,
                    f_curr=f_curr,
                    q_min=self.wegstein_q_min,
                    q_max=self.wegstein_q_max,
                )
                if self.print_flag:
                    print(f"  Wegstein [{guess_name}]: "
                          f"x={x_curr:.4g} f={f_curr:.4g} → {next_value:.4g}")
            else:
                next_value = f_curr

            self.set_cycle_guess(target=guess.target,
                                 **{guess.variable: next_value})

    # -----------------------------------------------------------------------
    # Newton-based residual function  (used by fsolve / root)
    # -----------------------------------------------------------------------

    def _build_guess_vector(self) -> np.ndarray:
        """Current guess values as a numpy vector (ordered as self.guesses)."""
        return np.array([g.value for g in self.guesses.values()], dtype=float)

    def _apply_guess_vector(self, x: np.ndarray):
        """Write a numpy vector back into the guess connectors."""
        for value, (guess_name, guess) in zip(x, self.guesses.items()):
            self.set_cycle_guess(target=guess.target,
                                 **{guess.variable: float(value)})

    def _residual_function(self, x: np.ndarray) -> np.ndarray:
        """
        F(x) = f(x) - x

        Apply x as guesses, run one sequential sweep, read back connector
        outputs.  Returns the mismatch vector.
        """
        try:
            self._apply_guess_vector(x)
            self._enforce_fixed_properties()
            self._sequential_sweep()

            f_x = np.array(
                [float(getattr(
                    getattr(self.components[gn.split(":")[0]].model,
                            gn.split(":")[1].split("-")[0]),
                    gn.split("-")[1]))
                 for gn in self.guesses],
                dtype=float)
            return f_x - x

        except Exception:
            return np.ones(len(self.guesses)) * 1e6

    # -----------------------------------------------------------------------
    # Convergence check helpers
    # -----------------------------------------------------------------------

    def _check_convergence(self, current_snapshot) -> bool:
        """Full convergence check: p/h connectors + iteration variables + all solved."""
        self.converged = True

        connector_converged, connector_messages = \
            self._check_connector_convergence(current_snapshot)
        if not connector_converged:
            self.converged = False
            self.messages.extend(connector_messages)

        for it_var in self.it_vars:
            if not it_var.converged:
                self.converged = False
                self.messages.append(
                    f"it_var not converged: {it_var.target}-{it_var.variable}")

        if not self.check_all_component_solved():
            self.converged = False
            self.messages.append("Not all components solved.")

        return self.converged

    # -----------------------------------------------------------------------
    # Main solve entry point
    # -----------------------------------------------------------------------

    def solve(self, max_iter: int = 30, method: str = 'fsolve',
              tol: float = 1.49012e-8):
        """
        Solve the circuit.

        Parameters
        ----------
        max_iter : int
            Maximum iterations (substitution-based methods only).
        method : str
            One of: 'successive_substitution', 'wegstein',
                    'fsolve', 'hybr', 'lm', 'broyden1', 'anderson',
                    'krylov', 'df-sane', ...
        tol : float
            Tolerance passed to scipy root/fsolve (Newton-based methods only).
        """
        if method not in _ALL_METHODS:
            raise ValueError(
                f"Unknown method '{method}'. "
                f"Valid options: {sorted(_ALL_METHODS)}")

        if self.print_flag:
            print(f"\nSolve method: {method}")

        # ------------------------------------------------------------------
        # 1. Setup pass — DFS to determine solving order
        # ------------------------------------------------------------------
        for component in self.components:
            if self.components[component].model.check_calculable():
                self.solve_start_components.append(component)

        for start in self.solve_start_components:
            self.recursive_solve(start)
            if self.check_all_component_solved():
                break

        if not self.check_all_component_solved():
            if self.print_flag:
                blocking = [c for c in self.components
                            if not self.components[c].model.check_calculable()]
                print(f"Setup failed. Blocking components: {blocking}")
            return

        if self.print_flag:
            print(f"Setup complete. Solving order: {self.solving_order}")

        # Initial snapshot
        self._prev_connector_states = self._snapshot_connector_states()
        self.guess_update = True

        # ------------------------------------------------------------------
        # 2a. Substitution-based solve loop
        # ------------------------------------------------------------------
        if method in _SUBST_METHODS:
            use_wegstein = (method == 'wegstein')

            for i in range(max_iter):
                self.messages = []
                if self.print_flag:
                    print(f"\n{'#'*30}\nIteration {i+1}\n{'#'*30}")

                self._enforce_iteration_variables()
                self._enforce_fixed_properties()
                self._substitution_step(use_wegstein)
                self._sequential_sweep()

                current_snapshot = self._snapshot_connector_states()
                if self._check_convergence(current_snapshot):
                    if self.print_flag:
                        print(f"Converged in {i+1} iteration(s).")
                    return
                self._prev_connector_states = current_snapshot

            if self.print_flag:
                print(f"Failed to converge in {max_iter} iterations.")

        # ------------------------------------------------------------------
        # 2b. Newton-based solve — residual built automatically from guesses
        # ------------------------------------------------------------------
        else:
            x0 = self._build_guess_vector()

            if method == 'fsolve':
                sol = fsolve(self._residual_function, x0,
                             full_output=False)
            else:
                result = root(
                    self._residual_function, x0,
                    method=method, tol=tol,
                    options={'maxiter': max_iter * 10})
                sol = result.x
                if self.print_flag and not result.success:
                    print(f"[{method}] did not converge: {result.message}")

            # Apply final solution and check
            self._apply_guess_vector(sol)
            self._enforce_fixed_properties()
            self._sequential_sweep()

            current_snapshot = self._snapshot_connector_states()
            self._check_convergence(current_snapshot)

            if self.print_flag:
                status = "Converged" if self.converged else "Did not converge"
                print(f"{status} [{method}].")

    # -----------------------------------------------------------------------
    # Print helpers
    # -----------------------------------------------------------------------

    def print_guesses(self):
        for name, g in self.guesses.items():
            print(f"{name}: {g.value}")

    def print_guess_history(self):
        for name, guess in self.guesses.items():
            print(f"\n--- {name} ---")
            print(f"  {'Iter':<6} {'Input (x)':<20} {'Output f(x)':<20}")
            n = max(len(guess.input_history), len(guess.output_history))
            for i in range(n):
                x  = guess.input_history[i]  if i < len(guess.input_history)  else "-"
                fx = guess.output_history[i] if i < len(guess.output_history) else "-"
                print(f"  {i:<6} "
                      f"{f'{x:.6g}' if isinstance(x, float) else x:<20} "
                      f"{f'{fx:.6g}' if isinstance(fx, float) else fx:<20}")

    def print_fixed_properties(self):
        for name, fp in self.fixed_properties.items():
            print(f"{name}: {fp.value}")
            
            