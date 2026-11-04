"""
Lumped pressure-loss + volume component based on a flow coefficient Kv.

Kv definition (IEC 60534 / VDI 2173): volumetric flow of water [m3/h]
(rho0 = 1000 kg/m3) passing through the component with a pressure drop of 1 bar.

    Q [m3/h] = Kv * sqrt( dp[bar] / (rho / rho0) )

Rewritten with mass flow (SI units):

    dp [Pa] = (1e5 / rho0) * (3600 * m_dot / Kv)**2 / rho  =  C * m_dot**2 / (rho * Kv**2)
    with C = 1e5 * 3600**2 / 1000 = 1.296e9

This is exactly an orifice law dp = m_dot^2 / (2 rho (Cd*A)^2) with
    (Cd*A)_eq = Kv / 50912   [m2]
"""
import warnings

import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import brentq

from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector

"Constant for Kv to pressure drop conversion (SI units)"
RHO0 = 1000.0          # reference density for Kv [kg/m3]
DP0 = 1e5              # reference pressure drop for Kv [Pa]
C_KV = DP0 * 3600.0**2 / RHO0   # 1.296e9
AV_PER_KV = np.sqrt(RHO0 / DP0) / 3600.0      # 2.778e-5 m2 per (m3/h)
# --------------------------------------------------------------------------
# Component
# --------------------------------------------------------------------------
class KvComponent(BaseComponent):
    """
    **Component**: Kv Component

    **Model**:


    **Description**:


    **Assumptions**:


    **Parameters**
    ----------
    Kv        : flow coefficient when open [m3/h]
    V_int     : internal volume [m3]
    opening   : 0..1, scales Kv (1 = fully open, 0 = closed). Allows ramps during
                mode transitions. Kv_eff = opening * Kv.          (default 1.0)
    x_T       : pressure-drop ratio factor (IEC 60534)            (default 0.63)
    F_L       : liquid pressure recovery factor (IEC 60534)       (default 0.9)
    x_lim     : dp/P_su above which vapour flow is reported as compressible (default 0.05)

    **Inputs**
    ----------
    fluid, P_su [Pa], h_su [J/kg], m_dot [kg/s]

    **Outputs**
    ----------
    dp [Pa], P_ex [Pa], h_ex [J/kg] (isenthalpic), charge [kg], regime [str], flags [dict]
    solve() also returns these results as a dict.

    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()  # Inlet
        self.ex = MassConnector()  # Outlet
        # Default parameters (can be overridden with set_parameters)
        self.set_parameters(opening=1.0, x_T=0.63, F_L=0.9, x_lim=0.05)

    def get_required_inputs(self):
        return ['P_su', 'h_su', 'm_dot', 'fluid']

    def get_required_parameters(self):
        return ['Kv', 'V_int', 'opening', 'x_T', 'F_L', 'x_lim']
    
    # ---------------------------------------------------------------- pure physics (no CoolProp: testable on their own)
    @staticmethod
    def dp_vapour(Av, rho, P_su, m_dot, Fxt):
        """Inverse of m_dot_vapour. Returns (dp, choked). If m_dot exceeds the choked
        flow, no outlet pressure can deliver it: returns dp = Fxt*P_su and choked=True."""
        # m = Av sqrt(rho P) * (s - s^3/(3 Fxt)), s = sqrt(x), increasing for s in [0, sqrt(Fxt)]
        c = m_dot / (Av * np.sqrt(rho * P_su))
        s_max = np.sqrt(Fxt)
        if c >= (2.0 / 3.0) * s_max:
            return Fxt * P_su, True
        s = brentq(lambda s: s - s**3 / (3 * Fxt) - c, 0.0, s_max)
        return s**2 * P_su, False

    @staticmethod
    def liquid_choked_dp(P_su, p_sat, p_crit, F_L):
        """Max effective dp before flashing at the vena contracta chokes the flow (IEC / ValveVaporizing)."""
        F_F = 0.96 - 0.28 * np.sqrt(p_sat / p_crit)
        return F_L**2 * (P_su - F_F * p_sat)


    # ---------------------------------------------------------------- state
    def _fluid_state(self, fluid, pressure, enthalpy):
        AS = CP.AbstractState("HEOS", fluid)
        AS.update(CP.HmassP_INPUTS, enthalpy, pressure)
        p_crit = AS.p_critical()
        rho = AS.rhomass()

        p_sat = None
        if pressure >= p_crit:
            phase = "supercritical"
        else:
            AS_sat = CP.AbstractState("HEOS", fluid)
            AS_sat.update(CP.PQ_INPUTS, pressure, 0.0)
            h_l = AS_sat.hmass()
            AS_sat.update(CP.PQ_INPUTS, pressure, 1.0)
            h_v = AS_sat.hmass()
            phase = "liquid" if enthalpy <= h_l else "vapour" if enthalpy >= h_v else "two-phase"
            if phase == "liquid":
                # saturation pressure at the (subcooled) liquid temperature
                AS_sat.update(CP.QT_INPUTS, 0.0, AS.T())
                p_sat = AS_sat.p()

        F_gamma = 1.0
        if phase in ("vapour", "supercritical"):
            F_gamma = AS.cpmass() / AS.cvmass() / 1.4
        return dict(rho=rho, phase=phase, p_sat=p_sat, p_crit=p_crit, F_gamma=F_gamma) 

    def _Av(self):
        return AV_PER_KV * self.params["Kv"] * self.params["opening"]

    # ---------------------------------------------------------------
    def solve(self):
        self.check_calculable()
        self.check_parametrized()
        if not (self.calculable and self.parametrized):
            self.solved = False
            return None

        st, Av = self._fluid_state(self.su.fluid, self.su.p, self.su.h), self._Av()
        P_su, m_dot = self.su.p, self.su.m_dot
        Fxt = st["F_gamma"] * self.params["x_T"]
        flags = {}

        if Av <= 0.0:
            if m_dot != 0:
                raise ValueError("Closed valve with non-zero mass flow.")
            dp = 0.0
        elif st["phase"] in ("vapour", "supercritical") and st["rho"] < 0.3 * RHO0:
            dp, flags["choked"] = self.dp_vapour(Av, st["rho"], P_su, m_dot, Fxt)
        else:
            dp = m_dot**2 / (st["rho"] * Av**2)
            if st["phase"] == "liquid":
                flags["choked"] = dp >= self.liquid_choked_dp(P_su, st["p_sat"], st["p_crit"], self.params["F_L"])

        self._finish(st, dp, m_dot, flags)
        self._compute_charge()
        self.solved = True
        return dict(dp=self.dp, P_ex=self.ex.p, h_ex=self.ex.h, m_dot=m_dot,
                    charge=self.charge, regime=self.regime, flags=self.flags)

    # ---------------------------------------------------------------- common post-processing
    def _finish(self, st, dp, m_dot, flags):
        fl, P_su, h_su = self.su.fluid, self.su.p, self.su.h
        x = dp / P_su
        P_ex, h_ex = P_su - dp, h_su
        self.ex.set_p(P_ex)
        self.ex.set_h(h_ex)
        self.ex.set_fluid(fl)
        self.ex.set_m_dot(m_dot)

        compressible = st["phase"] in ("vapour", "supercritical") and x > self.params["x_lim"]
        flags.setdefault("choked", False)
        flags["compressible"] = compressible
        flags["two_phase_inlet"] = st["phase"] == "two-phase"

        if flags["choked"]:
            regime = "choked (flashing)" if st["phase"] == "liquid" else "choked"
        elif flags["two_phase_inlet"]:
            regime = "two-phase (outside IEC 60534 validity)"
        elif compressible:
            regime = "compressible"
        else:
            regime = "incompressible"

        if flags["choked"] or flags["two_phase_inlet"]:
            warnings.warn(f"KvComponent: {regime} at P_su={P_su/1e5:.2f} bar, x={x:.3f}")

        self.dp, self.x, self.regime, self.flags, self.phase_su = dp, x, regime, flags, st["phase"]

    def _compute_charge(self):
        rho_su = self._fluid_state(self.su.fluid, self.su.p, self.su.h)['rho']
        rho_ex = self._fluid_state(self.ex.fluid, self.ex.p, self.ex.h)['rho']

        rho_mean = 0.5 * (rho_su + rho_ex)
        V_int = float(self.params.get('V_int'))
        self.charge = rho_mean * V_int







# """
# Lumped flow resistance + internal volume, parametrised by a flow coefficient Kv.
# The component detects the flow regime on its own (liquid / vapour / two-phase,
# incompressible / compressible / choked / flashing).

# Flow equations: IEC 60534-2-1, as implemented in Modelica.Fluid.Valves
# (ValveIncompressible, ValveCompressible, ValveVaporizing).

#     m_dot = Av * Y * sqrt(rho_su * dp_eff)        Av = Kv / 3600 * sqrt(rho0/dp0) = 2.778e-5 * Kv [m2]

#   - liquid  : Y = 1,  dp_eff = min(dp, F_L^2 (P_su - F_F p_sat))      (flashing / liquid choking)
#   - vapour  : Y = 1 - x/(3 Fxt), x = min(dp/P_su, Fxt), Fxt = F_gamma * x_T  (gas choking)
#               Y -> 1 when x -> 0, so the same equation covers the incompressible limit
#               (no switch, no discontinuity).
#   - 2-phase : homogeneous density, incompressible law, flagged as outside IEC 60534 validity.

# Two calling modes:
#   - solve_dp(m_dot)   : imposed mass flow  -> pressure drop   (steady-state cycle models)
#   - solve_m_dot(P_ex) : imposed pressures  -> mass flow       (transients; well-defined when choked)

# Defaults (to override when better data exist): x_T = 0.63 (Modelica Fxt_full = 0.5 for F_gamma ~ 0.79),
# F_L = 0.9 (Modelica Fl_nominal).

# NOTE: method names must be aligned with labothappy's BaseComponent/connectors before merging.
# """

# # =============================================================================
# # Pure physics (no CoolProp: testable on their own)
# # =============================================================================
# def liquid_choked_dp(P_su, p_sat, p_crit, F_L):
#     """Max effective dp before flashing at the vena contracta chokes the flow (IEC / ValveVaporizing)."""
#     F_F = 0.96 - 0.28 * np.sqrt(p_sat / p_crit)
#     return F_L**2 * (P_su - F_F * p_sat)


# def m_dot_liquid(Av, rho, P_su, P_ex, p_sat, p_crit, F_L):
#     dp = P_su - P_ex
#     dp_ch = liquid_choked_dp(P_su, p_sat, p_crit, F_L)
#     choked = dp >= dp_ch
#     return Av * np.sqrt(rho * min(dp, dp_ch)), choked


# def m_dot_vapour(Av, rho, P_su, P_ex, Fxt):
#     x = min((P_su - P_ex) / P_su, Fxt)
#     Y = 1.0 - x / (3.0 * Fxt)
#     return Av * Y * np.sqrt(rho * x * P_su), x >= Fxt, x, Y



# # =============================================================================
# # Component
# # =============================================================================
# class KvComponent:
#     """
#     Parameters: Kv [m3/h], V_int [m3], opening [0-1], x_T [-], F_L [-],
#                 x_lim [-] threshold on dp/P_su above which compressibility is *reported*
#                 (the Y correction is always applied to vapour anyway).
#     Inputs:     fluid, P_su [Pa], h_su [J/kg], and either m_dot [kg/s] or P_ex [Pa]
#     Outputs:    dp, P_ex, h_ex (isenthalpic), m_dot, M [kg], regime, flags
#     """

#     def __init__(self, Kv, V_int, opening=1.0, x_T=0.63, F_L=0.9, x_lim=0.05):
#         self.p = dict(Kv=Kv, V_int=V_int, opening=opening, x_T=x_T, F_L=F_L, x_lim=x_lim)
#         self.inputs, self.outputs = {}, {}

#     def set_parameters(self, **kw): self.p.update(kw)
#     def set_inputs(self, **kw):     self.inputs.update(kw)

#     # ---------------------------------------------------------------- mode 1: m_dot -> dp
#     def solve_dp(self):
#         st, Av = self._inlet_state(), self._Av()
#         P_su, m_dot = self.inputs["P_su"], self.inputs["m_dot"]
#         Fxt = st["F_gamma"] * self.p["x_T"]
#         flags = {}

#         if Av <= 0.0:
#             if m_dot != 0:
#                 raise ValueError("Closed valve with non-zero mass flow.")
#             dp = 0.0
#         elif st["phase"] in ("vapour", "supercritical") and st["rho"] < 0.3 * RHO0:
#             dp, flags["choked"] = dp_vapour(Av, st["rho"], P_su, m_dot, Fxt)
#         else:
#             dp = m_dot**2 / (st["rho"] * Av**2)
#             if st["phase"] == "liquid":
#                 flags["choked"] = dp >= liquid_choked_dp(P_su, st["p_sat"], st["p_crit"], self.p["F_L"])

#         return self._finish(st, dp, m_dot, flags)

#     # ---------------------------------------------------------------- mode 2: P_ex -> m_dot
#     def solve_m_dot(self):
#         st, Av = self._inlet_state(), self._Av()
#         P_su, P_ex = self.inputs["P_su"], self.inputs["P_ex"]
#         Fxt = st["F_gamma"] * self.p["x_T"]
#         flags = {}

#         if Av <= 0.0 or P_ex >= P_su:          # closed, or no driving dp (reverse flow not handled here)
#             m_dot = 0.0
#         elif st["phase"] in ("vapour", "supercritical") and st["rho"] < 0.3 * RHO0:
#             m_dot, flags["choked"], _, _ = m_dot_vapour(Av, st["rho"], P_su, P_ex, Fxt)
#         elif st["phase"] == "liquid":
#             m_dot, flags["choked"] = m_dot_liquid(Av, st["rho"], P_su, P_ex,
#                                                   st["p_sat"], st["p_crit"], self.p["F_L"])
#         else:
#             m_dot = Av * np.sqrt(st["rho"] * (P_su - P_ex))

#         return self._finish(st, P_su - P_ex, m_dot, flags)

#     # ---------------------------------------------------------------- common post-processing
#     def _finish(self, st, dp, m_dot, flags):
#         from CoolProp.CoolProp import PropsSI
#         fl, P_su, h_su = self.inputs["fluid"], self.inputs["P_su"], self.inputs["h_su"]
#         x = dp / P_su
#         P_ex, h_ex = P_su - dp, h_su
#         rho_ex = PropsSI("D", "P", P_ex, "H", h_ex, fl)

#         compressible = st["phase"] in ("vapour", "supercritical") and x > self.p["x_lim"]
#         flags.setdefault("choked", False)
#         flags["compressible"] = compressible
#         flags["two_phase_inlet"] = st["phase"] == "two-phase"

#         if flags["choked"]:
#             regime = "choked (flashing)" if st["phase"] == "liquid" else "choked"
#         elif flags["two_phase_inlet"]:
#             regime = "two-phase (outside IEC 60534 validity)"
#         elif compressible:
#             regime = "compressible"
#         else:
#             regime = "incompressible"

#         if flags["choked"] or flags["two_phase_inlet"]:
#             warnings.warn(f"KvComponent: {regime} at P_su={P_su/1e5:.2f} bar, x={x:.3f}")

#         rho_mean = 0.5 * (st["rho"] + rho_ex)
#         self.outputs = dict(dp=dp, P_ex=P_ex, h_ex=h_ex, m_dot=m_dot, x=x,
#                             phase_su=st["phase"], regime=regime, flags=flags,
#                             rho_su=st["rho"], rho_ex=rho_ex,
#                             M=rho_mean * self.p["V_int"])
#         return self.outputs


# if __name__ == "__main__":
#     from CoolProp.CoolProp import PropsSI
#     fl = "R1233ZDE"
#     evr = KvComponent(Kv=9.8, V_int=50e-6)      # V_int: placeholder, to measure
#     cases = [("liquid line",    10e5, PropsSI("T", "P", 10e5, "Q", 0, fl) - 5),
#              ("suction vapour",  2e5, PropsSI("T", "P", 2e5,  "Q", 1, fl) + 5)]
#     for label, P, T in cases:
#         evr.set_inputs(fluid=fl, P_su=P, h_su=PropsSI("H", "P", P, "T", T, fl), m_dot=0.3)
#         o = evr.solve_dp()
#         print(f"{label:15s} {o['phase_su']:8s} dp={o['dp']/1e5:6.3f} bar  x={o['x']:.3f}  "
#               f"regime={o['regime']}  M={o['M']*1e3:.2f} g")