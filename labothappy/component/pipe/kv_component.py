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
import numpy as np
import CoolProp.CoolProp as CP

from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.correlations.properties.two_phase import compute_two_phase_density, get_saturated_phase_properties
from labothappy.correlations.void_fraction.void_fraction import compute_void_fraction

"Constant for Kv to pressure drop conversion (SI units)"
RHO0 = 1000.0          # reference density for Kv [kg/m3]
DP0 = 1e5              # reference pressure drop for Kv [Pa]
C_KV = DP0 * 3600.0**2 / RHO0   # 1.296e9

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
                mode transitions. Kv_eff = opening * Kv.
    compressible : if True, applies the IEC 60534 expansion factor Y for vapour/gas.
    x_T       : pressure-drop ratio factor (IEC 60534), ~0.7 if unknown.

    **Inputs**
    ----------
    fluid, P_su [Pa], h_su [J/kg], m_dot [kg/s]

    **Outputs**
    ----------
    dp [Pa], P_ex [Pa], h_ex [J/kg] (isenthalpic), M [kg], rho_mean [kg/m3], choked [bool]

    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()  # Inlet
        self.ex = MassConnector()  # Outlet

    def get_required_inputs(self):
        return ['P_su', 'h_su', 'm_dot', 'fluid']

    def get_required_parameters(self):
        return ['Kv', 'V_int', 'opening', 'compressible', 'x_T']

    # ---------------------------------------------------------------
    def solve(self):
        # Get inputs and parameters
        fluid = self.su.fluid
        P_su, h_su, m_dot = self.su.p, self.su.h, self.su.m_dot
        opening = float(self.params.get("opening", 1.0))
        Kv = self.params["Kv"] * opening

        # Compute inlet density (separately for single-phase and two-phase)
        if self.su.is_two_phase():
            # Two-phase: compute saturated properties and mixture density
            rho_l, rho_g = get_saturated_phase_properties(self.su.p, fluid)
            x = self.su.x
            void_fraction = compute_void_fraction(x, rho_l, rho_g, self.params.get("void_fraction_model", "homogeneous"))
            rho_su = compute_two_phase_density(x, rho_l, rho_g, void_fraction)
        else:
            rho_su = self.su.D # Single-phase: use density from the connector

        if Kv <= 0.0:
            if abs(m_dot) > 0:
                raise ValueError("Closed valve (Kv=0) with non-zero mass flow.")
            dp, choked = 0.0, False
        else:
            dp = C_KV * m_dot * abs(m_dot) / (rho_su * Kv**2)  # Incompressible flow

        




#         if Kv <= 0.0:                    # closed valve
#             if abs(m_dot) > 0:
#                 raise ValueError("Closed valve (Kv=0) with non-zero mass flow.")
#             dp, choked = 0.0, False
#         else:
#             dp = C_KV * m_dot * abs(m_dot) / (rho_su * Kv**2)     # incompressible law
#             choked = False
#             if self.params["compressible"]:
#                 dp, choked = self._compressible_dp(fl, P_su, h_su, m_dot, rho_su, Kv)

#         P_ex = P_su - dp
#         h_ex = h_su                       # adiabatic, no work -> isenthalpic
#         rho_ex = PropsSI("D", "P", P_ex, "H", h_ex, fl)
#         rho_mean = 0.5 * (rho_su + rho_ex)

#         self.outputs = dict(dp=dp, P_ex=P_ex, h_ex=h_ex,
#                             M=rho_mean * self.params["V_int"],
#                             rho_su=rho_su, rho_ex=rho_ex, rho_mean=rho_mean,
#                             choked=choked)
#         return self.outputs






# # --------------------------------------------------------------------------
# # Conversions: any datasheet format -> Kv
# # --------------------------------------------------------------------------
# def kv_from_zeta(zeta, A):
#     """Loss coefficient zeta (dp = zeta*rho*v^2/2) on cross-section A [m2] -> Kv [m3/h]."""
#     return np.sqrt(2.0 * C_KV) * A / np.sqrt(zeta)          # = 50912 * A / sqrt(zeta)


# def kv_from_point(dp_nom, m_dot_nom, rho_nom):
#     """One measured/datasheet point (dp [Pa], m_dot [kg/s], rho [kg/m3]) -> Kv [m3/h]."""
#     return m_dot_nom * np.sqrt(C_KV / (rho_nom * dp_nom))


# def kv_from_cv(Cv):
#     """US flow coefficient Cv [gpm, psi] -> Kv."""
#     return Cv / 1.156


# def kv_series(*kvs):
#     """Components in series: 1/Kv^2 = sum(1/Kv_i^2)."""
#     return 1.0 / np.sqrt(sum(1.0 / k**2 for k in kvs))




#     # ---------------------------------------------------------------
#     def solve(self):
#         fl = self.inputs["fluid"]
#         P_su, h_su, m_dot = self.inputs["P_su"], self.inputs["h_su"], self.inputs["m_dot"]
#         Kv = self.params["Kv"] * self.params["opening"]

#         rho_su = PropsSI("D", "P", P_su, "H", h_su, fl)   # homogeneous density if two-phase

#         if Kv <= 0.0:                    # closed valve
#             if abs(m_dot) > 0:
#                 raise ValueError("Closed valve (Kv=0) with non-zero mass flow.")
#             dp, choked = 0.0, False
#         else:
#             dp = C_KV * m_dot * abs(m_dot) / (rho_su * Kv**2)     # incompressible law
#             choked = False
#             if self.params["compressible"]:
#                 dp, choked = self._compressible_dp(fl, P_su, h_su, m_dot, rho_su, Kv)

#         P_ex = P_su - dp
#         h_ex = h_su                       # adiabatic, no work -> isenthalpic
#         rho_ex = PropsSI("D", "P", P_ex, "H", h_ex, fl)
#         rho_mean = 0.5 * (rho_su + rho_ex)

#         self.outputs = dict(dp=dp, P_ex=P_ex, h_ex=h_ex,
#                             M=rho_mean * self.params["V_int"],
#                             rho_su=rho_su, rho_ex=rho_ex, rho_mean=rho_mean,
#                             choked=choked)
#         return self.outputs

#     # ---------------------------------------------------------------
#     def _compressible_dp(self, fl, P_su, h_su, m_dot, rho_su, Kv):
#         """
#         IEC 60534 gas sizing: m_dot = (Kv/3600) * Y * sqrt(rho0 * rho_su * dp / dp0) with
#         Y = 1 - x / (3 F_g x_T), x = dp/P_su, F_g = gamma/1.4, choked for x >= F_g x_T.
#         Solved for dp by fixed-point iteration.
#         """
#         cp = PropsSI("CPMASS", "P", P_su, "H", h_su, fl)
#         cv = PropsSI("CVMASS", "P", P_su, "H", h_su, fl)
#         F_g = (cp / cv) / 1.4
#         x_ch = F_g * self.params["x_T"]

#         # maximum flow the component can pass (choked)
#         Y_ch = 2.0 / 3.0
#         m_max = (Kv / 3600.0) * Y_ch * np.sqrt(RHO0 * rho_su * x_ch * P_su / DP0)
#         if m_dot >= m_max:
#             return x_ch * P_su, True

#         dp = C_KV * m_dot**2 / (rho_su * Kv**2)
#         for _ in range(50):
#             Y = 1.0 - (dp / P_su) / (3.0 * x_ch)
#             dp_new = C_KV * m_dot**2 / (rho_su * Kv**2 * Y**2)
#             if abs(dp_new - dp) < 1e-6 * P_su:
#                 break
#             dp = dp_new
#         return dp_new, False

#     # Inverse use: imposed P_su and P_ex -> m_dot (useful in a network solver)
#     def m_dot_from_dp(self, fl, P_su, h_su, P_ex):
#         Kv = self.params["Kv"] * self.params["opening"]
#         rho_su = PropsSI("D", "P", P_su, "H", h_su, fl)
#         dp = P_su - P_ex
#         return np.sign(dp) * Kv * np.sqrt(abs(dp) * rho_su / C_KV)


# if __name__ == "__main__":
#     # EVR 25, Kv = 9.8 m3/h; internal volume to be measured/estimated (here 50 cm3 as a placeholder)
#     evr = KvComponent(Kv=9.8, V_int=50e-6)

#     fl = "R1233ZDE"
#     cases = [("liquid line",    10e5, PropsSI("T", "P", 10e5, "Q", 0, fl) - 5),   # 5 K subcooled
#              ("suction vapour",  2e5, PropsSI("T", "P", 2e5,  "Q", 1, fl) + 5)]   # 5 K superheated
#     for label, P, T in cases:
#         h = PropsSI("H", "P", P, "T", T, fl)
#         evr.set_inputs(fluid=fl, P_su=P, h_su=h, m_dot=0.3)
#         out = evr.solve()
#         print(f"{label:15s} rho={out['rho_su']:7.1f} kg/m3  dp={out['dp']/1e5:6.3f} bar  M={out['M']*1e3:6.2f} g")