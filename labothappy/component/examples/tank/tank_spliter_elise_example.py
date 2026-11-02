# -*- coding: utf-8 -*-
"""
Examples of the TankSpliter model with pressure drop and charge (spliter_elise.py)

Cases:
    1) No geometry           : isobaric split (original behaviour)
    2) Liquid water          : symmetric split, 2 outlets
    3) Liquid water          : uneven split and different outlet tubes
    4) Two-phase R134a       : refrigerant distributor after an expansion valve, 4 outlets
    5) Charge from V_int only

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

from labothappy.component.tank.spliter_elise import TankSpliter


def print_results(spliter, title):
    print(f"\n=== {title} ===")
    for i in range(len(spliter.outlet_repartition)):
        outlet = getattr(spliter, f"ex_{i+1}")
        print(f"Outlet {i+1}: m_dot = {outlet.m_dot:.4f} kg/s, p = {outlet.p:.0f} Pa, h = {outlet.h:.0f} J/kg")
    if spliter.DP is not None:
        d = spliter.DP_detail
        print(f"Two-phase : {d['two_phase']},  Re_su = {d['Re_su']:.0f},  Re_ex = {d['Re_ex'].round(0)}")
        print(f"Transition: K = {d['K_transition']:.3f}, theta = {d['theta_deg']:.1f} deg, area ratio = {d['area_ratio']:.3f}")
        print(f"DP inlet body   : {d['inlet_body']:.1f} Pa")
        print(f"DP transition   : {d['transition'].round(1)} Pa")
        print(f"DP outlet tubes : {d['outlet_tubes'].round(1)} Pa")
        print(f"DP total        : {spliter.DP.round(1)} Pa")
    if spliter.charge is not None:
        print(f"Charge : {spliter.charge*1e3:.2f} g  ({spliter.charge_detail})")


# 1) No geometry: isobaric split ------------------------------------------------------------------

spliter = TankSpliter(outlet_repartition=[0.3, 0.4, 0.3])
spliter.set_inputs(T_su=10 + 273.15, m_dot=13.8, P_su=0.8e5, fluid='Cyclopentane')
spliter.solve()
print_results(spliter, "1) No geometry (isobaric)")

# 2) Liquid water, symmetric split ----------------------------------------------------------------

spliter = TankSpliter(outlet_repartition=[0.5, 0.5])
spliter.set_inputs(T_su=20 + 273.15, m_dot=0.5, P_su=3e5, fluid='Water')
spliter.set_parameters(D_su=0.025, L_su=0.10,     # Inlet body: DN25, 10 cm
                       D_ex=0.016, L_ex=0.20,     # Outlet tubes: 16 mm, 20 cm
                       L_trans=0.03)              # Transition length [m]
spliter.solve()
print_results(spliter, "2) Water, symmetric split")

# 3) Liquid water, uneven split and different outlet tubes ----------------------------------------

spliter = TankSpliter(outlet_repartition=[0.2, 0.3, 0.5])
spliter.set_inputs(T_su=20 + 273.15, m_dot=0.5, P_su=3e5, fluid='Water')
spliter.set_parameters(D_su=0.025, L_su=0.10,
                       D_ex=[0.010, 0.012, 0.016], L_ex=[0.30, 0.20, 0.20],
                       L_trans=0.03)
spliter.solve()
print_results(spliter, "3) Water, uneven split")

# 4) Two-phase R134a distributor (after an expansion valve) ---------------------------------------

P_evap = PropsSI('P', 'T', 0 + 273.15, 'Q', 0.5, 'R134a')
h_su = PropsSI('H', 'P', P_evap, 'Q', 0.25, 'R134a')

spliter = TankSpliter(outlet_repartition=[0.25, 0.25, 0.25, 0.25])
spliter.set_inputs(h_su=h_su, m_dot=0.05, P_su=P_evap, fluid='R134a')
spliter.set_parameters(D_su=0.010, L_su=0.03,      # Distributor body
                       D_ex=0.004, L_ex=0.30,      # Distributor tubes (capillary-like)
                       L_trans=0.015)
spliter.solve()
print_results(spliter, "4) R134a two-phase distributor")

# 5) Charge from a total internal volume only ------------------------------------------------------

spliter = TankSpliter(outlet_repartition=[0.5, 0.5])
spliter.set_inputs(h_su=h_su, m_dot=0.05, P_su=P_evap, fluid='R134a')
spliter.set_parameters(V_int=20e-6)                # 20 cm^3
spliter.solve()
print_results(spliter, "5) R134a, charge from V_int")
