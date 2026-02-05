from machine.circuit_it import IterativeCircuit
from connector.mass_connector import MassConnector
from component.compressor.compressor_csteff import CompressorCstEff
from component.heat_exchanger.hex_cstpinch import HexCstPinch
from component.valve.valve_isenthalpic import ValveIsenthalpic
from component.heat_exchanger.hex_csteff import HexCstEff

from CoolProp.CoolProp import PropsSI

# Instanciate Circuit
fluid = "Propane"
HP = IterativeCircuit(fluid)

# Create components
Compressor = CompressorCstEff()
Condenser = HexCstPinch()
ExpansionValve = ValveIsenthalpic()
Evaporator = HexCstPinch()
Recuperator = HexCstEff()

# Set component parameters
eta_is_cp = 0.8
Pinch_cd = 3  # K
SC_cd = 3  # K
Pinch_ev = 3  # K
SH_ev = 3  # K
eff_rec = 0.8 # -
Compressor.set_parameters(eta_is=0.8)
Condenser.set_parameters(Pinch=Pinch_cd, Delta_T_sh_sc=SC_cd, HX_type="condenser")
Evaporator.set_parameters(Pinch=Pinch_ev, Delta_T_sh_sc=SH_ev, HX_type="evaporator")
Recuperator.set_parameters(eta=eff_rec)

# Add components to circuit
HP.add_component(Compressor, "Compressor")
HP.add_component(Condenser, "Condenser")
HP.add_component(ExpansionValve, "ExpansionValve")
HP.add_component(Evaporator, "Evaporator")
HP.add_component(Recuperator, "Recuperator")

# Link components with mass connectors
HP.link_components("Compressor", "m-ex", "Condenser", "m-su_H")
HP.link_components("Condenser", "m-ex_H", "Recuperator", "m-su_H")
HP.link_components("Recuperator", "m-ex_H", "ExpansionValve", "m-su")
HP.link_components("ExpansionValve", "m-ex", "Evaporator", "m-su_C")
HP.link_components("Evaporator", "m-ex_C", "Recuperator", "m-su_C")
HP.link_components("Recuperator", "m-ex_C", "Compressor", "m-su")


# Add fluid sources
CD_source = MassConnector('Water')
T_su_w_cd = 40+273.15
P_su_w_cd = 3e5
m_dot_w_cd = 0.5  # kg/s
EV_source = MassConnector('Water')
T_su_w_ev = 20+273.15
P_su_w_ev = 1e5
m_dot_w_ev = 5  # kg/s

HP.add_source("CD_Water", CD_source, HP.components["Condenser"], "m-su_C")
HP.set_source_properties(T=T_su_w_cd, fluid='Water', P=P_su_w_cd, m_dot = m_dot_w_cd, target="CD_Water")
HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
HP.set_source_properties(T=T_su_w_ev, fluid='Water', P=P_su_w_ev, m_dot = m_dot_w_ev, target="EV_Water")

# Set cycle inputs
m_dot_ref = 0.2 # kg/s
SC_cd = 3  # K
SH_ev = 3  # K

HP.set_cycle_input(target="ExpansionValve:su", m_dot = m_dot_ref, SC=SC_cd)
HP.set_cycle_input(target="Evaporator:ex_C", SH=SH_ev)
HP.set_cycle_input(target="Condenser:ex_H", m_dot = m_dot_ref, SC=SC_cd)

# Set iteration variables
P_HP_guess = PropsSI("P", "T", T_su_w_cd+10, "Q", 0, fluid)
P_LP_guess = PropsSI("P", "T", T_su_w_ev-10, "Q", 1, fluid)
T_su_vlv_guess = PropsSI('T', 'P', P_HP_guess, 'Q', 0, fluid) - SC_cd

HP.set_iteration_variable(
    target=["ExpansionValve:su", "Compressor:ex", "Condenser:ex_H"],
    variable="p",
    guess=P_HP_guess,
    tolerance=1e-6
)

HP.set_iteration_variable(
    target="ExpansionValve:ex",
    variable="p",
    guess=P_LP_guess,
    tolerance=1e-6
)

HP.set_iteration_variable(
    target="ExpansionValve:su",
    variable="T",
    guess=T_su_vlv_guess,
    tolerance=1e-6
)

# Set residual variables
HP.set_residual_variable(
    pre_target="ExpansionValve:su",
    post_target="Condenser:ex_H",
    variable="p",
    tolerance=1e-3
)

HP.set_residual_variable(
    pre_target="ExpansionValve:ex",
    post_target="Evaporator:su_C",
    variable="p",
    tolerance=1e-3
)

HP.set_residual_variable(
    pre_target="ExpansionValve:su",
    post_target="ExpansionValve:su",
    variable="h",
    tolerance=1e-3
)

# Solve circuit
HP.solve()
print(f"Converged at P_HP = {Compressor.ex.p}, P_LP = {Compressor.su.p}")
