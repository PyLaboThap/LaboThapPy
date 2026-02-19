from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.valve.valve_isenthalpic import ValveIsenthalpic

from CoolProp.CoolProp import PropsSI

# Instanciate Circuit
fluid = "Propane"
HP = RecursiveCircuit(fluid)

# Create components
Compressor = CompressorCstEff()
Condenser = HexCstPinch()
ExpansionValve = ValveIsenthalpic()
Evaporator = HexCstPinch()

# Set component parameters
eta_is_cp = 0.8
Pinch_cd = 3  # K
SC_cd = 3  # K
Pinch_ev = 3  # K
SH_ev = 3  # K
Compressor.set_parameters(eta_is=0.8)
Condenser.set_parameters(Pinch=Pinch_cd, Delta_T_sh_sc=SC_cd, HX_type="condenser")
Evaporator.set_parameters(Pinch=Pinch_ev, Delta_T_sh_sc=SH_ev, HX_type="evaporator")

# Add components to circuit
HP.add_component(Compressor, "Compressor")
HP.add_component(Condenser, "Condenser")
HP.add_component(ExpansionValve, "ExpansionValve")
HP.add_component(Evaporator, "Evaporator")

# Link components with mass connectors
HP.link_components("Compressor", "m-ex", "Condenser", "m-su_H")
HP.link_components("Condenser", "m-ex_H", "ExpansionValve", "m-su")
HP.link_components("ExpansionValve", "m-ex", "Evaporator", "m-su_C")
HP.link_components("Evaporator", "m-ex_C", "Compressor", "m-su")

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

# Inputs
T_HS = 40 + 273.15  # K
P_HS = 3e5  # Pa
m_dot_HS = 0.5  # kg/s
T_CS = 20 + 273.15  # K
P_CS = 1e5  # Pa
m_dot_CS = 5 # kg/s
m_dot_ref = 0.2 # kg/s

# Cycle guess values
P_low = PropsSI("P", "T", T_CS-10, "Q", 1, fluid)
P_high = PropsSI("P", "T", T_HS+10, "Q", 0, fluid)

HP.set_cycle_guess(target="Compressor:su", m_dot = m_dot_ref, SH=SH_ev, p=P_low)
HP.set_cycle_guess(target="Compressor:ex", p=P_high)

HP.set_cycle_guess(target="ExpansionValve:su", p=P_high, m_dot = m_dot_ref, SC=SC_cd)
HP.set_cycle_guess(target="ExpansionValve:ex", p=P_low)
# Set residual variables
HP.set_residual_variable(target="Compressor:ex", variable="h", tolerance=1e3)
HP.set_residual_variable(target="Compressor:ex", variable="p", tolerance=1e3)

HP.set_residual_variable(target="Condenser:ex_H", variable="h", tolerance=1e-3)
HP.set_residual_variable(target="Condenser:ex_H", variable="p", tolerance=1e-3)

HP.set_residual_variable(target="Evaporator:ex_C", variable="h", tolerance=1e-3)
HP.set_residual_variable(target="Evaporator:ex_C", variable="p", tolerance=1e-3)

HP.set_residual_variable(target="ExpansionValve:ex", variable="h", tolerance=1e3)
HP.set_residual_variable(target="ExpansionValve:ex", variable="p", tolerance=1e3)

HP.solve()
print(f"Converged at P_HP = {Compressor.ex.p}, P_LP = {Compressor.su.p}")
HP.print_states()

