from labothappy.machine.circuit_rec import RecursiveCircuit
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.compressor.compressor_csteff import CompressorCstEff
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff
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
Recuperator = HexCstEff()

# Set component parameters
eta_is_cp = 0.8
Pinch_cd = 3  # K
SC_cd = 3  # K
Pinch_ev = 3  # K
SH_ev = 3  # K
eta_rec = 0.8
Compressor.set_parameters(eta_is=0.8)
Condenser.set_parameters(Pinch=Pinch_cd, Delta_T_sh_sc=SC_cd, HX_type="condenser")
Evaporator.set_parameters(Pinch=Pinch_ev, Delta_T_sh_sc=SH_ev, HX_type="evaporator")
Recuperator.set_parameters(eta=eta_rec)


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


# Inputs
T_HS = 40 + 273.15  # K
P_HS = 3e5  # Pa
m_dot_HS = 0.5  # kg/s
T_CS = 20 + 273.15  # K
P_CS = 1e5  # Pa
m_dot_CS = 5 # kg/s
m_dot_ref = 0.2 # kg/s

# Add fluid sources
CD_source = MassConnector('Water')
EV_source = MassConnector('Water')

HP.add_source("CD_Water", CD_source, HP.components["Condenser"], "m-su_C")
HP.set_source_properties(T=T_HS, P=P_HS, m_dot=m_dot_HS, fluid='Water', target="CD_Water")
HP.add_source("EV_Water", EV_source, HP.components["Evaporator"], "m-su_H")
HP.set_source_properties(T=T_CS, P=P_CS, m_dot=m_dot_CS, fluid='Water', target="EV_Water")

# Evaporator.print_setup()
# Cycle guess values
P_low = PropsSI("P", "T", T_CS-5, "Q", 1, fluid)
P_high = PropsSI("P", "T", T_HS+5, "Q", 0, fluid)

HP.set_cycle_guess(target="Compressor:su", m_dot = m_dot_ref, SH=SH_ev, p=P_low)
HP.set_cycle_guess(target="Compressor:ex", p=P_high)

HP.set_cycle_guess(target="Recuperator:su_C", p=P_low, m_dot=m_dot_ref, T=T_CS+5)

HP.set_cycle_guess(target="ExpansionValve:su", p=P_high)
HP.set_cycle_guess(target="ExpansionValve:ex", p=P_low)
# Compressor.print_setup()
# Set residual variables
HP.set_residual_variable(target="Compressor:ex", variable="h", tolerance=1e3)
HP.set_residual_variable(target="Compressor:ex", variable="p", tolerance=1e3)

HP.set_residual_variable(target="Condenser:ex_H", variable="h", tolerance=1e-3)
HP.set_residual_variable(target="Condenser:ex_H", variable="p", tolerance=1e-3)

HP.set_residual_variable(target="Evaporator:ex_C", variable="h", tolerance=1e-3)
HP.set_residual_variable(target="Condenser:ex_C", variable="p", tolerance=1e-3)

HP.set_residual_variable(target="ExpansionValve:ex", variable="h", tolerance=1e3)
HP.set_residual_variable(target="ExpansionValve:ex", variable="p", tolerance=1e3)


HP.solve()


