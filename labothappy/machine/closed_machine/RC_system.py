
import sys

from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
from labothappy.connector.heat_connector import HeatConnector

from component.heat_exchanger.steady_state.pinch_cst.simulation_model import HXPinchCst
from component.volumetric_machine.expander.steady_state.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.pump.steady_state.constant_efficiency.simulation_model import PumpCstEff

from component.heat_exchanger.steady_state.cst_efficiency.simulation_model import HXEffCst

from machine.circuit import Circuit

from machine.boundary_conditions.mass_sink import MassSink
from machine.boundary_conditions.mass_source import MassSource

from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize
from scipy.optimize import fsolve

from scipy.optimize import minimize
import numpy as np


# Example usage of the Rankine Cycle (RC)
if __name__ == "__main__":
    orc_cycle = Circuit(fluid='R1233zd(E)')

    # Add components
    PUMP = PumpCstEff()
    EVAP = HXPinchCst()
    EXP = ExpanderCstEff()
    COND = HXPinchCst()

    # Set component parameters
    PUMP.set_parameters(eta_is=0.6)
    EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
    COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
    EXP.set_parameters(eta_is=0.8)

    # Add components to the cycles
    orc_cycle.add_component(PUMP, "Pump")
    orc_cycle.add_component(EVAP, "Evaporator")
    orc_cycle.add_component(EXP, "Expander")
    orc_cycle.add_component(COND, "Condenser")

    # Link components
    orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
    orc_cycle.link_components("Evaporator", "m-ex_C", "Expander", "m-su")
    orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_H")
    orc_cycle.link_components("Condenser", "m-ex_H", "Pump", "m-su")

    # Set the cycle properties
    orc_cycle.set_cycle_properties(m_dot=0.4, target='Pump:su')

    orc_cycle.set_cycle_properties(T=15 + 273.15, fluid='Water', m_dot=2, target='Condenser:su_C', p = 4e5)
    orc_cycle.set_cycle_properties(cp=4186, target='Condenser:su_C')
    orc_cycle.set_cycle_properties(fluid='Water', target='Condenser:ex_C')

    orc_cycle.set_cycle_properties(T=90 + 273.15, fluid='Water', m_dot=1, target='Evaporator:su_H', p = 4e5)
    orc_cycle.set_cycle_properties(cp=4186, target='Evaporator:su_H')
    orc_cycle.set_cycle_properties(fluid='Water', target='Evaporator:ex_H')

    # Set parameters for the cycle
    SC_cd = 5
    SH_ev = 5
    orc_cycle.set_cycle_parameters(SC_cd=SC_cd, SH_ev=SH_ev)

    # Initial guesses for pressures
    T_ev_guess = 90 + 273.15
    P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R1233zd(E)')

    T_cd_guess = 15 + 273.15
    P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R1233zd(E)')

    # Define guesses and residuals
    guesses = {
        "Evaporator:ex_C-p": P_ev_guess,
        "Condenser:ex_H-p": P_cd_guess,
        "Pump:su-T": T_cd_guess - SC_cd,
        "Pump:ex-p": P_ev_guess,
        "Expander:ex-p": P_cd_guess, 
        "Evaporator:ex_C-T": T_ev_guess + SH_ev
    }

    residuals_var = [
        "Evaporator:ex_C-h",
        "Condenser:ex_H-h"
    ]

    orc_cycle.set_cycle_guesses_residuals(guesses, residuals_var)

    # Solve the cycle
    orc_cycle.solve(start_key="Pump")

    cold_water_loop = Circuit(fluid='Water')

    # Add components to the cycle
    cold_water_loop.add_component(COND, "Condenser")

    # Add sink and source components
    mass_sink = MassSink()
    mass_source = MassSource(properties = {'T': 10 + 273.15, 'fluid': 'Water', 'p': 1e5, 'm_dot': 2})

    # Add components to the cycle
    cold_water_loop.add_component(mass_source, "ColdSource")
    cold_water_loop.add_component(mass_sink, "ColdSink")


    # Link components
    cold_water_loop.link_components("ColdSource", "m-ex", "Condenser", "m-su_C")
    cold_water_loop.link_components("Condenser", "m-ex_C", "ColdSink", "m-su")

    # Set the cycle properties
    # cold_water_loop.set_cycle_properties(m_dot=2, target='Pump:su')
    # cold_water_loop.set_cycle_properties(T=20 + 273.15, fluid='Water', target='Dry Cooler:su_C', P = 1e5)

    # FAIRE AVEC MASS SINK ET MASS SOURCE
