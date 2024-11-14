from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import numpy as np

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from machine.circuit import Circuit
from machine.boundary_conditions.mass_source import MassSource
from machine.boundary_conditions.mass_sink import MassSink

from labothappy.component.heat_exchanger.epsilon_NTU.simulation_model import HXeNTU
from component.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from component.pump.constant_efficiency.simulation_model import PumpCstEff

class System:
    def __init__(self):
        self.cycle = Circuit()
        self.source = MassSource()
        self.sink = MassSink()


if __name__ == "__main__":
    # Create a cycle
    # POC = System()
    # POC.cycle.components["Pump"] = PumpCstEff()
    # POC.cycle.components["Evaporator"] = HXeNTU()
    # POC.cycle.components["Condenser"] = HXeNTU()
    # POC.cycle.components["Expander"] = ExpanderCstEff()

    # # Create a source
    # Tank = Source()
    
    # print(POC.cycle.components)

    # # Create a sink
    # River = Sink()
    # Example usage in the main script
    # Create a cycle
    ORC = Circuit()
    
    # Create components
    Pump = PumpCstEff()
    Evaporator = HXeNTU()
    Condenser = HXeNTU()
    Expander = ExpanderCstEff()

    # Heat transfer area  
    A_htx = 0.752 # m^2
    
    # HTX dimensions
    L = 0.393 # [m] : length
    w = 0.243 # [m] : width
    h = 0.0446 # [m] : height
    l_v = 0.324 # [m] : length between ports
    casing_t = 0.005 # [m] : casing thickness # !!! arbitrary
    
    # Number and thickness of plates
    n_plates = 10
    t_plates = 0.0008 # [m] # !!! arbitrary
    
    # Fooling factor
    fooling = 98.51/1000 # (m^2*K)/W

    # Number of canals
    C_n_canals = 5
    H_n_canals = 4
    
    # Plate values 
    plate_cond = 45 # [W/(m*K)] : plate conduction
    plate_pitch_co = 0.005 # 0.00745870973 # corrugated pitch # !!! arbitrary
    chevron_angle = 20*np.pi/180 # !!! arbitrary

    # Total volume of each part
    V_tot = 0.9*1e-3 # [m^3]

    # Canal thickness
    C_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*C_n_canals)
    H_canal_t = ((h-2*casing_t) - n_plates*t_plates)/(2*H_n_canals)

    # Total Canal Surface
    C_CS = C_canal_t*(w-2*casing_t)*C_n_canals
    H_CS = H_canal_t*(w-2*casing_t)*H_n_canals

    # Dh : hydraulic diameter
    C_Dh = (4*C_canal_t*w)/(2*C_canal_t+2*w)
    H_Dh = (4*H_canal_t*w)/(2*H_canal_t+2*w)

    # Set parameters for components
    Pump.set_parameters(eta_is=0.8)
    Evaporator.set_parameters(
        A_htx=0.752, L_HTX=0.393, V_HTX=0.9*1e-3, Flow_Type = 'CounterFlow',
        A_canal_h=H_CS, A_canal_c=C_CS, D_h=H_Dh, 
        k_plate=45, t_plate=0.0008, n_plates = 10,
        co_pitch=0.005, chevron_angle=20*np.pi/180, fouling=98.51/1000
    )
    Condenser.set_parameters(NTU=1.5, Cmin=1000, hot_fluid='R1233ZD', cold_fluid='R1233ZD')
    Expander.set_parameters(eta_is=0.8)

    # Add components to the cycle
    ORC.add_component(Pump, "Pump") # :!\ Est-ce que les noms des composants sont importants?
    ORC.add_component(Evaporator, "Evaporator")
    ORC.add_component(Condenser, "Condenser")
    ORC.add_component(Expander, "Expander")
    
    # Link components
    ORC.link_components("Pump", "ex", "Evaporator", "su_cold")
    # ORC.link_components("Evaporator", "ex", "Expander", "su")
    # ORC.link_components("Expander", "ex", "Condenser", "su")
    # ORC.link_components("Condenser", "ex", "Pump", "su")