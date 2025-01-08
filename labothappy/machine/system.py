from CoolProp.CoolProp import PropsSI

import matplotlib.pyplot as plt
import numpy as np

from labothappy.connector.mass_connector import MassConnector
from labothappy.connector.work_connector import WorkConnector
from labothappy.connector.heat_connector import HeatConnector

from labothappy.machine.circuit import Circuit
from machine.boundary_conditions.mass_source import MassSource
from machine.boundary_conditions.mass_sink import MassSink

from labothappy.component.heat_exchanger.steady_state.epsilon_NTU.simulation_model import HXeNTU
from labothappy.component.volumetric_machine.expander.steady_state.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from labothappy.component.pump.steady_state.constant_efficiency.simulation_model import PumpCstEff

class System:
    def __init__(self):
        self.cycle = Circuit()
        self.source = MassSource()
        self.sink = MassSink()

    def solve(self):
        return

