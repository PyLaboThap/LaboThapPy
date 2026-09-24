# LaboThApPy: An open-source tool for advanced thermodynamic cycle simulation, designed for researchers and engineers.
#
# Copyright (C) <2025> <Université catholique de Louvain (UCLouvain), Belgique
#                       Université de Liège (ULiège), Belgique
#                       Université de Mons (UMONS), Belgique>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# List of the contributors to the development of LaboThApPy: see AUTHORS file.
# Description and complete License: see NOTICE & LICENSE files.

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

