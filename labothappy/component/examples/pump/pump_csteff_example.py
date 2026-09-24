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

from labothappy.component.pump.pump_csteff import PumpCstEff
import numpy as np

# Example usage
PP = PumpCstEff()

# Set initial conditions
PP.su.set_properties(P=319296.56, T=331.03, fluid='R1233ZDE')
PP.su.set_m_dot(1.0)  
PP.ex.set_properties(P=606240.14, fluid='R1233ZDE')
PP.set_parameters(eta_is=0.9)
PP.solve()
PP.print_results()

