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

from labothappy.component.compressor.compressor_csteff import CompressorCstEff

# Example usage
CP = CompressorCstEff()
# CP.print_setup()

"If the inputs are not set directly BUT through the connectors"
# CP.su.set_properties(P=319296.5575177148, T=331.033964665788, fluid='R1233ZDE', m_dot = 0.1)
# CP.ex.set_properties(P=606240.1433176235)

"If the inputs are set directly"
CP.set_inputs(
    P_su=319296.5575177148,
    T_su=331.033964665788,
    P_ex=606240.1433176235,
    fluid='R1233ZDE',  # Make sure to include fluid information
    m_dot=0.1  # Mass flow rate
)
CP.set_parameters(eta_is=0.8)
# CP.print_setup()

CP.solve()
CP.print_results()

fig = CP.plot_Ts()
fig.show()
