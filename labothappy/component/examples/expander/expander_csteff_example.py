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

from labothappy.component.expander.expander_csteff import ExpanderCstEff

# Example usage
EXP = ExpanderCstEff()
EXP.print_setup()

# "If the inputs are not set directly BUT throught the connectors"
EXP.su.set_properties(P=955214.9, T=374.18, fluid='R134a', m_dot = 0.1)
EXP.ex.set_properties(P=293940.1)
EXP.set_parameters(eta_is=0.8)
EXP.print_setup()

EXP.solve()
EXP.print_results()

fig = EXP.plot_Ts()
fig.show()

