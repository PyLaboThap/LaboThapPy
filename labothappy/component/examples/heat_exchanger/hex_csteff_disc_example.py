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

from labothappy.component.heat_exchanger.hex_csteff_disc import HexCstEffDisc


"Simple test - CO2 HTHP GasCooler"

#Exo ORC M&S
HTX = HexCstEffDisc()

HTX.set_inputs(
    fluid_C = 'Water',
    T_su_C = 273.15 + 15,
    m_dot_C = 0.1,
    P_su_C = 10e5,

    fluid_H = 'CO2',
    T_su_H = 450,
    m_dot_H = 0.16,
    P_su_H = 140*1e5,
)

# HTX.set_inputs(
#     fluid_C = 'CO2',
#     T_su_C = 270.15,
#     m_dot_C = 0.16,
#     P_su_C = 2963161,

#     fluid_H = 'CO2',
#     T_su_H = 314.75,
#     m_dot_H = 0.16,
#     P_su_H = 120*1e5,
# )

HTX.set_parameters(**{
    'eta_max' : 0.95,
    'n_disc' : 100, 
    'Pinch_min' : 10,
    'DP_c' : 50*1e3,
    'DP_h' : 50*1e3,    
})

HTX.solve()
HTX.plot_disc()

fig = HTX.plot_Ts(choose_HX_side='H')
fig.show()
