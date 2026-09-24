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

"""
Created on Fri June 19 2024

@author: Elise Neven
@email: elise.neven@uliege.be

"""

from CoolProp.CoolProp import PropsSI
from connector.mass_connector import MassConnector


#------------------Example of a mass connector------------------#
"What you can do with a mass connector:"
# Create an instance of the MassConnector class
point = MassConnector()
point.set_properties(T=500, m_dot=0.5, fluid = 'INCOMP::DowQ', P=101325)

# Print the state of the connector
point.print_resume(unit_T='C', unit_p='bar')

# Reset the pressure to 200000 Pa
point.set_p(200000)

"What you cannot do with a mass connector:"
# Put three properties at the same time
point = MassConnector()
point.set_properties(T=500, m_dot=0.5, fluid = 'INCOMP::DowQ', P=101325, H=100000)


