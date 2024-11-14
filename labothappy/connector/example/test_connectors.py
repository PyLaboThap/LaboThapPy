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
# Put three properteis at the same time
point = MassConnector()
point.set_properties(T=500, m_dot=0.5, fluid = 'INCOMP::DowQ', P=101325, H=100000)


