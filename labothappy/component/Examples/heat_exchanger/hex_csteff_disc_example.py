import __init__

from component.heat_exchanger.hex_csteff_disc import HXEffCstDisc

import numpy as np

"Simple test - CO2 HTHP GasCooler"

#Exo ORC M&S
HTX = HXEffCstDisc()

HTX.set_inputs(
    Csu_fluid = 'Water',
    Csu_T = 273.15 + 15,
    Csu_m_dot = 0.1,
    Csu_p = 10e5,

    Hsu_fluid = 'CO2',
    Hsu_T = 450,
    Hsu_m_dot = 0.16,
    Hsu_p = 140*1e5,
)

# HTX.set_inputs(
#     Csu_fluid = 'CO2',
#     Csu_T = 270.15,
#     Csu_m_dot = 0.16,
#     Csu_p = 2963161,

#     Hsu_fluid = 'CO2',
#     Hsu_T = 314.75,
#     Hsu_m_dot = 0.16,
#     Hsu_p = 120*1e5,
# )

HTX.set_parameters(**{
    'eta': 0.9, 'n_disc' : 50, 'Pinch_min' : 5
})

HTX.solve()

HTX.plot_disc()
