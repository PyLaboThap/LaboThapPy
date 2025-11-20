# import __init__

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
