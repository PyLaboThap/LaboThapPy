
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch

# from simulation_model import HXPinchCst
import numpy as np

case_test = 'EVAP_C3'

"Evaporator test"

if case_test == 'EVAP_C5':

    # # # Exo ORC M&S
    EVAP = HexCstPinch()
    
    EVAP.set_inputs(
        fluid_C = 'Cyclopentane',
        T_su_C = 110+273.15,
        P_su_C = 831.8*1e3,
        m_dot_C = 51.03,
    
        fluid_H = 'Water', #Oil
        T_su_H = 145+273.15,
        P_su_H = 5*1e5,
        m_dot_H = 400,
    )
    
    EVAP.set_parameters(**{
        'Pinch': 4,
        'Delta_T_sh_sc': 10,
        'HX_type': 'evaporator'
    })
    
    EVAP.solve()
    EVAP.print_results()
    EVAP.print_states_connectors()
    EVAP.plot_disc()
    
    fig = EVAP.plot_Ts(choose_HX_side='C')
    fig.show()
    
    EVAP.equivalent_effectiveness()

elif case_test == 'EVAP_C3':

    # Exo ORC M&S
    EVAP = HexCstPinch()
    
    EVAP.set_inputs(
        fluid_C = 'Propane',
        T_su_C = 15+273.15,
        P_su_C = 793237,
        m_dot_C = 0.1,
    
        fluid_H = 'Water', #Oil
        T_su_H = 25+273.15,
        P_su_H = 3*1e5,
        m_dot_H = 2,
    )
    
    EVAP.set_parameters(**{
        'Pinch': 4,
        'Delta_T_sh_sc': 1,
        'HX_type': 'evaporator',
        'DP_c' : 10*1e3,
        'DP_h' : 10*1e3,
    })
    
    EVAP.solve()
    EVAP.print_results()
    EVAP.print_states_connectors()
    EVAP.plot_disc()
    
    # EVAP.equivalent_effectiveness()

elif case_test == 'COND_C5':

    "Condenser test"
    
    COND = HexCstPinch()
    
    COND.set_inputs(
        fluid_H = 'Cyclopentane',
        T_su_H = 41.2+273.15,
        P_su_H = 68.3*1e3,
        m_dot_H = 46.18,
        
        fluid_C = 'Air',
        T_su_C = 20+273.15,
        P_su_C = 1e5,
        m_dot_C = 1911
    )
    
    COND.set_parameters(**{
        'Pinch': 5,
        'Delta_T_sh_sc': 5,
        'HX_type': 'condenser'
    })
    
    COND = HexCstPinch()
    
    COND.set_inputs(
        fluid_H = 'CO2',
        T_su_H = 30+273.15,
        P_su_H = 55*1e5,
        m_dot_H = 30,
        
        fluid_C = 'Water',
        T_su_C = 15+273.15,
        P_su_C = 5*1e5,
        m_dot_C = 100
    )
    
    COND.set_parameters(**{
        'Pinch': 5,
        'Delta_T_sh_sc': 0.1,
        'HX_type': 'condenser',
        # 'DP_c' : 10*1e3,
        # 'DP_h' : 10*1e3,
    })
    
    
    COND.solve()
    
    COND.print_results()
    COND.print_states_connectors()
    COND.plot_disc()
    
    COND.equivalent_effectiveness()
    
