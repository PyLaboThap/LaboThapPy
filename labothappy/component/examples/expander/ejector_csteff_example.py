
from labothappy.component.expander.ejector_csteff import EjectorCstEff

if __name__ == "__main__":

    ej = EjectorCstEff()
    
    P_su_1 = 100
    P_su_2 = 39
    
    ej.set_inputs(
        T_su_1 = 40 + 273.15, # K
        P_su_1 = P_su_1*1e5, # Pa
        m_dot_su_1 = 0.18, # kg/s
        fluid_su_1 = 'CO2',
        
        T_su_2 = 15 + 273.15, # K
        P_su_2 = P_su_2*1e5, # Pa
        m_dot_su_2 = 0.07, # kg/s
        fluid_su_2 = 'CO2',
        )
    
    ej.set_parameters(
        eta_m = 0.9,
        eta_s = 0.8,
        eta_mix = 0.8,
        eta_d =  0.85,
        C_t = 0.7
        )

    ej.solve()
    