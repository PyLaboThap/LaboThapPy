
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

# Second case
PP.su.set_properties(P=400000, fluid='R1233ZDE')
PP.su.set_m_dot(1.0)  
PP.ex.set_properties(P=800000, fluid='R1233ZDE')
PP.solve()
PP.print_results()

# Third case (set T via inputs)
PP.set_inputs(su_T=400)
PP.su.set_m_dot(1.0)  
PP.solve()
PP.print_states_connectors()


