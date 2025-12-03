
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
