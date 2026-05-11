from labothappy.component.expander.expander_csteff import ExpanderCstEff

# # Example usage
# EXP = ExpanderCstEff()
# EXP.print_setup()

# # "If the inputs are not set directly BUT throught the connectors"
# EXP.su.set_properties(P=955214.9, T=374.18, fluid='R134a', m_dot = 0.1)
# EXP.ex.set_properties(P=293940.1)
# EXP.set_parameters(eta_is=0.8)
# EXP.print_setup()

# EXP.solve()
# EXP.print_results()

# fig = EXP.plot_Ts()
# fig.show()


# Example usage
EXP = ExpanderCstEff()
EXP.print_setup()

# "If the inputs are not set directly BUT throught the connectors"
EXP.su.set_properties(P=14258720.709528737, h=522325.64157148026, fluid='CO2', m_dot = 327.9722712943479)
EXP.ex.set_properties(P=4503029.698839702)
EXP.set_parameters(eta_is=0.892)
EXP.print_setup()

EXP.solve()
EXP.print_results()

fig = EXP.plot_Ts()
fig.show()


