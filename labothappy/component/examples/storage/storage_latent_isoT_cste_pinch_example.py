
import __init__
from component.storage.storage_latent_isoT_cste_pinch import StorageLatentIsothermalCstePinch

# from simulation_model import HXPinchCst
import numpy as np

"Ice storage case"

# # Exo ORC M&S
ICE_Storage = StorageLatentIsothermalCstePinch()

ICE_Storage.set_inputs(
    fluid = 'CO2',
    T_su = -15+273.15,
    P_su = 50*1e5,
    m_dot = 10,
    sto_fluid = 'Water'
)

ICE_Storage.set_parameters(**{
    'Pinch': 3,
    'Delta_T_sh_sc': 5
    })

ICE_Storage.solve()
ICE_Storage.print_results()
ICE_Storage.print_states_connectors()
ICE_Storage.plot_disc()

