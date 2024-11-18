"""
Author: Elise Neven
email: elise.neven@uliege.be
Date: 18/11/2024
"""

from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector


class PumpSE(BaseComponent):
    """
    Component: Pump

    Model: Based on the Thesis of Rémi Dickes. Semi-empirical model.

    **Description**:
        This model is used to simulate the performances and the behavior of pumps.
        The parameters need to be calibrated with experimental data.

    **Assumption**:
        Steady-state operation.

    **Connectors**:
        su (MassConnector): Mass connector for the suction side.

        ex (MassConnector): Mass connector for the exhaust side.

        W_pp (WorkConnector): Work connector to connect to a potential motor model.

    **Parameters**:
     


    
    """
    def __init__(self):
