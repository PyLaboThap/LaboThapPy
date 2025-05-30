Mass Connector
==============

The *MassConnector* class links components by transferring mass flow and thermodynamic properties between them. 
If two thermodynamic properties are provided, the class automatically calculates the remaining properties. 

The class has two key flags:

- **state_known**: This boolean flag indicates whether all thermodynamic properties are known.
- **state_completely_known**: This boolean flag indicates whether both the mass flow rate and all thermodynamic properties are known.

.. image:: ../../../../figures/connector/MassConnector.png
   :alt: Connectors description.
   :width: 100%

.. autoclass:: connector.mass_connector.MassConnector

Example of usage
----------------

You can see an interactive example in the following notebook:

.. toctree::
   :maxdepth: 1

   ../../../notebooks/mass_connector_demo


..    :maxdepth: 2
..    :caption: Example of usage:

..    ../../../notebooks/mass_connector_demo

