Mass Connector
==============

This section includes details on the *'MassConnector'* Class of the PyLaboThap library.

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

- **What you can do:**
  
  Create an instance of the *MassConnector* class:
  
  .. literalinclude:: ../../../../../library/connector/example/test_connectors.py
      :language: python
      :lines: 16-17

  The properties can be put in any order. The class will automatically calculate the remaining properties.

  Print the state of the MassConnector:
  
  .. literalinclude:: ../../../../../library/connector/example/test_connectors.py
      :language: python
      :lines: 20

  Reset the pressure:
  
  .. literalinclude:: ../../../../../library/connector/example/test_connectors.py
      :language: python
      :lines: 23

- **What you cannot do:**
  
  Set a third property:
  
  .. literalinclude:: ../../../../../library/connector/example/test_connectors.py
      :language: python
      :lines: 27-28

