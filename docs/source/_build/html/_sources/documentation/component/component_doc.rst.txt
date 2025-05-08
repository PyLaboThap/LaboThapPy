Components
==========

The components are the building blocks of the library. The LaboThapPy library provides a variety of component models that can be used to simulate different thermodynamic cycles by interconnecting the models using connectors.
These models are designed to represent various types of components commonly found in thermodynamic systems, such as compressors, expanders, heat exchangers, pumps, solar collectors, tanks, valves,...
The documentation for each component model is organized into separate sections, making it easy to find the information about a certain model..

Each component model is implemented as a class that inherits from the `BaseComponent` class.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   compressor/compressor_doc
   expander/expander_doc
   heat_exchanger/ha_heat_exchanger_doc
   heat_exchanger/heat_exchanger_doc
   pump/pump_doc
   solar/solar_doc
   tank/tank_doc
   valve/valve_doc


BaseComponent Class
^^^^^^^^^^^^^^^^^^^^
The `BaseComponent` class is the parent class for all components in the PyLaboThap library. 
All components inherit the properties and methods of the `BaseComponent` class.

.. autoclass:: component.base_component.BaseComponent

.. Example of a Steady-State Component Model


.. The following example demonstrates how to create a steady-state model using the PyLaboThap library nomenclature. 
.. The different methods to implement are described below.

.. This example is based on the semif-empirical model of a volumetric compressor (XXX). For more details on this model, 
.. refer to the documentation [here](lien).

.. Inputs for the model can be defined using either the connector approach or the input/output approach. 
.. The figures below illustrate these two approaches for the compressor model:

.. .. raw:: html

..    <div class="side-by-side">

.. .. image:: ../../../../../figures/component/compressor_connectors.png
..    :alt: Connectors approach for a compressor model.
..    :width: 100%

.. .. image:: ../../../../../figures/component/compressor_in_out.png
..    :alt: Input/Output approach for a compressor model.
..    :width: 100%

.. .. raw:: html

..    </div>

.. For each model, the following methods need to be implemented:

.. 1. **Creating the class**:
.. Name the model accoring to PyLaboThap's conventions. In this example, the model is named 'CompressorSE' (SE stands for Semi-Empirical). The model inherits from `BaseComponent` its methods and attributes.
.. In the '__init__(self)' method, define the connectors for the component. *Note: super().__init__() is used to call the parent class constructor*.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 21-27

.. 2. **Defining the required inputs**: 
.. In the 'get_required_inputs' method, specif the inputs necessary for the models.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 29-32

.. 3. **Synchronizing the inputs**:
.. In the `sync_inputs` method synchronizes the inputs dictionary with the connectors' states. If the inputs are provided through the connectors, this methods ensures that the model reads them correctly.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 34-47

.. 4. **Setting the inputs**:
.. In the `set_inputs` method, the inputs are set directly by the user. This method ensures that the connectors are updated automatically.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 49-65

.. 5. **Defining required parameters**: 
.. Use the `get_required_parameters` method to list all parameters necessary for the model.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 67-71

.. 6. **Printing setup information**:
.. The `print_setup` method prints details of the connectors and inputs needed to run the simulation.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 73-95

.. 7. **Solving the model**: 
.. The `solve` method checks if the model is ready to be calculated (i.e., if all inputs and parameters are set) before performing the necessary calculations and updating the connectors with the results.

..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 282-335

.. 8. **Updating the connectors**: 
.. After solving the model, the `update_connectors` method updates the state of the connectors based on the calculated results.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 337-346

.. 9. **Printing results**:
.. Use the `print_results` method to display the results of the model calculations.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 348-355

.. 10. **Printing the states of connectors**:
.. The `print_states_connectors` method prints the current state of the connectors after the simulation has been run.
  
..   .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
..      :language: python
..      :lines: 359-370

.. For further details on this model, refer to the [documentation](link).


Inputs and Parameters Setting
-------------------------------

Components' inputs can be set in two different ways: through connectors or directly within the components. The components is defined by its inputs and its connectors. 

Connector Approach
^^^^^^^^^^^^^^^^^^^^
In this approach, the inputs for components are set through connectors. This is particularly useful when multiple components are linked 
to form a larger system. Connectors facilitate the transfer of inputs and outputs between components, enabling the simulation of system-wide 
behavior.

.. image:: ../../../figures/component/Component_connectors.png
   :alt: Components connectors description.
   :width: 100%


Input/Output Approach
^^^^^^^^^^^^^^^^^^^^^
In this approach, inputs are set directly within each component. This is more suitable when using components individually, without linking 
them to other components.

.. image:: ../../../figures/component/Component_in_out.png
   :alt: Components inlet/outlet description.
   :width: 100%

.. toctree::
   :maxdepth: 2
   :caption: Example of usage:

   ../../../notebooks/semi_empirical_cp_example


