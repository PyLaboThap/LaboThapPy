Steady-state models
===================

This section provides an overview of the steady-state models for components in the PyLaboThap library.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   base_component/base_component_doc
   models/models_doc


Steady-state models serve as the fundamental building blocks for simulating systems. They are used to represent
individual components under steady-state conditions and can be interconnected using connectors to model more complex systems. Components 
can be defined in two different ways. The first approach is to use connectors to set the inputs for components, while the second approach
is to set the inputs directly within each component.

Connector Approach
-------------------
In this approach, the inputs for components are set through connectors. This is particularly useful when multiple components are linked 
to form a larger system. Connectors facilitate the transfer of inputs and outputs between components, enabling the simulation of system-wide 
behavior.

.. image:: ../../../../figures/component/Component_connectors.png
   :alt: Components connectors description.
   :width: 100%


Input/Output Approach
---------------------
In this approach, inputs are set directly within each component. This is more suitable when using components individually, without linking 
them to other components.

.. image:: ../../../../figures/component/Component_in_out.png
   :alt: Components inlet/outlet description.
   :width: 100%

.. toctree::
   :maxdepth: 2
   :caption: Example of usage:

   ../../../notebooks/semi_empirical_cp_example

.. toctree::
   :maxdepth: 2
   :caption: Example of setup:

   base_component/exampleofsetup_doc

