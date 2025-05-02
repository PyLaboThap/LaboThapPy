Documentation
=============

LaboThApPy is a Python library for the simulation of thermodynamic cycles. It is designed to be modular 
and easy to use, allowing users to create and simulate different thermodynamic cycles with minimal effort. 
The library is built using an object-oriented programming approach, which makes it highly adaptable and flexible.

At the core of LaboThApPy are three primary modules:
- **Connector**: Allows the linking of components through energy transfers that occur between them.
- **Component**: Serves as the fundamental building blocks, representing individual parts of a system, 
such as compressors, turbines, and heat exchangers.
- **Machine**: Acts as a container for each simulation, representing entire systems and organizing the assembly of 
components and connections into cohesive models.

In this documentation, the different modules and their functionalities are described in detail.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   connector/connector_doc
   component/component_doc