Connector
=========

This section includes details on the connectors of the PyLaboThap library.

The purpose of the **Connector class** is to link components between each other. There are three different types of connectors that represent different
types of energy transfer within the system:

- **Mass Connector**: energy transfer in the form of mass flow (:math:`\dot{m} \cdot h`).
- **Heat Connector**: energy transfer in the form of heat flow (:math:`\dot{Q}`).
- **Work Connector**: energy transfer in the form of work flow (:math:`\dot{W}`).

.. image:: ../../../figures/connector/connectors.png
   :alt: Connectors description.
   :width: 100%

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   mass_connector/mass_connector_doc
   heat_connector/heat_connector_doc
   work_connector/work_connector_doc

Example
-------

Consider a compressor connected to a motor. 
The compressor is connected with connectors:

.. raw:: html

   <div class="side-by-side">

.. image:: ../../../figures/connector/compressor_example.png
   :alt: Connectors for a compressor model.
   :width: 70%

.. raw:: html

   <ul>
        <li><strong>Two mass connectors</strong>: supply and exhaust mass flow that can be connected to an evaporator and a condenser, for example.</li>
        <li><strong>One heat connector</strong>: heat flow is connected to the ambient air.</li>
        <li><strong>One work connector</strong>: work flow is connected to a motor.</li>
   </ul>

.. raw:: html

   </div>
