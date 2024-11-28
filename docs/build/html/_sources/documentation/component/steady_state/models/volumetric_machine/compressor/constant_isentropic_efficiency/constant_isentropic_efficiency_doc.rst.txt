Constant Isentropic Efficiency Compressor model
===============================================

.. raw:: html

   <div class="side-by-side">

.. image:: ../../../../../../../../figures/component/constant_isentropic_efficiency_compressor_in_out.png
   :alt: Description of the plot
   :width: 100%

.. image:: ../../../../../../../../figures/component/constant_isentropic_efficiency_compressor_connectors.png
   :alt: Description of the plot
   :width: 100%

.. raw:: html

   </div>




Model description
-----------------

The constant isentropic efficiency model is a simple model based on the assumption that the isentropic efficiency stays constant.

.. math::

   \epsilon_{is} = \frac{h_{su} - h_{ex, is}}{h_{su} - h_{ex}}

where :math:`\epsilon_{is}` is the isentropic efficiency, :math:`h_{su}` is the specific enthalpy at the supply inlet,
:math:`h_{ex, is}` is the specific enthalpy at the exhaust outlet in isentropic conditions and :math:`h_{ex}` is the specific enthalpy at 
the exhaust outlet.

Based on the isentropic efficiency definition, the specific enthalpy at the exhaust outlet can be calculated and thus also the exhaust temperature.

References
----------


