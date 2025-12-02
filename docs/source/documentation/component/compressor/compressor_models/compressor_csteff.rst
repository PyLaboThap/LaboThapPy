Compressor - Constant Isentropic Efficiency Model
================================================

Model description
-----------------

The constant isentropic efficiency model is a simple model based on the assumption that the isentropic efficiency stays constant.

.. math::

   \varepsilon_{is} = \frac{h_{su} - h_{ex, is}}{h_{su} - h_{ex}}

where :math:`\epsilon_{is}` is the isentropic efficiency, :math:`h_{su}` is the specific enthalpy at the supply inlet,
:math:`h_{ex, is}` is the specific enthalpy at the exhaust outlet in isentropic conditions and :math:`h_{ex}` is the specific enthalpy at 
the exhaust outlet.

Based on the isentropic efficiency definition, the exhaust specific enthalpy can be calculated and thus also the exhaust temperature.

Class description
-----------------
.. autoclass:: component.compressor.compressor_csteff.CompressorCstEff

Example of use
-----------------
.. literalinclude:: ../../../../../../labothappy/component/examples/compressor/compressor_csteff_example.py

References
----------
/

