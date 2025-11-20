Pump - Constant Isentropic Efficiency Model
================================================

Model description
-----------------

The constant isentropic efficiency model is a simple model based on the assumption that the isentropic efficiency stays constant.

.. math::

   \epsilon_{is} = \frac{h_{su} - h_{ex, is}}{h_{su} - h_{ex}}

where :math:`\epsilon_{is}` is the isentropic efficiency, :math:`h_{su}` is the supply specific enthalpy,
:math:`h_{ex, is}` is the isentropic exhaust specific enthalpy and :math:`h_{ex}` is the exhaust specific enthalpy.

Based on the isentropic efficiency definition, the exhaust specific enthalpy can be calculated and thus also the exhaust temperature.

.. autoclass:: component.pump.pump_csteff.PumpCstEff

References
----------
/

