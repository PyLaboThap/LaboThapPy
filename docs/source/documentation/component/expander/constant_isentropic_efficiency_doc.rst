Constant isentropic efficiency compressor model
===============================================

Model description
-----------------

.. autoclass:: component.steady_state.volumetric_machine.expander.constant_isentropic_efficiency.simulation_model.ExpanderCstEff

The constant isentropic efficiency model is a simple model based on the assumption that the isentropic efficiency stays constant.

.. math::

   \epsilon_{is} = \frac{h_{ex} - h_{su}}{h_{ex, is} - h_{su}}

where :math:`\epsilon_{is}` is the isentropic efficiency, :math:`h_{su}` is the specific enthalpy at the supply inlet,
:math:`h_{ex, is}` is the specific enthalpy at the exhaust outlet in isentropic conditions and :math:`h_{ex}` is the specific enthalpy at 
the exhaust outlet.

Based on the isentropic efficiency definition, the specific enthalpy at the exhaust outlet can be calculated and thus also the exhaust temperature.

References
----------


