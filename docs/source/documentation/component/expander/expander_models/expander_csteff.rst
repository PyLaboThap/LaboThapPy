Expander - Constant Isentropic Efficiency Model
===============================================

Model description
-----------------

The constant isentropic efficiency model is a simple model based on the assumption that the isentropic efficiency stays constant.

.. math::

   \varepsilon_{is} = \frac{h_{ex} - h_{su}}{h_{ex, is} - h_{su}}

where :math:`\varepsilon_{is}` is the isentropic efficiency, :math:`h_{su}` is the supply specific enthalpy,
:math:`h_{ex, is}` is the exhaust specific isentropic enthalpy and :math:`h_{ex}` is the exhaust specific enthalpy.

Based on the isentropic efficiency definition, the specific enthalpy at the exhaust outlet can be calculated and thus also the exhaust temperature.

Class description
-----------------
.. autoclass:: component.expander.expander_csteff.ExpanderCstEff

Example of use
-----------------
.. literalinclude:: ../../examples/expander_csteff_example.py
   :language: python

References
----------
/

