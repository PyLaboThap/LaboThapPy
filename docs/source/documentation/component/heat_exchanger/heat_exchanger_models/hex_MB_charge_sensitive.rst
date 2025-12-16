Heat Exchanger - Moving Boundary Charge Sensitive Model
==========================================

Model description
-----------------

The model implemented is based on the moving boundary approach as described in [1]. This model can be used for the
following geometris of heat exchangers: plate, tube and fins, shell and tube, and PCHE. An example of use for 
each geometry is provided in the examples section.

Class description
-----------------
.. autoclass:: component.heat_exchanger.hex_MB_charge_sensitive.HexMBChargeSensitive

Example of use
-----------------
**Plate heat exchanger:**
.. literalinclude:: ../../examples/hex_plate_MB_charge_sensitive_example.py
    :language: python

**Tube and fins heat exchanger:**
.. literalinclude:: ../../examples/hex_T&F_MB_charge_sensitive_example.py
    :language: python

**Shell and tube heat exchanger:**
.. literalinclude:: ../../examples/hex_S&T_MB_charge_sensitive_example.py
    :language: python

**PCHE heat exchanger:**
.. literalinclude:: ../../examples/hex_PCHE_MB_charge_sensitive_example.py
    :language: python

References
----------
[1] I. H. Bell et al., ‘A generalized moving-boundary algorithm to predict the heat transfer rate of counterflow heat exchangers for any phase configuration’, Applied Thermal Engineering, vol. 79, pp. 192–201, Mar. 2015, doi: 10.1016/j.applthermaleng.2014.12.028.
