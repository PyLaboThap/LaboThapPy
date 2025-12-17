Crossflow Fin and Tube Heat Exchanger - Finite Volume Model
==================================================================

Model description
-----------------

The model geometrical implementation is similar to the model with 8 circuits from the reference [1] (see Figure 17). 

.. image:: ../../../../_static/figures/hex_crossflowfintube_finitevolume/finite_volume_model.png
   :alt: Scheme of the heat exchanger discretization.
   :width: 100%

Class description
-----------------
.. autoclass:: component.heat_exchanger.hex_crossflowfintube_finitevolume.HexCrossFlowTubeAndFinsFiniteVolume

Example of use
-----------------
.. literalinclude:: ../../../../../../labothappy/component/examples/heat_exchanger/hex_crossflowfintube_finitevolume_example.py
   :language: python

References
----------
[1] S. Macchitella, G. Colangelo, and G. Starace, ‘Performance Prediction of Plate-Finned Tube Heat Exchangers for Refrigeration: A Review on Modeling and Optimization Methods’, Energies, vol. 16, no. 4, p. 1948, Feb. 2023, doi: 10.3390/en16041948.

