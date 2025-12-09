Pump - Curves Similarity Law Model
===================================

Model description
-----------------


The :code:`PumpCurveSimilarity` component models a pump using manufacturer
performance curves together with the classical pump affinity (similarity) laws.
It allows computing the pump operating point in three possible modes:

- ``P_N`` : Given suction and discharge pressures and rotational speed,
  compute flow rate and shaft power.
- ``P_M`` : Given suction and discharge pressures and mass flow rate,
  compute the rotational speed required to reach this operating point.
- ``M_N`` : Given suction conditions, mass flow rate and speed,
  compute the discharge pressure and power.

The model uses four characteristic curves at rated speed:

- Head rise :math:`\Delta H(Q)`
- Isentropic efficiency :math:`\eta_{\text{is}}(Q)`
- Required NPSH :math:`\text{NPSH}_r(Q)`
- (Power curve is derived from head and efficiency)

These curves are supplied at a *rated* speed :math:`N_\text{rated}` and are
scaled to other rotational speeds using similarity laws.
The model description is available in the `PDF document <../../../../_static/pdf_files/PumpCurveSimLaw.pdf>`_.

Class description
-----------------

.. autoclass:: component.pump.pump_xurve_similarity.PumpCurveSimilarity

References
----------
- `PDF document <../../../../_static/pdf_files/PumpCurveSimilarity.pdf>`_