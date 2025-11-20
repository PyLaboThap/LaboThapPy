Heat Exchanger - Constant Efficiency Model
==========================================

Model description
-----------------

This document describes a counterflow heat exchanger model that uses a
**constant effectiveness** approach. The user provides the heat
exchanger effectiveness, denoted :math:`\varepsilon` (epsilon), as a
model parameter.

The heat transfer rate is computed as:

.. math::

   \dot{Q} = \varepsilon \cdot \dot{Q}_{\max}

where :math:`\dot{Q}` is the actual heat transfer rate and
:math:`\dot{Q}_{\max}` is the **maximum possible heat transfer rate**
(based on the inlet conditions and an ideal exchanger).

The maximum possible heat transfer rate is the smaller of the maximum
energy that the hot stream can give up and the maximum energy that the
cold stream can accept:

.. math::

   \dot{Q}_{\max} = \min \left(
     \dot{m}_H \cdot (h_{su,H} - h_{ex,id,H}),
     \; \dot{m}_C \cdot (h_{ex,id,C} - h_{su,C})
   \right)

**Definitions:**

- :math:`\dot{m}_H`, :math:`\dot{m}_C`  
  Mass flow rates of the hot and cold streams, respectively (kg/s).

- :math:`h_{su,H}`, :math:`h_{su,C}`  
  Specific enthalpies at the **inlets** of the hot and cold streams
  (J/kg).

- :math:`h_{ex,id,H}`, :math:`h_{ex,id,C}`  
  Specific enthalpies at the **outlets** of the hot and cold streams
  in an **ideal** heat exchanger (i.e., infinite area / infinite NTU).

**Procedure to compute heat transfer rate:**

1. Compute or obtain the inlet temperatures and specific enthalpies  
   :math:`h_{su,H}` and :math:`h_{su,C}`.

2. Determine the *ideal* outlet temperatures for a perfect counterflow exchanger:

   - The ideal cold-stream outlet temperature equals the hot-stream inlet temperature.
   - The ideal hot-stream outlet temperature equals the cold-stream inlet temperature.

   Convert these ideal temperatures to enthalpies:
   :math:`h_{ex,id,C}` and :math:`h_{ex,id,H}`.

3. Compute the two candidate maximum heat transfer rates:

   - Hot-limited: :math:`\dot{Q}_{H,\max} = \dot{m}_H \cdot (h_{su,H} - h_{ex,id,H})`
   - Cold-limited: :math:`\dot{Q}_{C,\max} = \dot{m}_C \cdot (h_{ex,id,C} - h_{su,C})`

4. Compute:

   .. math::

      \dot{Q}_{\max} = \min(\dot{Q}_{H,\max}, \dot{Q}_{C,\max})

5. The actual heat transfer rate is then:

   .. math::

      \dot{Q} = \varepsilon \cdot \dot{Q}_{\max}

Outlet enthalpies are computed as:

.. math::

   h_{out,H} = h_{su,H} - \frac{\dot{Q}}{\dot{m}_H}

.. math::

   h_{out,C} = h_{su,C} + \frac{\dot{Q}}{\dot{m}_C}

**Assumptions:**

- Counterflow geometry.
- User-specified constant effectiveness :math:`\varepsilon`,  
  with :math:`0 < \varepsilon \le 1`.
- Ideal outlet enthalpies computed using the ideal counterflow limit.
- No heat losses to the environment.
- No pressure drop inside the heat exchanger.

Class description
-----------------
.. autoclass:: component.heat_exchanger.hex_csteff.HexCstEff

References
----------
/