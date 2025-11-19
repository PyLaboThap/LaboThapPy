Heat Exchanger - Constant Efficiency Model
==================================================================

Model description
-----------------
Heat Exchanger — Constant Effectiveness Model
============================================

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

   \dot{Q}_{\max} = \min\!\left(
     \dot{m}_H \cdot \bigl(h_{su,H} - h_{ex,id,H}\bigr),
     \; \dot{m}_C \cdot \bigl(h_{ex,id,C} - h_{su,C}\bigr)
   \right)

**Definitions**:

- :math:`\dot{m}_H`, :math:`\dot{m}_C`  
  Mass flow rates of the hot and cold streams, respectively (units: kg/s).

- :math:`h_{su,H}`, :math:`h_{su,C}`  
  Specific enthalpies at the **inlets** of the hot and cold streams,
  respectively (units: J/kg).

- :math:`h_{ex,id,H}`, :math:`h_{ex,id,C}`  
  Specific enthalpies at the **outlets** of the hot and cold streams
  in an **ideal** heat exchanger (i.e., infinite area / infinite NTU).
  These are computed from the inlet temperatures assuming complete
  counterflow exchange (see "Ideal outlet conditions" below).

**Procedure to compute heat transfer rate**:
1. From the inlet conditions compute or obtain the inlet temperatures
   and specific enthalpies :math:`h_{su,H}` and :math:`h_{su,C}`.

2. Determine the *ideal* outlet temperatures for a perfect
   counterflow exchanger:
   - The ideal cold-stream outlet temperature equals the hot-stream
     inlet temperature.
   - The ideal hot-stream outlet temperature equals the cold-stream
     inlet temperature.

   Convert those ideal outlet temperatures to enthalpies to obtain
   :math:`h_{ex,id,C}` and :math:`h_{ex,id,H}` (using the fluid property
   relations appropriate for each stream).

3. Compute the two candidate heat transfer rates:
   - Hot-limited: :math:`\dot{Q}_{H,\max} = \dot{m}_H \cdot (h_{su,H} - h_{ex,id,H})`
   - Cold-limited: :math:`\dot{Q}_{C,\max} = \dot{m}_C \cdot (h_{ex,id,C} - h_{su,C})`

4. Take :math:`\dot{Q}_{\max} = \min(\dot{Q}_{H,\max}, \dot{Q}_{C,\max})`.

5. Multiply by the specified effectiveness to get the actual heat
   transfer rate:
   :math:`\dot{Q} = \varepsilon \cdot \dot{Q}_{\max}`.

.. math::

   h_{out,H} = h_{su,H} - \frac{\dot{Q}}{\dot{m}_H}

.. math::

   h_{out,C} = h_{su,C} + \frac{\dot{Q}}{\dot{m}_C}

**Assumptions:**

- The model assumes **counterflow** geometry.
- Effectiveness :math:`\varepsilon` is supplied by the user and is
  treated as constant (0 < :math:`\varepsilon` \le 1).
- The ideal outlet enthalpies are computed by assuming that the
  maximum possible temperature difference is the inlet temperature of
  the opposing stream (i.e., ideal counterflow limit).
- No heat losses to the environment are considered.
- No pressure drops are considered in the heat exchanger.


.. autoclass:: component.heat_exchanger.hex_csteff.HXCstEff

References
----------
/