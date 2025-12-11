Heat Exchanger - Discretized Constant Efficiency Model
==========================================

Model description
-----------------

This document describes a discretized counterflow heat exchanger model that both a
**constant effectiveness** and a **Minimum Pinch** approach. The user provides the heat
exchanger effectiveness, denoted :math:`\varepsilon` (epsilon), as a
model parameter.

The heat transfer rate is computed as:

.. math::

   \dot{Q} = \varepsilon \cdot \dot{Q}_{\max}

where :math:`\dot{Q}` is the actual heat transfer rate and
:math:`\dot{Q}_{\max}` is the **maximum possible heat transfer rate**
(based on an external and internal pinching analysis and an ideal exchanger 
assumption).

The maximum possible heat transfer rate is the one that corresponds to a 
temperature pinch (\Delta T_{pp}) equal to zero. It is determined before 
solving by varying a fictive heat exchanger effectiveness from 1 to 0 to by 
increment of 1% until a pinch larger than zero is obtained.

.. math::

   \dot{Q}_{\max} = \dot{Q} such that \Delta T_{pp} = 0

Then, the heat exchanger is solved by decrementing its efficiency (starting 
from the user-defined value acting as a maximum efficiency : \varepsilon) 
until the imposed minimum temperature pinch (\Delta T_{pp,min}) is satisfied.

Pressure drops can be imposed on both sides of the heat exchanger, these are
equally distributed over the discretizations.

**Definitions**

- :math:`\dot{m}_H`, :math:`\dot{m}_C`  
  Mass flow rates of the hot and cold streams, respectively (kg/s).

- :math:`h_{su,H}`, :math:`h_{su,C}`  
  Specific enthalpies at the **inlets** of the hot and cold streams
  (J/kg).

- :math:`h_{ex,id,H}`, :math:`h_{ex,id,C}`  
  Specific enthalpies at the **outlets** of the hot and cold streams
  in an **ideal** heat exchanger (i.e., infinite area / infinite NTU).

**Procedure to compute heat transfer rate**


1. Compute or obtain the inlet temperatures and specific enthalpies  
   :math:`h_{su,H}` and :math:`h_{su,C}`.

2. Impose and distribute the pressure drops.  

3. Determine :math:`\dot{Q}_{max}` by decrementing a fictive heat exchanger 
effectiveness (starting from 1 and down 1% per iteration). :math:`\dot{Q}_{max}` is the
heat rate value when the minimum pinch (at every discretization) becomes larger 
or equal to 0.

4. Compute :math:`\dot{Q}` from :math:`\dot{Q}_{max}` in a similar way: by decrementing the 
heat exchanger effectiveness (starting from the user-defined maximum value) 
until the imposed minimum temperature pinch is satisfied at every 
discretization.

.. math::

      \dot{Q} = \varepsilon \cdot \dot{Q}_{\max}

5. Outlet enthalpies are computed as:

.. math::

   h_{out,H} = h_{su,H} - \frac{\dot{Q}}{\dot{m}_H}

.. math::

   h_{out,C} = h_{su,C} + \frac{\dot{Q}}{\dot{m}_C}

**Assumptions**

- Counterflow geometry.
- User-specified maximum effectiveness :math:`\varepsilon`,  
  with :math:`0 < \varepsilon \le 1`.
- Ideal outlet enthalpies computed using an iterative method for internal 
pinching (a solving method could be employed but proved to be slower for 
efficient heat exchangers).
- No heat losses to the environment.
- Pressure drop equally distributed along the heat exchanger discretizations.

Class description
-----------------

.. autoclass:: component.heat_exchanger.hex_csteff_disc.HexCstEffDisc

Example of use
-----------------
.. literalinclude:: ../../../../../../labothappy/component/examples/heat_exchanger/hex_csteff_disc_example.py
   :language: python


References
----------
/