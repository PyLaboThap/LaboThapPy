Heat Exchanger - Epsilon NTU Model
==========================================

Model description
-----------------
The HexeNTU component models a steady-state heat exchanger using the effectiveness–NTU (ε-NTU) method. 
This approach is well suited for situations where outlet temperatures are unknown and must be determined from inlet conditions, 
flow rates, and exchanger geometry.

The `HexeNTU` component models a steady-state heat exchanger using the **effectiveness–NTU (ε-NTU)** method. This method computes the 
heat transfer rate directly from inlet conditions, without requiring outlet temperatures a priori.

**Heat Capacity Rates**

For each fluid, the heat capacity rate is defined as:

.. math::

    C = \dot{m} \cdot c_p


where (:math:`\dot{m}`) is the mass flow rate and (:math:`c_p`) is the specific heat capacity evaluated at inlet conditions.
The minimum and maximum heat capacity rates are:

.. math::

   C_{\min} = \min(C_H, C_C), \qquad
    C_{\max} = \max(C_H, C_C)

The heat capacity ratio is:

.. math::

   C_r = \frac{C_{\min}}{C_{\max}}


**Overall Heat Transfer Coefficient (AU)**

The overall heat transfer coefficient multiplied by area, (AU), is computed from a series of thermal resistances:

.. math::

   \frac{1}{AU} = \frac{1}{h_H A} \cdot \frac{1}{h_C A} \cdot \frac{t_{plate}}{k_{plate} A} \cdot \frac{R_{fouling}}{A}


The convective heat transfer coefficients (:math:`h_H`) and (:math:`h_C`) are evaluated using the **Gnielinski correlation** for internal turbulent flow.


**Number of Transfer Units (NTU)**

The number of transfer units is defined as:

.. math::

   NTU = \frac{AU}{C_{\min}}


NTU represents the relative size of the heat exchanger compared to the ability of the fluids to store thermal energy.

**Heat Exchanger Effectiveness**

The heat exchanger effectiveness (:math:`\varepsilon`) is computed using a standard **ε-NTU correlation**, which depends on:

* the flow configuration (counterflow, parallel flow, crossflow, etc.),
* the number of transfer units (NTU),
* the heat capacity ratio (:math:`C_r`).

.. math::

   \varepsilon = f(NTU, C_r, \text{Flow Type})


**Heat Transfer Rate**

The maximum possible heat transfer is defined as:

.. math::

   Q_{\max} = C_{\min} \left( T_{H,in} - T_{C,in} \right)

The actual heat transfer rate is then:

.. math::

   \dot{Q} = \varepsilon , Q_{\max}


**Assumptions**

* Steady-state operation
* No heat loss to the environment
* No pressure drop
* Single-phase flow
* Thermophysical properties evaluated at inlet conditions


Class description
-----------------
.. autoclass:: component.heat_exchanger.hex_eNTU.HexeNTU

Example of use
-----------------
.. literalinclude:: ../../../../../../labothappy/component/examples/heat_exchanger/hex_eNTU_example.py
   :language: python

References
----------
/