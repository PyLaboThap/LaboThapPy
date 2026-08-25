from labothappy.component.pipe.pipe import Pipe
from CoolProp.CoolProp import PropsSI

# Simple: chain fluently
pipe = Pipe()
pipe.add_straight(D=0.01, L=0.5, theta = 0).add_curved_elbow(D=0.01, delta=90, R0=0.5).add_straight(D=0.01, L=0.5, theta = 0.0) # 10mm diameter, 1.5m long

pipe.su.set_p(101325)
pipe.su.set_h(100000)
pipe.su.set_m_dot(0.5)
pipe.su.set_fluid('Water')

pipe.solve()
print(f"Mass inventory: {pipe.m_charge:.2f} kg")


# Test two phase
# Test two-phase with realistic heat pump conditions
pipe = Pipe()
pipe.add_curved_elbow(D=0.01, delta=90, R0=0.5).add_straight(D=0.010, L=2.0, theta = -90) # 10mm diameter, 2m long
# add_straight(D=0.010, L=2.0).
# R410A at evaporator outlet: low pressure, high quality
# Typical evaporator exit: P ≈ 400 kPa, x ≈ 0.9 (superheated vapor)
P = 400000  # Pa
x_target = 0.9
h = PropsSI('H', 'P', P, 'Q', x_target, 'R410A')

pipe.su.set_p(P)
pipe.su.set_h(h)
pipe.su.set_m_dot(0.05)  # 50 g/s (typical for small heat pump)
pipe.su.set_fluid('R410A')

pipe.solve()
print(f"Inlet pressure: {pipe.su.p/1000:.1f} kPa")
print(f"Outlet pressure: {pipe.ex.p/1000:.1f} kPa")
print(f"Pressure drop: {(pipe.su.p - pipe.ex.p)/1000:.2f} kPa")
print(f"Mass inventory: {pipe.m_charge:.6f} kg")