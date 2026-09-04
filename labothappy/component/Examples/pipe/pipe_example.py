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
print(f"Mass inventory: {pipe.charge:.2f} kg")


# Test two phase
# Test two-phase with realistic heat pump conditions
pipe = Pipe()
pipe.add_straight(D=0.010, L=2.0).add_curved_elbow(D=0.01, delta=90, R0=0.5).add_straight(D=0.010, L=2.0, theta = -90) # 10mm diameter, 2m long
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
print(f"Pressure drop (Friedel): {(pipe.su.p - pipe.ex.p)/1000:.2f} kPa")
print(f"Mass inventory: {pipe.charge:.6f} kg")


# Test two-phase with the MSH (Muller-Steinhagen & Heck) correlation instead
# of the default Friedel correlation, same conditions as above.
pipe_msh = Pipe()
pipe_msh.add_straight(D=0.010, L=2.0, two_phase_correlation='MSH').add_curved_elbow(D=0.01, delta=90, R0=0.5).add_straight(D=0.010, L=2.0, theta=-90, two_phase_correlation='MSH')

pipe_msh.su.set_p(P)
pipe_msh.su.set_h(h)
pipe_msh.su.set_m_dot(0.05)
pipe_msh.su.set_fluid('R410A')

pipe_msh.solve()
print(f"\nInlet pressure: {pipe_msh.su.p/1000:.1f} kPa")
print(f"Outlet pressure: {pipe_msh.ex.p/1000:.1f} kPa")
print(f"Pressure drop (MSH): {(pipe_msh.su.p - pipe_msh.ex.p)/1000:.2f} kPa")
print(f"Mass inventory: {pipe_msh.charge:.6f} kg")


# Test two-phase charge inventory with different void fraction models,
# same conditions as above (Friedel pressure drop kept fixed).
void_fraction_models = [
    'Homogeneous', 'Zivi', 'Fauske', 'Premoli', 'Lockhart-Martinelli',
    'Hughmark', 'Armand-Treschev', 'Bankoff', 'Rouhani-Axelsson', 'DIX',
    'Woldesemayat-Ghajar', 'Cioncolini-Thome',
]

print(f"\n{'Void fraction model':<22} {'alpha [-]':<12} {'Mass inventory [kg]':<20}")
print("-" * 56)
for vf_model in void_fraction_models:
    pipe_vf = Pipe()
    pipe_vf.add_straight(D=0.010, L=2.0, void_fraction_model=vf_model).add_curved_elbow(D=0.01, delta=90, R0=0.5).add_straight(D=0.010, L=2.0, theta=-90, void_fraction_model=vf_model)

    pipe_vf.su.set_p(P)
    pipe_vf.su.set_h(h)
    pipe_vf.su.set_m_dot(0.05)
    pipe_vf.su.set_fluid('R410A')

    pipe_vf.solve()
    alpha = pipe_vf.components[0].void_fraction
    print(f"{vf_model:<22} {alpha:<12.4f} {pipe_vf.charge:<20.6f}")