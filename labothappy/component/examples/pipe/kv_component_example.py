
from CoolProp.CoolProp import PropsSI

from labothappy.component.pipe.kv_component import KvComponent

fl = "R1233ZDE"
evr = KvComponent()      # V_int: placeholder, to measure
evr.set_parameters(Kv=9.8, V_int=50e-6)  # V_int: placeholder, to measure
cases = [("liquid line",    10e5, PropsSI("T", "P", 10e5, "Q", 0, fl) - 5),
            ("suction vapour",  2e5, PropsSI("T", "P", 2e5,  "Q", 1, fl) + 5)]
for label, P, T in cases:
    evr.set_inputs(fluid=fl, P_su=P, h_su=PropsSI("H", "P", P, "T", T, fl), m_dot=0.3)
    o = evr.solve()
    print(f"{label:15s} dp={o['dp']/1e5:6.3f} bar"
            f"M={o['charge']*1e3:.2f} g")


# I NEED TO CHANGE TO BE ABLE TO TAKE ZIVI INTO ACCOUNT!!!
# ALRIGHT I THINK I AM LOOSING MY MIND OVER NOTHING AHAHAH
# OKAY LEAVE IT FOR NOW, I WILL COME BACK TO IT LATER