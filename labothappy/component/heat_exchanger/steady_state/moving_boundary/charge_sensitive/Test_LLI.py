
from CoolProp.CoolProp import AbstractState
import CoolProp as CP

state = AbstractState("HEOS", "Cyclopentane")
state1 = CP.AbstractState("SRK", "Cyclopentane")

print(state)

h = 100000
p = 100000

state.update(CP.HmassP_INPUTS, h, p)

print(state)

state1.update(CP.QT_INPUTS, 1, 250)

T = state.T()
T1 = state1.T()

print(T)
print(T1)
