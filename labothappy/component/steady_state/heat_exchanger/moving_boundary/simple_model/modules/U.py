
from CoolProp.CoolProp import PropsSI

from math import pi
import numpy as np

from correlations.properties.thermal_conductivity import conducticity_R1233zd

def U_Gnielinski_calibrated(m, D, fluid, P):
    T = PropsSI('T', 'P', P, 'Q', 0, fluid)
    mu = PropsSI('V', 'Q', 0, 'P', P, fluid)
    cp = PropsSI('CPMASS', 'Q', 0, 'P', P, fluid)
    if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
        k = conducticity_R1233zd(T, P)
    else:
        k = PropsSI('L', 'Q', 0, 'P', P, fluid)
    Re = m * 4 / (pi * mu * D)
    Pr = mu * cp / k
    f = (0.79 * np.log(Re) - 1.64) ** (-2)
    c = 4.616163048309070 # calibrated for old model
    # c = 6
    Nu = c * (f / 8 * (Re - 1000) * Pr) / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
    U = k * Nu / D
    return U

def U_DittusBoelter(m, heat_transfer_direction, D, fluid, P, T):

    if heat_transfer_direction == 'heated':
        a = 0.4
    elif heat_transfer_direction == 'cooled':
        a = 0.5

    if T is not None:
        mu = PropsSI('V', 'T', T, 'P', P, fluid)
        cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
        if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
            k = conducticity_R1233zd(T, P)
        else:
            k = PropsSI('L', 'T', T, 'P', P, fluid)
        Re = m * 4 / (pi * mu * D)
        Pr = mu * cp / k
        U = 0.023 * k * Re ** 0.8 * Pr ** a / D
    elif T is None:
        mu = (PropsSI('V', 'P', P, 'Q', 0, fluid) + PropsSI('V', 'P', P, 'Q', 1, fluid)) / 2
        cp = (PropsSI('C', 'P', P, 'Q', 0, fluid) + PropsSI('C', 'P', P, 'Q', 1, fluid)) / 2
        k = (PropsSI('L', 'P', P, 'Q', 0, fluid) + PropsSI('L', 'P', P, 'Q', 1, fluid)) / 2
        Re = m * 4 / (pi * mu * D)
        Pr = mu * cp / k
        U = 0.023 * k * Re ** 0.8 * Pr ** a / D
    return U

def U_Cooper_calibrater(Q, HX_A, P, fluid):
    pr = P/PropsSI('Pcrit', '', 0, ' ', 0, fluid)
    M = PropsSI('MOLARMASS', fluid)
    c = [0.978552748683140, 1.07700234277466] # calibrated for old HX model
    # c = [3, 1.5]
    U = c[0] * ((Q / HX_A)**(0.67 * c[1])) * 55 * (pr**(0.12 - 0.2 * np.log(pr))) * ((-np.log(pr))**(-0.55 * c[1])) * (M**(-0.5))
    return U

def U_Thonon(m, D, fluid, P, T):
    mu = PropsSI('V', 'T', T, 'P', P, fluid)
    cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
    if fluid == 'R1233zd(E)' or fluid == 'R1233ZDE':
        k = conducticity_R1233zd(T, P)
    else:
        k = PropsSI('L', 'T', T, 'P', P, fluid)
    Re = m * 4 / (pi * mu * D)
    Pr = mu * cp / k
    U = 0.2946 * k * Re ** 0.7 * Pr ** (1 / 3) / D
    return U