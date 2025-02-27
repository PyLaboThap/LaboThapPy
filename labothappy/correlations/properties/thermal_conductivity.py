
from CoolProp.CoolProp import PropsSI

def conducticity_R1233zd(T, P):
    """ Richard A. Perkins and Marcia L. Huber

        "Measurement and Correlation of the Thermal Conductivity of trans-1-Chloro-3,3,3-trifluoropropene (R1233zd(E))"
        J. Chem. Eng. Data 2017, 62, 2659-2665
        Note that the critical enhancement contribution is not implemented.

        Same source as EES
       """
    # Dilute gas thermal conductivity
    T_c = 382.52 #[K]
    A_0 = -0.0103589
    A_1 = 0.0308929
    A_2 = 0.000230348

    k_0 = A_0 + A_1*(T/T_c) + A_2*(T/T_c)**2

    # Residual gas thermal conductivity
    try:
        Rho = PropsSI('D', 'T', T, 'P', P, 'R1233zd(E)')
    except:
        Rho = PropsSI('D', 'T', T, 'Q', 0, 'R1233zd(E)')
    Rho_c = 489.24 #[kg/m^3]

    B_11 = -0.0428296
    B_12 = 0.0434288
    B_21 = 0.0927099
    B_22 = -0.0605844
    B_31 = -0.0702107
    B_32 = 0.0440187
    B_41 = 0.0249708
    B_42 = -0.0155082
    B_51 = -0.00301838
    B_52 = 0.0021019

    Delta_kr = (B_11 + B_12*(T/T_c))*Rho/Rho_c + (B_21 + B_22*(T/T_c))*(Rho/Rho_c)**2 + (B_31 + B_32*(T/T_c))*(Rho/Rho_c)**3 + (B_41 + B_42*(T/T_c))*(Rho/Rho_c)**4 + (B_51 + B_52*(T/T_c))*(Rho/Rho_c)**5

    k = k_0 + Delta_kr
    return k