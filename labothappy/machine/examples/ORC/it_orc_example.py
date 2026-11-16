import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

from labothappy.machine.circuit_it import IterativeCircuit
from labothappy.connector.mass_connector import MassConnector

from labothappy.component.expander.expander_semi_empirical import ExpanderSE
from labothappy.component.expander.expander_csteff import ExpanderCstEff

from labothappy.component.heat_exchanger.hex_MB_charge_sensitive import HexMBChargeSensitive
from labothappy.component.heat_exchanger.hex_cstpinch import HexCstPinch
from labothappy.component.heat_exchanger.hex_csteff import HexCstEff

from labothappy.component.pump.pump_curve_similarity import PumpCurveSimilarity
from labothappy.component.pump.pump_csteff import PumpCstEff

# Geometries
from labothappy.toolbox.geometries.heat_exchanger.geometry_plate_hx_swep import PlateGeomSWEP

#%%

def preheated_rec_orc(fluid, CSource, HSource, PREsource, P_high, P_low, m_dot, eta_pp, eta_exp, PP_cd, SC_cd, PP_ev, SH_ev, eff_rec):

    # 1) Instanciate Circuit 
    orc = IterativeCircuit(fluid)
    # orc.mute_print()
    
    # 2) Create components 
    Pump = PumpCstEff()
    Condenser = HexCstPinch()
    Expander = ExpanderCstEff()
    Evaporator = HexCstPinch()
    Recuperator = HexCstEff()
    Preheater = HexCstEff()
    
    # 3) Components parameters
    # Pump
    Pump.set_parameters(eta_is=eta_pp)
    
    # Condenser
    Condenser.set_parameters(**{
        'Pinch': PP_cd,
        'Delta_T_sh_sc': SC_cd,
        'HX_type': 'condenser'
    })
    
    # Evaporator
    Evaporator.set_parameters(**{
        'Pinch': PP_ev,
        'Delta_T_sh_sc': SH_ev,
        'HX_type': 'evaporator'
    })

    # Expander
    Expander.set_parameters(eta_is=eta_exp)

    # Recuperator 
    Recuperator.set_parameters(
        eta=eff_rec)
    
    # Preheater 
    Preheater.set_parameters(
        eta=eff_rec)
    
    # 4) Add and link components
    # Add components to circuit
    orc.add_component(Pump, "Pump")
    orc.add_component(Condenser, "Condenser")
    orc.add_component(Expander, "Expander")
    orc.add_component(Evaporator, "Evaporator")
    orc.add_component(Recuperator, "Recuperator")
    orc.add_component(Preheater, "Preheater")
    
    # Link components
    orc.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
    orc.link_components("Recuperator", "m-ex_C", "Preheater", "m-su_C")
    orc.link_components("Preheater", "m-ex_C", "Evaporator", "m-su_C")
    orc.link_components("Evaporator", "m-ex_C", "Expander", "m-su")
    orc.link_components("Expander", "m-ex", "Recuperator", "m-su_H")
    orc.link_components("Recuperator", "m-ex_H", "Condenser", "m-su_H")
    orc.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    # 5) Add sources
    orc.add_source("CD_Water", CSource, orc.components["Condenser"], "m-su_C")
    orc.set_source_properties(T=CSource.T, fluid='Water', P=CSource.p, m_dot = CSource.m_dot, target="Condenser:su_C")
    
    orc.add_source("EV_Water", HSource, orc.components["Evaporator"], "m-su_H")
    orc.set_source_properties(T=HSource.T, fluid='Water', P=HSource.p, m_dot = HSource.m_dot, target="Evaporator:su_H")
    
    orc.add_source("PRE_Water", PREsource, orc.components["Preheater"], "m-su_H")
    orc.set_source_properties(T=PREsource.T, fluid='Water', P=PREsource.p, m_dot = PREsource.m_dot, target="Preheater:su_H")
    
    # 6) Set cycle guesses
    h_pp_guess = PropsSI("H", "T", CSource.T-SC_cd, "P", P_low, fluid)
    h_exp_guess = PropsSI("H", "T", HSource.T+SH_ev, "P", P_high, fluid)

    orc.set_cycle_guess(target="Pump:su", m_dot = m_dot, h=h_pp_guess, p=P_low)
    orc.set_cycle_guess(target="Pump:ex", p=P_high)
        
    orc.set_cycle_guess(target="Expander:su", p=P_high, m_dot = m_dot, SH=SH_ev*2)
    orc.set_cycle_guess(target="Expander:ex", p=P_low)

    # 7) Set iteration variables
    orc.set_iteration_variable(
        target  = 'Expander:ex',
        variable = 'p',
        guess = P_low,
        tolerance=1e-6
        )
                    
    orc.set_iteration_variable(
        target  = ['Pump:ex', 'Expander:su'],
        variable = 'p',
        guess=P_high,
        tolerance=1e-6
        )
    
    orc.set_iteration_variable(
        target="Expander:su",
        variable="h",
        guess=h_exp_guess,
        tolerance=1e-6
        )
    
    # 8) Set residual variables
    orc.set_residual_variable(
        target="Expander:su-p",
        tolerance=1e-3
    )
    
    orc.set_residual_variable(
        target="Pump:su-p",
        tolerance=1e-3
    )
    
    orc.set_residual_variable(
        target="Expander:su-h",
        tolerance=1e-3
    )
    
    return orc

#%%

def off_design_orc(fluid, CSource, HSource, guesses, inputs, exp_params, evap_SWEP_model, cond_SWEP_model, pump_params):
    
    # 1) Instanciate Circuit 
    orc = IterativeCircuit(fluid)
    
    # 2) Create components 
    Expander = ExpanderSE()
    Condenser = HexMBChargeSensitive('Plate')
    Pump = PumpCurveSimilarity()
    Evaporator = HexMBChargeSensitive('Plate')

    # 3) Set component parameters 
    # Expander
    Expander.set_parameters(**exp_params)

    # Evaporator
    evaporator_geom = PlateGeomSWEP()
    evaporator_geom.set_parameters(evap_SWEP_model)

    Evaporator.set_parameters(
        # Set the geometry of the evaporator
        A_c=evaporator_geom.A_c, A_h=evaporator_geom.A_h, h=evaporator_geom.h, l=evaporator_geom.l, l_v=evaporator_geom.l_v, w_v=evaporator_geom.w_v,
        C_CS=evaporator_geom.C_CS, C_Dh=evaporator_geom.C_Dh, C_V_tot=evaporator_geom.C_V_tot, C_canal_t=evaporator_geom.C_canal_t, C_n_canals=evaporator_geom.C_n_canals, 
        H_CS=evaporator_geom.H_CS, H_Dh=evaporator_geom.H_Dh, H_V_tot=evaporator_geom.H_V_tot, H_canal_t=evaporator_geom.H_canal_t, H_n_canals=evaporator_geom.H_n_canals,
        casing_t=evaporator_geom.casing_t, chevron_angle=evaporator_geom.chevron_angle, fooling=evaporator_geom.fooling, 
        n_plates=evaporator_geom.n_plates, plate_cond=evaporator_geom.plate_cond, plate_pitch_co=evaporator_geom.plate_pitch_co, t_plates=evaporator_geom.t_plates, w=evaporator_geom.w, 
        amplitude=evaporator_geom.amplitude, phi=evaporator_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0,
    )

    Evaporator.set_htc(htc_type="Correlation", Corr_H={"1P": "water_plate_HTC"}, Corr_C={"1P": "martin_holger_plate_HTC", "2P": "amalfi_plate_HTC"}) # 'User-Defined' or 'Correlation' # 28
    Evaporator.set_DP()

    # Condenser
    condenser_geom = PlateGeomSWEP()
    condenser_geom.set_parameters(cond_SWEP_model)

    Condenser.set_parameters(
        # Set the geometry of the condenser
        A_c=condenser_geom.A_c, A_h=condenser_geom.A_h, h=condenser_geom.h, l=condenser_geom.l, l_v=condenser_geom.l_v, w_v=condenser_geom.w_v,
        C_CS=condenser_geom.C_CS, C_Dh=condenser_geom.C_Dh, C_V_tot=condenser_geom.C_V_tot, C_canal_t=condenser_geom.C_canal_t, C_n_canals=condenser_geom.C_n_canals,
        H_CS=condenser_geom.H_CS, H_Dh=condenser_geom.H_Dh, H_V_tot=condenser_geom.H_V_tot, H_canal_t=condenser_geom.H_canal_t, H_n_canals=condenser_geom.H_n_canals,
        casing_t=condenser_geom.casing_t, chevron_angle=condenser_geom.chevron_angle, fooling=condenser_geom.fooling,
        n_plates=condenser_geom.n_plates, plate_cond=condenser_geom.plate_cond, plate_pitch_co=condenser_geom.plate_pitch_co, t_plates=condenser_geom.t_plates, w=condenser_geom.w,
        amplitude=condenser_geom.amplitude, phi=condenser_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=10,
    )

    Condenser.set_htc(htc_type="Correlation", Corr_H={"1P": "martin_holger_plate_HTC", "2P": "shah_condensation_plate_HTC"}, Corr_C={"1P": "water_plate_HTC"}) # 'User-Defined' or 'Correlation' # 28
    Condenser.set_DP()

    # Pump
    Pump.set_parameters(**pump_params)

    # 4) Add components to circuit
    orc.add_component(Expander, "Expander")
    orc.add_component(Condenser, "Condenser")
    orc.add_component(Pump, "Pump")
    orc.add_component(Evaporator, "Evaporator")

    # 5) Link components
    orc.link_components("Expander", "m-ex", "Condenser", "m-su_H")
    orc.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    orc.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
    orc.link_components("Evaporator", "m-ex_C", "Expander", "m-su")
    
    # 6) Add fluid sources
    orc.add_source("CD_Water", CSource, orc.components["Condenser"], "m-su_C")
    orc.set_source_properties(T=CSource.T, fluid=CSource.fluid, P=CSource.p, m_dot = CSource.m_dot, target="Condenser:su_C")
    orc.add_source("EV_Water", HSource, orc.components["Evaporator"], "m-su_H")
    orc.set_source_properties(T=HSource.T, fluid=HSource.fluid, P=HSource.p, m_dot = HSource.m_dot, target="Evaporator:su_H")
        
    # 7) Set inputs and guesses
    P_low = guesses["P_low"]
    P_high = guesses["P_high"]
    
    SC_cd = inputs["SC_cd"]
    m_dot_ref = inputs["m_dot_ref"]
    T_amb = inputs["T_amb"]
    N_exp = inputs["N_exp"]
    
    T_sat_LP_guess = PropsSI("T", "P", P_low, "Q", 0.5, fluid)

    h_SC_guess = PropsSI("H", "P", P_low, "T", T_sat_LP_guess - SC_cd, fluid)

    orc.set_cycle_input(target="Pump:su", m_dot = m_dot_ref, SC=SC_cd, P=P_high)
    orc.set_cycle_input(target="Expander:su", p = P_low)
    orc.set_cycle_input(target="Expander:W", N_rot = N_exp)
    orc.set_cycle_input(target="Expander:Q_amb", T_amb=T_amb)

    # 8) Set iteration variables
    orc.set_iteration_variable(
        target=["Expander:ex"],
        variable="p",
        guess=P_low,
        tolerance=1e-6
    )
    
    orc.set_iteration_variable(
        target="Pump:ex",
        variable="p",
        guess=P_high,
        tolerance=1e-6
    )
    
    # -------- 9) Set residual variables --------
    orc.set_residual_variable(target="Condenser:ex_H-SC",  target_value=SC_cd, tolerance=1e-3)
    orc.set_residual_variable(target="Expander:W-N_rot" ,  target_value=N_exp, scale=N_exp, tolerance=1e-3)
    
    return orc

#%%

if __name__ == "__main__":
    
    study_case = "off_design"
    
    if study_case == "off_design":
        
        fluid = 'R1233zd(E)'
        
        "Component parameters"
        
        # Expander
        exp_params = {
            'AU_amb'      : 9.3, 
            'AU_su_n'     : 4.75, 
            'AU_ex_n'     : 17.7, 
            'd_su1'       : 6.48e-3, 
            'm_dot_n'     : 0.1, 
            'A_leak'      : 9.99e-06, 
            'W_dot_loss_0': 2.37e+1 , 
            'alpha'       : 1.16e-1, 
            'C_loss'      : 1.13, 
            'rv_in'       : 1.7, 
            'V_s'         : 3*0.0000712, 
            'mode'        :'P_M',
            'T_amb'       : 293, 
            }
        
        # Evaporator
        evap_SWEP_model = "P200THx140/1P_Evaporator"
        
        # Condenser 
        cond_SWEP_model = "P200THx140/1P_Condenser"
        
        # Pump
        pump_params = {
            'V_dot_curve'   : np.array([20, 30, 40, 50, 60, 70, 80]),   # m3/h
            'Delta_H_curve' : np.array([57, 55, 52, 49, 45, 42, 36]),  # m (head falls with flow)
            'eta_is_curve'  : np.array([0.45, 0.59, 0.69, 0.75, 0.79, 0.79, 0.75]),  # eff peaks near mid-flow
            'NPSH_r_curve'  : np.array([1.1, 1.1, 1.4, 1.8, 2, 3, 4.7]),  # m, increases again near max flow
            'N_rot_rated'   : 2900, # RPM
            'mode'          : "P_M"
        }
        
        "Sources"
        CSource = MassConnector('Water')
        CSource.set_properties(fluid = 'Water', T = 10 + 273.15, p = 2e5, m_dot = 2)
    
        HSource = MassConnector('Water')
        HSource.set_properties(fluid = 'Water', T = 70+273.15, p = 2e5, m_dot = 3)
        
        "Cycle inputs & guesses"
        inputs = {
            'm_dot_ref' : 0.45, # kg/s
            'SC_cd'     : 7,  # K
            'N_exp'     : 6000, # RPM
            'T_amb'     : 293, # K
            'damping_factor_SC' : 0.3
        }
        
        guesses = {
            'P_low'  : PropsSI("P", "T", CSource.T+10, "Q", 0, fluid),
            'P_high' : PropsSI("P", "T", HSource.T-10, "Q", 1, fluid)
        }
        
        orc = off_design_orc(fluid = fluid, CSource = CSource, HSource = HSource, guesses = guesses, inputs = inputs, exp_params=exp_params, evap_SWEP_model = evap_SWEP_model, cond_SWEP_model = cond_SWEP_model, pump_params=pump_params)
        orc.solve()
        
        orc.plot_cycle_Ts()
    
    elif study_case == "Zorlu":
        
        fluid = "Cyclopentane"
        
        # Hot Source
        T_HS = 141+273.15 # K
        p_HS = 10e5 # bar
        fluid_HS = 'Water'
        m_dot_HS = 10000 # kg/s : emulates PCM
    
        # Cold Source
        T_CS = 24 + 273.15
        fluid_CS = 'Water'
        p_CS = 3e5
        m_dot_CS = 900
    
        # Preheater Source
        T_PRE = 113.1 + 273.15
        fluid_PRE = 'Water'
        p_PRE = 2e5
        m_dot_PRE = 60
    
        # Pressure Guesses
        P_high = PropsSI('P', 'T', T_HS, 'Q', 0.5, fluid)
        P_low  = PropsSI('P', 'T', T_CS, 'Q', 0.5, fluid)
        
        m_dot = 34.51      

        HSource = MassConnector()
        HSource.set_properties(fluid = 'Water', T = T_HS, p = p_HS, m_dot = m_dot_HS)
        
        CSource = MassConnector()
        CSource.set_properties(fluid = fluid_CS, T = T_CS, p = p_CS, m_dot = m_dot_CS)
        
        PREsource = MassConnector('Water')
        PREsource.set_properties(fluid = fluid_CS, T = T_CS, p = p_CS, m_dot = m_dot_CS)

        orc = preheated_rec_orc(fluid=fluid, HSource=HSource, CSource=CSource, PREsource = PREsource, P_high=P_high, P_low=P_low, m_dot=m_dot, eta_pp=0.7, eta_exp=0.8, PP_cd = 5, SC_cd = 3, PP_ev = 5, SH_ev = 3, eff_rec = 0.8)
        orc.solve(method="newton")
        
        orc.plot_cycle_Ts()
