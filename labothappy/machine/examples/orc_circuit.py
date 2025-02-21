from machine.circuit import Circuit
from CoolProp.CoolProp import PropsSI

from connector.mass_connector import MassConnector

from component.heat_exchanger.steady_state.moving_boundary.charge_sensitive.simulation_model_elise_2 import HeatExchangerMB
from component.heat_exchanger.steady_state.moving_boundary.charge_sensitive.modules.geometry_plate_hx_swep import PlateGeomSWEP
from component.volumetric_machine.expander.steady_state.semi_empirical.simulation_model import ExpanderSE
from component.volumetric_machine.pump.steady_state.polynomial_efficiency.simulation_model_Pout import PumpPolyEff
from component.tank.spliter.simulation_model import Spliter
from component.tank.mixer.simulation_model import Mixer


#%% DEFINE MODELS

if __name__ == "__main__":

    case_study = 'Recuperator' # 'Recuperator' or 'Simple'    

    if case_study == 'Simple':
        ORC = Circuit('R1233zd(E)')
        
        # Create components
        Expander_1 = ExpanderSE()
        Expander_2 = ExpanderSE()
        Expander_3 = ExpanderSE()
        
        Spliter_expander = Spliter(outlet_repartition = [1/3,1/3,1/3])
        Mixer_expander = Mixer(n_inlets = 3)
        
        Condenser = HeatExchangerMB('Plate')
        Pump = PumpPolyEff()
        Evaporator = HeatExchangerMB('Plate')
    
    #%% CONDENSER PARAMETERS
    
        # Condenser
        condenser_geom = PlateGeomSWEP()
        condenser_geom.set_parameters("P200THx140/1P_Condenser")
    
        Condenser.set_parameters(
            # Set the geometry of the condenser
            A_c=condenser_geom.A_c, A_h=condenser_geom.A_h, h=condenser_geom.h, l=condenser_geom.l, l_v=condenser_geom.l_v, w_v=condenser_geom.w_v,
            C_CS=condenser_geom.C_CS, C_Dh=condenser_geom.C_Dh, C_V_tot=condenser_geom.C_V_tot, C_canal_t=condenser_geom.C_canal_t, C_n_canals=condenser_geom.C_n_canals,
            H_CS=condenser_geom.H_CS, H_Dh=condenser_geom.H_Dh, H_V_tot=condenser_geom.H_V_tot, H_canal_t=condenser_geom.H_canal_t, H_n_canals=condenser_geom.H_n_canals,
            casing_t=condenser_geom.casing_t, chevron_angle=condenser_geom.chevron_angle, fooling=condenser_geom.fooling,
            n_plates=condenser_geom.n_plates, plate_cond=condenser_geom.plate_cond, plate_pitch_co=condenser_geom.plate_pitch_co, t_plates=condenser_geom.t_plates, w=condenser_geom.w,
            amplitude=condenser_geom.amplitude, phi=condenser_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)
    
        Corr_H = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_cond_BPHEX"}
        Corr_C = {"1P" : "water_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}
    
        # Set the pressure drop correlations of the condenser
        Condenser.set_DP()
        
        # Set the heat transfer coefficients correlations of the condenser           
        Condenser.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)
        
    #%% EVAPORATOR PARAMETERS
    
        # Evaporator
        evaporator_geom = PlateGeomSWEP()
        evaporator_geom.set_parameters("P200THx140/1P_Evaporator")
    
        Evaporator.set_parameters(
            # Set the geometry of the evaporator
            A_c=evaporator_geom.A_c, A_h=evaporator_geom.A_h, h=evaporator_geom.h, l=evaporator_geom.l, l_v=evaporator_geom.l_v, w_v=evaporator_geom.w_v,
            C_CS=evaporator_geom.C_CS, C_Dh=evaporator_geom.C_Dh, C_V_tot=evaporator_geom.C_V_tot, C_canal_t=evaporator_geom.C_canal_t, C_n_canals=evaporator_geom.C_n_canals, 
            H_CS=evaporator_geom.H_CS, H_Dh=evaporator_geom.H_Dh, H_V_tot=evaporator_geom.H_V_tot, H_canal_t=evaporator_geom.H_canal_t, H_n_canals=evaporator_geom.H_n_canals,
            casing_t=evaporator_geom.casing_t, chevron_angle=evaporator_geom.chevron_angle, fooling=evaporator_geom.fooling, 
            n_plates=evaporator_geom.n_plates, plate_cond=evaporator_geom.plate_cond, plate_pitch_co=evaporator_geom.plate_pitch_co, t_plates=evaporator_geom.t_plates, w=evaporator_geom.w, 
            amplitude=evaporator_geom.amplitude, phi=evaporator_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)
        
        Corr_H = {"1P" : "water_plate_HTC", "2P" : "Han_cond_BPHEX"}
        Corr_C = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}
    
        # Set the pressure drop correlations of the condenser
        Evaporator.set_DP()
        
        # Set the heat transfer coefficients correlations of the condenser           
        Evaporator.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)
    
    #%% EXPANDER PARAMETERS AND INPUTS
    
        N_exp = 6000
        T_amb = 293
    
        Expander_1.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)
    
        Expander_1.set_inputs(N_rot=N_exp, T_amb=T_amb)
    
        Expander_2.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)
    
        Expander_2.set_inputs(N_rot=N_exp, T_amb=T_amb)
        
        Expander_3.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)
    
        Expander_3.set_inputs(N_rot=N_exp, T_amb=T_amb)
    
    #%% PUMP INPUTS
    
        Pump.set_inputs(N_pp = 2500)
    
    #%% ADD AND LINK COMPONENTS
        ORC.add_component(Expander_1, "Expander_1")
        ORC.add_component(Expander_2, "Expander_2")
        ORC.add_component(Expander_3, "Expander_3")
    
        ORC.add_component(Mixer_expander, "Mixer")
        
        ORC.add_component(Condenser, "Condenser")
        ORC.add_component(Pump, "Pump") # :!\ Est-ce que les noms des composants sont importants?
        ORC.add_component(Evaporator, "Evaporator")
    
        ORC.add_component(Spliter_expander, "Spliter")
    
        # Link components+
        ORC.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
        ORC.link_components("Evaporator", "m-ex_C", "Spliter", "m-su")
    
        ORC.link_components("Spliter", "m-ex_1", "Expander_1", "m-su")
        ORC.link_components("Spliter", "m-ex_2", "Expander_2", "m-su")
        ORC.link_components("Spliter", "m-ex_3", "Expander_3", "m-su")
    
        ORC.link_components("Expander_1", "m-ex", "Mixer", "m-su_1")
        ORC.link_components("Expander_2", "m-ex", "Mixer", "m-su_2")
        ORC.link_components("Expander_3", "m-ex", "Mixer", "m-su_3")
    
        ORC.link_components("Mixer", "m-ex", "Condenser", "m-su_H")
        ORC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    #%% SOURCES AND SINKS
    
        Hot_oil_source = MassConnector()
        ORC.add_source("EV_Water", Hot_oil_source, ORC.components["Evaporator"], "m-su_H")
        ORC.set_source_properties(T=90 + 273.15, fluid='Water', m_dot=3, target='EV_Water', P = 2e5)
        
        CD_water_source = MassConnector()
        ORC.add_source("CD_Water", CD_water_source, ORC.components["Condenser"], "m-su_C")
        ORC.set_source_properties(T=15 + 273.15, fluid='Water', m_dot=2.6, target='CD_Water', P = 2e5)
    
    #%% CYCLE GUESSES
        
        P_low = 2*1e5 # 119079.39300547051
    
        ORC.set_cycle_guess(target='Pump:su', m_dot = 0.3, SC = 5, p = P_low)
        
        ORC.set_cycle_guess(target='Spliter:su', T = 340, m_dot = 0.3, p = 4*1e5)
        
        ORC.set_cycle_guess(target='Expander_1:ex', p = P_low)
    
        ORC.set_cycle_guess(target='Expander_2:ex', p = P_low)
        
        ORC.set_cycle_guess(target='Expander_3:ex', p = P_low)
        
    #%% CYCLE FIXED PROPERTIES
        
        ORC.set_fixed_properties(target='Pump:su', SC = 5)
    
    #%% CYCLE ITERATION VARIABLES
    
        ORC.set_iteration_variable(target=['Expander_1:ex','Expander_2:ex','Expander_3:ex'], variable='p', objective = 'Pump:su-SC', tol = 1e-2, rel = 1)
    
    #%% CYCLE RESIDUAL VARIABLES
        ORC.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
        ORC.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)
        ORC.set_residual_variable(target='Evaporator:ex_C', variable='m_dot', tolerance= 1e-3)
        ORC.set_residual_variable(target='Condenser:ex_H', variable='m_dot', tolerance= 1e-3)
    
        ORC.solve()
    
    if case_study == 'Recuperator':
        ORC = Circuit('R1233zd(E)')
        
        # Create components
        Expander_1 = ExpanderSE()
        Expander_2 = ExpanderSE()
        Expander_3 = ExpanderSE()
        
        Spliter_expander = Spliter(outlet_repartition = [1/3,1/3,1/3])
        Mixer_expander = Mixer(n_inlets = 3)
        
        Condenser = HeatExchangerMB('Plate')
        Recuperator = HeatExchangerMB('Plate')
        Pump = PumpPolyEff()
        Evaporator = HeatExchangerMB('Plate')
    
    #%% CONDENSER PARAMETERS
    
        # Condenser
        condenser_geom = PlateGeomSWEP()
        condenser_geom.set_parameters("P200THx140/1P_Condenser")
    
        Condenser.set_parameters(
            # Set the geometry of the condenser
            A_c=condenser_geom.A_c, A_h=condenser_geom.A_h, h=condenser_geom.h, l=condenser_geom.l, l_v=condenser_geom.l_v, w_v=condenser_geom.w_v,
            C_CS=condenser_geom.C_CS, C_Dh=condenser_geom.C_Dh, C_V_tot=condenser_geom.C_V_tot, C_canal_t=condenser_geom.C_canal_t, C_n_canals=condenser_geom.C_n_canals,
            H_CS=condenser_geom.H_CS, H_Dh=condenser_geom.H_Dh, H_V_tot=condenser_geom.H_V_tot, H_canal_t=condenser_geom.H_canal_t, H_n_canals=condenser_geom.H_n_canals,
            casing_t=condenser_geom.casing_t, chevron_angle=condenser_geom.chevron_angle, fooling=condenser_geom.fooling,
            n_plates=condenser_geom.n_plates, plate_cond=condenser_geom.plate_cond, plate_pitch_co=condenser_geom.plate_pitch_co, t_plates=condenser_geom.t_plates, w=condenser_geom.w,
            amplitude=condenser_geom.amplitude, phi=condenser_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)
    
        Corr_H = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_cond_BPHEX"}
        Corr_C = {"1P" : "water_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}
    
        # Set the pressure drop correlations of the condenser
        Condenser.set_DP()
        
        # Set the heat transfer coefficients correlations of the condenser           
        Condenser.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)
        
    #%% EVAPORATOR PARAMETERS
    
        # Evaporator
        evaporator_geom = PlateGeomSWEP()
        evaporator_geom.set_parameters("P200THx140/1P_Evaporator")
    
        Evaporator.set_parameters(
            # Set the geometry of the evaporator
            A_c=evaporator_geom.A_c, A_h=evaporator_geom.A_h, h=evaporator_geom.h, l=evaporator_geom.l, l_v=evaporator_geom.l_v, w_v=evaporator_geom.w_v,
            C_CS=evaporator_geom.C_CS, C_Dh=evaporator_geom.C_Dh, C_V_tot=evaporator_geom.C_V_tot, C_canal_t=evaporator_geom.C_canal_t, C_n_canals=evaporator_geom.C_n_canals, 
            H_CS=evaporator_geom.H_CS, H_Dh=evaporator_geom.H_Dh, H_V_tot=evaporator_geom.H_V_tot, H_canal_t=evaporator_geom.H_canal_t, H_n_canals=evaporator_geom.H_n_canals,
            casing_t=evaporator_geom.casing_t, chevron_angle=evaporator_geom.chevron_angle, fooling=evaporator_geom.fooling, 
            n_plates=evaporator_geom.n_plates, plate_cond=evaporator_geom.plate_cond, plate_pitch_co=evaporator_geom.plate_pitch_co, t_plates=evaporator_geom.t_plates, w=evaporator_geom.w, 
            amplitude=evaporator_geom.amplitude, phi=evaporator_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)
        
        Corr_H = {"1P" : "water_plate_HTC", "2P" : "Han_cond_BPHEX"}
        Corr_C = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}
    
        # Set the pressure drop correlations of the condenser
        Evaporator.set_DP()
        
        # Set the heat transfer coefficients correlations of the condenser           
        Evaporator.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)
    
    #%% EVAPORATOR PARAMETERS
    
        # Evaporator
        recuperator_geom = PlateGeomSWEP()
        recuperator_geom.set_parameters("P200THx140/1P_Evaporator")
    
        Recuperator.set_parameters(
            # Set the geometry of the evaporator
            A_c=recuperator_geom.A_c, A_h=recuperator_geom.A_h, h=recuperator_geom.h, l=recuperator_geom.l, l_v=recuperator_geom.l_v, w_v=recuperator_geom.w_v,
            C_CS=recuperator_geom.C_CS, C_Dh=recuperator_geom.C_Dh, C_V_tot=recuperator_geom.C_V_tot, C_canal_t=recuperator_geom.C_canal_t, C_n_canals=recuperator_geom.C_n_canals, 
            H_CS=recuperator_geom.H_CS, H_Dh=recuperator_geom.H_Dh, H_V_tot=recuperator_geom.H_V_tot, H_canal_t=recuperator_geom.H_canal_t, H_n_canals=recuperator_geom.H_n_canals,
            casing_t=recuperator_geom.casing_t, chevron_angle=recuperator_geom.chevron_angle, fooling=recuperator_geom.fooling, 
            n_plates=recuperator_geom.n_plates, plate_cond=recuperator_geom.plate_cond, plate_pitch_co=recuperator_geom.plate_pitch_co, t_plates=recuperator_geom.t_plates, w=recuperator_geom.w, 
            amplitude=recuperator_geom.amplitude, phi=recuperator_geom.phi, Flow_Type='CounterFlow', H_DP_ON=True, C_DP_ON=True, n_disc=0)
        
        Corr_H = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_cond_BPHEX"}
        Corr_C = {"1P" : "martin_holger_plate_HTC", "2P" : "Han_Boiling_BPHEX_HTC"}
    
        # Set the pressure drop correlations of the condenser
        Recuperator.set_DP()
        
        # Set the heat transfer coefficients correlations of the condenser           
        Recuperator.set_htc(htc_type = 'Correlation', Corr_H = Corr_H, Corr_C = Corr_C)
    
    #%% EXPANDER PARAMETERS AND INPUTS
    
        N_exp = 6000
        T_amb = 293
    
        Expander_1.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)
    
        Expander_1.set_inputs(N_rot=N_exp, T_amb=T_amb)
    
        Expander_2.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)
    
        Expander_2.set_inputs(N_rot=N_exp, T_amb=T_amb)
        
        Expander_3.set_parameters(AU_amb=8.33758799e+00, AU_su_n=6.67152053e-01, AU_ex_n=3.21181352e+01, d_su1=6.31789061e-03, m_dot_n=0.1, 
                    A_leak=1.00000000e-10, W_dot_loss_0=8.19123951e-01, alpha= 7.79756524e-02, C_loss=4.68294054e-01, rv_in=1.7, V_s=0.0000712)
    
        Expander_3.set_inputs(N_rot=N_exp, T_amb=T_amb)
    
    #%% PUMP INPUTS
    
        Pump.set_inputs(N_pp = 2500)
    
    #%% ADD AND LINK COMPONENTS
        ORC.add_component(Expander_1, "Expander_1")
        ORC.add_component(Expander_2, "Expander_2")
        ORC.add_component(Expander_3, "Expander_3")
    
        ORC.add_component(Mixer_expander, "Mixer")
        
        ORC.add_component(Condenser, "Condenser")
        ORC.add_component(Recuperator, "Recuperator")
        ORC.add_component(Pump, "Pump") # :!\ Est-ce que les noms des composants sont importants?
        ORC.add_component(Evaporator, "Evaporator")
    
        ORC.add_component(Spliter_expander, "Spliter")
    
        # Link components+
        ORC.link_components("Pump", "m-ex", "Recuperator", "m-su_C")
        ORC.link_components("Recuperator", "m-ex_C", "Evaporator", "m-su_C")
        ORC.link_components("Evaporator", "m-ex_C", "Spliter", "m-su")
    
        ORC.link_components("Spliter", "m-ex_1", "Expander_1", "m-su")
        ORC.link_components("Spliter", "m-ex_2", "Expander_2", "m-su")
        ORC.link_components("Spliter", "m-ex_3", "Expander_3", "m-su")
    
        ORC.link_components("Expander_1", "m-ex", "Mixer", "m-su_1")
        ORC.link_components("Expander_2", "m-ex", "Mixer", "m-su_2")
        ORC.link_components("Expander_3", "m-ex", "Mixer", "m-su_3")
    
        ORC.link_components("Mixer", "m-ex", "Recuperator", "m-su_H")
        ORC.link_components("Recuperator", "m-ex_H", "Condenser", "m-su_H")
        ORC.link_components("Condenser", "m-ex_H", "Pump", "m-su")
    
    #%% SOURCES AND SINKS
    
        Hot_oil_source = MassConnector()
        ORC.add_source("EV_Water", Hot_oil_source, ORC.components["Evaporator"], "m-su_H")
        ORC.set_source_properties(T=90 + 273.15, fluid='Water', m_dot=3, target='EV_Water', P = 2e5)
        
        CD_water_source = MassConnector()
        ORC.add_source("CD_Water", CD_water_source, ORC.components["Condenser"], "m-su_C")
        ORC.set_source_properties(T=15 + 273.15, fluid='Water', m_dot=2.6, target='CD_Water', P = 2e5)
    
    #%% CYCLE GUESSES
    
        P_low = 2*1e5
    
        ORC.set_cycle_guess(target='Pump:su', m_dot = 0.3, SC = 3, p = P_low)
    
        ORC.set_cycle_guess(target='Spliter:su', T = 340, m_dot = 0.3, p = 4*1e5)
        
        ORC.set_cycle_guess(target='Expander_1:ex', p = P_low)
    
        ORC.set_cycle_guess(target='Expander_2:ex', p = P_low)
        
        ORC.set_cycle_guess(target='Expander_3:ex', p = P_low)

    #%% CYCLE FIXED VARIABLES AND ITERATION VARIABLE
    
        ORC.set_fixed_properties(target='Pump:su', SC = 5)    
        ORC.set_iteration_variable(target=['Expander_1:ex','Expander_2:ex','Expander_3:ex'], variable='p', objective = 'Pump:su-SC', tol = 1e-2, rel = 1, damping_factor = 0.1)
    
    #%% CYCLE RESIDUAL VARIABLES
        ORC.set_residual_variable(target='Evaporator:ex_C', variable='h', tolerance= 1e-3)
        ORC.set_residual_variable(target='Condenser:ex_H', variable='h', tolerance= 1e-3)
        ORC.set_residual_variable(target='Evaporator:ex_C', variable='m_dot', tolerance= 1e-3)
        ORC.set_residual_variable(target='Condenser:ex_H', variable='m_dot', tolerance= 1e-3)
    
        ORC.solve()

    # ORC.print_states()
