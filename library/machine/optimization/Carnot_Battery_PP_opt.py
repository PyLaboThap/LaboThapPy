from CoolProp.CoolProp import PropsSI

"""
ZORLU GEOTHERMAL POWER PLANT
"""

"0) Global paremeters"

E_density_PCM = 90 # [kWh/m^3] Energy density of the phase change storage
T_PCM = 150 + 273.15 # [K] PCM Temperature
wf = 'Cyclopentane' # Working fluid  
 
"Zorlu Requierments : Proposal information"

Q_dot_cd_HP = 7000 # [kW]
W_dot_exp_orc = 2000 # [kW]
PCM_cap = 29000 # [kWh]
t_charge = 4 # [h]
t_discharge = 2 # [h]

#%%
"1) High Temperature Heat Pump"
  
"### Parameters ####"
DELTAT_sh_HP = 1 # [K] 
DELTAT_sc_HP = 3 # [K]

epsilon_s_cp = 0.65 # [-]
epsilon_v_cp = 0.95 # [-]
N_cp = 50 # [Hz] 

" Evaporator secondary fluid, several options n°8, n°58, for the n°8 T_res = 85°C "
sf = 'Water'
T_sf_su_ev = 113.1 + 273.15 # [K]
p_sf_su_ev = 159*1e3 # [Pa]
m_dot_sf_su_ev = (3088.6*1000/3600) # [kg/s]
 
"! #### CONDENSER - HP #### "
" Connections "
T_r_su_cd_HP = T_r_ex_cp_HP
 
"Refrigerant phase-change "
T_r_cd_x1_HP = temperature(wf_fluid$; x=1; P=p_sat_cd_HP)
T_r_cd_x0_HP = temperature(wf_fluid$; x=0; P=p_sat_cd_HP)
T_r_sat_cd_HP = (T_r_cd_x1_HP + T_r_cd_x0_HP)/2
 
"Refrigerants properties"
T_r_ex_cd_HP = T_r_cd_x0_HP - DELTAT_sc_HP
h_r_ex_cd_HP = enthalpy(wf_fluid$; T =T_r_ex_cd_HP; P = p_sat_cd_HP)
h_r_x0_cd_HP = enthalpy(wf_fluid$; x = 0; P = p_sat_cd_HP)
h_r_x1_cd_HP = enthalpy(wf_fluid$; x =1; P = p_sat_cd_HP)
h_r_su_cd_HP = enthalpy(wf_fluid$; T =T_r_su_cd_HP; P = p_sat_cd_HP) 
 
" Pinch-point "
DELTAT_pp_cd_HP = T_r_ex_cd_HP - T_PCM    " For this case the pintch point is in the exhaust of the condenser "
{p_sat_cd_HP = 1200}
 
" Energy Balance "
Q_dot_cd_HP = m_dot_wf_hp*(h_r_su_cd_HP-h_r_ex_cd_HP)
Q_dot_sh_cd_HP = m_dot_wf_hp*(h_r_su_cd_HP-h_r_x1_cd_HP)
Q_dot_tp_cd_HP = m_dot_wf_hp*(h_r_x1_cd_HP- h_r_x0_cd_HP)
 
"! #### COMPRESSOR - HP #### "
"Connections: " 
T_r_su_cp_HP = T_r_ex_cold_ihx_hp
p_r_su_cp_HP = p_sat_ev_HP
p_r_ex_cp_HP = p_sat_cd_HP
 
"Refrigerant properties"
h_r_su_cp_HP = enthalpy(wf_fluid$; T = T_r_su_cp_HP; P=p_r_su_cp_HP)
s_r_su_cp_HP =  entropy(wf_fluid$; T = T_r_su_cp_HP; P=p_r_su_cp_HP)
v_r_su_cp_HP =  volume(wf_fluid$; T = T_r_su_cp_HP; P=p_r_su_cp_HP)
T_r_ex_cp_HP = temperature(wf_fluid$; P = p_sat_cd_HP; h = h_r_ex_cp_HP)
  
"Issentropic compression"
h_r_ex_s_cp_HP = enthalpy(wf_fluid$; P = p_sat_cd_HP; s = s_r_su_cp_HP)
epsilon_s_cp = (h_r_ex_s_cp_HP - h_r_su_cp_HP)/(h_r_ex_cp_HP - h_r_su_cp_HP)
epsilon_v_cp = m_dot_wf_hp*v_r_su_cp_HP/V_dot_s_cp
V_dot_s_cp = V_s_cp*N_cp 
 
"Energy balance"
W_dot_cp_HP= m_dot_wf_hp*(h_r_ex_cp_HP - h_r_su_cp_HP)
 
"! #### INTERNAL -HX  #### " 
"epsilon_ihx_hp = 0,8"
 
" Connections "
h_r_su_hot_ihx_hp = h_r_ex_cd_HP
p_r_su_hot_ihx_hp = p_sat_cd_HP
T_r_su_hot_ihx_hp = T_r_ex_cd_HP
 
h_r_su_cold_ihx_hp = h_r_ex_ev_HP
p_r_su_cold_ihx_hp = p_sat_ev_HP
T_r_su_cold_ihx_hp = T_r_ex_ev_HP
 
" Thermodynamic refrigerants" 
h_r_ex_hot_ihx_hp = enthalpy(wf_fluid$; T =T_r_ex_hot_ihx_hp; P = p_r_su_hot_ihx_hp) 
h_r_ex_cold_ihx_hp = enthalpy(wf_fluid$; T =T_r_ex_cold_ihx_hp; P = p_r_su_cold_ihx_hp) 
x_r_ex_hot_ihx_hp = quality(wf_fluid$;h=h_r_ex_hot_ihx_hp; P=p_r_su_hot_ihx_hp)
 
cp_hot_ihx_hp = cp(wf_fluid$; h = h_r_su_hot_ihx_hp; P=p_r_su_hot_ihx_hp)
cp_cold_ihx_hp = cp(wf_fluid$; h = h_r_su_cold_ihx_hp; P=p_r_su_cold_ihx_hp)
C_dot_hot_ihx_hp = cp_hot_ihx_hp*m_dot_wf_hp
C_dot_cold_ihx_hp = cp_cold_ihx_hp*m_dot_wf_hp
C_dot_min_ihx_hp = min(C_dot_hot_ihx_hp; C_dot_cold_ihx_hp)
 
" Energy balance"
Q_dot_ihx = m_dot_wf_hp*(h_r_su_hot_ihx_hp - h_r_ex_hot_ihx_hp)
Q_dot_ihx = m_dot_wf_hp*(h_r_ex_cold_ihx_hp - h_r_su_cold_ihx_hp)
 
" Energy Heat Transfer"
Q_dot_ihx = C_dot_min_ihx_hp*epsilon_ihx_hp*(T_r_su_hot_ihx_hp - T_r_su_cold_ihx_hp)
 
"! #### EVAPORATOR - HP #### " 
" Connections: "
h_r_su_ev_HP = h_r_ex_hot_ihx_hp
 
" Saturation conditions"
x_r_su_ev_HP = quality(wf_fluid$;h=h_r_su_ev_HP; P=p_sat_ev_HP)
T_r_su_ev_HP = temperature(wf_fluid$;h=h_r_su_ev_HP; P=p_sat_ev_HP)
 
{x_r_su_ev_HP = 0.04}
T_r_ev_x1_HP =  temperature(wf_fluid$; x=1; P=p_sat_ev_HP)
h_r_ev_x1_HP = enthalpy(wf_fluid$; x=1; P=p_sat_ev_HP)
T_r_ev_xsuev_HP =  temperature(wf_fluid$; x=x_r_su_ev_HP; P=p_sat_ev_HP)
T_r_sat_ev_HP = (T_r_ev_x1_HP + T_r_ev_xsuev_HP)/2
 
"Refrigerant conditions "
T_r_ex_ev_HP = T_r_ev_x1_HP + DELTAT_sh_HP
h_r_ex_ev_HP = enthalpy(wf_fluid$; T =T_r_ex_ev_HP; P = p_sat_ev_HP)
 
"Secondary fluid "
T_sf_sat_ev_HP = t_sat(sf_fluid$; P = p_sf_su_ev)
h_sf_su_ev_HP = enthalpy(sf_fluid$; T = T_sf_su_ev; P = p_sf_su_ev)
v_sf_su_ev_HP =volume(sf_fluid$; T = T_sf_su_ev; P = p_sf_su_ev)
h_sf_ex_ev_HP = enthalpy(sf_fluid$; T = T_sf_ex_ev_HP; P = p_sf_su_ev)
 
"Pinch-point"
DELTAT_pp_ev_HP = T_sf_ex_ev_HP - T_r_ev_xsuev_HP
{p_sat_ev_HP = 160 }
 
"Energy balance"
Q_dot_ev_HP = m_dot_wf_hp*(h_r_ex_ev_HP - h_r_su_ev_HP)
Q_dot_ev_HP = m_dot_sf_su_ev*(h_sf_su_ev_HP - h_sf_ex_ev_HP)
Q_dot_hf_ev_HP =  m_dot_wf_hp*(h_r_ev_x1_HP - h_r_su_ev_HP)
 
"! ########## Global coefficients ######### "
r_p_HP = p_sat_cd_HP/p_sat_ev_HP
COP_HP = Q_dot_cd_HP/W_dot_cp_HP 
 
"########################################"
"!2.-                 Phace Change Storage                        "
 "########################################"
Q_PCM = Q_dot_cd_HP
E_PCM = Q_PCM*t_charge
E_loss = 0,05*E_PCM
E_PCM_f = E_PCM-E_loss
Vol_PCM = E_PCM/E_density_PCM
Q_loss_PCM = (E_PCM*0,05)/(t_charge + t_discharge)                  "Heat losses 5% of the total energy per day"
Q_dis_PCM = E_PCM_f/t_discharge 
 
 
"########################################"
"!3.-                 Organic Rankine Cycle                        "
 "########################################"
"Secondary fluid for ORC could be water as well as air"
sf_fluid_orc$ = 'Water'
T_sf_su_cd_orc = 24 [C]
p_sf_su_cd_orc = 101,3 [kPa]
m_dot_sf_su_cd = 2000 [kg/s]
 
"### Parameters ####"
DELTAT_sh_orc = 3   [K]
DELTAT_sc_orc = 3  [K]
{DELTAT_pp_cd_orc = 3 [K]
DELTAT_pp_ev_orc= 3 [K]}
epsilon_s_exp = 0,65
epsilon_v_exp = 0,95
N_exp = 50 [Hz]
epsilon_s_pp = 0,75
 
"! #### EVAPORATOR - ORC #### " 
"Connection:"
T_r_su_ev_orc =  T_r_ex_cold_rec
Q_dot_ev_orc = Q_dis_PCM 
 
"Refrigerant phase-change "
T_r_ev_orc_x1 = temperature(wf_fluid$; x=1; P=p_sat_ev_orc)
T_r_ev_orc_x0 = temperature(wf_fluid$; x=0; P=p_sat_ev_orc)
T_r_sat_ev_orc = (T_r_ev_orc_x1 + T_r_ev_orc_x0)/2
 
"Refrigerants properties"
T_r_ex_ev_orc = T_r_ev_orc_x1  +  DELTAT_sh_orc
h_r_ex_ev_orc = enthalpy(wf_fluid$; T =T_r_ex_ev_orc; P = p_sat_ev_orc)
h_r_x0_ev_orc = enthalpy(wf_fluid$; x = 0; P = p_sat_ev_orc)
h_r_x1_ev_orc = enthalpy(wf_fluid$; x =1; P = p_sat_ev_orc)
h_r_su_ev_orc = enthalpy(wf_fluid$; T =T_r_su_ev_orc; P = p_sat_ev_orc) 
 
" Pinch-point "
DELTAT_pp_ev_orc = T_PCM - T_r_ex_ev_orc
 
" Energy Balance "
Q_dot_ev_orc = m_dot_wf_orc*(h_r_ex_ev_orc - h_r_su_ev_orc) 
 
"! ####           PUMP - ORC        #### " 
 "Connections:"
T_r_su_pp = T_r_ex_cd_orc
p_r_su_pp = p_sat_cd_orc
p_r_ex_pp = p_sat_ev_orc
 
"Refrigerant properties"
h_r_su_pp = enthalpy(wf_fluid$; T=T_r_su_pp; P=p_r_su_pp)
v_r_su_pp = volume(wf_fluid$; T=T_r_su_pp; P=p_r_su_pp)
h_r_ex_pp = enthalpy(wf_fluid$; T=T_r_ex_pp; P=p_r_ex_pp)
 
"Refrigerant properties"
h_r_ex_pp = h_r_su_pp + v_r_su_pp*(p_r_ex_pp - p_r_su_pp)/epsilon_s_pp
 
"Energy balances"
W_dot_pp = m_dot_wf_orc*(h_r_ex_pp - h_r_su_pp)
 
"! ####   CONDENSER - ORC  #### " 
"Connections:"
T_r_su_cd_orc =T_r_ex_hot_rec
 
" Saturation conditions"
T_r_cd_orc_x1 =  temperature(wf_fluid$; x=1; P=p_sat_cd_orc)
T_r_cd_orc_x0 =  temperature(wf_fluid$; x=0; P=p_sat_cd_orc)
T_r_sat_cd_orc = (T_r_cd_orc_x1 + T_r_cd_orc_x0)/2
h_r_x0_cd_orc = enthalpy(wf_fluid$; x = 0; P = p_sat_cd_orc)
h_r_x1_cd_orc = enthalpy(wf_fluid$; x =1; P = p_sat_cd_orc)
 
 "Refrigerant conditions "
T_r_ex_cd_orc = T_r_cd_orc_x0 - DELTAT_sc_orc
h_r_ex_cd_orc = enthalpy(wf_fluid$; T =T_r_ex_cd_orc; P = p_sat_cd_orc)
h_r_su_cd_orc = enthalpy(wf_fluid$; T = T_r_su_cd_orc; P = p_sat_cd_orc)
x_r_su_cd_orc = quality(wf_fluid$;h=h_r_su_cd_orc; P=p_sat_cd_orc)
 
"Secondary fluid"
h_sf_su_cd_orc = enthalpy(sf_fluid_orc$; T = T_sf_su_cd_orc; P = p_sf_su_cd_orc) 
h_sf_pp_cd_orc = enthalpy(sf_fluid_orc$; T = T_sf_pp_cd_orc; P = p_sf_su_cd_orc)
h_sf_ex_cd_orc = enthalpy(sf_fluid_orc$; T = T_sf_ex_cd_orc; P = p_sf_su_cd_orc)
 
" Pinch-point "
DELTAT_pp_cd_orc = T_r_cd_orc_x1 - T_sf_pp_cd_orc 
Q_dot_pp_cd_orc = m_dot_wf_orc*( h_r_x1_cd_orc- h_r_ex_cd_orc ) 
Q_dot_pp_cd_orc = m_dot_sf_su_cd*( h_sf_pp_cd_orc - h_sf_su_cd_orc ) 
 
" Energy balance"
Q_dot_cd_orc = m_dot_wf_orc*(h_r_su_cd_orc - h_r_ex_cd_orc) 
Q_dot_cd_orc = m_dot_sf_su_cd*(h_sf_ex_cd_orc - h_sf_su_cd_orc ) 
 
"! ####       TURBINE - ORC     #### " 
"Connections:"
T_r_su_exp = T_r_ex_ev_orc
p_r_su_exp = p_sat_ev_orc
p_r_ex_exp = p_sat_cd_orc 
 
"Refrigerant properties"
h_r_su_exp = enthalpy(wf_fluid$; P=p_r_su_exp; T=T_r_su_exp)
s_r_su_exp = entropy(wf_fluid$; P=p_r_su_exp; T=T_r_su_exp)
v_r_su_exp = volume(wf_fluid$; P=p_r_su_exp; T=T_r_su_exp)
h_r_ex_s_exp = enthalpy(wf_fluid$; P=p_r_ex_exp; s=s_r_su_exp)
h_r_ex_exp = enthalpy(wf_fluid$; P=p_r_ex_exp; T=T_r_ex_exp) 
x_r_ex_exp = quality(wf_fluid$;h=h_r_ex_exp; P=p_r_ex_exp)
 
" Efficacy" 
epsilon_s_exp =  (h_r_su_exp - h_r_ex_exp)/(h_r_su_exp - h_r_ex_s_exp)
epsilon_v_exp = m_dot_wf_orc*v_r_su_exp/V_dot_s_exp
V_dot_s_exp = V_s_exp*N_exp
 
" Energy balance"
W_dot_exp = m_dot_wf_orc*(h_r_su_exp - h_r_ex_exp)
 
 
"! #### RECUPERATOR - ORC #### " 
 "epsilon_rec = 0,8"
 
"Connection:"
T_r_su_hot_rec = T_r_ex_exp
T_r_su_cold_rec = T_r_ex_pp
p_r_cold_rec = p_r_ex_pp
p_r_hot_rec = p_r_ex_exp
h_r_su_cold_rec =  h_r_ex_pp
h_r_su_hot_rec = h_r_ex_exp
 
"Refrigerant properties"
T_cold_bar_rec = (T_r_su_cold_rec + T_r_ex_cold_rec)/2
T_hot_bar_rec = (T_r_su_hot_rec + T_r_ex_hot_rec)/2
cp_cold_rec =  cp(wf_fluid$;T=T_cold_bar_rec;P=p_r_cold_rec)
cp_hot_rec =  cp(wf_fluid$;T=T_hot_bar_rec;P=p_r_hot_rec)
C_dot_cold_rec = cp_cold_rec*m_dot_wf_orc
C_dot_hot_rec = cp_hot_rec*m_dot_wf_orc
C_dot_min = min(C_dot_cold_rec; C_dot_hot_rec) 
h_r_ex_cold_rec = enthalpy(wf_fluid$; T=T_r_ex_cold_rec; P=p_r_cold_rec)
h_r_ex_hot_rec = enthalpy(wf_fluid$; T=T_r_ex_hot_rec; P=p_r_hot_rec)
x_r_ex_cold_rec = quality(wf_fluid$;h=h_r_ex_cold_rec; P=p_r_cold_rec)
 
"Energy balance"
Q_dot_rec = m_dot_wf_orc*(h_r_ex_cold_rec - h_r_su_cold_rec)
Q_dot_rec = m_dot_wf_orc*(h_r_su_hot_rec - h_r_ex_hot_rec)
 
"Heat transfer equation" 
Q_dot_rec = C_dot_min*epsilon_rec*(T_r_su_hot_rec - T_r_su_cold_rec)
 
 
"! ########## Global coefficients ######### "
r_p_orc = p_sat_ev_orc/p_sat_cd_orc
eta_orc = (W_dot_exp-W_dot_pp)/(Q_dot_ev_orc)
RTE = eta_orc*COP_HP
 
"! ########## Heat Exchanger Areas Computation ######### "
 
"HP EVAP"
 
DTA_ev_HP = (T_sf_su_ev - T_r_ex_ev_HP)
DTB_ev_HP = (T_sf_ex_ev_HP - T_r_sat_ev_HP)
 
LMTD_ev_HP = (DTA_ev_HP - DTB_ev_HP)/ln(DTA_ev_HP/DTB_ev_HP)
 
UA_ev_HP = 1000*Q_dot_ev_HP/LMTD_ev_HP
U_ev_HP = 1500
 
A_ev_HP = UA_ev_HP/U_ev_HP
 
"HP COND"
 
DTA_cd_HP = (T_r_su_cd_HP - T_PCM)
DTB_cd_HP = (T_r_ex_cd_HP - T_PCM)
 
LMTD_cd_HP = (DTA_cd_HP - DTB_cd_HP)/ln(DTA_cd_HP/DTB_cd_HP)
 
UA_cd_HP = 1000*Q_dot_cd_HP/LMTD_cd_HP
U_cd_HP = 2500
 
A_cd_HP = UA_cd_HP/U_cd_HP
 
"HP IHX"
 
DTA_ihx = (T_r_su_hot_ihx_hp - T_r_ex_cold_ihx_hp) 
DTB_ihx = (T_r_ex_hot_ihx_hp - T_r_su_cold_ihx_hp)
 
LMTD_ihx_HP = (DTA_ihx - DTB_ihx)/ln(DTA_ihx/DTB_ihx)
 
UA_ihx_HP = 1000*Q_dot_ihx/LMTD_ihx_HP
U_ihx_HP = 600
 
A_ihx_HP = UA_ihx_HP/U_ihx_HP
 
"ORC EVAP"
 
DTA_ev_orc = (T_PCM - T_r_ex_ev_orc)
DTB_ev_orc = (T_PCM - T_r_su_ev_orc)
 
LMTD_ev_orc = (DTA_ev_orc - DTB_ev_orc)/ln(DTA_ev_orc/DTB_ev_orc)
 
UA_ev_orc = 1000*Q_dot_ev_orc/LMTD_ev_orc
U_ev_orc = 1500
 
A_ev_orc = UA_ev_orc/U_ev_orc
 
"ORC COND"
 
DTA_cd_orc = (T_r_su_cd_orc - T_sf_ex_cd_orc)
DTB_cd_orc = (T_r_ex_cd_orc - T_sf_su_cd_orc)
 
LMTD_cd_orc = (DTA_cd_orc - DTB_cd_orc)/ln(DTA_cd_orc/DTB_cd_orc)
 
UA_cd_orc = 1000*Q_dot_cd_orc/LMTD_cd_orc
U_cd_orc = 2500
 
A_cd_orc = UA_cd_orc/U_cd_orc
 
"ORC REC"
 
DTA_rec = (T_r_su_hot_rec - T_r_ex_cold_rec)
DTB_rec = (T_r_ex_hot_rec - T_r_su_cold_rec)
 
LMTD_rec_orc = (DTA_rec - DTB_rec)/ln(DTA_rec/DTB_rec)
 
UA_rec_orc = 1000*Q_dot_rec/LMTD_rec_orc
U_rec_orc = 600
 
A_rec_orc = UA_rec_orc/U_rec_orc
 
"Total HTX A"
 
A_tot_HTX = A_ev_HP + A_cd_HP + A_ihx_HP + A_ev_orc + A_cd_orc + A_rec_orc
