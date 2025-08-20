# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 13:52:57 2025

@author: Basile
"""

import numpy as np
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize

from pyswarms.single.global_best import GlobalBestPSO
from tqdm import tqdm

#%%

class PressuzedWaterStorageDesign(object):

    def __init__(self, fluid):
        # Inputs
        self.inputs = {}
        
        # Params
        self.params = {}

        # Abstract State
        self.fluid = fluid
        self.AS = CP.AbstractState('HEOS', fluid)
        
        # Gravity Constant
        self.g = 9.81 # m/s^2
        
        return

    # ---------------- Data Handling ----------------------------------------------------------------------
    
    def set_inputs(self, **parameters):
        for key, value in parameters.items():
            self.inputs[key] = value
            
    def set_parameters(self, **parameters):
            for key, value in parameters.items():
                self.params[key] = value

    # ---------------- Computation Related Methods ----------------------------------------------------------------------

    def req_thickness(self, P_des, R_i, sigma_max):
        """
        Inputs:
        -------------
        P_des     : Design pressure [Pa]
        sigma_max : Maximum Allowable tensile stress of the material [Pa]
        R_i         : For pipe, the inside radius R is determined by the nominal outside radius minus the nominal wall thickness [m]
        
        Reference 
        -------------
        BVPC 2007 
        """
                
        E_joint = 1 #0.85 # Joint efficiency : Butt joints welded from both sides, or from one side only without a permanent
                       # backing strip, that are verified as having achieved full penetration and fusion as
                       # required by UW-35. Case b 
                
        t_req_circ = (P_des*R_i)/(2*sigma_max*E_joint - 1.2*P_des) # Circumferential Stress (Longitudinal Joints)
        
        # print(f"t_req_P : {max(t_req_circ)}")
        # print(f"t_req_P : {t_req_circ}")
        # print(f"t_req_min : {self.t_req_min}")
        
        # return max(max(t_req_circ, t_req_long, t_head)/QF, self.t_req_min)
        return  max(t_req_circ, self.t_req_min)
    
    def material_properties(self):
        """
        CHEMICAL ENGINEERING DESIGN PRINCIPLES, PRACTICE AND ECONOMICS
        OF PLANT AND PROCESS DESIGN : THIRD EDITION - Table 14.2
        
        Gavin Towler Ray Sinnott
        """
        def linear_interp_extrap(x_data, y_data, x):
            # If x is less than the first point, extrapolate using the first two points
            if x < x_data[0]:
                x0, x1 = x_data[0], x_data[1]
                y0, y1 = y_data[0], y_data[1]
            # If x is greater than the last point, extrapolate using the last two points
            elif x > x_data[-1]:
                x0, x1 = x_data[-2], x_data[-1]
                y0, y1 = y_data[-2], y_data[-1]
            else:
                # Interpolation
                for i in range(len(x_data) - 1):
                    x0, x1 = x_data[i], x_data[i + 1]
                    y0, y1 = y_data[i], y_data[i + 1]
                    if x0 <= x <= x1:
                        break

            # Linear formula (works for both interpolation and extrapolation)
            y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            return y
        
        if self.params['Material'] == 'Carbon_Steel': # Grade : A285 Gr A
            sigma_max_vec = np.array([88945500, 88945500, 88945500, 79292500, 40680500]) # Maximum Allowable stress 
            T_vec = np.array([37.78, 148.89, 260, 371.11, 482.22]) + 273.15
            
            self.sigma_max = linear_interp_extrap(T_vec, sigma_max_vec, self.T_des)
            self.rho_mat = 7300 # kg/m^3
            self.t_req_min = 0.6*1e-3 # m
            
            if self.params['Arrang'] == 'Vertical':
                self.a = 11600
                self.b = 34
                self.n = 0.85
            
            elif self.params['Arrang'] == 'Horizontal':
                self.a = 10200
                self.b = 31
                self.n = 0.85
                
            else:
                print("Parameter 'Arrang' shall be either 'Vertical' or 'Horizontal'")
                
        elif self.params['Material'] == '304ss':
            sigma_max_vec = np.array([1.379e+8, 1.034e+8, 8.8942e+7, 8.0669e+7, 74463379]) # Maximum Allowable stress 
            T_vec = np.array([37.78, 148.89, 260, 371.11, 482.22]) + 273.15
            
            self.sigma_max = linear_interp_extrap(T_vec, sigma_max_vec, self.T_des)
            self.rho_mat = 7850 # kg/m^3
            self.t_req_min = 0.6*1e-3
        
            if self.params['Arrang'] == 'Vertical':
                self.a = 17400
                self.b = 79
                self.n = 0.85
            
            elif self.params['Arrang'] == 'Horizontal':
                self.a = 12800
                self.b = 73
                self.n = 0.85
                
            else:
                print("Parameter 'Arrang' shall be either 'Vertical' or 'Horizontal'")
        
        else: 
            raise ValueError(f"Unknown material: {self.params['Material']}")
        
        return    
    
    def capital_cost_est(self):
        """
        CHEMICAL ENGINEERING DESIGN PRINCIPLES, PRACTICE AND ECONOMICS
        OF PLANT AND PROCESS DESIGN : THIRD EDITION - Table 7.2
        
        Gavin Towler Ray Sinnott
        
        Cost in 2022 $
        """
        
        if self.params['Arrang'] == 'Vertical':
            if 250000 > self.m_mat > 10: # >160
                self.CAPEX = self.a + self.b*self.m_mat**self.n
                # print("Tank mass is not in the cost correlation bounds")
                
        elif self.params['Arrang'] == 'Horizontal':
            if 50000 > self.m_mat > 10: # >160
                self.CAPEX = self.a + self.b*self.m_mat**self.n
                # print("Tank mass is not in the cost correlation bounds")
        return
    
    def create_CAPEX_correlation(self):
        """
        From optimization of CAPEX per stored MWh : Iterating on D with fixed V_sto and T_sto for vertical tanks
        For Water Storage
        """
        # Optimization data : [T, V] couples
        points = np.array([[100.0, 1.0], [105.0, 1.0], [110.0, 1.0], [115.0, 1.0], [120.0, 1.0], [125.0, 1.0], [130.0, 1.0], [135.0, 1.0], [140.0, 1.0], [145.0, 1.0], [150.0, 1.0], [155.0, 1.0], [160.0, 1.0], [165.0, 1.0], [170.0, 1.0], [175.0, 1.0], [180.0, 1.0], [185.0, 1.0], [190.0, 1.0], [195.0, 1.0], [200.0, 1.0], [100.0, 11.642857142857142], [105.0, 11.642857142857142], [110.0, 11.642857142857142], [115.0, 11.642857142857142], [120.0, 11.642857142857142], [125.0, 11.642857142857142], [130.0, 11.642857142857142], [135.0, 11.642857142857142], [140.0, 11.642857142857142], [145.0, 11.642857142857142], [150.0, 11.642857142857142], [155.0, 11.642857142857142], [160.0, 11.642857142857142], [165.0, 11.642857142857142], [170.0, 11.642857142857142], [175.0, 11.642857142857142], [180.0, 11.642857142857142], [185.0, 11.642857142857142], [190.0, 11.642857142857142], [195.0, 11.642857142857142], [200.0, 11.642857142857142], [100.0, 22.285714285714285], [105.0, 22.285714285714285], [110.0, 22.285714285714285], [115.0, 22.285714285714285], [120.0, 22.285714285714285], [125.0, 22.285714285714285], [130.0, 22.285714285714285], [135.0, 22.285714285714285], [140.0, 22.285714285714285], [145.0, 22.285714285714285], [150.0, 22.285714285714285], [155.0, 22.285714285714285], [160.0, 22.285714285714285], [165.0, 22.285714285714285], [170.0, 22.285714285714285], [175.0, 22.285714285714285], [180.0, 22.285714285714285], [185.0, 22.285714285714285], [190.0, 22.285714285714285], [195.0, 22.285714285714285], [200.0, 22.285714285714285], [100.0, 32.92857142857143], [105.0, 32.92857142857143], [110.0, 32.92857142857143], [115.0, 32.92857142857143], [120.0, 32.92857142857143], [125.0, 32.92857142857143], [130.0, 32.92857142857143], [135.0, 32.92857142857143], [140.0, 32.92857142857143], [145.0, 32.92857142857143], [150.0, 32.92857142857143], [155.0, 32.92857142857143], [160.0, 32.92857142857143], [165.0, 32.92857142857143], [170.0, 32.92857142857143], [175.0, 32.92857142857143], [180.0, 32.92857142857143], [185.0, 32.92857142857143], [190.0, 32.92857142857143], [195.0, 32.92857142857143], [200.0, 32.92857142857143], [100.0, 43.57142857142857], [105.0, 43.57142857142857], [110.0, 43.57142857142857], [115.0, 43.57142857142857], [120.0, 43.57142857142857], [125.0, 43.57142857142857], [130.0, 43.57142857142857], [135.0, 43.57142857142857], [140.0, 43.57142857142857], [145.0, 43.57142857142857], [150.0, 43.57142857142857], [155.0, 43.57142857142857], [160.0, 43.57142857142857], [165.0, 43.57142857142857], [170.0, 43.57142857142857], [175.0, 43.57142857142857], [180.0, 43.57142857142857], [185.0, 43.57142857142857], [190.0, 43.57142857142857], [195.0, 43.57142857142857], [200.0, 43.57142857142857], [100.0, 54.21428571428571], [105.0, 54.21428571428571], [110.0, 54.21428571428571], [115.0, 54.21428571428571], [120.0, 54.21428571428571], [125.0, 54.21428571428571], [130.0, 54.21428571428571], [135.0, 54.21428571428571], [140.0, 54.21428571428571], [145.0, 54.21428571428571], [150.0, 54.21428571428571], [155.0, 54.21428571428571], [160.0, 54.21428571428571], [165.0, 54.21428571428571], [170.0, 54.21428571428571], [175.0, 54.21428571428571], [180.0, 54.21428571428571], [185.0, 54.21428571428571], [190.0, 54.21428571428571], [195.0, 54.21428571428571], [200.0, 54.21428571428571], [100.0, 64.85714285714286], [105.0, 64.85714285714286], [110.0, 64.85714285714286], [115.0, 64.85714285714286], [120.0, 64.85714285714286], [125.0, 64.85714285714286], [130.0, 64.85714285714286], [135.0, 64.85714285714286], [140.0, 64.85714285714286], [145.0, 64.85714285714286], [150.0, 64.85714285714286], [155.0, 64.85714285714286], [160.0, 64.85714285714286], [165.0, 64.85714285714286], [170.0, 64.85714285714286], [175.0, 64.85714285714286], [180.0, 64.85714285714286], [185.0, 64.85714285714286], [190.0, 64.85714285714286], [195.0, 64.85714285714286], [200.0, 64.85714285714286], [100.0, 75.5], [105.0, 75.5], [110.0, 75.5], [115.0, 75.5], [120.0, 75.5], [125.0, 75.5], [130.0, 75.5], [135.0, 75.5], [140.0, 75.5], [145.0, 75.5], [150.0, 75.5], [155.0, 75.5], [160.0, 75.5], [165.0, 75.5], [170.0, 75.5], [175.0, 75.5], [180.0, 75.5], [185.0, 75.5], [190.0, 75.5], [195.0, 75.5], [200.0, 75.5], [100.0, 86.14285714285714], [105.0, 86.14285714285714], [110.0, 86.14285714285714], [115.0, 86.14285714285714], [120.0, 86.14285714285714], [125.0, 86.14285714285714], [130.0, 86.14285714285714], [135.0, 86.14285714285714], [140.0, 86.14285714285714], [145.0, 86.14285714285714], [150.0, 86.14285714285714], [155.0, 86.14285714285714], [160.0, 86.14285714285714], [165.0, 86.14285714285714], [170.0, 86.14285714285714], [175.0, 86.14285714285714], [180.0, 86.14285714285714], [185.0, 86.14285714285714], [190.0, 86.14285714285714], [195.0, 86.14285714285714], [200.0, 86.14285714285714], [100.0, 96.78571428571428], [105.0, 96.78571428571428], [110.0, 96.78571428571428], [115.0, 96.78571428571428], [120.0, 96.78571428571428], [125.0, 96.78571428571428], [130.0, 96.78571428571428], [135.0, 96.78571428571428], [140.0, 96.78571428571428], [145.0, 96.78571428571428], [150.0, 96.78571428571428], [155.0, 96.78571428571428], [160.0, 96.78571428571428], [165.0, 96.78571428571428], [170.0, 96.78571428571428], [175.0, 96.78571428571428], [180.0, 96.78571428571428], [185.0, 96.78571428571428], [190.0, 96.78571428571428], [195.0, 96.78571428571428], [200.0, 96.78571428571428], [100.0, 107.42857142857142], [105.0, 107.42857142857142], [110.0, 107.42857142857142], [115.0, 107.42857142857142], [120.0, 107.42857142857142], [125.0, 107.42857142857142], [130.0, 107.42857142857142], [135.0, 107.42857142857142], [140.0, 107.42857142857142], [145.0, 107.42857142857142], [150.0, 107.42857142857142], [155.0, 107.42857142857142], [160.0, 107.42857142857142], [165.0, 107.42857142857142], [170.0, 107.42857142857142], [175.0, 107.42857142857142], [180.0, 107.42857142857142], [185.0, 107.42857142857142], [190.0, 107.42857142857142], [195.0, 107.42857142857142], [200.0, 107.42857142857142], [100.0, 118.07142857142857], [105.0, 118.07142857142857], [110.0, 118.07142857142857], [115.0, 118.07142857142857], [120.0, 118.07142857142857], [125.0, 118.07142857142857], [130.0, 118.07142857142857], [135.0, 118.07142857142857], [140.0, 118.07142857142857], [145.0, 118.07142857142857], [150.0, 118.07142857142857], [155.0, 118.07142857142857], [160.0, 118.07142857142857], [165.0, 118.07142857142857], [170.0, 118.07142857142857], [175.0, 118.07142857142857], [180.0, 118.07142857142857], [185.0, 118.07142857142857], [190.0, 118.07142857142857], [195.0, 118.07142857142857], [200.0, 118.07142857142857], [100.0, 128.71428571428572], [105.0, 128.71428571428572], [110.0, 128.71428571428572], [115.0, 128.71428571428572], [120.0, 128.71428571428572], [125.0, 128.71428571428572], [130.0, 128.71428571428572], [135.0, 128.71428571428572], [140.0, 128.71428571428572], [145.0, 128.71428571428572], [150.0, 128.71428571428572], [155.0, 128.71428571428572], [160.0, 128.71428571428572], [165.0, 128.71428571428572], [170.0, 128.71428571428572], [175.0, 128.71428571428572], [180.0, 128.71428571428572], [185.0, 128.71428571428572], [190.0, 128.71428571428572], [195.0, 128.71428571428572], [200.0, 128.71428571428572], [100.0, 139.35714285714286], [105.0, 139.35714285714286], [110.0, 139.35714285714286], [115.0, 139.35714285714286], [120.0, 139.35714285714286], [125.0, 139.35714285714286], [130.0, 139.35714285714286], [135.0, 139.35714285714286], [140.0, 139.35714285714286], [145.0, 139.35714285714286], [150.0, 139.35714285714286], [155.0, 139.35714285714286], [160.0, 139.35714285714286], [165.0, 139.35714285714286], [170.0, 139.35714285714286], [175.0, 139.35714285714286], [180.0, 139.35714285714286], [185.0, 139.35714285714286], [190.0, 139.35714285714286], [195.0, 139.35714285714286], [200.0, 139.35714285714286], [100.0, 150.0], [105.0, 150.0], [110.0, 150.0], [115.0, 150.0], [120.0, 150.0], [125.0, 150.0], [130.0, 150.0], [135.0, 150.0], [140.0, 150.0], [145.0, 150.0], [150.0, 150.0], [155.0, 150.0], [160.0, 150.0], [165.0, 150.0], [170.0, 150.0], [175.0, 150.0], [180.0, 150.0], [185.0, 150.0], [190.0, 150.0], [195.0, 150.0], [200.0, 150.0]])
        
        # CAPEX corresponding data
        CAPEX_vec = np.array([18392.267668261557, 18296.25823894743, 18707.894179236922, 19424.004343875113, 18500.05619359577, 19405.940074013823, 18644.358552157093, 18878.98020549728, 18242.72292386092, 20260.639886635487, 20737.240582636754, 18752.494028035697, 18936.716640798077, 20215.037741297027, 21422.36029000573, 19058.48259634056, 18765.62782090656, 23999.556351209394, 25639.38664777278, 26434.23804097734, 18527.630386910052, 29405.74174099577, 29324.502300657085, 29099.732935060943, 29955.467824286723, 31491.34848304245, 33219.4860415188, 35137.61641885675, 37330.25565969137, 39802.795354567395, 42540.48165924593, 45603.977731919855, 48830.94278305167, 52360.326283708986, 56220.31935326285, 60525.48938850523, 65195.18598039591, 70327.03592666995, 75736.32129458207, 81870.72873701747, 88355.29843782721, 95635.91324918423, 37732.54974384038, 37573.58501325523, 37118.3826545419, 37679.049439250084, 40235.796619212466, 43145.52549352378, 46398.80341990883, 50217.191773490056, 54286.377883886096, 58958.02472059048, 64108.39521538311, 69541.41517908828, 75673.2571432345, 82343.63613947983, 89633.33819261304, 97528.7839454286, 106267.63247397667, 115930.76284433393, 126230.2373392406, 137775.985219391, 149815.63640142733, 45364.548257385184, 45084.70307311279, 44309.32895071161, 44235.77526182787, 47720.59022908528, 51670.2580232679, 56157.301820021814, 61182.87569667855, 66869.1830550865, 73218.99456295895, 80367.08910465869, 87969.15056838168, 96205.12624114341, 105396.51015282754, 115457.51521427587, 126181.67737041079, 138473.02822022105, 151782.05147374052, 166099.64213249818, 181907.8074859281, 199022.37688858874, 52503.04094906923, 52063.00298353686, 50942.408627065975, 50067.769972490954, 54471.834863783035, 59317.29670775778, 64913.60564361377, 71268.83288335905, 78352.29318871592, 86382.98504845687, 95279.62016241222, 104788.30113422261, 115113.83874501473, 126606.95086121044, 139298.02727739676, 153413.58112076626, 168558.6887917496, 185346.3830379233, 203628.62607801147, 223433.6159537519, 245017.83788349773, 59303.091611225944, 58637.46182872077, 57232.5343125029, 55442.80095490066, 60565.67028863193, 66350.50121099943, 73025.02698206308, 80502.72127540223, 89023.58450409412, 98581.43727075479, 109181.31045684913, 120516.80591250856, 132992.81326439712, 146959.6168174259, 162078.64588958668, 179233.85043452503, 197077.79621037657, 217212.16329151046, 239215.19745479702, 262794.84229963284, 288741.8964761115, 65809.90189655963, 64992.2602069468, 63096.8423418961, 60368.79132107711, 66278.38145300862, 73076.02996199076, 80585.27187756557, 89212.92422584337, 99081.38396568019, 110049.45581048033, 122320.82289866127, 135249.1885255766, 150028.3781482254, 166036.66142563694, 183692.467111848, 203105.5785626544, 224370.7648285722, 247747.75008802526, 273171.38506775105, 300883.6410364323, 331081.49131340435, 72077.25083237662, 71062.61808265864, 68804.30143817601, 64952.08620578276, 71486.40342736754, 79138.82428795223, 87741.44294708232, 97319.65134517415, 108605.07967151179, 121018.77717824373, 134824.59678005156, 149675.58951202518, 166253.42779826146, 184287.9511864904, 204430.21782549596, 226527.96237136834, 250578.37940414902, 277254.31300221826, 306136.0384739281, 337642.53383523243, 371489.69474096515, 78153.94274939218, 76952.21149102591, 74279.85345270843, 69399.77022255809, 76814.40080303999, 85067.23996139657, 94550.3139439409, 105377.02473683222, 117662.87947507077, 131316.9093954696, 146972.243379638, 163509.52795332135, 181835.14347533145, 202193.22723944485, 224631.63737093235, 249091.70829311726, 276130.5909321933, 305738.13615193445, 338104.0994058333, 373231.1863557489, 411637.9295355743, 84024.05532433577, 82582.62346224123, 79598.77479929489, 73537.00159877198, 81493.30673150196, 90819.7282427578, 101195.30615462647, 112957.31403954396, 126462.77746684957, 141667.87868822674, 158610.6028257952, 176841.88360721496, 197056.2455928804, 219312.36045660957, 244089.9944395824, 271163.343314335, 301023.07559613854, 333581.3819436542, 369400.3860151866, 408135.6466204421, 450319.50871888804, 89830.21134599061, 88210.9716303698, 84715.89825445713, 77604.94986285796, 86146.60995077598, 96146.13885181151, 107324.44206288228, 120236.39651129591, 135042.84022576438, 151525.22738118516, 169938.96533631856, 189755.80502208782, 211728.18932027413, 236027.65181939368, 263002.7577366219, 292586.15634228423, 325089.4671445853, 360768.55712052283, 399569.39512118377, 442020.33182616666, 488287.69447110716, 95438.13499711816, 93586.67285909729, 89726.3739569308, 81353.94914661856, 90643.36097078671, 101326.85144569461, 113498.97660357888, 127350.02630018877, 143142.3164673673, 160978.15386990758, 180948.8399278423, 202248.5252213173, 226166.02524471094, 252472.16286579557, 281546.36835144606, 313573.21262908564, 348847.41231290094, 387358.00829829986, 429543.4170547492, 475486.74538864597, 525543.046121406, 100956.69001317178, 98909.34431269954, 94646.15292171558, 85124.108994894, 94971.74324120942, 106361.19247085703, 119351.63117678725, 134191.43326731215, 151013.07761340885, 170329.9440253746, 191829.931912252, 214910.30908296854, 240236.02299413554, 268475.9083410715, 299782.4730902847, 334482.646053479, 372167.3424384517, 413497.3244142854, 458633.5458494597, 508241.4970666799, 561975.9048484723, 106363.77297809818, 104083.506772689, 99393.27349514692, 89146.25085171398, 99138.1172379357, 111289.39327364763, 125013.63752184468, 140865.06734485386, 158864.93203861619, 179325.76145251663, 202139.9622041541, 226676.01513761442, 253893.63580767278, 284307.2377790243, 317572.89929793024, 354334.0411423184, 394823.0589352352, 439067.4743976819, 487572.3217459507, 540531.6523003523, 597727.6190665967, 111589.05041051836, 109156.53227512873, 104094.17097231765, 92952.70871140448, 103378.26605521196, 115947.07929072778, 130625.73843142879, 147306.22017598973, 166374.13405738666, 188093.15468554993, 212383.43914101884, 238446.71768272825, 267479.7405469859, 299564.57805819134, 335108.27070429176, 374251.64518034784, 417248.7868811624, 464341.60617713956, 515931.2057622455, 572183.4967351712, 633520.2680455706])
        
        from scipy.interpolate import LinearNDInterpolator
        
        self.CAPEX_T_V_model = LinearNDInterpolator(points, CAPEX_vec)
        
        return 
    
    # ---------------- Design Methods ----------------------------------------------------------------------
    
    def design_system(self, x):
        # self.T_des = x[0]
        self.D_i   = x[0]

        self.H = (self.inputs['V'] - (4/3)*np.pi*(self.D_i/2)**3)/(np.pi*(self.D_i/2)**2)
        
        # print(f"D_i : {self.D_i}")
        # print(f"H : {self.H}")
        
        try:
            if self.H <= 0 or self.H >= 20*self.D_i:
                raise ValueError()
            
            self.T_des = self.inputs['T_des']
        
            self.material_properties()
        
            if self.T_des >= 100+273.15:
                (P_sat, rho_l) = PropsSI(('P','D'), 'T', self.T_des, 'Q', 0, self.fluid)
            else:
                P_sat = 101325
                rho_l = PropsSI('D', 'T', self.T_des, 'P', P_sat, self.fluid)
        
            P_stat = P_sat + rho_l*self.g*self.H
            self.P_des = max(P_stat*(1+self.params['P_OD']), 2.5*1e5)
        
            self.t_req = self.req_thickness(self.P_des, self.D_i, self.sigma_max)
        
            self.lateral_plate_vol = np.pi*(self.D_i/2 + self.t_req)**2 * self.t_req * 2
            self.shell_vol = (np.pi*(self.D_i/2 + self.t_req)**2 - np.pi*(self.D_i/2)**2)*self.H
            self.metal_vol = self.shell_vol #+ self.lateral_plate_vol
            
            self.m_mat = self.rho_mat*self.metal_vol  
        
            self.fluid_vol = np.pi*(self.D_i/2)**2 * self.H
            h_cap = (PropsSI('H', 'T', self.T_des, 'P', self.P_des, self.fluid) - PropsSI('H', 'T', self.inputs['T_in'], 'P', self.P_des, self.fluid))
            E_cap = h_cap*self.fluid_vol*rho_l
            
            self.E_cap_MWh = E_cap/(3600*1e6)
            
            self.capital_cost_est()
            
            return self.CAPEX/self.E_cap_MWh # self.Cost_cap_dollar/self.E_cap_MWh
        except:
            return 100000

    # def opt_design(self):

    #     # Bounds: [T_des, D_i, H]
    #     bounds = [
    #         # tuple(self.params['bounds_T_des']),
    #         tuple(self.params['bounds_D_i']),  # inner diameter in meters
    #         tuple(self.params['bounds_H']),
    #     ]
        
    #     # Objective function to minimize
    #     def objective(x):
    #         try:
    #             return self.design_system(x)
    #         except Exception as e:
    #             print(f"Error at x={x}: {e}")
    #             return np.inf
        
    #     # x0 = [self.params['bounds_T_des'][-1], self.params['bounds_D_i'][-1], self.params['bounds_H'][-1]]
    #     x0 = [self.params['bounds_D_i'][-1], self.params['bounds_H'][-1]]
        
    #     result = minimize(
    #         fun=objective,
    #         x0=x0,
    #         method='L-BFGS-B',
    #         bounds=bounds,
    #         options={'disp': False}
    #     )
        
    #     return result

    def opt_design(self):
        # Bounds in format ([lower_bounds], [upper_bounds]) as required by pyswarms
        lower_bounds = np.array([self.params['bounds_D_i'][0]])
        upper_bounds = np.array([self.params['bounds_D_i'][1]])
                                 
        bounds = (lower_bounds, upper_bounds)
    
        # Objective wrapper for pyswarms — input is a 2D array: shape (n_particles, n_dimensions)
        def objective_wrapper(X):
            results = []
            for x in X:
                try:
                    res = self.design_system(x)
                except Exception as e:
                    print(f"Error at x={x}: {e}")
                    res = 1e10  # Large penalty
                results.append(res)
            return np.array(results)
    
        # Initialize PSO optimizer
        self.optimizer = GlobalBestPSO(
            n_particles=200,
            dimensions=1,
            options={'c1': 1.5, 'c2': 2.0, 'w': 0.7},
            bounds=bounds
        )
    
        best_cost = np.inf
        no_improve_counter = 0
        patience = 5
        tol = 1e-3
        max_iter = 30
    
        for i in tqdm(range(max_iter), desc="Optimizing", ncols=80):
            self.optimizer.optimize(objective_wrapper, iters=1, verbose=False)
            current_best = self.optimizer.swarm.best_cost
    
            if current_best < best_cost - tol:
                best_cost = current_best
                no_improve_counter = 0
            else:
                no_improve_counter += 1
    
            # print(f"[{i+1:03}] Best cost: {best_cost:.6f}")
            if no_improve_counter >= patience:
                # print("Stopping early due to stagnation.")
                break
    
        # Extract best solution
        best_position = self.optimizer.swarm.best_pos
    
        # Recalculate final tank state using best solution
        self.design_system(best_position)
    
        # Return result in a similar format to scipy.optimize
        from types import SimpleNamespace
        result = SimpleNamespace()
        result.x = best_position
        result.fun = best_cost
        result.success = True
        return result


#%% 

if __name__ == "__main__":
    
    T_list = []
    V_list = []
    
    T_des_vec = np.linspace(100,200,21) + 273.15
    V_obj_vec = np.linspace(1, 150, 15)
    
    D_opt_vec = []

    # Storage
    fluid = 'Water'
    fun_vec = []
    cap_vec = []
    CAPEX_vec = []
    
    for V in V_obj_vec:
        for T_des in T_des_vec:
            Tank = PressuzedWaterStorageDesign(fluid)
            
            T_in = 15 + 273.15 # K
            P_OD = 0.15
            
            Tank.set_inputs(
                T_in = T_in, # [K]
                T_des = T_des, # [K]
                V = V 
            )
            
            Tank.set_parameters(
                  P_OD = P_OD,
                  Material = '304ss', # 'Carbon_Steel' or '304ss'
                  Arrang = 'Vertical',
                  bounds_D_i = [0.1, 20],
            )
    
            opt_result = Tank.opt_design()
            
            D_opt_vec.append(opt_result.x[0])
            T_list.append(T_des - 273.15)
            V_list.append(V)
            fun_vec.append(opt_result.fun) # $/MWh
            cap_vec.append(Tank.E_cap_MWh * 1e3 / Tank.fluid_vol)  # kWh/m³
            CAPEX_vec.append(Tank.CAPEX)  # $
        
    T_grid = np.unique(T_list)
    V_grid = np.unique(V_list)
    
    fun_array = np.array(fun_vec).reshape(len(V_grid), len(T_grid))
    cap_array = np.array(cap_vec).reshape(len(V_grid), len(T_grid))
    CAPEX_array = np.array(CAPEX_vec).reshape(len(V_grid), len(T_grid))
        
    T_mesh, V_mesh = np.meshgrid(T_grid, V_grid)
    
    # Plot: Capital Cost per energy stored [$/MWh]
    plt.figure(figsize=(8, 6))
    cp = plt.contourf(T_mesh, V_mesh, fun_array, cmap='viridis', levels=20)
    plt.colorbar(cp, label='Capital Cost per energy stored [$/MWh]')
    plt.xlabel('Temperature [°C]')
    plt.ylabel('Storage Volume [m³]')
    plt.title('Capital Cost per energy stored vs Volume and Temperature')
    plt.tight_layout()
    plt.grid(True)
    
    # Plot: Capital Cost [$]
    plt.figure(figsize=(8, 6))
    cp = plt.contourf(T_mesh, V_mesh, CAPEX_array, cmap='viridis', levels=20)
    plt.colorbar(cp, label='Capital Cost [$]')
    plt.xlabel('Temperature [°C]')
    plt.ylabel('Storage Volume [m³]')
    plt.title('Capital Cost vs Volume and Temperature')
    plt.tight_layout()
    plt.grid(True)
    
    # Plot: Energy Density [kWh/m³]
    plt.figure(figsize=(8, 6))
    cp = plt.contourf(T_mesh, V_mesh, cap_array, cmap='plasma', levels=20)
    plt.colorbar(cp, label='Energy Density [kWh/m³]')
    plt.xlabel('Temperature [°C]')
    plt.ylabel('Storage Volume [m³]')
    plt.title('Energy Density vs Volume and Temperature')
    plt.tight_layout()
    plt.grid(True)
    
    plt.show()


    if opt_result.success:
        print("\n✅ Optimisation réussie !")
        print(f"D_i optimal   : {opt_result.x[0]:.3f} m")
        print(f"Coût [$/MWh] : {opt_result.fun:.2f} $/MWh")
    else:
        print("\n❌ Optimisation échouée.")
        print(opt_result.message)

    from scipy.interpolate import RegularGridInterpolator

    T_grid = np.unique(T_list)  # Temperature in °C
    V_grid = np.unique(V_list)  # Volume in m³
    
    # Already reshaped from earlier:
    fun_array = np.array(fun_vec).reshape(len(V_grid), len(T_grid))  # $/MWh
    
    # Create the interpolator
    cost_interp = RegularGridInterpolator(
        points=(V_grid, T_grid),
        values=fun_array,
        bounds_error=False,
        fill_value=None  # or np.nan
    )
    
    # Redefine everything from scratch to ensure consistency

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import LinearNDInterpolator

    T_vec = np.array(T_list)
    V_vec = np.array(V_list)

    # Mock dataset for ΔT and Volume (replace with your real data)
    np.random.seed(0)
    T_grid, V_grid = np.meshgrid(T_vec, V_vec)
    T_vals = T_grid.flatten()
    V_vals = V_grid.flatten()

    # Create interpolator
    points = np.column_stack((T_vec, V_vec))
    interpolator = LinearNDInterpolator(points, CAPEX_vec)

    # Predict using interpolation
    CAPEX_pred_interp = interpolator(T_vec, V_vec)

    # Compute relative error
    relative_error_interp = (CAPEX_pred_interp - CAPEX_vec) / CAPEX_vec

    # Plot relative error
    plt.figure(figsize=(8, 5))
    plt.scatter(range(len(relative_error_interp)), relative_error_interp * 100, color='green', alpha=0.7)
    plt.axhline(0, color='black', linestyle='--')
    plt.xlabel("Data Point Index")
    plt.ylabel("Relative Error (%)")
    plt.title("Relative Error of CAPEX Prediction using Interpolation")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Function to compute your cost model
    def cost_model_custom(V, T):
        return (
            # 0.00002 * T**3.838 * V**0.785 + 17400
            0.00003 * (T)**3.725 * V**0.786 + 17400
            # 13.23602 * (T - 100)**1.464 * V**0.786 + 17400
            # 1.19e-04 * T**3.493 * V**0.748 + 17400
        )

    # 2009 reference cost model
    def cost_model_2009(V):
        return 4500*1.45 * V ** 0.6 + 14500

    def cost_model_2017(V):
        return 2945.5*1.25 * V ** 0.667

    # Volume range
    V_range = np.linspace(1, 150, 200)  # avoid V=0 for power laws

    # Temperatures from 100°C to 150°C
    Temp = np.arange(100, 151, 10)

    # Plotting
    plt.figure(figsize=(10, 6))

    for T in Temp:
        cost_custom = interpolator(T, V_range)
        plt.plot(V_range, cost_custom, label=f"Your Model (T={T}°C)", linewidth=2)

    # Add reference model (single curve)
    cost_ref = cost_model_2009(V_range)
    cost_ref_2017 = cost_model_2017(V_range)

    plt.plot(V_range, cost_ref, label="2009 Model", linestyle="--", color='black', linewidth=2)
    plt.plot(V_range, cost_ref_2017, label="2017 Model", linestyle="--", color='red', linewidth=2)

    # Labels and layout
    plt.xlabel("Volume [m³]")
    plt.ylabel("Cost [USD]")
    plt.title("Cost vs Volume for Various Temperatures")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    
    # 2009 reference cost model
    def cost_model_2009(V):
        return 4500*1.45 * V ** 0.6 + 14500

    def cost_model_2017(V):
        return 2945.5*1.25 * V ** 0.667

    # Volume range
    V_range = np.linspace(1, 150, 200)  # avoid V=0 for power laws

    # Temperatures from 100°C to 150°C
    Temp = np.arange(100, 151, 10)

    # Plotting
    plt.figure(figsize=(10, 6))

    for T in Temp:
        cost_custom = Tank.interpolator(T, V_range)
        plt.plot(V_range, cost_custom, label=f"Your Model (T={T}°C)", linewidth=2)

    # Add reference model (single curve)
    cost_ref = cost_model_2009(V_range)
    cost_ref_2017 = cost_model_2017(V_range)

    plt.plot(V_range, cost_ref, label="2009 Model", linestyle="--", color='black', linewidth=2)
    plt.plot(V_range, cost_ref_2017, label="2017 Model", linestyle="--", color='red', linewidth=2)

    # Labels and layout
    plt.xlabel("Volume [m³]")
    plt.ylabel("Cost [USD]")
    plt.title("Cost vs Volume for Various Temperatures")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()