
import __init__
import numpy as np
import matplotlib.pyplot as plt
from library.component.sizing.heat_exchanger.basic_sizing import find_UA

case_1 = "Zorlu ht-hp-pcm-orc-rec"
case_2 = "Zorlu ht-hp-pcm-orc-rec-pr"
case_3 = "Zorlu ht Sizing"

flow_config = "CounterFlow"

"""
CounterFlow
CrossFlow
Shell&Tube
ParallelFlow
"""

params = {'Flow_Type':flow_config}

if flow_config == "Shell&Tube":
    params['n_series'] = 2

case = case_3

if case == case_1:
    data = {
            "Evap_HP" : {
                        "T_c_i" : 36.65,
                        "T_c_o" : 109.5,
                        "T_h_i" : 113.1,
                        "T_h_o" : 111.5,
                        "Q_dot" : 5722*1e3},

            "IHX_HP" : {
                        "T_c_i" : 109.5,
                        "T_c_o" : 143.1,
                        "T_h_i" : 153,
                        "T_h_o" : 129.1,
                        "Q_dot" : 1177*1e3},

            "COND_HP" : {
                        "T_c_i" : 150,
                        "T_c_o" : 150,
                        "T_h_i" : 186.7,
                        "T_h_o" : 153,
                        "Q_dot" : 7000*1e3},

            "EVAP_ORC" : {
                        "T_c_i" : 57.07,
                        "T_c_o" : 147,
                        "T_h_i" : 150,
                        "T_h_o" : 150,      
                        "Q_dot" : 13300*1e3},
    
            "RECUP_ORC" : {
                        "T_c_i" : 25.85,
                        "T_c_o" : 57.05,
                        "T_h_i" : 79.86,
                        "T_h_o" : 36.65,
                        "Q_dot" : 1573*1e3},     

            "COND_ORC" : {
                        "T_c_i" : 24,
                        "T_c_o" : 25.34,
                        "T_h_i" : 36.65,
                        "T_h_o" : 25.3,
                        "Q_dot" : 11132*1e3}                   
            }

if case == case_2:
    data = {
            "Evap_HP" : {
                        "T_c_i" : 36.65,
                        "T_c_o" : 109.5,
                        "T_h_i" : 113.1,
                        "T_h_o" : 111.5,
                        "Q_dot" : 5722*1e3},

            "IHX_HP" : {
                        "T_c_i" : 109.5,
                        "T_c_o" : 143.1,
                        "T_h_i" : 153,
                        "T_h_o" : 129.1,
                        "Q_dot" : 1177*1e3},

            "COND_HP" : {
                        "T_c_i" : 150,
                        "T_c_o" : 150,
                        "T_h_i" : 186.7,
                        "T_h_o" : 153,
                        "Q_dot" : 7000*1e3},

            "EVAP_ORC" : {
                        "T_c_i" : 102.2,
                        "T_c_o" : 147,
                        "T_h_i" : 150,
                        "T_h_o" : 150,      
                        "Q_dot" : 13300*1e3},

            "PR_ORC" : {
                        "T_c_i" : 58.64,
                        "T_c_o" : 102.2,
                        "T_h_i" : 113.1,
                        "T_h_o" : 112.3,      
                        "Q_dot" : 3017*1e3},

            "RECUP_ORC" : {
                        "T_c_i" : 27.75,
                        "T_c_o" : 58.64,
                        "T_h_i" : 81.15,
                        "T_h_o" : 38.43,
                        "Q_dot" : 1931*1e3},     

            "COND_ORC" : {
                        "T_c_i" : 24,
                        "T_c_o" : 27.29,
                        "T_h_i" : 38.43,
                        "T_h_o" : 27.21,
                        "Q_dot" : 13757*1e3}                   
            }

import warnings
warnings.filterwarnings("ignore") 

if case == case_1 or case == case_2:
    print('\n')

    if flow_config == 'Shell&Tube':
        print(f"Case : {case}, \nFlow configuration : {flow_config}, \nNumber of passes per shell : {params['n_series']}. \n")
    else:
        print(f"Case : {case}, \nFlow configuration : {flow_config}. \n")

    UA_tot = 0

    for HTX in data:
        data[HTX]['UA'] = find_UA(data[HTX]['Q_dot'], data[HTX]['T_c_i'], data[HTX]['T_c_o'], data[HTX]['T_h_i'], data[HTX]['T_h_o'], flow = flow_config, params = params)
        print(f"Required UA for {HTX}: {data[HTX]['UA']} (W/K)")
        
        if not np.isnan(data[HTX]['UA']):
            UA_tot += data[HTX]['UA']

    print(f"Total required UA : {UA_tot} [W/K]")

    print('\n')

#%%

if case == case_3:
        data = {
            "1" : {
                        "DT_PP" : 1,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 36.57,
                        "T_h_o" : 25.2,
                        "Q_dot" : 13779*1e3,
                        "W_dot_exp" : 2675,
                        "RTE" : 0.8741},

            "2" : {
                        "DT_PP" : 2,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 37.5,
                        "T_h_o" : 26.2,
                        "Q_dot" : 13768*1e3,
                        "W_dot_exp" : 2646,
                        "RTE" : 0.8668},

            "3" : {
                        "DT_PP" : 3,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 38.43,
                        "T_h_o" : 27.2,
                        "Q_dot" : 13757*1e3,
                        "W_dot_exp" : 2618,
                        "RTE" : 0.8595},

            "4" : {
                        "DT_PP" : 4,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 39.37,
                        "T_h_o" : 28.2,
                        "Q_dot" : 13746*1e3,
                        "W_dot_exp" : 2591,
                        "RTE" : 0.8521},

            "5" : {
                        "DT_PP" : 5,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 40.3,
                        "T_h_o" : 29.2,
                        "Q_dot" : 13735*1e3,
                        "W_dot_exp" : 2563,
                        "RTE" : 0.8448},

            "6" : {
                        "DT_PP" : 6,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 41.23,
                        "T_h_o" : 30.2,
                        "Q_dot" : 13723*1e3,
                        "W_dot_exp" : 2535,
                        "RTE" : 0.8374},   

            "7" : {
                        "DT_PP" : 7,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 42.16,
                        "T_h_o" : 31.2,
                        "Q_dot" : 13712*1e3,
                        "W_dot_exp" : 2507,
                        "RTE" : 0.8301},    

                        
            "8" : {
                        "DT_PP" : 8,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 43.09,
                        "T_h_o" : 32.2,
                        "Q_dot" : 13700*1e3,
                        "W_dot_exp" : 2480,
                        "RTE" : 0.8228},

                        
            "9" : {
                        "DT_PP" : 9,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 44.01,
                        "T_h_o" : 33.2,
                        "Q_dot" : 13688*1e3,
                        "W_dot_exp" : 2452,
                        "RTE" : 0.8154},       
            
            "10" : {
                        "DT_PP" : 10,
                        "T_c_i" : 24,
                        "T_c_o" : 27.3,
                        "T_h_i" : 44.94,
                        "T_h_o" : 34.2,
                        "Q_dot" : 13676*1e3,
                        "W_dot_exp" : 2425,
                        "RTE" : 0.808},         
            }


if case == case_3:
    UA_val = []
    RTE_val = []

    print('\n')

    if flow_config == 'Shell&Tube':
        print(f"Case : {case}, \nFlow configuration : {flow_config}, \nNumber of passes per shell : {params['n_series']}. \n")
    else:
        print(f"Case : {case}, \nFlow configuration : {flow_config}. \n")

    for pp_val in data:
        data[pp_val]['UA'] = find_UA(data[pp_val]['Q_dot'], data[pp_val]['T_c_i'], data[pp_val]['T_c_o'], data[pp_val]['T_h_i'], data[pp_val]['T_h_o'], flow = flow_config, params = params)
        print(f"Required UA for pinch of {pp_val} K: {data[pp_val]['UA']} (W/K)")
        
        UA_val.append(data[pp_val]['UA'])
        RTE_val.append(data[pp_val]['RTE'])

    print('\n')

    # Example data
    DT_pinch = np.linspace(1, 10, 10)

    # Create the first subplot
    fig, ax1 = plt.subplots()

    # Plot the first curve on the left y-axis
    ax1.plot(DT_pinch, UA_val, 'b-', label='sin(x)')
    ax1.set_xlabel('$\Delta$T$_{pinch,cd}$ [K]')  # Set the x-axis label
    ax1.set_ylabel('UA$_{req}$ [W/K]', color='b')  # Set the y-axis label for the left side
    ax1.tick_params(axis='y', labelcolor='b')

    # Create the second y-axis (right side)
    ax2 = ax1.twinx()

    # Plot the second curve on the right y-axis
    ax2.plot(DT_pinch, RTE_val, 'r-', label='RTE')
    ax2.set_ylabel('RTE [-]', color='r')  # Set the y-axis label for the right side
    ax2.tick_params(axis='y', labelcolor='r')

    # Axis of the plot
    ax1.set_xlim(1,10)
    ax1.set_ylim(0,3.6*1e6)
    ax2.set_ylim(0.8, 0.92)

    # Show the plot
    ax1.grid(True, axis='x')
    plt.show()

