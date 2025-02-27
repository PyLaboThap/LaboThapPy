# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:35:25 2023

@author: Basile
"""

import pandas as pd
import numpy as np
import os
import fnmatch
import re
import matplotlib.pyplot as plt
import matplotlib

def res_encoding(Matrix, Name, Sheet_names, Column_names, Row_names):
    # Replace nan values with 0 (4D matrix)
    for i in range(len(Matrix)):
        for j in range(len(Matrix[i])):
            for k in range(len(Matrix[i][j])):
                if np.isnan(Matrix[i][j][k]):
                    Matrix[i][j][k] = 0
    
    # Construct a column name list with information from Column_names
    cnames = [f'L_HTX = {round(col,1)}' for col in Column_names]
    
    # Construct a row name list with information from Row_names
    rnames = [f'L_pt = {round(row,1)}' for row in Row_names]

    
    # Create an Excel writer
    with pd.ExcelWriter('Results/' + Name + '.xlsx', engine='openpyxl') as writer:
        for i, two_d_array in enumerate(Matrix):
            # Construct a sheet name with information from Sheet_names
            sname = f'D_o = {Sheet_names[i]}'
             
            # Create a DataFrame from each 2D array with headers
            df = pd.DataFrame(two_d_array, columns=cnames)
            
            # Add row headers
            df.index = rnames  
            
            # Write the DataFrame to a separate sheet in the Excel file
            df.to_excel(writer, sheet_name=sname, index=True)
    
    print("Encoded 4D matrix saved to Excel with separate sheets, headers, and row headers.")

def get_data(parent_directory, search_strings):
    files_data = {}

    for search_string in search_strings:
        search_data = {}

        for foldername, subfolders, filenames in os.walk(parent_directory):
            for filename in fnmatch.filter(filenames, '*.xlsx'):
                # Skip temporary Excel files
                if filename.startswith('~$'):
                    continue
                
                if search_string in filename:
                    file_path = os.path.join(foldername, filename)
                    xls = pd.ExcelFile(file_path)

                    file_data = {}

                    for sheet_name in xls.sheet_names:
                        df = xls.parse(sheet_name, index_col=0)
                        file_data[sheet_name] = df

                    search_data[filename] = file_data

        files_data[search_string] = search_data

    return files_data

def extract_numerical_value_from_filename(file_name):
    match = re.search(r'\d+(\_\d+)?', file_name)

    if match:
        return float(match.group().replace('_', '.'))
    else:
        return None

def extract_variable_and_value(name_string):
    parts = name_string.split('=')
    variable_name = parts[0].strip()
    numerical_value = float(parts[1].strip())
    return variable_name, numerical_value

def collect_data(parent_directory, search_strings):
    files_data = get_data(parent_directory, search_strings)

    result_dict = {}

    for search_string, search_data in files_data.items():
        matching_results = []

        for file_name, file_data in search_data.items():
            numerical_value = extract_numerical_value_from_filename(file_name)

            for sheet_name, df in file_data.items():
                for index, row in df.iterrows():
                    for col, value in row.items():
                        matching_results.append({
                            search_string: value,
                            'W_HTX': numerical_value,
                            'D_o': extract_variable_and_value(sheet_name)[1],
                            'L_pt': extract_variable_and_value(index)[1],
                            'L_HTX': extract_variable_and_value(col)[1],
                        })

        result_df = pd.DataFrame(matching_results)
        result_dict[search_string] = result_df

    return result_dict

def concatenate_dataframes(dataframes):
    # Use the first DataFrame as the base
    result_df = dataframes[list(dataframes.keys())[0]]

    # Iterate over the rest of the DataFrames and merge on common columns
    for df_name, df in dataframes.items():
        if df_name != list(dataframes.keys())[0]:
            result_df = pd.merge(result_df, df, on=['W_HTX', 'D_o', 'L_pt', 'L_HTX'], how='outer', suffixes=('', f'_{df_name}'))

    return result_df

opt = 0

if opt == 1:
    
    t = 3.5/1000 # Thickness
    F_r = 0.6
    
    rho_steel = 8000
    rho_water = 1000
    
    # Example usage:
    parent_directory = 'Results\Results_2'  # Adjust this to your actual parent directory
    search_strings = ['Q_dot', 'DP_air', 'DP_oil','coef','N_T','P_in_max','P_max_adm','R_c','R_e','R_tube','T_in_max']  # Adjust this to your specific search strings
    
    result_dict = collect_data(parent_directory, search_strings)
    
    # Concatenate the DataFrames in the dictionary
    concatenated_df = concatenate_dataframes(result_dict)
    
    concatenated_df['Volume'] = concatenated_df['L_pt'] * concatenated_df['L_HTX'] * concatenated_df['W_HTX']
    Tube_volume = concatenated_df['L_pt'] * concatenated_df['D_o'] * concatenated_df['D_o'] * (np.pi/4)
    Water_volume = concatenated_df['L_pt'] * (concatenated_df['D_o'] - t) * (concatenated_df['D_o'] - t) * (np.pi/4)

    concatenated_df['T_Mass'] = ( (Tube_volume-Water_volume)*rho_steel + Water_volume*rho_water*F_r ) * concatenated_df['N_T']
    
    # Assuming 'Q_dot' is one of the columns in the concatenated DataFrame
    concatenated_df_filtered_Q_dot = concatenated_df[concatenated_df['Q_dot'] >= 7.9e6].reset_index(drop=True)
    concatenated_df_filtered = concatenated_df_filtered_Q_dot[concatenated_df_filtered_Q_dot['DP_air'] <= 1e4].reset_index(drop=True)
      
    #%% Find min rows
    
    # Assuming 'chosen_column' is the column on which you want to minimize
    chosen_column = 'Q_dot'
    
    # Find the index corresponding to the minimum value in 'chosen_column'
    min_index_vol = concatenated_df_filtered['Volume'].idxmin()
    min_index_mass = concatenated_df_filtered['T_Mass'].idxmin()
    
    # Extract the row with the minimum value in 'chosen_column'
    min_row_vol = concatenated_df_filtered.loc[min_index_vol]
    min_row_mass = concatenated_df_filtered.loc[min_index_mass]
    
    # Find best compromise
    vol_values = concatenated_df_filtered['Volume'].values
    mass_values = concatenated_df_filtered['T_Mass'].values
    
    vol_values_unit = vol_values/min(vol_values)
    mass_values_unit = mass_values/min(mass_values)
    
    score_vect = vol_values_unit + mass_values_unit
    min_score_index = np.argmin(score_vect)
    
    min_score_row = concatenated_df_filtered.loc[min_score_index]
    
    # Now 'min_row' contains the row with the minimum value in 'chosen_column'
    print("Optimal volume configuration")
    print(min_row_vol)
    print("------------------ \n")
    
    print("Optimal tube mass configuration")
    print(min_row_mass)
    print("------------------ \n")
    
    print("Optimal compromise")
    print(min_score_row)
    print("------------------ \n")
    
    #%% Plot Results
    
    "1) All Possible Points"
    
    "1.1) Power VS Volume"
    
    # Set the font size
    font_size = 25
    label_font_size = 35
    
    # Create the scatter plot
    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.grid()
    ax1.set_ylabel('$Volume \, [m^3]$', fontsize=label_font_size)
    ax1.set_xlabel(r'$\dot{Q}_{tot} \, [MW]$', fontsize=label_font_size)
    
    # Scatter plot
    ax1.scatter(concatenated_df['Q_dot'] * 1e-6, concatenated_df['Volume'])
    plt.title(label = "All Possible Points", fontsize=font_size)
    plt.axis([0,10,0,250])
    
    plt.show()
    
    "1.2) Power VS Tube Mass"
    
    # Set the font size
    font_size = 25
    label_font_size = 35
    
    # Create the scatter plot
    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.grid()
    ax1.set_ylabel('$Total Tube Mass \, [t]$', fontsize=label_font_size)
    ax1.set_xlabel(r'$\dot{Q}_{tot} \, [MW]$', fontsize=label_font_size)
    
    # Scatter plot
    ax1.scatter(concatenated_df['Q_dot'] * 1e-6, concatenated_df['T_Mass']*1e-3)
    plt.title(label = "All Possible Points", fontsize=font_size)
    plt.axis([0,10,0,700])
    
    plt.show()
    
    "2) All Q_dot-wise valid Points"
    
    "2.1) Power VS Volume"
    
    # Set the font size
    font_size = 25
    label_font_size = 35
    
    # Create the scatter plot
    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.grid()
    ax1.set_ylabel('$Volume \, [m^3]$', fontsize=label_font_size)
    ax1.set_xlabel(r'$\dot{Q}_{tot} \, [MW]$', fontsize=label_font_size)
    
    # Scatter plot
    ax1.scatter(concatenated_df_filtered_Q_dot['Q_dot'] * 1e-6, concatenated_df_filtered_Q_dot['Volume'])
    plt.title(label = "Q_dot valid Points", fontsize=font_size)
    plt.axis([0,10,0,250])
    
    plt.show()
    
    "2.2) Power VS Tube Mass"
    
    # Set the font size
    font_size = 25
    label_font_size = 35
    
    # Create the scatter plot
    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.grid()
    ax1.set_ylabel('$Total Tube Mass \, [t]$', fontsize=label_font_size)
    ax1.set_xlabel(r'$\dot{Q}_{tot} \, [MW]$', fontsize=label_font_size)
    
    # Scatter plot
    ax1.scatter(concatenated_df_filtered_Q_dot['Q_dot'] * 1e-6, concatenated_df_filtered_Q_dot['T_Mass']*1e-3)
    plt.title(label = "Q_dot valid Points", fontsize=font_size)
    plt.axis([0,10,0,700])
    
    plt.show()
    
    "3) All Q_dot-wise and DP-wise valid Points"
    
    "3.1) Power VS Volume"
    
    # Set the font size
    font_size = 25
    label_font_size = 35
    
    # Create the scatter plot
    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.grid()
    ax1.set_ylabel('$Volume \, [m^3]$', fontsize=label_font_size)
    ax1.set_xlabel(r'$\dot{Q}_{tot} \, [MW]$', fontsize=label_font_size)
    
    # Scatter plot
    ax1.scatter(concatenated_df_filtered['Q_dot'] * 1e-6, concatenated_df_filtered['Volume'])
    plt.title(label = "Q_dot and DP valid Points", fontsize=font_size)
    plt.axis([0,10,0,250])
    
    plt.show()
    
    "3.2) Power VS Tube Mass"
    
    # Set the font size
    font_size = 25
    label_font_size = 35
    
    # Create the scatter plot
    fig, ax1 = plt.subplots(figsize=(20, 10))
    ax1.grid()
    ax1.set_ylabel('$Total Tube Mass \, [t]$', fontsize=label_font_size)
    ax1.set_xlabel(r'$\dot{Q}_{tot} \, [MW]$', fontsize=label_font_size)
    
    # Scatter plot
    ax1.scatter(concatenated_df_filtered['Q_dot'] * 1e-6, concatenated_df_filtered['T_Mass']*1e-3)
    plt.title(label = "Q_dot and DP valid Points", fontsize=font_size)
    plt.axis([0,10,0,700])
    
    plt.show()




