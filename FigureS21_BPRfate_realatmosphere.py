# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 14:41:07 2025

@author: Willis Lab
tucker wrote this im adapting it
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

Temp = 298

#these are concentrations for the two plots in format (lowrange, highrange, numberofpoints)
HO2_1= np.linspace(0.1, 100, 100)  #ppt
NO_1 = np.linspace(0.1, 10, len(HO2_1)) #ppb

HO2_2= np.linspace(0.1, 100, 100) #ppt
NO_2 = np.linspace(50, 750, len(HO2_2))  #ppb


HO2_3 = np.linspace(0.1, 500, 100) #ppt
NO_3 = np.linspace(1, 10, len(HO2_2))  #ppb




#this uses the temp-dependent rates from GECKO (no autoox), calculates bimolec lifetime at 298K
#NO*2.7e10 converts to 
def calculate_bimolecular_lifetime_ortho(NO, HO2):
    t_bi= 1 / ((2.59E-12 * np.exp(360.00 / Temp) * NO*2.7E10) + (2.58E-13 * np.exp(1300.00 / Temp) * HO2*2.7E7))
    return t_bi

common_NO_HO2_pairs = []
def looper_ortho(NO_val, HO2_val, df):
    total_iterations = len(NO_val) * len(HO2_val)
    with tqdm(total=total_iterations, desc="Processing") as pbar:
        for NO in NO_val:
            for HO2 in HO2_val:
                t_bi = calculate_bimolecular_lifetime_ortho(NO, HO2)
                df.loc[len(df)] = {"NO": NO, "HO2": HO2, "t_bi": t_bi}
                pbar.update(1)  # Update the progress bar by one step




def calculate_bimolecular_lifetime_ipso(NO, HO2):
    t_bi= 1 / ((2.44E-12 * np.exp(360.00 / Temp) * NO*2.7E10) + (2.58E-13 * np.exp(1300.00 / Temp) * HO2*2.7E7))
    return t_bi

common_NO_HO2_pairs = []
def looper_ipso(NO_val, HO2_val, df):
    total_iterations = len(NO_val) * len(HO2_val)
    with tqdm(total=total_iterations, desc="Processing") as pbar:
        for NO in NO_val:
            for HO2 in HO2_val:
                t_bi = calculate_bimolecular_lifetime_ipso(NO, HO2)
                df.loc[len(df)] = {"NO": NO, "HO2": HO2, "t_bi": t_bi}
                pbar.update(1)  # Update the progress bar by one step



#%%
result_df_1 = pd.DataFrame(columns=["NO", "HO2", "t_bi"])          
looper_ipso(NO_1,HO2_1,result_df_1)

result_df_1_pivot = result_df_1.pivot(index='HO2', columns='NO', values='t_bi')
NO_grid, HO2_grid = np.meshgrid(result_df_1_pivot.columns, result_df_1_pivot.index)
#%%

contour_style = {
    0.75: '#9F9F9F', #lighter
    1: '#858585',
    2: '#6A6A6A',
    4: '#353535',
    8: '#1B1B1B',
    60: '#000000'} #darker

#grab the levels and colors from the dictionary
custom_levels = list(contour_style.keys())
contour_colors = list(contour_style.values())

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 9), dpi=500)

#plot the contour lines with specified colors
cp = ax.contour(NO_grid, HO2_grid, result_df_1_pivot.values,
    levels=custom_levels, colors=contour_colors, linewidths=4)

ax.clabel(cp, inline=True, fontsize=17) #add labels on contour lines

#adds dashed line where tuni = tbi for BPR
ax.contour(NO_grid, HO2_grid, result_df_1_pivot.values, levels=[1 / 0.95], 
    colors='#498FFF', linewidths=4, linestyles='dashed')
#shades to the left of dashed line (where tuni < tbi)
ax.contourf(NO_grid, HO2_grid, result_df_1_pivot.values, 
    levels=[1/0.95, result_df_1_pivot.values.max()], colors=['#4AC6FF'], alpha=0.3) 

#axes formatting
ax.set_xlabel('NO Concentration (ppb$_v$)', fontsize=20)
ax.set_ylabel('HO$_2$ Concentration (ppt$_v$)', fontsize=20)
ax.set_xticks(np.sort(np.append(ax.get_xticks(), 0.1)))
ax.set_yticks(np.sort(np.append(ax.get_yticks(), 0)))
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlim(left=0.1)
ax.set_ylim(bottom = 0.1)

plt.show()
