# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 10:57:19 2026

@author: audlaw
"""

import pandas as pd
import matplotlib.pyplot as plt


#this plots figure S6

carb_dist_time = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_gecko_carb_dist_time.csv"
carb_aer = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_gecko_carbon_chain_aer_ppbC_time.csv"
carb_gas = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_gecko_carbon_chain_gas_ppbC_time.csv"
oh_exposure = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_gecko_oh_exposure.csv"

total_carb = pd.read_csv(carb_dist_time, header=0).dropna()
C3_aer = pd.read_csv(carb_aer, header=0).dropna()
C3_gas = pd.read_csv(carb_gas, header=0).dropna()
ohexposure = pd.read_csv(oh_exposure, header=0).dropna()

#%%

time_col = [c for c in C3_gas.columns if 'Time' in c][0]

carbon_cols = ['            1', '            2', '            3', '            4',
       '            5', '            6', '            7']


time_col = [c for c in total_carb.columns if 'Time' in c][0]

cols = [' Precursor', ' Gas Phase', ' Aerosol Phase', ' CO+CO2', ' Aqueous Phase', ' Wall']

#%%


time_col_gas = [c for c in C3_gas.columns if 'Time' in c][0]
time_col_total = [c for c in total_carb.columns if 'Time' in c][0]

oh_exp = ohexposure[' cumulative OH exposure [molec cm-3 s]                                                                                         ']


carbon_cols = [
    '            1',
    '            2',
    '            3',
    '            4',
    '            5',
    '            6',
    '            7'
]

#order the columns
layers = [
    ('Benzyl alcohol', total_carb[' Precursor']),
    ('Gas-phase C$_1$', C3_gas[carbon_cols[0]]),
    ('Gas-phase C$_2$', C3_gas[carbon_cols[1]]),
    ('Gas-phase C$_3$', C3_gas[carbon_cols[2]]),
    ('Gas-phase C$_4$', C3_gas[carbon_cols[3]]),
    ('Gas-phase C$_5$', C3_gas[carbon_cols[4]]),
    ('Gas-phase C$_6$', C3_gas[carbon_cols[5]]),
    ('Gas-phase C$_7$', C3_gas[carbon_cols[6]]),
    ('CO+CO$_2$', total_carb[' CO+CO2']),
    ('Aerosol Phase', total_carb[' Aerosol Phase']),
    ('Aqueous Phase', total_carb[' Aqueous Phase']),
    ('Wall', total_carb[' Wall'])
]


#a platypus???
colors = [
    'k',       # Precursor
    '#E6FFFF',        # C1
    '#ADE0E0',     # C2
    '#90D0D0',       # C3
    '#73C0C0',      # C4
    '#57B0B0',       # C5
    '#3AA0A0',       # C6
    '#008080',     # C7
    '#F58700',       # CO + CO2
    '#FFBC85',      # aerosol
    '#F8ECB4',       # Aqueous
    '#633C2C'       # Wall
]

fig, ax = plt.subplots(figsize=(11, 7))

ax.stackplot(
    oh_exp,
    *[values for name, values in layers],
    labels=[name for name, values in layers],
    colors=colors
)



ax.set_xlim(0, 1.11e11)
ax.set_ylim(0,265)
ax.set_xlabel('OH Exposure (molecules cm$^{-3}$ s)')
ax.set_ylabel('Concentration of Carbon (ppbC)')
#ax.set_title('Fate of carbon predicted by GECKO-A')

handles, labels = ax.get_legend_handles_labels()

ax.legend(
    handles[::-1],
    labels[::-1],
    loc='center left',
    bbox_to_anchor=(1.02, 0, 0.3, 1),
    bbox_transform=ax.transAxes,
    mode='expand',
    frameon=True,
    fontsize = 15
)

fig.subplots_adjust(right=0.75)

OH_atm = 1.5e6  # molec cm^-3

secax = ax.secondary_xaxis(
    'top',
    functions=(
        lambda OHexp: OHexp / (OH_atm * 3600),   
        lambda age_hr: age_hr * OH_atm * 3600
    )
)

secax.set_xlabel('Photochemical age (hours)')

plt.tight_layout()
plt.show()
