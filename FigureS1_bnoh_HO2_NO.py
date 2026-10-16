# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 12:06:04 2026

@author: audlaw
"""

#this plots figure S1
#adapted from tucker melles' script to plot HO2 and NOx

import pandas as pd
import re
from pathlib import Path
import matplotlib.pyplot as plt
import glob
import os
from datetime import datetime



#mean OH concentrations (molec/cm3) for each exp, derived from the OH_exp slope in huiyings' soa data
OH_values = {
    "20220328": 6.9e6,
    "20220425": 6.4e6,
    "20220706": 9.1e6,
    "20220708": 1.1e7,
    "20220722": 7.7e6
}




folder_path = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/NOx_and_HO2_finder/NOx_allexps"

csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

def convert_datetime(datetime_str):
    # Parse the string to a datetime object with timezone info
    dt = datetime.fromisoformat(datetime_str)
    return dt.strftime('%Y/%m/%d %H:%M')

# Dictionary to store DataFrames with date as key
dataframes = {}

# Read each CSV file into a DataFrame and store it in the dictionary
for file in csv_files:
    filename = os.path.basename(file)
    date_part = filename.split('_')[0]
    
    df = pd.read_csv(file)
    
    # Combine existing Date and Time into datetime
    df['actual_datetime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'], format='%Y/%m/%d %H:%M')
    
    dataframes[date_part] = df
    
def calculate_elapsed_time_from_first(df, datetime_column):
    df_copy = df.copy()
    df_copy[datetime_column] = pd.to_datetime(df_copy[datetime_column])
    first_time = df_copy[datetime_column].iloc[0]
    elapsed_time = df_copy[datetime_column] - first_time
    elapsed_time_seconds = elapsed_time.dt.total_seconds()
    
    return elapsed_time_seconds/60


exp_labels = {
    "20220328": "Exp. 1",
    "20220425": "Exp. 2",
    "20220706": "Exp. 3",
    "20220708": "Exp. 4",
    "20220722": "Exp. 5"
}

exp_colors = {
    "20220328": "#DBE5FA",
    "20220425": "#94AEF0",
    "20220706": "#5F83E7",
    "20220708": "#1D49C3",
    "20220722": "#0F2563"
}

times = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/NOx_and_HO2_finder/bnoh_experiment_timing.csv")

times['datetime_on'] = pd.to_datetime(times['date'].astype(str) + ' ' + times['timeon'])
times['datetime_off'] = pd.to_datetime(times['date'].astype(str) + ' ' + times['timeoff'])

NO_info = pd.DataFrame({
    'exp date': pd.Series(dtype='int'),
    'mean NO': pd.Series(dtype='float'),
    'min NO': pd.Series(dtype='float'),
    'max NO': pd.Series(dtype='float'),
    'avg NO (0–1hr)': pd.Series(dtype='float'),
    'avg NO (0–2hr)': pd.Series(dtype='float')
})





#%%

#HO2

file_path = Path("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/NOx_and_HO2_finder/HO2_BnOH_GECKO_11142025_7GEN.csv")


def resample_HO2_data(HO2_df, time_col, value_col, frequency='min'):
    HO2_df[time_col] = pd.to_timedelta(HO2_df[time_col], unit='s')
    HO2_df.set_index(time_col, inplace=True)
    HO2_resampled = HO2_df.resample(frequency).mean()
    HO2_resampled[value_col] = HO2_resampled[value_col].interpolate(method='linear')
    return HO2_resampled


def HO2_printer(file_path):
    HO2_df = pd.read_csv(file_path)

    # identify HO2 + time columns
    ho2_cols = [c for c in HO2_df.columns if "HO2" in c]
    time_cols = [c for c in HO2_df.columns if "Time" in c]

    out = []

    for ho2_col in ho2_cols:
        date_match = re.findall(r"\d{8}", ho2_col)
        if not date_match:
            continue

        date_str = date_match[0]
        time_col = next((t for t in time_cols if date_str in t), None)
        if time_col is None:
            continue

        tmp = HO2_df[[time_col, ho2_col]].copy()
        tmp.columns = ["time_sec", "HO2"]

        tmp["exp_date"] = date_str
        out.append(tmp)

    return pd.concat(out, ignore_index=True)
    
 
    HO2_info = pd.DataFrame(columns=[
        'exp date', 'mean HO2', 'min HO2', 'max HO2',
        'resampled mean HO2', 'HO2 at 1hr', 'HO2 at 2hr', 'HO2 at 3hr',
        'avg HO2 (0–1hr)', 'avg HO2 (0–2hr)', 'avg HO2 (0–3hr)',
        'HO2 at 50% completion', 'OH exposure (50%)', 'Photochem age (days, 50%)'
    ])

ho2_data = HO2_printer(file_path)
ho2_data['exp_date'] = ho2_data['exp_date'].apply(
    lambda x: x[4:] + x[:4]   # MMDDYYYY → YYYYMMDD
)


#%%
#this version plots and then saves a combined df of all the exps to csv
all_exp_data = []

fig, axs = plt.subplots(2, 2, figsize=(14, 10), dpi=400, sharex=True)

ax_tl, ax_tr = axs[0]
ax_bl, ax_br = axs[1]

order = ["20220328", "20220425", "20220706", "20220708", "20220722"]

for exp, on, off in zip(times['date'], times['datetime_on'], times['datetime_off']):

    exp = str(exp)

    exp_data = dataframes[exp].copy()
    exp_data['actual_datetime'] = pd.to_datetime(exp_data['actual_datetime'])

    after_lights_on = exp_data[exp_data['actual_datetime'] > on]
    only_lights_on = after_lights_on[
        after_lights_on['actual_datetime'] < off
    ].copy()

    only_lights_on['elapsed_sec'] = (
        only_lights_on['actual_datetime'] - on
    ).dt.total_seconds()

    OH_conc = OH_values[exp]
    only_lights_on['OH_exp'] = only_lights_on['elapsed_sec'] * OH_conc

    only_lights_on['NO/NO2'] = (
        only_lights_on['NO Calibrated [ppbv]'] /
        only_lights_on['NO2 Calibrated [ppbv]'])

    
    only_lights_on['experiment'] = exp
    all_exp_data.append(only_lights_on)

    #NO
    ax_tr.scatter(
        only_lights_on['OH_exp'],
        only_lights_on['NO Calibrated [ppbv]'],
        label=exp_labels[exp],
        color=exp_colors[exp])

    #NO2
    ax_br.scatter(
        only_lights_on['OH_exp'],
        only_lights_on['NO2 Calibrated [ppbv]'],
        color=exp_colors[exp])

    #NO/NO2
    only_lights_on['NO/NO2'] = (
        only_lights_on['NO Calibrated [ppbv]'] /
        only_lights_on['NO2 Calibrated [ppbv]'])

    ax_bl.scatter(
        only_lights_on['OH_exp'],
        only_lights_on['NO/NO2'],
        color=exp_colors[exp])

    #HO2
    HO2_exp = ho2_data[ho2_data['exp_date'] == exp].copy()
    
    HO2_OH_exp = HO2_exp['time_sec'] * OH_values[exp]
    
    ax_tl.scatter(
        HO2_OH_exp,
        HO2_exp['HO2'],
        color=exp_colors[exp])

#axis limits for all
for ax in [ax_tl, ax_tr, ax_bl, ax_br]:
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

#labellingggggg
ax_bl.set_xlabel("OH exposure (molec cm$^{-3}$ s)", fontsize=18)
ax_br.set_xlabel("OH exposure (molec cm$^{-3}$ s)", fontsize=18)

ax_tl.set_ylabel("[HO$_2$] (ppb$_v$)", fontsize=18)
ax_tr.set_ylabel("[NO] (ppb$_v$)", fontsize=18)
ax_bl.set_ylabel("NO/NO$_2$", fontsize=18)
ax_br.set_ylabel("[NO$_2$] (ppb$_v$)", fontsize=18)

#throw pAge in there
OH_atm = 1.5e6

def exp_to_age(x):
    return x / (OH_atm * 3600)

def age_to_exp(x):
    return x * (OH_atm * 3600)

axes = [ax_tl, ax_tr, ax_bl, ax_br]

secaxes = [
    ax.secondary_xaxis('top', functions=(exp_to_age, age_to_exp))
    for ax in axes]

secaxes[0].tick_params(labeltop=True)
secaxes[1].tick_params(labeltop=True)
secaxes[2].tick_params(labeltop=False)
secaxes[3].tick_params(labeltop=False)

secaxes[0].set_xlabel("Photochemical age (hrs)", fontsize=18)
secaxes[1].set_xlabel("Photochemical age (hrs)")

for secax in secaxes[2:]:
    secax.set_xlabel("")

#org labels
handles, labels = ax_tr.get_legend_handles_labels()
map_h = dict(zip(labels, handles))

ordered_handles = [map_h[exp_labels[e]] for e in order]
ordered_labels = [exp_labels[e] for e in order]

fig.legend(
    ordered_handles,
    ordered_labels,
    loc="upper center",
    ncol=5,
    fontsize=18,
    frameon=False)

plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.show()  


    
    
#this is where they save to csvs 
combined_df = pd.concat(all_exp_data, ignore_index=True)

#combined_df.to_csv("all_experiments_NO_NO2.csv", index=False)

all_ho2 = []

for exp in order:

    HO2_exp = ho2_data[
        ho2_data['exp_date'].astype(str) == exp
    ].copy()

    HO2_exp['OH_exp'] = (
        HO2_exp['time_sec'] * OH_values[exp]
    )

    HO2_exp['experiment'] = exp
    HO2_exp['experiment_label'] = exp_labels[exp]

    all_ho2.append(HO2_exp)

combined_ho2 = pd.concat(all_ho2, ignore_index=True)

#combined_ho2.to_csv("all_experiments_HO2.csv", index=False)