# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 13:24:07 2025

@author: audlaw
"""
#plotting EESI data! includes SOA, NH4, and PTR (NH4 and PTR quantified)
#plots figure 3, plots for SOA yield/conc maxes, signal maxes, and SI figure for C7H8O6 from BnAld

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import re
import matplotlib.cm as cm
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

mpl.rcParams.update({'font.size': 20})

titles = "0425"

NH4_file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_NH4_ppb_r2.csv"
PTR_file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_PTR_ppb.csv"

#EESI_file_data =  "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_EESI_CES_PWLcorr.csv"
EESI_file_data ='C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/20220425_BnOH_EESI_AVL_PWLcorr.csv'
EESI_file_time =  "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_EESI_time.csv"
EESI_file_RO2 = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/workups_CES/20220425_BnOH_EESI_AVL.csv"
EESI_zeros = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/workups_CES/0425_EESI_zeros.csv"

OH_exp_file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SOAamssum.csv"
SOA_data = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SOAamssum.csv"
SOA_df = pd.read_csv(SOA_data, header = 0)

lighton = '2022-04-25 10:55:00'
lightoff = '2022-04-25 15:55:00'
#%%
def e_time_finder(time_list):#this function finds the elapsed time from datetime
    e_time_list = []
    for i in range(len(time_list)):
        seconds_elapsed = (time_list.iloc[i]-time_list.iloc[0]).total_seconds()
        e_time_list.append(seconds_elapsed/60)
    return e_time_list

def rolling_avg_2(reg_data):
    reg_data_1min = reg_data.resample('1min', on='time').mean()
    reg_data_2min = reg_data.resample('2min', on='time').mean()
    reg_data_5min = reg_data.resample('5min', on='time').mean()
    return(reg_data_1min,reg_data_2min,reg_data_5min)

class CustomScalarFormatter(ScalarFormatter):
    def __call__(self, x, pos=None):
        if x == 0.0:
            return "0"
        return super().__call__(x, pos)
 
def format_y_ticks(value, pos):
    if value == 0:
        return "0"
    return f"{value:.1f}"  # two decimal places for all other ticks



bnoh_ppb_data_source = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/20220425_benzyl_alcohol_50_200_bg_sub_r4.csv"
bnoh_ppb = pd.read_csv(bnoh_ppb_data_source)
bnoh_ppb = bnoh_ppb.loc[:, ['time (MDT)', 'C7H8O_ppb']]

lighton_dt = datetime.strptime(lighton, "%Y-%m-%d %H:%M:%S")
lightoff_dt = datetime.strptime(lightoff, "%Y-%m-%d %H:%M:%S")

bnoh_ppb['time (MDT)'] = pd.to_datetime(bnoh_ppb['time (MDT)'])

bnoh_ppb = bnoh_ppb.loc[bnoh_ppb['time (MDT)'] <=lightoff_dt]#removes data after lights off
bnoh_ppb = bnoh_ppb.loc[bnoh_ppb['time (MDT)'] >=lighton_dt]

bnoh_ppb.rename(columns={'time (MDT)': 'time'}, inplace=True)

bnoh_ppb_1min, bnoh_ppb_2min, bnoh_ppb_5min = rolling_avg_2(bnoh_ppb)
bnoh_ppb_1min.reset_index(inplace=True)

bnoh_ppb_1min['e_time'] = e_time_finder(bnoh_ppb_1min['time'])

def func(x,m,b):
    return(x*m+b)

OH_exp_data = pd.read_csv(OH_exp_file)
popt, pcov =curve_fit(func, OH_exp_data.dropna()['Time_Lon_min'], OH_exp_data.dropna()['photohr'], bounds = (0,2e11))
OH_exp_slope = popt[0]*3600*1.5E6#this is pulls OH exposure from jathar group 

NH4_ppb_df = pd.read_csv(NH4_file)
PTR_ppb_df = pd.read_csv(PTR_file)



NH4_ppb_df['date_time'] = pd.to_datetime(  NH4_ppb_df['date_time'])
PTR_ppb_df['t_start'] = pd.to_datetime(PTR_ppb_df['t_start'])
NH4_ppb_df['date_time'] = NH4_ppb_df['date_time'] - pd.to_timedelta(6, unit='h')#convert to local time
PTR_ppb_df['t_start'] = PTR_ppb_df['t_start'] - pd.to_timedelta(6, unit='h')#convert to local time
lighton_dt = datetime.strptime(lighton, "%Y-%m-%d %H:%M:%S") #converting from string to dt so I can subtract 5 min to get pre lighton levels
lighton_dt_minus_5 = lighton_dt - timedelta(minutes=5)
NH4_ppb_df_pre_light_on = NH4_ppb_df.loc[NH4_ppb_df['date_time'] >=lighton_dt_minus_5]
NH4_ppb_df_pre_light_on = NH4_ppb_df_pre_light_on.loc[NH4_ppb_df_pre_light_on['date_time'] <=lighton_dt]
PTR_ppb_df_pre_light_on = PTR_ppb_df.loc[PTR_ppb_df['t_start'] >=lighton_dt_minus_5]
PTR_ppb_df_pre_light_on = PTR_ppb_df_pre_light_on.loc[PTR_ppb_df_pre_light_on['t_start'] <=lighton_dt]

NH4_ppb_df = NH4_ppb_df.loc[NH4_ppb_df['date_time'] <=lightoff]#removes data after lights off
NH4_ppb_df = NH4_ppb_df.loc[NH4_ppb_df['date_time'] >=lighton]
    
PTR_ppb_df = PTR_ppb_df.loc[PTR_ppb_df['t_start'] <=lightoff]#removes data after lights off
PTR_ppb_df = PTR_ppb_df.loc[PTR_ppb_df['t_start'] >=lighton]
    
NH4_ppb_df.rename(columns={'date_time': 'time'}, inplace=True)
PTR_ppb_df.rename(columns={'t_start': 'time'}, inplace=True)


NH4_ppb_df.set_index('time', inplace=True)
PTR_ppb_df.set_index('time', inplace=True)
#NH4_ppb_df = NH4_ppb_df.resample('min').mean()


vapor_pressure_set = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/BnOH_vapor_pressures.csv"
vapor_pressures = pd.read_csv(vapor_pressure_set)

vapor_pressures['mw'] = np.nan
for index, row in vapor_pressures.iterrows():
    formula = row['species']
    formula_clean = re.sub(r'\(NH4\)', '', formula)
    atomic_weights = {
        'H': 1.0079,
        'C': 12.0107,
        'O': 15.9994,
        'N': 14.0067
    }

    matches = re.findall(r'([CHON])(\d+)', formula_clean)
    element_counts = {el: int(count) for el, count in matches}
    
    molecular_weight = sum(atomic_weights[el] * count for el, count in element_counts.items())
    vapor_pressures.at[index, 'mw'] = molecular_weight


def VWL_corrector_2(data_df, vp_df):
    for VOC in data_df.columns:
        vp_row = vp_df[vp_df['species'] == VOC]

        if vp_row.empty:
            print(f"[Skipping] '{VOC}' not found in vapor pressure table.")
            continue

        Pvap = vp_row['vp_Pa'].iloc[0]
        mw = vp_row['mw'].iloc[0]
        cstar = (Pvap * mw * 1e6) / (8.314 * 298)  # R = m3 Pa K-1 mol-1
 

        if pd.isna(Pvap):
            print(f'missing vapor pressure for {VOC}')
        else:
            # Estimate wall-phase concentration c_wall (mg/m³)
            if cstar >= 1e4:
                cwall = 10
            elif cstar >= 1:
                cwall = 16 * cstar**0.6 / 1000
            else:
                cwall = 16 / 1000

            # Rate constants
            k_vapon = 1.28e-3  # 1/s
            k_vapoff = k_vapon * cstar / (cwall * 1e3)

            wall_col = VOC + '_wlc'
            data_df[wall_col] = 0.0

            voc_series = data_df[VOC]/1.0E9*101325/9.314/298*mw*1000
            time_seconds = (data_df.index - data_df.index[0]).total_seconds().to_numpy()
            valid_mask = ~pd.isna(voc_series)
            
            # Filter voc_series and time_seconds to keep only valid entries
            voc_series_valid = voc_series[valid_mask].reset_index(drop=True)
            time_seconds_valid = time_seconds[valid_mask.to_numpy()]
            
            # Initialize the wall column with zeros for the valid indices only
            wall_col = VOC + '_wlc'
            data_df[wall_col] = 0.0

            valid_indices = data_df.index[valid_mask]

            for i in range(1, len(voc_series_valid)):
                dt = time_seconds_valid[i] - time_seconds_valid[i - 1]
    
                flux_to_wall = voc_series_valid.iat[i - 1] * k_vapon * dt
                flux_from_wall = data_df.at[valid_indices[i - 1], wall_col] * k_vapoff * dt
                new_wall_amount = data_df.at[valid_indices[i - 1], wall_col] + flux_to_wall - flux_from_wall
                data_df.at[valid_indices[i], wall_col] = max(0, new_wall_amount)
            data_df[wall_col] = data_df[wall_col]/1000/mw*298*8.314/101325*1.0E9
    return data_df



def e_time_finder_index(time_list):
    start_time = time_list[0]
    return [(t - start_time).total_seconds() / 60 for t in time_list]

correct_ptr = VWL_corrector_2(PTR_ppb_df, vapor_pressures)
correct_data = VWL_corrector_2(NH4_ppb_df, vapor_pressures)
time_e = e_time_finder_index(correct_data.index)
time_e_ptr = e_time_finder_index(correct_ptr.index)


for col in correct_data.columns:
    if col.endswith('_wlc'):
        base_col = col.replace('_wlc', '')
        if base_col in correct_data.columns:
            # Create the new total column
            total_col = base_col + '_total'
            correct_data[total_col] = correct_data[base_col] + correct_data[col]

for col in correct_ptr.columns:
    if col.endswith('_wlc'):
        base_col = col.replace('_wlc', '')
        if base_col in correct_ptr.columns:
            # Create the new total column
            total_col = base_col + '_total'
            correct_ptr[total_col] = correct_ptr[base_col] + correct_ptr[col]

        
corr_df = correct_data.copy()
corr_ptr = correct_ptr.copy()
corr_df = corr_df.resample('min').mean()
corr_ptr = corr_ptr.resample('min').mean()


NH4_ppb_df = NH4_ppb_df.resample('min').mean()#resamples to get the same time index for NH4 and ptr data
PTR_ppb_df = PTR_ppb_df.resample('min').mean()

PTR_ppb_df_pre_light_on = PTR_ppb_df_pre_light_on.drop('t_start', axis=1)#droping the time column so I can find the mean
NH4_ppb_df_pre_light_on = NH4_ppb_df_pre_light_on.drop('date_time', axis=1)

PTR_ppb_df_pre_light_on_mean = PTR_ppb_df_pre_light_on.mean()
NH4_ppb_df_pre_light_on_mean = NH4_ppb_df_pre_light_on.mean() 
      
PTR_ppb_df = PTR_ppb_df-PTR_ppb_df_pre_light_on_mean
NH4_ppb_df = NH4_ppb_df-NH4_ppb_df_pre_light_on_mean 

bnoh_ppb_1min = bnoh_ppb_1min.set_index('time')
NH4_ppb_df['C7H8O_ppb'] = bnoh_ppb_1min['C7H8O_ppb']

PTR_ppb_df = PTR_ppb_df.clip(lower=0)
NH4_ppb_df = NH4_ppb_df.clip(lower=0)
    

merged_df = pd.merge(NH4_ppb_df, PTR_ppb_df, left_index=True, right_index=True, how='outer')#combines NH4 ppb data and PTR ppb dataso stacked plots can be made
merged_df = merged_df.reset_index()
merged_df_tets = merged_df
e_time = e_time_finder(merged_df['time'])#finds elapsed time in order to get at OH exposure
OH_exp = np.array(e_time, dtype=float)*OH_exp_slope#calculates OH exposure 
#merged_df = merged_df.drop('time', axis=1)#gets rid of time axis
merged_df = merged_df.dropna(axis=1, how='all')#drops any columns that were all nan (these are species that we did not quantify)
merged_df.loc[merged_df.isna().any(axis=1)] = np.nan

corr_df_pre_light_on_mean = NH4_ppb_df_pre_light_on.mean() 

exclude_col = 'C7H8O_ppb_total'
corr_df['142.086255;(NH4)C7H8O2+_ppb_total'] = corr_df['142.086255;(NH4)C7H8O2+_ppb_total'] - 2.167282 #manually background subtracting HbnOH


PTR_ppb_df = PTR_ppb_df.clip(lower=0)
corr_ptr = corr_ptr.clip(lower = 0)
corr_df = corr_df.clip(lower=0)

merged_df_corr = pd.merge(corr_df, PTR_ppb_df, left_index=True, right_index=True, how='outer')#combines NH4 ppb data and PTR ppb dataso stacked plots can be made
merged_df_corr = merged_df_corr.reset_index()
e_time = e_time_finder(merged_df_corr['time'])#finds elapsed time in order to get at OH exposure
OH_exp = np.array(e_time, dtype=float)*OH_exp_slope#calculates OH exposure 
#merged_df = merged_df.drop('time', axis=1)#gets rid of time axis
merged_df_corr = merged_df_corr.dropna(axis=1, how='all')#drops any columns that were all nan (these are species that we did not quantify)
merged_df_corr.loc[merged_df_corr.isna().any(axis=1)] = np.nan

#%%

VWLdf = merged_df_corr.copy()

def is_valid_species(col):
    # formulas with format "mass;(...formula...)+_ppb_total"
    NH4_pattern = r"^\d+\.\d+;.*_ppb_total$"
    # formulas with format "mXX_formula_ppb"
    PTR_pattern = r"^m\d+_.+_ppb$"
    special_pattern = r"^m\d+_.+_ppb_matt$"
    bnoh_pattern = "C7H8O_ppb_total"
    
    return re.match(NH4_pattern, col) or re.match(PTR_pattern, col) or re.match(special_pattern, col) or re.match(bnoh_pattern, col)

exclude_species = ["124.075690;(NH4)C7H6O+_ppb_total",  "m107_C7H6OH_ppb", 'm91_C7H6H_ppb'] #taking these out because theyre duplicates of the BnAld signal



#species_cols = [c for c in VWLdf.columns if c != "time"]  #drop any non-time columns from grouping
species_cols = [c for c in merged_df_corr.columns if is_valid_species(c)]
species_cols = ["time"] + species_cols
VWLdf_filtered_all = merged_df_corr[species_cols].copy()

drop_cols = [
    "124.075690;(NH4)C7H6O+_ppb_total",
    "126.091340;(NH4)C7H8O+_ppb_total", 
    "m107_C7H6OH_ppb",
    "m33_CH4OH_ppb", 
    "m42_C2H3NH_ppb",
    "m123_C7H6O2H_ppb"]

VWLdf_filtered = VWLdf_filtered_all.drop(columns=drop_cols, errors="ignore").copy()


def normalize_species_name(col):
    # replace possible endings with "_ppb_VWL"
    col = re.sub(r"_ppb_total$", "_ppb_VWL", col)
    col = re.sub(r"_ppb_matt$", "_ppb_VWL", col)
    col = re.sub(r"_ppb$", "_ppb_VWL", col)
    return col

# apply to dataframe
VWLdf_filtered.rename(columns=normalize_species_name, inplace=True)

#%%

def get_carbon_number_vwl(col_name):
    """
    Extract number of carbons from VWL species columns.
    Handles:
    - NH4 format: "80.070605;(NH4)C2H6O2+_ppb_VWL"
    - PTR format: "m140_C6H5NO3H_ppb_VWL"
    """
    # Remove prefix if present (mass or mXX_)
    col_clean = col_name.split(";")[-1]  # for NH4 style
    col_clean = col_clean.split("_")[-2] # for PTR style, get formula part
    
    # Look for "C<number>"
    match = re.search(r"C(\d+)", col_name)
    if match:
        return int(match.group(1))
    else:
        return 0

#*********
#THIS PUTS THE VWL CORRECTED DF INTO PPBC
#*********
VWLdf_ppbC = VWLdf_filtered.copy()
for col in VWLdf_ppbC.columns:
    if col != "time":
        nC = get_carbon_number_vwl(col)
        VWLdf_ppbC[col] = VWLdf_ppbC[col] * nC

print(VWLdf_filtered)
print(VWLdf_ppbC)

#%%

#okay want to merge EESI csv into the corrected csv (GP data in ppbC)

eesi_data = pd.read_csv(EESI_file_data)
eesi_ro2 = pd.read_csv(EESI_file_RO2, header = 0)
eesi_time = pd.read_csv(EESI_file_time)
eesi_zeros = pd.read_csv(EESI_zeros)


# Attach and convert time column
eesi_time['time'] = pd.to_datetime(eesi_data['t_series'])
eesi_data['time'] = eesi_time['time']
eesi_ro2['time'] = eesi_time['time']


eesi_ro2['zeros'] = eesi_zeros['USER_DO0 monitor []']
eesi_ro2 = eesi_ro2[eesi_ro2["zeros"] == 0]

#eesi_data['zeros'] = eesi_zeros['USER_DO0 monitor []']
#eesi_data = eesi_data[eesi_data["zeros"] == 0]

plt.scatter(eesi_data['time'], eesi_data['197.042044;C7H10O5Na+_pwlcorr'])
plt.show()


# Restrict to light-on/off window
lighton_dt = datetime.strptime(lighton, "%Y-%m-%d %H:%M:%S")
lightoff_dt = datetime.strptime(lightoff, "%Y-%m-%d %H:%M:%S")

eesi_data = eesi_data.loc[
    (eesi_data['time'] >= lighton_dt) &
    (eesi_data['time'] <= lightoff_dt)]

eesi_ro2 = eesi_ro2.loc[
    (eesi_ro2['time'] >= lighton_dt) &
    (eesi_ro2['time'] <= lightoff_dt)]

plt.scatter(eesi_data['time'], eesi_data['197.042044;C7H10O5Na+_pwlcorr'])
plt.show()


# Keep only columns matching desired format
pattern1 = r"^m\d+_[A-Za-z0-9]+_cps_bg_corr_pwlcorr$"
pattern2 = r"^\d+\.\d+;C\d+[A-Za-z0-9]*\+_pwlcorr$"

eesi_data = eesi_data[['time'] + [c for c in eesi_data.columns if re.match(pattern1, c) or re.match(pattern2, c)]]

pattern_ro2 = r"^\d+\.\d+;[A-Za-z0-9]+\+$"
eesi_ro2 = eesi_ro2[['time'] + [c for c in eesi_ro2.columns if re.match(pattern_ro2, c)]]

#this is to smooth the data for the only EESI plots, I just have to do it here before the time column gets too wonky
eesi_smoothed_1min, eesi_smoothed_2min, eesi_smoothed_5min = rolling_avg_2(eesi_data)
eesi_smoothed_ro2_1min, eesi_smoothed_ro2_2min, eesi_smoothed_ro2_5min = rolling_avg_2(eesi_ro2)

# Set index and resample
eesi_data.set_index('time', inplace=True)
eesi_data = eesi_data.resample('min').mean()

eesi_ro2.set_index('time', inplace=True)
eesi_ro2 = eesi_ro2.resample('min').mean()

#plt.scatter(eesi_smoothed_2min.index, eesi_smoothed_2min['197.042044;C7H10O5Na+_pwlcorr'])
plt.scatter(eesi_smoothed_2min.index, eesi_smoothed_2min['243.011138;C7H8O8Na+_pwlcorr'], label = "C7H8O8")
#plt.scatter(eesi_smoothed_2min.index, eesi_smoothed_2min['211.021309;C7H8O6Na+_pwlcorr'])

plt.legend()
plt.xlabel("time")
plt.ylabel("EESI cps")
plt.show()

#%%
c3_cols = ['175.021309;C4H8NaO6+', '203.016223;C5H8NaO7+', '245.026788;C7H10O8Na+']
c4_cols = ['187.021309;C5H8O6Na+', '215.016223;C6H8O7Na+', '245.026788;C7H10O8Na+']
c5_cols = ['217.031873;C6H10NaO7+', '245.026788;C7H10O8Na+', '287.037353;C9H12NaO9+']
c7_cols = ['243.047523;C8H12O7Na+']

groups = {
    "C3": c3_cols,
    "C4": c4_cols,
    "C5": c5_cols,
    "C7": c7_cols
}

fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=False)
axes = axes.flatten()

for ax, (group_name, cols) in zip(axes, groups.items()):
    
    for col in cols:
        if col in eesi_smoothed_ro2_2min.columns:
            y = eesi_smoothed_ro2_2min[col]
            
            # Normalize (avoid divide by zero)
            y_norm = y / y.max() if y.max() != 0 else y
            
            # Extract formula for label (after semicolon)
            label = col.split(";")[1]
            
            ax.scatter(eesi_smoothed_ro2_2min.index, y_norm, label=label)
    
    ax.set_title(group_name)
    ax.set_ylabel("Normalized Intensity")
    ax.legend(fontsize=7)

for ax in axes:
    ax.set_xlabel("Time")

plt.tight_layout()
plt.show()
#%%

# Keep merged_df also as datetime index
VWLdf_ppbC['time'] = pd.to_datetime(VWLdf_ppbC['time'])
VWLdf_ppbC = VWLdf_ppbC.set_index('time')
VWLdf_ppbC = VWLdf_ppbC.resample('min').mean()


# Now merge on aligned datetime indexes
merged_all = pd.merge(VWLdf_ppbC, eesi_data, left_index=True, right_index=True, how='outer')

# Only reset index if you really need a time column
merged_all = merged_all.reset_index()

#%%


#Calculating OH_exp
merged_all['time'] = pd.to_datetime(merged_all['time'])
merged_all = merged_all.set_index('time')
e_time = (merged_all.index - merged_all.index[0]).total_seconds() / 60
OH_exp = e_time * OH_exp_slope
merged_all['OH_exp'] = OH_exp


#%%
#throwing SOA into that merged df

SOA_df['Time'] = pd.to_datetime(SOA_df['Time'])

#merged_df_aerosol = pd.merge(SOA_df, ams_df, left_on='Time', right_on='datetime', how='inner')#merged the SMSP (SOA) df with the AMS df so you can find the amount of C
merged_df_aerosol = SOA_df
merged_df_aerosol['organic_carbon ug_m3'] = (merged_df_aerosol['SOA_AVGpwl']/2.6) #finds amount of organic carbon with OM:OC ratio and mass of organic matter (OM/OC is 2.6 based off of scents 2 data)

moles_air_L3 = (0.83*1000)/(0.08206*298)
merged_df_aerosol['OC (ppb)'] = merged_df_aerosol['organic_carbon ug_m3']/(1.0E6*12)/moles_air_L3*1.0E9 #converting mass of organic carbon to ppb
to_merge_aerosol = merged_df_aerosol[['OC (ppb)','Time', 'SOAAMSAVGpwl_vocobs_yield']]
to_merge_aerosol = to_merge_aerosol.set_index('Time')

#******************************************************************************
#merging the VWL corrected df
#******************************************************************************

to_merge_aerosol = to_merge_aerosol.copy()
to_merge_aerosol.index = pd.to_datetime(to_merge_aerosol.index)
to_merge_aerosol_shifted = to_merge_aerosol.copy()
to_merge_aerosol_shifted.index = to_merge_aerosol_shifted.index + pd.Timedelta(minutes=10)

merged_all = merged_all.copy()
merged_all.index = pd.to_datetime(merged_all.index)

# Merge the OC (ppb) column into merged_all
merged_all_with_OC = pd.merge(
    merged_all, 
    to_merge_aerosol_shifted, 
    left_index=True, 
    right_index=True, 
    how="outer"   # or "inner" if you only want exact matches
)

merged_all_with_OC = merged_all_with_OC.loc[lighton:lightoff]

merged_all_with_OC["OC (ppb)"] = merged_all_with_OC["OC (ppb)"].interpolate() #merging aerosol data with gas phase data

# finding the OH exposure from the aerosol + gas phase data separately since it has additional zeros 
df_for_aerosol_OH = pd.DataFrame()
df_for_aerosol_OH['Time'] = merged_all_with_OC.index 
OH_exp_for_aerosol = np.array(e_time_finder(df_for_aerosol_OH['Time']), dtype=float) * OH_exp_slope 


plt.scatter(OH_exp_for_aerosol, merged_all_with_OC['197.042044;C7H10O5Na+_pwlcorr'])
#%%
#okay EESI and SOA

#for subscripts for labels
def convert_to_latex(formula):#this is used to generate subscripts for formulas in legend labels 
    subscript_dict = {
        '0': r'$_0$', '1': r'$_1$', '2': r'$_2$', '3': r'$_3$', '4': r'$_4$',
        '5': r'$_5$', '6': r'$_6$', '7': r'$_7$', '8': r'$_8$', '9': r'$_9$',  '+': r'$^+$'}  
    pattern = r'(\d+|\+)'
    formula_with_latex = re.sub(pattern, lambda match: ''.join(subscript_dict[digit] for digit in match.group()), formula)
    return formula_with_latex


# secondary x-axis: Photochemical Age
def OHexp_to_age(OH_exp):
    return OH_exp/1.5e6/3600

def age_to_OHexp(age):
    return age * 1.5e6 * 3600



def normalize(df_to_normalize):#this functin normalizes the data
    normalized_df = (df_to_normalize-df_to_normalize.min())/(df_to_normalize.max()-df_to_normalize.min())
    return normalized_df

normalized_eesi= normalize(eesi_smoothed_2min)
normalized_eesi_with_OC = normalized_eesi.merge(to_merge_aerosol_shifted, 
                             left_index=True, right_index=True, how='inner')
'''
normalized_eesi_with_OC = normalized_eesi_with_OC.loc[lighton:lightoff]
normalized_eesi_with_OC["OC (ppb)"] = normalized_eesi_with_OC["OC (ppb)"].interpolate() #merging aerosol data with gas phase data
'''

# --- interpolate SOA onto EESI time grid ---

# Make sure indices are datetime
normalized_eesi.index = pd.to_datetime(normalized_eesi.index)
to_merge_aerosol_shifted.index = pd.to_datetime(to_merge_aerosol_shifted.index)

# Reindex aerosol data onto EESI timestamps
soa_on_eesi_time = (
    to_merge_aerosol_shifted
    .reindex(normalized_eesi.index)
    .interpolate(method="time"))

# Merge interpolated SOA back into EESI dataframe
normalized_eesi_with_OC = normalized_eesi.copy()
normalized_eesi_with_OC["OC (ppb)"] = soa_on_eesi_time["OC (ppb)"]
normalized_eesi_with_OC["SOAAMSAVGpwl_vocobs_yield"] = soa_on_eesi_time["SOAAMSAVGpwl_vocobs_yield"]

# Apply lights-on window
normalized_eesi_with_OC = normalized_eesi_with_OC.loc[lighton:lightoff]

e_time = (normalized_eesi_with_OC.index - normalized_eesi_with_OC.index[0]).total_seconds() / 60
OHexp_normalizedEESI = e_time * OH_exp_slope



plt.scatter(normalized_eesi_with_OC.index, normalized_eesi_with_OC['197.042044;C7H10O5Na+_pwlcorr'])
#%%
plt.scatter(OHexp_normalizedEESI, normalized_eesi_with_OC['197.042044;C7H10O5Na+_pwlcorr'])
plt.show

#%%

c3_cols = ['175.021309;C4H8NaO6+_pwlcorr', '203.016223;C5H8NaO7+_pwlcorr', '245.026788;C7H10O8Na+_pwlcorr']
c4_cols = ['187.021309;C5H8O6Na+_pwlcorr', '199.021309;C6H8O6Na+_pwlcorr', '215.016223;C6H8O7Na+_pwlcorr', '245.026788;C7H10O8Na+_pwlcorr']
c5_cols = ['201.036959;C6H10O6Na+_pwlcorr', '217.031873;C6H10NaO7+_pwlcorr', '245.026788;C7H10O8Na+_pwlcorr']
c7_cols = ['243.047523;C8H12O7Na+_pwlcorr']
SOA = ['SOAAMSAVGpwl_vocobs_yield']

col_labels_ro2 = {
    '175.021309;C4H8NaO6+_pwlcorr': convert_to_latex("C4H8O6"),
    '203.016223;C5H8NaO7+_pwlcorr': convert_to_latex("C5H8O7"), 
    '245.026788;C7H10O8Na+_pwlcorr': convert_to_latex("C7H10O8"),
    '187.021309;C5H8O6Na+_pwlcorr': convert_to_latex("C5H8O6"), 
    '199.021309;C6H8O6Na+_pwlcorr': convert_to_latex("C6H8O6"),
    '201.036959;C6H10O6Na+_pwlcorr': convert_to_latex("C6H10O6"),
    '215.016223;C6H8O7Na+_pwlcorr': convert_to_latex("C6H8O7"),  
    '217.031873;C6H10NaO7+_pwlcorr': convert_to_latex("C6H10O7"),
    '287.037353;C9H12NaO9+_pwlcorr': convert_to_latex("C9H12O9"), 
    '243.047523;C8H12O7Na+_pwlcorr': convert_to_latex("C8H12O7"),
    'm85_C4H4O2H_ppb_VWL': convert_to_latex("C4H4O2"),
    '132.065520;(NH4)C5H6O3+_ppb_VWL': convert_to_latex("C5H6O3"),
    }

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, dpi = 300)
axes = axes.flatten()

groups = {
    "C3": c3_cols,
    "C4": c4_cols,
    "C5": c5_cols,
    "C7": c7_cols
}

for ax, (group_name, cols) in zip(axes, groups.items()):
    if "SOAAMSAVGpwl_vocobs_yield" in normalized_eesi_with_OC.columns:
        soa = normalized_eesi_with_OC["SOAAMSAVGpwl_vocobs_yield"]
        soa_norm = soa / soa.max() if soa.max() != 0 else soa
        
        soa_roll = soa_norm.rolling(window=10, center=True).mean()
        ax.plot(OHexp_normalizedEESI, soa_roll, color='gray', linewidth=6, alpha = 0.8)
        
        ax.scatter(OHexp_normalizedEESI, soa_norm,
                color='grey', linestyle='--', linewidth=0,
                label="SOA", marker="D", edgecolor="k", s=50)
        y = soa_norm.values
        x = OHexp_normalizedEESI
        
        # 40% relative error
        y_err = 0.40 * y
        y_low = y - y_err
        y_high = y + y_err
        
        # Shaded uncertainty
        ax.fill_between(x, y_low, y_high,
                             alpha=0.16, color="gray")

    for col in cols:
        if col in normalized_eesi_with_OC.columns:
            y = normalized_eesi_with_OC[col]
            
            y_norm = y / y.max() if y.max() != 0 else y
            
            label=col_labels_ro2.get(col, col)
            
            y_roll = y_norm.rolling(window=10, center=True).mean()
            ax.plot(OHexp_normalizedEESI, y_roll, linewidth=6, alpha = 0.8)

            ax.scatter(OHexp_normalizedEESI, y_norm, label=label)
            
    if group_name == "C4":
        col = "m85_C4H4O2H_ppb_VWL"  # change to your actual column name
        if col in merged_all.columns:
            y = merged_all[col]
            y_norm = y / y.max() if y.max() != 0 else y
            
            label=col_labels_ro2.get(col, col)
            
            ax.scatter(merged_all['OH_exp'], y_norm,
                    color='w', edgecolor = "#BFC0DE", linewidth=2, label=label)
            y_roll = y_norm.rolling(window=10, center=True).mean()
            ax.plot(merged_all['OH_exp'], y_roll, color="#BFC0DE", linewidth=6, alpha = 0.8)
            
            y_roll_csv_C4 = y_norm.to_frame(name="m85_C4H4O2H_ppb_VWL")
            y_roll_csv_C4["OH Exp"] = merged_all['OH_exp'].values
            #y_roll_csv_C4.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS19/20220425_C4H4O2.csv")
    
    
    if group_name == "C5":
        col = "132.065520;(NH4)C5H6O3+_ppb_VWL"  # change to your actual column name
        if col in merged_all.columns:
            y = merged_all[col]
            y_norm = y / y.max() if y.max() != 0 else y
            
            label=col_labels_ro2.get(col, col)
            
            ax.scatter(merged_all['OH_exp'], y_norm,
                    color='w', edgecolor = "#6F5BA8", linewidth=2, label=label)
            y_roll = y_norm.rolling(window=10, center=True).mean()
            ax.plot(merged_all['OH_exp'], y_roll, color="#6F5BA8", linewidth=6, alpha = 0.8)
            
            y_roll_csv_C5 = y_norm.to_frame(name="132.065520;(NH4)C5H6O3+_ppb_VWL")
            y_roll_csv_C5["OH Exp"] = merged_all["OH_exp"].values
            #y_roll_csv_C5.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS19/20220425_C5H6O3.csv",index=False)
    
    #ax.set_title(group_name)
    ax.set_ylabel("Normalized Intensity (c/s)", fontsize = 18)
    ax.legend(fontsize=12)
    ax.set_ylim(0, 1.1)
for ax in axes:
    ax.set_xlabel("OH exposure (molec cm$^{-3}$ s)", fontsize = 16)


plt.tight_layout(h_pad = 0, w_pad = 0.15)
plt.show()

normalized_eesi_with_OC_csv = normalized_eesi_with_OC.reset_index()
normalized_eesi_with_OC_csv["OH Exp"] = OHexp_normalizedEESI
#normalized_eesi_with_OC_csv.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS19/20220425_particle_phase.csv")
#merged_all.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS19/GP.csv")
#OHexp_normalizedEESI.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/MainText/Figure3/20220425_particle_phase_ohexp.csv")
#%%

tl_HO2_H = ["211.021309;C7H8O6Na+_pwlcorr", "213.036959;C7H10NaO6+_pwlcorr"]
tl_HO2_H_closed = ['227.016223;C7H8O7Na+_pwlcorr', '229.031873;C7H10O7Na+_pwlcorr']
tr_RO2_OH_ether = ['195.026394;C7H8NaO5+_pwlcorr', '197.042044;C7H10O5Na+_pwlcorr']
tr_RO2_OH_ether_closed = ['171.026394;C5H8NaO5+_pwlcorr', '187.021309;C5H8O6Na+_pwlcorr', '199.021309;C6H8O6Na+_pwlcorr', '201.036959;C6H10O6Na+_pwlcorr']
bl_NO_decomp_sq = ["154.995094;C4H4NaO5+_pwlcorr", "157.010744;C4H6NaO5+_pwlcorr"]
bl_NO_decomp = ['169.010744;C5H6NaO5+_pwlcorr', '171.026394;C5H8NaO5+_pwlcorr', "185.005659;C5H6NaO6+_pwlcorr", "187.021309;C5H8O6Na+_pwlcorr", "203.016223;C5H8NaO7+_pwlcorr"]
br_RO2_perox = ['175.021309;C4H8NaO6+_pwlcorr', "187.021309;C5H8O6Na+_pwlcorr", '203.016223;C5H8NaO7+_pwlcorr', '215.016223;C6H8O7Na+_pwlcorr', '217.031873;C6H10NaO7+_pwlcorr', '245.026788;C7H10O8Na+_pwlcorr', '243.047523;C8H12O7Na+_pwlcorr']
SOA_col = ["SOAAMSAVGpwl_vocobs_yield"]

col_labels = {
    "211.021309;C7H8O6Na+_pwlcorr": convert_to_latex("C7H8O6"), 
    "213.036959;C7H10NaO6+_pwlcorr": convert_to_latex("C7H10O6"), 
    '227.016223;C7H8O7Na+_pwlcorr': convert_to_latex("C7H8O7"), 
    '229.031873;C7H10O7Na+_pwlcorr': convert_to_latex("C7H10O7"),
    
    '195.026394;C7H8NaO5+_pwlcorr': convert_to_latex("C7H8O5"), 
    '197.042044;C7H10O5Na+_pwlcorr': convert_to_latex("C7H10O5"),
    '171.026394;C5H8NaO5+_pwlcorr': convert_to_latex("C5H8O5"), 
    '187.021309;C5H8O6Na+_pwlcorr': convert_to_latex("C5H8O6"),
    '199.021309;C6H8O6Na+_pwlcorr': convert_to_latex("C6H8O6"), 
    '201.036959;C6H10O6Na+_pwlcorr': convert_to_latex("C6H10O6"),
    
    "154.995094;C4H4NaO5+_pwlcorr": convert_to_latex("C4H4O5"),
    "157.010744;C4H6NaO5+_pwlcorr": convert_to_latex("C4H6O5"),
    '169.010744;C5H6NaO5+_pwlcorr': convert_to_latex("C5H6O5"), 
    #'171.026394;C5H8NaO5+_pwlcorr': convert_to_latex("C5H8O5"), 
    "185.005659;C5H6NaO6+_pwlcorr": convert_to_latex("C5H6O6"), 
    #"187.021309;C5H8O6Na+_pwlcorr": convert_to_latex("C5H8O6"), 
    "203.016223;C5H8NaO7+_pwlcorr": convert_to_latex("C5H8O7"),
    
    '175.021309;C4H8NaO6+_pwlcorr': convert_to_latex("C4H8O6"), 
    #"187.021309;C5H8O6Na+_pwlcorr": convert_to_latex("C5H8O6"), 
    #'203.016223;C5H8NaO7+_pwlcorr': convert_to_latex("C5H8O7"), 
    '215.016223;C6H8O7Na+_pwlcorr': convert_to_latex("C6H8O7"), 
    '217.031873;C6H10NaO7+_pwlcorr': convert_to_latex("C6H10O7"),
    '245.026788;C7H10O8Na+_pwlcorr': convert_to_latex("C7H10O8"),
    '243.047523;C8H12O7Na+_pwlcorr': convert_to_latex("C8H12O7"),
    "SOAAMSAVGpwl_vocobs_yield": convert_to_latex("SOA")
    }

colors_tl = {
    "211.021309;C7H8O6Na+_pwlcorr": "#EDAD0D", 
    "213.036959;C7H10NaO6+_pwlcorr": "#E8770E", 
    '227.016223;C7H8O7Na+_pwlcorr': "#4677DB", 
    '229.031873;C7H10O7Na+_pwlcorr': "#1B2C96",
    "SOAAMSAVGpwl_vocobs_yield": "#EBEBEB"
    }

colors_tr = {
    
    '195.026394;C7H8NaO5+_pwlcorr': "#52001A", 
    '197.042044;C7H10O5Na+_pwlcorr': "#8F002D",
    '171.026394;C5H8NaO5+_pwlcorr': "#B82207", 
    '187.021309;C5H8O6Na+_pwlcorr': "#DB522C",
    '199.021309;C6H8O6Na+_pwlcorr': "#E88564", 
    '201.036959;C6H10O6Na+_pwlcorr': "#FCCFC5",
    "SOAAMSAVGpwl_vocobs_yield": "#EBEBEB"
    }

colors_bl = {
    "154.995094;C4H4NaO5+_pwlcorr": "#510B6F",
    "157.010744;C4H6NaO5+_pwlcorr": "#4612A5",
    '169.010744;C5H6NaO5+_pwlcorr': "m", 
    '171.026394;C5H8NaO5+_pwlcorr': "#6723E7", 
    "185.005659;C5H6NaO6+_pwlcorr": "#875EFF", 
    "187.021309;C5H8O6Na+_pwlcorr": "#8712BA",
    "203.016223;C5H8NaO7+_pwlcorr": "#AD20E9",
    "SOAAMSAVGpwl_vocobs_yield": "#EBEBEB"
}

colors_br = {
    '175.021309;C4H8NaO6+_pwlcorr': "#142306", 
    "187.021309;C5H8O6Na+_pwlcorr": "#445A0C", 
    '203.016223;C5H8NaO7+_pwlcorr': "#87A314", 
    '215.016223;C6H8O7Na+_pwlcorr': "#B5D819", 
    '217.031873;C6H10NaO7+_pwlcorr': "#207410", 
    '245.026788;C7H10O8Na+_pwlcorr': "#3AA428",
    '243.047523;C8H12O7Na+_pwlcorr': "#72D161",
    "SOAAMSAVGpwl_vocobs_yield": "#EBEBEB"
    }


tl_desired_order = [
    convert_to_latex("C7H8O6"),
    convert_to_latex("C7H10O6"),
    convert_to_latex("C7H8O7"),
    convert_to_latex("C7H10O7"),
    convert_to_latex("SOA")]

tr_desired_order = [
    convert_to_latex("C7H8O5"),
    convert_to_latex("C7H10O5"),
    convert_to_latex("C5H8O5"),
    convert_to_latex("C5H8O6"),
    convert_to_latex("C6H8O6"),
    convert_to_latex("C6H10O6"),
    convert_to_latex("SOA")]


groups_tl_2 = [
    {"cols": tl_HO2_H, "marker": "o", "fill": "open"}]
groups_tl_1 = [
    {"cols": tl_HO2_H_closed, "marker": "o", "fill": "filled"}]

groups_tl = groups_tl_1 + groups_tl_2

groups_tr_2 = [
    {"cols": tr_RO2_OH_ether, "marker": "h", "fill": "open"}]
groups_tr_1 = [
    {"cols": tr_RO2_OH_ether_closed, "marker": "p", "fill": "filled"}]

groups_tr = groups_tr_1 + groups_tr_2

groups_bl = [
    {"cols": bl_NO_decomp_sq, "marker": "s"},
    {"cols": bl_NO_decomp, "marker": "^"}]
groups_br = [
    {"cols": br_RO2_perox, "marker": "v"}]



def plot_panel(ax, x, df_main, groups, df_soa, soa_cols, colors, labels):

    ax2 = ax.twinx()
    ax2.set_zorder(1)
    ax.set_zorder(2)
    ax.patch.set_visible(False)

    # -------------------
    # SOA (diamonds)
    # -------------------
    for col in soa_cols:
        y = df_soa[col]
        y_vals = y.values

        ax2.scatter(x, y,
                    color=colors.get(col, 'grey'),
                    marker="D", edgecolor="k", s=60, zorder=1,
                    label=labels.get(col, col))

        y_roll = y.rolling(window=10, center=True, min_periods=1).mean()
        ax2.plot(x, y_roll, color="k", linewidth=4, alpha=0.4, zorder=3)

        y_err = 0.40 * y_vals
        ax2.fill_between(x, y_vals - y_err, y_vals + y_err,
                         alpha=0.16, color="gray", zorder=2)

    # -------------------
    # MAIN SPECIES
    # -------------------
    for group in groups:
        cols = group["cols"]
        marker = group["marker"]
        fill = group.get("fill", "filled")  # default filled

        for col in cols:
            y = df_main[col]
            c = colors.get(col, "grey")

            # fill behavior
            if fill == "open":
                face = "white"
                edge = c
                lw = 2.2
            else:
                face = c
                edge = c
                lw = 1.5

            ax.scatter(x, y,
                       label=labels.get(col, col),
                       facecolor=face,
                       edgecolor=edge,
                       marker=marker,
                       linewidths=lw,
                       s=50,
                       zorder=2)

            # rolling
            y_roll = y.rolling(window=10, center=True, min_periods=1).mean()
            ax.plot(x, y_roll,
                    color=c,
                    linewidth=4,
                    alpha=0.8,
                    zorder=2)

    return ax2

fig, axs = plt.subplots(2, 2, figsize=(17, 10), dpi=300, sharex=True)
(ax_tl, ax_tr), (ax_bl, ax_br) = axs

ax_tl2 = plot_panel(ax_tl, OHexp_normalizedEESI,
                    normalized_eesi_with_OC,
                    groups_tl,
                    normalized_eesi_with_OC,
                    SOA_col,
                    colors_tl, col_labels)

ax_tr2 = plot_panel(ax_tr, OHexp_normalizedEESI,
                    normalized_eesi_with_OC,
                    groups_tr,
                    normalized_eesi_with_OC,
                    SOA_col,
                    colors_tr, col_labels)


ax_bl2 = plot_panel(ax_bl, OHexp_normalizedEESI,
                    normalized_eesi_with_OC,
                    groups_bl,
                    normalized_eesi_with_OC,
                    SOA_col,
                    colors_bl, col_labels)

ax_br2 = plot_panel(ax_br, OHexp_normalizedEESI,
                    normalized_eesi_with_OC,
                    groups_br,
                    normalized_eesi_with_OC,
                    SOA_col,
                    colors_br, col_labels)



sec_axes = []
for ax in [ax_tl, ax_tr, ax_bl, ax_br]:
    secax = ax.secondary_xaxis('top',
                              functions=(OHexp_to_age, age_to_OHexp))
    secax.tick_params(labelsize=18)
    sec_axes.append(secax)
sec_axes[2].tick_params(top=True, labeltop=False, length = 10)
sec_axes[3].tick_params(top=True, labeltop=False, length = 10)


#AXIS LIMITS

for ax in [ax_tl2, ax_tr2, ax_bl2, ax_br2]:
    ax.set_ylim(0, 1.05)
for ax in [ax_tl2, ax_tr2, ax_bl2, ax_br2]:
    ax.set_xlim(0, 1.14e11)
for ax2 in [ax_tl2, ax_tr2, ax_bl2, ax_br2]:
    ax2.set_ylim(0, 0.5)
for ax2 in [ax_tl2, ax_tr2, ax_bl2, ax_br2]:
    ax2.set_xlim(left=0)

handles = []
labels = []

for ax, ax2 in zip([ax_tl, ax_tr, ax_bl, ax_br],
                   [ax_tl2, ax_tr2, ax_bl2, ax_br2]):
    
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    
    handles = h1 + h2
    labels = l1 + l2

    if ax is ax_tr:
        order_map = {label: i for i, label in enumerate(tr_desired_order)}
        sorted_pairs = sorted(zip(handles, labels),
                              key=lambda x: order_map.get(x[1], 999))
        
        if sorted_pairs:
            handles, labels = zip(*sorted_pairs)
        else:
            handles, labels = [], []
            
    if ax is ax_tl:
        order_map = {label: i for i, label in enumerate(tl_desired_order)}
        sorted_pairs = sorted(zip(handles, labels),
                              key=lambda x: order_map.get(x[1], 999))
        
        if sorted_pairs:
            handles, labels = zip(*sorted_pairs)
        else:
            handles, labels = [], []

    ax.legend(handles, labels,
              loc="lower right",
              fontsize=15,
              ncol=3,
              frameon=True,
              handletextpad=0.01,
              labelspacing=0.3,
              borderpad=0.3)

ax_tr.tick_params(labelleft=False, length = 10)
ax_br.tick_params(labelleft=False, length = 10)
ax_tl2.tick_params(labelright=False, length = 10)
ax_tl.tick_params(labelright=False, length = 10)
ax_bl2.tick_params(labelright=False, length = 10)

formatter_sci = CustomScalarFormatter(useMathText=True)
formatter_sci.set_scientific(True)
formatter_sci.set_powerlimits((0, 0))  # force sci notation

for ax in [ax_bl, ax_br]:
    ax.xaxis.set_major_formatter(formatter_sci)

def clean_int(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:.0f}"  # no decimals

for secax in sec_axes:
    secax.xaxis.set_major_formatter(FuncFormatter(clean_int))
        
fig.text(0.5, 0.06, "OH exposure (molec cm$^{-3}$ s)", ha = "center", fontsize=19)
fig.text(0.05, 0.5, "Normalized EESI-ToFMS Signal", va='center', rotation='vertical', fontsize=19)
fig.text(0.94, 0.5, "SOA mass yield",
         va='center', rotation='vertical', fontsize=19)
fig.text(0.5, 0.93, "Photochemical age (hr)",
         ha='center', fontsize=19)


labels = ['(a)', '(b)', '(c)', '(d)']
axes = [ax_tl, ax_tr, ax_bl, ax_br]

for ax, lab in zip(axes, labels):
    ax.text(
        0.01, 0.97, lab,          # position (left, up)
        transform=ax.transAxes,    # axis-relative coords
        fontsize=22,
        fontweight='bold',
        va='top',
        ha='left'
    )

plt.tight_layout(rect=[0.06, 0.06, 0.94, 0.94])
plt.show

#%%
window = 5
# === this was for me to look at specific product maxes ===

peak_results = []

def plot_group(ax, groups, colors_dict, title):
    for group in groups:
        for col in group["cols"]:
            y_raw = normalized_eesi_with_OC[col]
            smooth = y_raw.rolling(window=window, center=True, min_periods=1).mean()

            # --- find peak ---
            max_idx = smooth.idxmax()
            max_val = smooth[max_idx]
            max_oh = OHexp_normalizedEESI[smooth.index.get_loc(max_idx)]
            max_age = OHexp_to_age(max_oh)

            label = col_labels.get(col, col)
            peak_results.append({
                "compound": label,
                "column": col,
                "group": title,
                "peak norm intensity": max_val,
                "OH_exposure at peak": max_oh,
                "age_hr at peak": max_age
            })


            # --- plot ---
            ax.plot(OHexp_normalizedEESI, smooth,
                    color=colors_dict.get(col, "gray"),
                    label=col_labels.get(col, col),
                    marker=group.get("marker", None),
                    markersize=4,
                    linewidth=1.5)

            # --- highlight peak ---
            ax.scatter(max_oh, max_val,
                       color=colors_dict.get(col, "gray"),
                       edgecolor='gold',
                       s=400,
                       zorder=5,
                       marker='*',
                       linewidth = 1.5)
    
    #ax.set_title(title)
    ax.legend(ncol=3)
    
fig, axs = plt.subplots(2, 2, figsize=(14, 12))

ax_tl, ax_tr = axs[0]
ax_bl, ax_br = axs[1]

plot_group(ax_tl, groups_tl, colors_tl, "TL: HO2 + H")
plot_group(ax_tr, groups_tr, colors_tr, "TR: RO2 / Ether")
plot_group(ax_bl, groups_bl, colors_bl, "BL: NO Decomposition")
plot_group(ax_br, groups_br, colors_br, "BR: RO2 Perox")

# labels
ax_bl.set_xlabel("OH Exposure (molec cm$^{-3}$ s)")
ax_br.set_xlabel("OH Exposure (molec cm$^{-3}$ s)")
ax_tl.set_ylabel("Normalized intensity")
ax_bl.set_ylabel("Normalized intensity")

plt.tight_layout()
plt.show()


peak_df = pd.DataFrame(peak_results)
print(peak_df)

peak_df.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_EESI_peak_contribs.csv")
#%%

tl = ["213.036959;C7H10NaO6+_pwlcorr"]
tl_closed = ['183.026394;C6H8O5Na+_pwlcorr', '245.026788;C7H10O8Na+_pwlcorr', '215.016223;C6H8O7Na+_pwlcorr']
SOA_col = ["SOAAMSAVGpwl_vocobs_yield"]

col_labels = {
    "183.026394;C6H8O5Na+_pwlcorr": convert_to_latex("C6H8O5"), 
    "213.036959;C7H10NaO6+_pwlcorr": convert_to_latex("C7H10O6"), 
    '215.016223;C6H8O7Na+_pwlcorr': convert_to_latex("C6H8O7"), 
    '245.026788;C7H10O8Na+_pwlcorr': convert_to_latex("C7H10O8"),
    
    
    "SOAAMSAVGpwl_vocobs_yield": convert_to_latex("SOA")
    }

colors_tl = {
    "183.026394;C6H8O5Na+_pwlcorr": "#B61A0F", 
    "213.036959;C7H10NaO6+_pwlcorr": "#E8770E",
    '215.016223;C6H8O7Na+_pwlcorr': "#FFD900",
    '245.026788;C7H10O8Na+_pwlcorr': "#ED620C", 
    
    "SOAAMSAVGpwl_vocobs_yield": "#EBEBEB"}




tl_desired_order = [
    
    convert_to_latex("C7H10O6"),
    convert_to_latex("C6H8O5"),
    convert_to_latex("C7H10O8"),
    convert_to_latex("SOA")]



groups_tl_2 = [
    {"cols": tl, "marker": "o", "fill": "open"}]
groups_tl_1 = [
    {"cols": tl_closed, "marker": "o", "fill": "filled"}]

groups_tl_hv = groups_tl_2 + groups_tl_1


def plot_panel(ax, x, df_main, groups, df_soa, soa_cols, colors, labels):

    ax2 = ax.twinx()
    ax2.set_zorder(1)
    ax.set_zorder(2)
    ax.patch.set_visible(False)

    # -------------------
    # SOA (diamonds)
    # -------------------
    for col in soa_cols:
        y = df_soa[col]
        y_vals = y.values

        ax2.scatter(x, y,
                    color=colors.get(col, 'grey'),
                    marker="D", edgecolor="k", s=60, zorder=1,
                    label=labels.get(col, col))

        y_roll = y.rolling(window=10, center=True, min_periods=1).mean()
        ax2.plot(x, y_roll, color="k", linewidth=4, alpha=0.4, zorder=3)

        y_err = 0.40 * y_vals
        ax2.fill_between(x, y_vals - y_err, y_vals + y_err,
                         alpha=0.16, color="gray", zorder=2)

    # -------------------
    # MAIN SPECIES
    # -------------------
    for group in groups:
        cols = group["cols"]
        marker = group["marker"]
        fill = group.get("fill", "filled")  # default filled

        for col in cols:
            y = df_main[col]
            c = colors.get(col, "grey")

            # fill behavior
            if fill == "open":
                face = "white"
                edge = c
                lw = 2.2
            else:
                face = c
                edge = c
                lw = 1.5

            ax.scatter(x, y,
                       label=labels.get(col, col),
                       facecolor=face,
                       edgecolor=edge,
                       marker=marker,
                       linewidths=lw,
                       s=50,
                       zorder=2)

            # rolling
            y_roll = y.rolling(window=10, center=True, min_periods=1).mean()
            ax.plot(x, y_roll,
                    color=c,
                    linewidth=4,
                    alpha=0.8,
                    zorder=2)

    return ax2

fig, axs = plt.subplots(1,1, figsize=(8, 5), dpi=300, sharex=True)
(ax_tl) = axs

ax_tl2 = plot_panel(ax_tl, OHexp_normalizedEESI,
                    normalized_eesi_with_OC,
                    groups_tl_hv,
                    normalized_eesi_with_OC,
                    SOA_col,
                    colors_tl, col_labels)



sec_axes = []
for ax in [ax_tl, ax_tr, ax_bl, ax_br]:
    secax = ax.secondary_xaxis('top',
                              functions=(OHexp_to_age, age_to_OHexp))
    secax.tick_params(labelsize=18)
    sec_axes.append(secax)
sec_axes[2].tick_params(top=True, labeltop=False, length = 10)
sec_axes[3].tick_params(top=True, labeltop=False, length = 10)


#AXIS LIMITS

for ax in [ax_tl2]:
    ax.set_ylim(0, 1.05)
for ax in [ax_tl2]:
    ax.set_xlim(0, 1.14e11)
for ax2 in [ax_tl2]:
    ax2.set_ylim(0, 0.5)
for ax2 in [ax_tl2]:
    ax2.set_xlim(left=0)

handles = []
labels = []

for ax, ax2 in zip([ax_tl],
                   [ax_tl2]):
    
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    
    handles = h1 + h2
    labels = l1 + l2


            
    if ax is ax_tl:
        order_map = {label: i for i, label in enumerate(tl_desired_order)}
        sorted_pairs = sorted(zip(handles, labels),
                              key=lambda x: order_map.get(x[1], 999))
        
        if sorted_pairs:
            handles, labels = zip(*sorted_pairs)
        else:
            handles, labels = [], []

    ax.legend(handles, labels,
              loc="lower right",
              fontsize=15,
              ncol=3,
              frameon=True,
              handletextpad=0.01,
              labelspacing=0.3,
              borderpad=0.3)



formatter_sci = CustomScalarFormatter(useMathText=True)
formatter_sci.set_scientific(True)
formatter_sci.set_powerlimits((0, 0))  # force sci notation

for ax in [ax_tl]:
    ax.xaxis.set_major_formatter(formatter_sci)

def clean_int(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:.0f}"  # no decimals

for secax in sec_axes:
    secax.xaxis.set_major_formatter(FuncFormatter(clean_int))
        
fig.text(0.5, 0.04, "OH exposure (molec cm$^{-3}$ s)", ha = "center", fontsize=19)
fig.text(0.001, 0.5, "Normalized EESI-ToFMS Signal", va='center', rotation='vertical', fontsize=19)
fig.text(0.97, 0.5, "SOA mass yield",
         va='center', rotation='vertical', fontsize=19)
fig.text(0.5, 0.97, "Photochemical age (hr)",
         ha='center', fontsize=19)


plt.tight_layout()
plt.show
#%%
# === Smoothing and max detection ===
window = 5  # rolling window size

# --- Dataset 1: SOA yield ---
soa_smooth_yield = (
    normalized_eesi_with_OC['SOAAMSAVGpwl_vocobs_yield']
    .rolling(window=window, center=True, min_periods=1)
    .mean()
)
OH_exp_yield = OH_exp_for_aerosol

max_idx_yield = soa_smooth_yield.idxmax()
max_soa_yield = soa_smooth_yield[max_idx_yield]
max_oh_yield = OH_exp_yield[soa_smooth_yield.index.get_loc(max_idx_yield)]
max_age_yield = OHexp_to_age(max_oh_yield)

end_soa_yield = soa_smooth_yield.iloc[-5]
end_oh_yield = OH_exp_yield[-5]
end_age_yield = OHexp_to_age(end_oh_yield)

print(f"Max SOA yield: {max_soa_yield:.4f} at OHexp={max_oh_yield:.3e}, age={max_age_yield:.2f} h")
print(f"End SOA yield: {end_soa_yield:.4f} at OHexp={end_oh_yield:.3e}, age={end_age_yield:.2f} h")


# --- Dataset 2: SOA concentration ---
soa_smooth_conc = (
    SOA_df['SOA_AVGpwl']
    .rolling(window=window, center=True, min_periods=1)
    .mean()
)
OH_exp_conc = SOA_df['OHexp']

max_idx_conc = soa_smooth_conc.idxmax()
max_soa_conc = soa_smooth_conc[max_idx_conc]
max_oh_conc = OH_exp_conc[soa_smooth_conc.index.get_loc(max_idx_conc)]
max_age_conc = OHexp_to_age(max_oh_conc)
'''
N = 10

max_idx = soa_smooth_conc.idxmax()
center_loc = soa_smooth_conc.index.get_loc(max_idx)

start = max(center_loc - N//2, 0)
end = min(center_loc + N//2 + 1, len(soa_smooth_conc))

window_idx = soa_smooth_conc.index[start:end]

max_soa_conc = soa_smooth_conc.loc[window_idx].mean()
max_oh_conc = OH_exp_conc.loc[window_idx].mean()
max_age_conc = OHexp_to_age(max_oh_conc)
'''
end_soa_conc = soa_smooth_conc.iloc[-16]
end_oh_conc = OH_exp_conc.iloc[-16]
end_age_conc = OHexp_to_age(end_oh_conc)

print("  ")
print(f"Max SOA conc: {max_soa_conc:.4f} at OHexp={max_oh_conc:.3e}, age={max_age_conc:.2f} h")
print(f"End SOA conc: {end_soa_conc:.4f} at OHexp={end_oh_conc:.3e}, age={end_age_conc:.2f} h")


# === Plot both datasets ===
fig, axes = plt.subplots(2, 1, figsize=(8, 10), sharex=False)

# --- (a) SOA yield ---
ax1 = axes[0]
ax1.scatter(OHexp_normalizedEESI, normalized_eesi_with_OC['SOAAMSAVGpwl_vocobs_yield'],
            label='Original SOA Yield', alpha=0.5)
ax1.scatter(OHexp_normalizedEESI, soa_smooth_yield, color='orange', label='Smoothed Yield')
ax1.scatter(max_oh_yield, max_soa_yield, color='red', zorder=5, label='Max SOA Yield')

ax1.set_ylabel("SOA Mass Yield")
ax1.legend()

ax1_top = ax1.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
ax1_top.set_xlabel("Photochemical Age (hours)")
ax1.set_xlabel("OH Exposure (molec cm$^{-3}$ s)")

# --- (b) SOA concentration ---
ax2 = axes[1]
ax2.scatter(OH_exp_conc, SOA_df['SOA_AVGpwl'], label='Original SOA Conc', alpha=0.5)
ax2.scatter(OH_exp_conc, soa_smooth_conc, color='orange', label='Smoothed Conc')
ax2.scatter(max_oh_conc, max_soa_conc, color='red', zorder=5, label='Max SOA Conc')
'''ax2.scatter(
    OH_exp_conc.loc[window_idx],
    soa_smooth_conc.loc[window_idx],
    color="none",
    edgecolor="red",
    linewidth=2,
    s=50,
    zorder=6,
    label="Peak window (avg)"
)'''

ax2.set_xlabel("OH Exposure (molec cm$^{-3}$ s)")
ax2.set_ylabel("SOA Concentration (µg m$^{-3}$)")
ax2.legend()

ax2_top = ax2.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
ax2_top.set_xlabel("Photochemical Age (hours)")

plt.tight_layout()
plt.show()


#%%
#this is for that SI figure that looks at time dependence of C7H8O6 formation
fig, ax1 = plt.subplots(figsize=(10,5), dpi = 300)

nh4_ptr_cols = [#"C7H8O_ppb_VWL", 
                #"m107_C7H6OH_ppb_VWL", 
                "142.086255;(NH4)C7H8O2+_ppb_VWL"]

eesi_cols    = [#"211.021309;C7H8O6Na+_pwlcorr", 
                "213.036959;C7H10NaO6+_pwlcorr",
                "245.026788;C7H10O8Na+_pwlcorr",
                "247.042438;C7H12O8Na+_pwlcorr",
                "217.031873;C6H10NaO7+_pwlcorr"]

col_labels = {
    #"C7H8O_ppb_VWL": "BnOH (NH4)", 
    #"m107_C7H6OH_ppb_VWL": "BnAld (PTR)",
    "142.086255;(NH4)C7H8O2+_ppb_VWL": "HBnOH (NH4)",
    "211.021309;C7H8O6Na+_pwlcorr": "EESI C7H8O6",
    "m213_C7H10NaO6_cps_bg_corr_pwlcorr": "EESI C7H10O6",
    "m245_C7H10O8Na_cps_bg_corr_pwlcorr": "EESI C7H8O8", 
    "247.042438;C7H12O8Na+_pwlcorr": "EESI C7H12O8",
    "217.031873;C6H10NaO7+_pwlcorr": "EESI C6H10O7"}

color_names = {
    #"C7H8O_ppb_VWL": "k", 
    #"m107_C7H6OH_ppb_VWL": "#B83030",
    "142.086255;(NH4)C7H8O2+_ppb_VWL": "#52A7F1",
    #"211.021309;C7H8O6Na+_pwlcorr": "#FFD500",
    "213.036959;C7H10NaO6+_pwlcorr": "darkorange",
    "245.026788;C7H10O8Na+_pwlcorr": "c", 
    "247.042438;C7H12O8Na+_pwlcorr": "m",
    "217.031873;C6H10NaO7+_pwlcorr": "#3153B6"}

# Plot NH4/PTR
for col in nh4_ptr_cols:
    ax1.scatter(merged_all['OH_exp'], merged_all[col], label=col_labels.get(col, col), color=color_names.get(col, 'grey'), s = 50)

# Plot EESI on secondary y-axis
ax2 = ax1.twinx()
for col in eesi_cols:
    ax2.scatter(merged_all['OH_exp'], merged_all[col], label=col_labels.get(col, col), color=color_names.get(col, 'grey'), marker="+", s = 50)

# Labels
ax1.set_xlabel("OH exposure (molec cm$^{-3}$ s)", fontsize=20)
ax1.set_ylabel("Concentration (ppbC)", fontsize=20)
ax2.set_ylabel("EESI intenstity (cps)", fontsize=20)

# Legends outside the plot
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

plt.legend(lines1 + lines2, labels1 + labels2, loc='center left', bbox_to_anchor=(0.16, 0.2), frameon = False, fontsize = 15)


ax1.tick_params(axis='both', labelsize=23)
ax2.tick_params(axis='both', labelsize=23)

plt.xlim(0, 1.13e11)
plt.ylim(0)
plt.tight_layout()  # makes room for legend
plt.show()


#%%

#this is for that SI figure that looks at time dependence of C7H8O6 formation
fig, ax1 = plt.subplots(figsize=(10,7), dpi = 300)

nh4_ptr_cols = [#"C7H8O_ppb_VWL", 
                "m107_C7H6OH_ppb_VWL", 
                "142.086255;(NH4)C7H8O2+_ppb_VWL"]

eesi_cols    = ["211.021309;C7H8O6Na+_pwlcorr", 
                #"213.036959;C7H10NaO6+_pwlcorr",
                #"245.026788;C7H10O8Na+_pwlcorr",
                "229.031873;C7H10O7Na+_pwlcorr"]

col_labels = {
    #"C7H8O_ppb_VWL": "BnOH (NH4)", 
    "m107_C7H6OH_ppb_VWL": convert_to_latex("C7H6O (PTR)"),
    "142.086255;(NH4)C7H8O2+_ppb_VWL": convert_to_latex("C7H8O2 (NH4)"),
    "211.021309;C7H8O6Na+_pwlcorr": convert_to_latex("C7H8O6 (EESI)"),
    #"m213_C7H10NaO6_cps_bg_corr_pwlcorr": "EESI C7H10O6",
    #"m245_C7H10O8Na_cps_bg_corr_pwlcorr": "EESI C7H8O8", 
    "229.031873;C7H10O7Na+_pwlcorr": convert_to_latex("C7H10O7 (EESI)")}

color_names = {
    #"C7H8O_ppb_VWL": "k", 
    "m107_C7H6OH_ppb_VWL": "#B83030",
    "142.086255;(NH4)C7H8O2+_ppb_VWL": "#52A7F1",
    "211.021309;C7H8O6Na+_pwlcorr": "#FFD500",
    #"213.036959;C7H10NaO6+_pwlcorr": "darkorange",
    #"245.026788;C7H10O8Na+_pwlcorr": "darkorange", 
    "229.031873;C7H10O7Na+_pwlcorr": "#3153B6"}

# Plot NH4/PTR
for col in nh4_ptr_cols:
    ax1.scatter(merged_all['OH_exp'], merged_all[col], label=col_labels.get(col, col), color=color_names.get(col, 'grey'), s = 70)

# Plot EESI on secondary y-axis
ax2 = ax1.twinx()
for col in eesi_cols:
    ax2.scatter(merged_all['OH_exp'], merged_all[col], label=col_labels.get(col, col), color=color_names.get(col, 'grey'), marker="+", s = 70)

# Labels
ax1.set_xlabel("OH exposure (molec cm$^{-3}$ s)", fontsize=20)
ax1.set_ylabel("Gas-Phase Concentration (ppbC)", fontsize=20)
ax2.set_ylabel("EESI-ToF-MS intenstity (cps)", fontsize=20)

# Legends outside the plot
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

plt.legend(lines1 + lines2, labels1 + labels2, loc='center left', bbox_to_anchor=(0.16, 0.2), frameon = True, fontsize = 19)


ax1.tick_params(axis='both', labelsize=23)
ax2.tick_params(axis='both', labelsize=23)

plt.xlim(0, 1.13e11)
plt.ylim(0)
plt.tight_layout()  # makes room for legend
plt.show()

#merged_all.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS17/20220425_BnAld_Hshift.csv")