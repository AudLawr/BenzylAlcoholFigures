# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 09:30:40 2025

@author: audlaw
"""

#VERSION SEVEN BAYBEEEEEEEEEEE
#version eight i should do scents1/2 iterations
#going to go through and comment and clean this sucker so it makes sense to someone else
#and lowkey also to me

#move things around to work for new fig grouping
#fix zero BS and VWL and errorr on at least C# stacked plot, probs on pathway plot too

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import re
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter

mpl.rcParams.update({'font.size': 20})
titles = "0425" #this is to distinguish experiments

#load in files. could do this with a json?
NH4_file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_NH4_ppb_r3.csv" #NH4 quantified with vscans
PTR_file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_PTR_ppb.csv" #PTR quantified with kptr
OH_exp_file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SOAamssum.csv"
SOA_data = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SOAamssum.csv"
bnoh_ppb_data_source = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/20220425_benzyl_alcohol_50_200_bg_sub_r4.csv" #pVOC quantified with liquid cal
vapor_pressure_set = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/BnOH_vapor_pressures.csv"
relative_errors = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_BnOH_linear_errors.csv"

#make dataframes
NH4_ppb_df = pd.read_csv(NH4_file)
PTR_ppb_df = pd.read_csv(PTR_file)
OH_exp_data = pd.read_csv(OH_exp_file)
SOA_df = pd.read_csv(SOA_data, header = 0)
bnoh_ppb = pd.read_csv(bnoh_ppb_data_source)
vapor_pressures = pd.read_csv(vapor_pressure_set)
rel_err_df = pd.read_csv(relative_errors)

#load in lightson/off then convert to datetime
lighton = '2022-04-25 10:50:00'
lightoff = '2022-04-25 15:55:00'
lighton_dt = datetime.strptime(lighton, "%Y-%m-%d %H:%M:%S") 
lightoff_dt = datetime.strptime(lightoff, "%Y-%m-%d %H:%M:%S")

#calculate prelighton average
lighton_dt_minus_5 = lighton_dt - timedelta(minutes=5)

#find elapsed time from datetime for OH exposure and VWL estimates
def elapsed_minutes(time_list):
    # if it's a pandas Series, convert to a plain list of datetimes
    if isinstance(time_list, pd.Series):
        time_list = time_list.tolist()
    
    start_time = time_list[0]
    return [(t - start_time).total_seconds() / 60 for t in time_list]

#calculates a rolling average on concentration data to make it visually cleaner (1min, 2min, and 5min average)
def rolling_avg_2(reg_data):
    reg_data_1min = reg_data.resample('1min', on='time').mean()
    reg_data_2min = reg_data.resample('2min', on='time').mean()
    reg_data_5min = reg_data.resample('5min', on='time').mean()
    return(reg_data_1min,reg_data_2min,reg_data_5min)

#linear function to use for OH exposure calculations
def func(x,m,b):
    return(x*m+b)


#for subscripts for labels
def convert_to_latex(formula):#this is used to generate subscripts for formulas in legend labels 
    subscript_dict = {
        '0': r'$_0$', '1': r'$_1$', '2': r'$_2$', '3': r'$_3$', '4': r'$_4$',
        '5': r'$_5$', '6': r'$_6$', '7': r'$_7$', '8': r'$_8$', '9': r'$_9$',  '+': r'$^+$'}  
    pattern = r'(\d+|\+)'
    formula_with_latex = re.sub(pattern, lambda match: ''.join(subscript_dict[digit] for digit in match.group()), formula)
    return formula_with_latex

class CustomScalarFormatter(ScalarFormatter):
    def __call__(self, x, pos=None):
        if x == 0.0:
            return "0"
        return super().__call__(x, pos)
 
def format_y_ticks(value, pos):
    if value == 0:
        return "0"
    return f"{value:.1f}"  # two decimal places for all other ticks


#calculate OH exposure slope from Huiying's SOA data file
OH_exp_data = pd.read_csv(OH_exp_file)
popt, pcov =curve_fit(func, OH_exp_data.dropna()['Time_Lon_min'], OH_exp_data.dropna()['photohr'], bounds = (0,2e11))
OH_exp_slope = popt[0]*3600*1.5E6#this is pulls OH exposure from jathar group 

#********************************
#load in PTR and NH4 dataframes and cutoff data outside lights on
NH4_ppb_df['date_time'] = pd.to_datetime(  NH4_ppb_df['date_time'])
NH4_ppb_df['date_time'] = NH4_ppb_df['date_time'] - pd.to_timedelta(6, unit='h')#convert to local time

NH4_ppb_df_pre_light_on = NH4_ppb_df.loc[NH4_ppb_df['date_time'] >=lighton_dt_minus_5] #cuts off data before 5min before LO to end of df
NH4_ppb_df_pre_light_on = NH4_ppb_df_pre_light_on.loc[NH4_ppb_df_pre_light_on['date_time'] <=lighton_dt] #cuts off data after lights on

NH4_ppb_df = NH4_ppb_df.loc[NH4_ppb_df['date_time'] <=lightoff]#removes data after lights off
NH4_ppb_df = NH4_ppb_df.loc[NH4_ppb_df['date_time'] >=lighton]

NH4_ppb_df.rename(columns={'date_time': 'time'}, inplace=True)
NH4_ppb_df.set_index('time', inplace=True)

#resample on 1min, calculate the mean for prelighton, then subtract from dataframe to account for starting concentrations before oxidation
NH4_ppb_df = NH4_ppb_df.resample('T').mean()   #resamples to get the same time index for NH4 and ptr data
NH4_ppb_df_pre_light_on = NH4_ppb_df_pre_light_on.drop('date_time', axis=1)  #droping the time column so I can find the mean
NH4_ppb_df_pre_light_on_mean = NH4_ppb_df_pre_light_on.mean() 
NH4_ppb_df_deltaC = NH4_ppb_df-NH4_ppb_df_pre_light_on_mean   #making dataframe change in carbon relative to lightson (deltaC)


#***************************************
#all the same steps but for PTR
PTR_ppb_df['t_start'] = pd.to_datetime(PTR_ppb_df['t_start'])
PTR_ppb_df['t_start'] = PTR_ppb_df['t_start'] - pd.to_timedelta(6, unit='h')#convert to local time
PTR_ppb_df_pre_light_on = PTR_ppb_df.loc[PTR_ppb_df['t_start'] >=lighton_dt_minus_5]
PTR_ppb_df_pre_light_on = PTR_ppb_df_pre_light_on.loc[PTR_ppb_df_pre_light_on['t_start'] <=lighton_dt]
PTR_ppb_df = PTR_ppb_df.loc[PTR_ppb_df['t_start'] <=lightoff]#removes data after lights off
PTR_ppb_df = PTR_ppb_df.loc[PTR_ppb_df['t_start'] >=lighton]
PTR_ppb_df.rename(columns={'t_start': 'time'}, inplace=True)
PTR_ppb_df.set_index('time', inplace=True)

PTR_ppb_df = PTR_ppb_df.resample('T').mean()
PTR_ppb_df_pre_light_on = PTR_ppb_df_pre_light_on.drop('t_start', axis=1)
PTR_ppb_df_pre_light_on_mean = PTR_ppb_df_pre_light_on.mean()
PTR_ppb_df_deltaC = PTR_ppb_df-PTR_ppb_df_pre_light_on_mean 

#vapor wall loss estimations to correct dataframes
#this bit calculates molecular weights based off formulas in the vapor pressure csv
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

#this function takes the vapor pressures and molecular weights to estimate loss to wall based off of volatility
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



#find how much is lost to walls (_wlc column) and find elapsed time
loss_nh4 = VWL_corrector_2(NH4_ppb_df_deltaC, vapor_pressures)
time_e_nh4 = elapsed_minutes(loss_nh4.index)
loss_ptr = PTR_ppb_df_deltaC #vapor pressures are high enough for all ptr species they shouldnt need correcting
time_e_ptr = elapsed_minutes(loss_ptr.index)

#this takes the _wlc column and then adds whats in that column to the existing concentration to make a corrected column (_total)
for col in loss_nh4.columns:
    if col.endswith('_wlc'):
        base_col = col.replace('_wlc', '')
        if base_col in loss_nh4.columns:
            # Create the new total column
            total_col = base_col + '_total'
            loss_nh4[total_col] = loss_nh4[base_col] + loss_nh4[col]
loss_nh4 = loss_nh4[[c for c in loss_nh4.columns if c.endswith('_total') or c == 'time']]


#*****************************************************************************************
#make a copy so I know which corrections have been applied
#at this point, it should be VWL corrected change in concentration
#also just making sure the timestamps are on the minute so we can merge
corrected_nh4 = loss_nh4.copy()
corrected_nh4 = corrected_nh4.resample('min').mean()
corrected_ptr = loss_ptr.copy()
corrected_ptr = corrected_ptr.resample('min').mean()

#this merged the NH4 and PTR corrected dfs
merged_df = pd.merge(corrected_nh4, corrected_ptr, left_index=True, right_index=True, how='outer')#combines NH4 ppb data and PTR ppb dataso stacked plots can be made
merged_df = merged_df.reset_index()
merged_df = merged_df.dropna(axis=1, how='all')#drops any columns that were all nan (these are species that we did not quantify)
data_cols = [c for c in merged_df.columns if c != "time"] # set only the data columns to NaN if that row has *any* NaN
merged_df.loc[merged_df[data_cols].isna().any(axis=1), data_cols] = np.nan #drop the nans

#this reads in the direct calibrated BnOH, cuts off data outside of lights on, then averages on 1min to line up with other quantified data
bnoh_ppb = bnoh_ppb.loc[:, ['time (MDT)', 'C7H8O_ppb']] #read in concentration and time columns
bnoh_ppb['time (MDT)'] = pd.to_datetime(bnoh_ppb['time (MDT)'])

bnoh_ppb = bnoh_ppb.loc[bnoh_ppb['time (MDT)'] <=lightoff_dt]#removes data after lights off
bnoh_ppb = bnoh_ppb.loc[bnoh_ppb['time (MDT)'] >=lighton_dt]

bnoh_ppb.rename(columns={'time (MDT)': 'time'}, inplace=True) #find 1 min average
bnoh_ppb_1min, bnoh_ppb_2min, bnoh_ppb_5min = rolling_avg_2(bnoh_ppb)
bnoh_ppb_1min.reset_index(inplace=True)
bnoh_ppb_1min['e_time'] = elapsed_minutes(bnoh_ppb_1min['time'])

#now this merged dataframe holds all gas phase species and we can calculate the OH exposure for it
merged_df_corr = pd.merge(merged_df, bnoh_ppb_1min[['time', 'C7H8O_ppb']], on='time', how='outer')
e_time = elapsed_minutes(merged_df_corr['time'])#finds elapsed time in order to get at OH exposure
OH_exp = np.array(e_time, dtype=float)*OH_exp_slope#calculates OH exposure 


#%%
#********************************************************
#want to export the VWL corrected dataframe to share with Amel
#this function keeps corrected signals and removes duplicates (ie BnAld not corrected by matt's method)
#********************************************************
VWLdf = merged_df_corr.copy()

def is_valid_species(col):
    # formulas with format "mass;(...formula...)+_ppb_total"
    NH4_pattern = r"^\d+\.\d+;.*_ppb_total$"
    # formulas with format "mXX_formula_ppb"
    PTR_pattern = r"^m\d+_.+_ppb$"
    special_pattern = r"^m\d+_.+_ppb_matt$"
    bnoh_pattern = "C7H8O_ppb"
    
    return re.match(NH4_pattern, col) or re.match(PTR_pattern, col) or re.match(special_pattern, col) or re.match(bnoh_pattern, col)

exclude_species = ["124.075690;(NH4)C7H6O+_ppb_total",  "m107_C7H6OH_ppb", 'm91_C7H6H_ppb'] #taking these out because theyre duplicates of the BnAld signal


species_cols = [c for c in merged_df_corr.columns if is_valid_species(c)]
species_cols = ["time"] + species_cols
VWLdf_filtered_all = merged_df_corr[species_cols].copy()

drop_cols = ["124.075690;(NH4)C7H6O+_ppb_total", "126.091340;(NH4)C7H8O+_ppb_total", 
             "m107_C7H6OH_ppb", "m33_CH4OH_ppb", "m42_C2H3NH_ppb", "m123_C7H6O2H_ppb"]

VWLdf_filtered = VWLdf_filtered_all.drop(columns=drop_cols, errors="ignore").copy()

#this is so all the names end with the same format
def normalize_species_name(col):
    col = re.sub(r"_ppb_total$", "_ppb_VWL", col)
    col = re.sub(r"_ppb_matt$", "_ppb_VWL", col)
    col = re.sub(r"_ppb$", "_ppb_VWL", col)
    return col

VWLdf_filtered.rename(columns=normalize_species_name, inplace=True)
VWLdf_filtered.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/quantified_signals_4Amel/20220425_BnOH_quantified_VWLcorr.csv", index=False)


#%%
#************************************************************
#okay figure 1 (C# stacked plot and mass defect)
#************************************************************

#adding soa info to get pp carbon
def c_num_finder(df, list_):#function that finds the number of carbon based on formula 
    for species in df.columns:
        match = re.search(r'C(\d|[A-Z]|[a-z]?)', species)
        if match:
            next_char = match.group(1)  
                
            if next_char.isdigit():  
                list_.append(int(next_char))
            elif next_char.isupper():  
                list_.append(1)
            else:  
                list_.append(0)
        else:
            list_.append(0)  
    return list_

#functions for calculating photochemical age
def OHexp_to_age(OH_exp):
    return OH_exp/1.5e6/3600

def age_to_OHexp(age):
    return age * 1.5e6 * 3600


'''' this is for scents2 scripts
ams_df = ams_df.loc[~ams_df.isin([-9999]).any(axis=1)].copy() #removes rows that have the LDO flag 
ams_df['datetime'] = pd.to_datetime(ams_df[' AMS_MidPointTime Datetime'])
ams_df.set_index('datetime', inplace=True)
ams_df = ams_df.resample('T').mean()
ams_df = ams_df.reset_index()
'''

SOA_df['Time'] = pd.to_datetime(SOA_df['Time'])

#merged_df_aerosol = pd.merge(SOA_df, ams_df, left_on='Time', right_on='datetime', how='inner')#merged the SMSP (SOA) df with the AMS df so you can find the amount of C
merged_df_aerosol = SOA_df
merged_df_aerosol['organic_carbon ug_m3'] = (merged_df_aerosol['SOA_AVGpwl']/2.6) #finds amount of organic carbon with OM:OC ratio and mass of organic matter (OM/OC is 2.6 based off of scents 2 data)

moles_air_L3 = (0.83*1000)/(0.08206*298)
merged_df_aerosol['OC (ppb)'] = merged_df_aerosol['organic_carbon ug_m3']/(1.0E6*12)/moles_air_L3*1.0E9 #converting mass of organic carbon to ppb
to_merge_aerosol = merged_df_aerosol[['OC (ppb)','Time']]
to_merge_aerosol = to_merge_aerosol.set_index('Time')

#******************************************************************************
#i dont know what this does, but if i take it out, the whole thing breaks
#******************************************************************************

merged_df_corr = merged_df_corr.reset_index()
e_time = elapsed_minutes(merged_df_corr['time'])#finds elapsed time in order to get at OH exposure
merged_df_corr = merged_df_corr.set_index('time')  # sets time as index
OH_exp = np.array(e_time, dtype=float)*OH_exp_slope#calculates OH exposure 

merged_df_corr = merged_df_corr.dropna(axis=1, how='all')#drops any columns that were all nan (these are species that we did not quantify)
merged_df_corr.loc[merged_df_corr.isna().any(axis=1)] = np.nan
merged_index_corr = merged_df_corr.index
merged_df_corr = merged_df_corr.reset_index()
merged_df_corr = merged_df_corr.drop('time', axis=1)#gets rid of time axis


#%%
c_num_list = []
c_num_list = c_num_finder(merged_df_corr, c_num_list)#uses function to find the number of carbon for each 
c_num_list = [float(item) for item in c_num_list]
ppb_C_df = merged_df_corr.multiply(c_num_list, axis=1)#multiplying columns for specific species based on the carbon number to get it in terms of ppbC

ppb_C_df_for_aerosol = ppb_C_df.copy() #making a separate df for the aerosol data since it has additional zero periods
ppb_C_df_for_aerosol.index = merged_index_corr
ppb_C_df_for_aerosol['C7H8O_ppb_C'] = ppb_C_df_for_aerosol['C7H8O_ppb'] #adding calibration bnoh signal to aerosol df 

ppb_C_df_with_aerosol = pd.merge(ppb_C_df_for_aerosol, to_merge_aerosol, left_index=True, right_index=True, how='inner') #merging aerosol data with gas phase data

# finding the OH exposure from the aerosol + gas phase data separately since it has additional zeros 
df_for_aerosol_OH = pd.DataFrame()
df_for_aerosol_OH['Time'] = ppb_C_df_with_aerosol.index 
OH_exp_for_aerosol = np.array(elapsed_minutes(df_for_aerosol_OH['Time']), dtype=float) * OH_exp_slope 

def get_carbon_number(col_name):
    # Look for "C#" in the formula part
    match = re.search(r"C(\d+)", col_name)
    if match:
        return int(match.group(1))
    return None  # if no C number found

# Make a copy to avoid modifying original
df = merged_df_corr.copy()


def is_valid_species(col):
    # formulas with format "mass;(...formula...)+_ppb_total"
    NH4_pattern = r"^\d+\.\d+;.*_ppb_total$"
    # formulas with format "mXX_formula_ppb"
    PTR_pattern = r"^m\d+_.+_ppb$"
    special_pattern = r"^m\d+_.+_ppb_matt$"
    bnoh_pattern = "C7H8O_ppb"
    
    return re.match(NH4_pattern, col) or re.match(PTR_pattern, col) or re.match(special_pattern, col) or re.match(bnoh_pattern, col)

exclude_species = ["124.075690;(NH4)C7H6O+_ppb_total", "126.091340;(NH4)C7H8O+_ppb_total",  "m107_C7H6OH_ppb", 'm91_C7H6H_ppb', '140.070605;(NH4)C7H6O2+_ppb_total', 'm47_C2H6OH_ppb'] #taking these out because theyre duplicates or bad



species_cols = [c for c in df.columns if c != "time"]  #drop any non-time columns from grouping
species_cols = [c for c in merged_df_corr.columns if is_valid_species(c)]
species_cols = [c for c in species_cols if c not in exclude_species]
df_filtered = merged_df_corr[species_cols].copy()


carbon_groups = {}
primary = {"C7H8O_ppb": "Benzyl alcohol"}
for col in species_cols:
    cnum = get_carbon_number(col)
    if cnum is not None:
        if cnum not in carbon_groups:
            carbon_groups[cnum] = []
        carbon_groups[cnum].append(col)

# Sum within each C-number group
carbon_grouped_df = pd.DataFrame(index=df_filtered.index)
# Add special species separately
for species, label in primary.items():
    if species in df_filtered.columns:
        carbon_grouped_df[label] = df_filtered[species]

#add carbon number groups (excluding specials)
for cnum, cols in sorted(carbon_groups.items()):
    # drop species already handled as special
    cols = [c for c in cols if c not in primary]
    if cols:  # only if something is left
        carbon_grouped_df[f"C{cnum}"] = df_filtered[cols].sum(axis=1)
    
#make new df for plotting fig 1a
carbon_grouped_df = pd.DataFrame(index=df_filtered.index)

#add BnOH to grouped df
if "C7H8O_ppb" in df_filtered.columns:
    carbon_grouped_df["Benzyl alcohol"] = df_filtered["C7H8O_ppb"]

#add grouped by C# to df
for cnum, cols in sorted(carbon_groups.items()):
    cols = [c for c in cols if c != "C7H8O_ppb"]
    if cols:
        carbon_grouped_df[f"C{cnum}"] = df_filtered[cols].sum(axis=1)

    
#build a dictionary of whats in which group so later you can propagate errors correctly

group_to_species = {}

# Handle primary species first
for species, label in primary.items():
    if species in df_filtered.columns:
        group_to_species[label] = [species]

# Handle generic carbon-number groups
for cnum, cols in sorted(carbon_groups.items()):
    # drop already-handled specials
    cols = [c for c in cols if c not in primary]
    if cols:
        label = f"C{cnum}"
        group_to_species[label] = cols


#merged the GP with the PP carbon
carbon_grouped_with_SOA = pd.merge(carbon_grouped_df, to_merge_aerosol,
    left_index=True, right_index=True, how="inner")

#renames the OC (ppb) to SOA for the plot label
carbon_grouped_with_SOA = carbon_grouped_with_SOA.rename(columns={"OC (ppb)": "SOA"})

#grab labels from the column names and drop nans
labels = carbon_grouped_with_SOA.columns
carbon_grouped_df = carbon_grouped_df.dropna(axis=1, how="all")
carbon_grouped_df = carbon_grouped_df.loc[:, (carbon_grouped_df != 0).any(axis=0)]

#This section takes the carbon grouped df and multiplies each section by carbon # to get ppbC

#BnOH multiplication
for species, label in primary.items():
    if species in df_filtered.columns:
        carbon_grouped_df[label] = df_filtered[species] * 7  

#C# group multiplication
for cnum, cols in sorted(carbon_groups.items()):
    cols = [c for c in cols if c not in primary]  # avoid double counting
    if cols:
        summed = df_filtered[cols].sum(axis=1)
        carbon_grouped_df[f"C{cnum}"] = summed * cnum  # convert to ppbC

#this makes sure the OH exposure is aligned correctly
carbon_grouped_df = carbon_grouped_df.copy()
carbon_grouped_df["OH_exp"] = OH_exp
carbon_grouped_df = carbon_grouped_df.set_index("OH_exp")

#***************************************************
#this makes sure the SOA is aligned correctly
#***************************************************
# Make a temporary df with datetime + SOA
soa_tmp = ppb_C_df_with_aerosol[["OC (ppb)"]].copy()
soa_tmp = soa_tmp.reset_index().rename(columns={"index": "Time"})

# Make a datetime column for gas-phase species
gas_tmp = carbon_grouped_df.copy()
gas_tmp = gas_tmp.reset_index()  # OH_exp is index, datetime is lost if not kept
gas_tmp["Time"] = merged_index_corr.values  # datetime index from merged_df_corr

# Merge using nearest minute timestamp
merged = pd.merge_asof(gas_tmp.sort_values("Time"), soa_tmp.sort_values("Time"),
    on="Time", direction="nearest", tolerance=pd.Timedelta("1min"))

# Put SOA back into carbon_grouped_df aligned to OH_exp
carbon_grouped_df["SOA"] = merged["OC (ppb)"].values

#this bit drops gas phase data wherever theres a GP zero, but keeps as much SOA coverage as possible
gas_phase_cols = [c for c in carbon_grouped_df.columns if c != "SOA"]
mask = (carbon_grouped_df[gas_phase_cols].isna() | (carbon_grouped_df[gas_phase_cols] == 0)).all(axis=1)
carbon_grouped_df.loc[mask & carbon_grouped_df["SOA"].notna(), "SOA"] = np.nan
carbon_grouped_df = carbon_grouped_df.iloc[:-1]




# Replace zeros with NaN, then interpolate only those NaNs
carbon_grouped_df["SOA"] = carbon_grouped_df["SOA"].replace(0, np.nan)
carbon_grouped_df["SOA_i"] = carbon_grouped_df["SOA"].interpolate(method='linear')

# Identify gas-phase columns (exclude SOA and SOA_i)
gas_phase_cols = [c for c in carbon_grouped_df.columns if c not in ["SOA", "SOA_i"]]

# Mask both SOA and SOA_i when all gas-phase values are NaN or zero
mask = (carbon_grouped_df[gas_phase_cols].isna() | (carbon_grouped_df[gas_phase_cols] == 0)).all(axis=1)
carbon_grouped_df.loc[mask, ["SOA", "SOA_i"]] = np.nan

# Reorder for clarity
cols = [c for c in carbon_grouped_df.columns if c not in ["SOA", "SOA_i"]] + ["SOA", "SOA_i"]
carbon_grouped_df = carbon_grouped_df[cols]

# Fill NaNs in gas-phase for stackplot use
stackplot_df = carbon_grouped_df.copy()
stackplot_df[gas_phase_cols] = stackplot_df[gas_phase_cols].fillna(0)

# Diagnostics
print("Interpolated SOA_i count:", carbon_grouped_df['SOA_i'].notna().sum())
print("Masked rows:", mask.sum())



# --- Diagnostics (helpful to print after running) ---
print("Total rows:", len(carbon_grouped_df))
print("SOA original: NaNs =", carbon_grouped_df["SOA"].isna().sum())
print("SOA_i (interpolated): NaNs =", carbon_grouped_df["SOA_i"].isna().sum())
print("Rows where all gas-phase are NaN/zero (mask.sum):", mask.sum())
print("Example rows (masked rows):")
print(carbon_grouped_df.loc[mask, ["SOA", "SOA_i"]].head())



stackplot_df = carbon_grouped_df.copy()
stackplot_df[gas_phase_cols] = stackplot_df[gas_phase_cols].fillna(0) # Fill gas-phase NaNs with 0 for stacking
stackplot_df["SOA_i"] = stackplot_df["SOA_i"].fillna(method='ffill').fillna(method='bfill')# Create plotting-friendly SOA_i (interpolated) to avoid NaN gaps
plot_columns = gas_phase_cols + ["SOA"] # define columns to plot: gas-phase + SOA_i only

# --- Mask SOA_i_plot where all gas-phase species are zero ---
mask = (stackplot_df[gas_phase_cols] == 0).all(axis=1)  # True where all gas-phase = 0
stackplot_df.loc[mask, "SOA_i"] = 0  # or np.nan, depending on whether you want stackplot to ignore it



#%%
#section for loading in the relative errors
#makes a dictionary of errors then propagates
rel_err_dict = rel_err_df.set_index("compound")["rel_error"].to_dict()



#add OH exposure as index to df_filtered beacuse apparently its still on datetime
ppb_C_df_with_aerosol_with_OH = ppb_C_df_with_aerosol.copy()
ppb_C_df_with_aerosol_with_OH["OH_exp"] = OH_exp_for_aerosol
ppb_C_df_with_aerosol_with_OH = ppb_C_df_with_aerosol_with_OH.set_index("OH_exp")

#propagate errors at the last OH exposure point
last_oh_exp = ppb_C_df_with_aerosol_with_OH.index[-1]  # <- use the DataFrame with OH index


group_errors = {}      # absolute errors
group_err_rel = {}     # relative errors (fractions)
group_contribs = {}    # individual species abs errors

for group, species_list in group_to_species.items():
    variances = []
    contribs = {}  
    group_total = 0.0  
    
    for sp in species_list:
        if sp in ppb_C_df_with_aerosol.columns and sp in rel_err_dict:
            conc_val = ppb_C_df_with_aerosol_with_OH.loc[last_oh_exp, sp]
            group_total += conc_val
            rel_err = rel_err_dict[sp]
            abs_err = conc_val * rel_err
            variances.append(abs_err**2)
            contribs[sp] = abs_err
    
    if variances:
        group_abs_error = np.sqrt(sum(variances))
        group_errors[group] = group_abs_error
        group_contribs[group] = contribs
        
        # relative error (fraction)
        group_rel_error = group_abs_error / group_total if group_total != 0 else np.nan
        group_err_rel[group] = group_rel_error

# Example: inspect
for group in group_contribs:
    print(f"Group: {group}")
    for sp, err in group_contribs[group].items():
        print(f"  {sp}: abs error = {err:.4f}")
    print(f"  => Total abs error: {group_errors[group]:.4f}")
    print(f"  => Relative error: {group_err_rel[group]:.4f}\n")

#manually adding in SOA error based off of what tucker said shantanu said the error was (40%)
if "SOA" in carbon_grouped_df.columns:
    last_soa = carbon_grouped_df.loc[carbon_grouped_df.index[-1], "SOA"]
    group_errors["SOA"] = last_soa * 0.4
    
#manually adding in SOA error based off of what tucker said shantanu said the error was (40%)
if "SOA_i" in carbon_grouped_df.columns:
    last_soa = carbon_grouped_df.loc[carbon_grouped_df.index[-1], "SOA_i"]
    group_errors["SOA_i"] = last_soa * 0.4
    
'''           
#this just shows what the errors are for each group
for group in carbon_grouped_df.columns:
    if group in group_errors:
        val = carbon_grouped_df.loc[last_idx, group]
        err = group_errors[group]
        print(f"{group}: {val:.3f} ± {err:.3f}")
'''
#%%


#make a color dictionary for carbon# stackedplot
''' #green-blue scheme
colors_number_stacked = {
    "Benzyl alcohol": "k",
    "C2": "#1D3557", 
    "C3": "#234153", 
    "C4": "#284C4E", 
    "C5": "#336244", 
    "C6": "#3E783B",
    "C7": "#498E31",
    "SOA": "gray" }
'''
''' #blue-orange scheme
colors_number_stacked = {
    "Benzyl alcohol": "k",
    "SOA": "gray",
    "C2": "#1458B6", 
    "C3": "#485E8F", 
    "C4": "#7C6467", 
    "C5": "#8D665B", 
    "C6": "#AD6A42",
    "C7": "#CE6719" }
'''
colors_number_stacked = {
    "Benzyl alcohol": "k",
    "SOA": "gray",
    "C2": "#2B83FF", 
    "C3": "#5678C0", 
    "C4": "#806D80", 
    "C5": "#A06E60", 
    "C6": "#C97137",
    "C7": "#F1740D" }

#color dictionary if you want to color by assigned pathway
colors_pathway = {
    "Benzyl alcohol": "#000000",
    "BnAld_VWL": "#B83030",
    "HBnOH": "#52A7F1",
    "BPR Frag": "#5817AC",
    "Phen_VWL": "#B3D90B",
    "Unassigned GP Decomp Products": "#D4D4D4", 
    "BPR Ring-retaining": "#E28113"}

#dictionary for coloring by p hase
phase_colors = {"gas phase": '#AD6A42',
                "Benzyl alcohol": "k",
                "particle phase": 'grey'}

#************************************************************
#this is for the mass defect plot (fig 1b)
#load in all species names (theyre grouped here for my sanity, but there a line that just turns the dictionary into a list)
#************************************************************

#these are the gas phase groups
groups_gp = {
    "BnAld_VWL": ['m107_C7H6OH_ppb_matt', '130.049870;(NH4)C5H4O3+_ppb_total', 'm123_C7H6O2H_ppb', '156.065520;(NH4)C7H6O3+_ppb_total'],
    "HBnOH": ["142.086255;(NH4)C7H8O2+_ppb_total", 'm123_C7H6O2H_ppb', '156.065520;(NH4)C7H6O3+_ppb_total', '158.081170;(NH4)C7H8O3+_ppb_total', 'm170_C7H7NO4H_ppb', '174.076084;(NH4)C7H8O4+_ppb_total'],
    "Phen_VWL": ['m95_C6H6OH_ppb', "128.070605;(NH4)C6H6O2+_ppb_total", "157.060769;(NH4)(NO2)C6H5O+_ppb_total"],# '173.055683;(NH4)C6H5NO2O2+_ppb'],
    "BPR Frag": ['106.049870;(NH4)C3H4O3+_ppb_total', 'm85_C4H4O2H_ppb', '118.049870;(NH4)C4H4O3+_ppb_total', '122.044784;(NH4)C3H4O4+_ppb_total', '134.044784;(NH4)C4H4O4+_ppb_total', 
                   '136.060434;(NH4)C4H6O4+_ppb_total', '168.050263;(NH4)C4H6O6+_ppb_total', '132.065520;(NH4)C5H6O3+_ppb_total',
                   '148.060434;(NH4)C5H6O4+_ppb_total', '150.076084;(NH4)C5H8O4+_ppb_total', '174.076084;(NH4)C7H8O4+_ppb_total'],
    "Unassigned GP Decomp Products": ['m45_C2H4OH_ppb', 'm61_C2H4O2H_ppb', 'm59_C3H6OH_ppb', 'm73_C3H4O2H_ppb', #'80.070605;(NH4)C2H6O2+_ppb_total', 'm71_C4H6OH_ppb',
                   '120.065520;(NH4)C4H6O3+_ppb_total', '134.081170;(NH4)C5H8O3+_ppb_total']}#,#, 'm47_C2H6OH_ppb'], # 'm33_CH4OH_ppb'],}

#particle phase groups
groups_eesi = {"HBnOH": ["m227_C7H8O7Na", 'm229_C7H10O7Na'],
               "BPR Ring-retaining": ['m211_C7H8O6Na', "m213_C7H10NaO6", 'm195_C7H8NaO5', 'm197_C7H10O5Na'],
               "BPR RO2 accretion": ['m171_C5H8NaO5', 'm175_C4H8O6', 'm185_C5H6NaO6', "m183_C6H8O5Na", "m187_C5H8O6", 
                                 'm199_C6H8O6', 'm201_C6H10O6', 'm203_C5H8O7', 'm215_C6H8O7Na', 'm217_C6H10O7',
                                 'm243_C8H12O7', 'm245_C7H10O8'],
               "BPR Frag": ['m155_C4H4NaO5', 'm157_C4H6NaO5', 'm169_C5H6NaO5', 'm171_C5H8NaO5', 'm185_C5H6NaO6', 'm187_C5H8O6Na', 'm203_C5H8NaO7']}




#this goes in its own group to highlight it on the plot
BnOH = {"Benzyl alcohol": ['C7H8O_ppb']}

#for the mass defect plot to distinguish GP and PP by marker
markers_md = {"Benzyl alcohol": "X", "gas phase": "o", "particle phase": "^"} 
sizes = {"gas phase": 150, "Benzyl alcohol": 250, "particle phase": 150}

#this just flattens the groups into a list and makes a mega list of all detected species
gp_species = [sp for species_list in groups_gp.values() for sp in species_list]
eesi_species = [sp for species_list in groups_eesi.values() for sp in species_list]
bnoh = [sp for species_list in BnOH.values() for sp in species_list]
all_species = gp_species + eesi_species + bnoh

#now we gotta calculate the molecular mass from the formula
atomic_weights = {'C': 12.0107, 'H': 1.0079, 'O': 15.9994, 'N': 14.0067}
def get_ratios_and_mass(col_header):
    """
    Parses a chemical formula from a column header, including specific groups like NO2,
    ignoring NH4 and Na adducts, and calculates O/C, H/C, and approximate mass.
    """
    counts = {'C':0, 'H':0, 'O':0, 'N':0, 'Na':0}
    mass = 0

    # Step 1: handle adducts/groups in parentheses
    adduct_matches = re.findall(r'\(([A-Za-z0-9]+)\)', col_header)
    for adduct in adduct_matches:
        adduct_upper = adduct.upper()
        if adduct_upper in ['NH4', 'NA']:
            continue  # skip ammonium and sodium
        # parse other adducts like NO2
        elems = re.findall(r'([A-Z][a-z]*)(\d*)', adduct)
        for el, num in elems:
            num = int(num) if num else 1
            counts[el] += num

    # Step 2: find the main formula starting with C
    match = re.search(r'C[A-Za-z0-9]*', col_header)
    if not match:
        return None, None, None

    formula = re.sub(r'[^A-Za-z0-9]', '', match.group(0))

    # Step 3: parse main formula
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    seen = set()
    for i, (el, num) in enumerate(matches):
        # skip H in OH groups
        if el == 'H' and i > 0 and matches[i-1][0] == 'O':
            continue
        if el in seen:
            continue
        num = int(num) if num else 1
        counts[el] += num
        seen.add(el)

    # Step 4: calculate mass (ignore Na)
    for el, num in counts.items():
        if el != 'Na':
            mass += atomic_weights.get(el,0) * num

    # Step 5: calculate ratios
    C = counts['C']
    H = counts['H']
    O = counts['O']
    if C == 0:
        return None, None, mass

    oc = round(O / C, 3)
    hc = round(H / C, 3)
    return oc, hc, mass

#this fuction is to extract the formula for the label
def extract_formula(col_header):
    # Find the first substring starting with C
    adduct_matches = re.findall(r'\(([A-Za-z0-9]+)\)', col_header)
    pieces = []

    for adduct in adduct_matches:
        adduct_upper = adduct.upper()
        if adduct_upper in ["NH4", "NA"]:
            continue  # skip these
        pieces.append(adduct)

    # Step 2: find the main formula starting with C
    match = re.search(r'C[A-Za-z0-9]*', col_header)
    if match:
        pieces.append(match.group(0))

    # If nothing matched, return the raw header
    if not pieces:
        return col_header

    # Step 3: combine into one formula string (without symbols like +/_)
    formula = "".join(pieces)
    formula = re.sub(r'[^A-Za-z0-9]', '', formula)

    # Step 4: parse into elements
    element_matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    valid_elements = ['C','H','O','N']
    clean_formula = ""

    for i, (el, num) in enumerate(element_matches):
        if el not in valid_elements:
            continue
        # Ignore trailing H from OH
        if el == 'H' and i > 0 and element_matches[i-1][0] == 'O':
            continue
        # Default count is empty string if 1
        clean_formula += f"{el}{num}"

    return clean_formula

#make a df to store the info for plotting the mass defect plot
vk_df = pd.DataFrame(columns=['col_header','group','source','O:C','H:C','Mass','OSC'])

#avg ox state = 2(O/C - H/C)
#calculate avg ox state and add GP species to the df
for group, species_list in groups_gp.items():
    for sp in species_list:
        oc, hc, mass = get_ratios_and_mass(sp)
        osc = 2*oc - hc if oc is not None and hc is not None else None
        vk_df.loc[len(vk_df)] = [sp, group, 'gas phase', oc, hc, mass, osc]

#calculate avg ox state and add PP species to the df
for group, species_list in groups_eesi.items():
    for sp in species_list:
        oc, hc, mass = get_ratios_and_mass(sp)
        osc = 2*oc - hc if oc is not None and hc is not None else None
        vk_df.loc[len(vk_df)] = [sp, group, 'particle phase', oc, hc, mass, osc]

for group, species_list in BnOH.items():
    for sp in species_list:
        oc, hc, mass = get_ratios_and_mass(sp)
        osc = 2*oc - hc if oc is not None and hc is not None else None
        vk_df.loc[len(vk_df)] = [sp, group, 'Benzyl alcohol', oc, hc, mass, osc]


#rop NaNs
vk_df_clean = vk_df.dropna(subset=['Mass','OSC'])




def get_carbon_number_from_formula(col_header):
    clean_formula = extract_formula(col_header)
    match = re.search(r'C(\d+)', clean_formula)
    if match:
        return f"C{match.group(1)}"  # matches your color dict keys
    return None

vk_df_clean["Cnum"] = vk_df_clean["col_header"].apply(get_carbon_number_from_formula)

#*********************************************************************
#setup figure 1 with 1a: C# stackedplot and 1b: mass defect plot
#*********************************************************************
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(16, 12), dpi=300, gridspec_kw={"height_ratios": [1.3, 1.35]})

#top plot: carbon# stackedplot
ordered_labels = []
if "Benzyl alcohol" in carbon_grouped_df.columns:
    ordered_labels.append("Benzyl alcohol")

ordered_labels += [
    c for c in carbon_grouped_df.columns
    if c not in ["Benzyl alcohol", "SOA", "SOA_i"]
]

if "SOA_i" in carbon_grouped_df.columns:
    ordered_labels.append("SOA_i")
elif "SOA" in carbon_grouped_df.columns:
    ordered_labels.append("SOA")


ax1.stackplot(carbon_grouped_df.index, *[carbon_grouped_df[col] for col in ordered_labels],
    labels=ordered_labels, colors=[colors_number_stacked.get(lbl.split('_')[0], 'lightgray') for lbl in ordered_labels])


#add error bars!!
last_row = carbon_grouped_df.iloc[-1]
stack_heights = np.cumsum(last_row.values)
total_height = last_row[ordered_labels].sum()
errors = np.array([group_errors.get(lbl, 0.5) for lbl in ordered_labels])  # fallback 0.5 if missing
total_error = np.sqrt(np.sum(errors**2))

#positions error bars (outside right edge)
xmax = carbon_grouped_df.index[-1]
base_outside_x = xmax
x_offset_step = 0.0065 * xmax

stack_colors = [colors_number_stacked.get(lbl, "lightgray") for lbl in labels]

for i, (height, err, color) in enumerate(zip(stack_heights, errors, stack_colors)):
    lower_err = min(err, height)
    staggered_x = base_outside_x + i * x_offset_step
    overall_x = base_outside_x + len(stack_heights) * x_offset_step + x_offset_step
    ax1.errorbar(staggered_x, height, yerr=[[lower_err], [err]], fmt='none',
                 ecolor=color, capsize=1, elinewidth=5, clip_on=False)
    
overall_x = base_outside_x + len(stack_heights) * x_offset_step + x_offset_step
ax1.errorbar(overall_x, total_height, yerr=total_error, fmt='o',
             color='#D4D4D4', ecolor='#474747', capsize=1, elinewidth=7, markersize=3.5, clip_on=False)


ax1.set_xlabel("OH Exposure (molecules cm⁻³ s)", fontsize = 18)
ax1.set_ylabel("Concentration (ppb C)", fontsize = 18)
ax1.set_xlim(0, 1.15e11)
ax1.set_ylim(0, 310)

ax1.tick_params(axis='both', labelsize=18)


#this changes the 1e11 to x10^11
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
ax1.xaxis.set_major_formatter(formatter)



secax = ax1.secondary_xaxis("top", functions=(OHexp_to_age, age_to_OHexp))
secax.set_xlabel("Photochemical Age (hr)", fontsize = 18)
secax.tick_params(axis='both', labelsize=18)

formatter_x = CustomScalarFormatter(useMathText=True)
formatter_x.set_scientific(True)   # prevent sci notation
secax.xaxis.set_major_formatter(formatter_x)
#secax.xaxis.set_major_formatter(formatter)
#secax.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))

legend_labels = []
for lbl in ordered_labels:
    if lbl in ["SOA_i", "SOA"]:
        legend_labels.append("SOA")
    else:
        legend_labels.append(convert_to_latex(lbl))
        
#ax1.legend(legend_labels, ncol=2, loc="upper right", fontsize = 18)

from matplotlib.patches import Patch

left_col = [
    "Benzyl alcohol",
    "SOA",
]

right_col = [
    "C7",
    "C6",
    "C5",
    "C4",
    "C3",
    "C2"
]

# pad shorter column
max_len = max(len(left_col), len(right_col))
left_col += [""] * (max_len - len(left_col))

# Build entries COLUMN-WISE
legend_items = left_col + right_col

handles = []
labels = []

for item in legend_items:

    if item == "":
        handles.append(Patch(facecolor='none', edgecolor='none'))
        labels.append("")

    else:
        handles.append(
            Patch(facecolor=colors_number_stacked.get(item, "lightgray"))
        )

        if item in ["SOA", "SOA_i"]:
            labels.append("SOA")
        else:
            labels.append(convert_to_latex(item))

ax1.legend(
    handles,
    labels,
    ncol=2,
    fontsize=18,
    loc="upper right",
    columnspacing=0.8,
    handletextpad=0.3,
    labelspacing=0.2,
    frameon=True
)



# --------------------------
# Bottom subplot = mass plot
# --------------------------

#dictionary of manual label offsets
manual_positions = {
    "C7H10O7": (-3, -0.22),
    "C7H8O7" : (0.5, 0.01),
    "C5H8O4": (-5, 0.13),
    "C5H4O3": (-4.0, 0.35),
    "C7H6O3": (1.4, -0.05),
    "C6H6O2": (-18, 0),
    "C5H8O3": (3, -0.47),
    "NO2C6H5O": (1.6, 0.01),
    "C2H6O2": (-2, 0.15),
    "C7H6O": (1.4, -0.01),
    "C4H6O6": (1.8, 0.01),
    "C7H7NO4": (0, -0.45),
    "C4H6O3": (0.5, 0.05),
    "C8H12O7": (-3, -0.23),
    "C7H10O8": (-5, -0.23),
    "C7H8O5": (-5, -1),
    "C6H10O6": (5, -1),
    "C7H8O6": (10, -1),
    "C6H8O5": (0.5, 0.2),
    "C2H4O2": (0.5, 0.1),
    "C3H4O2": (0.5, 0.1)
    
    }#,   # (x_offset, y_offset)
'''  "C7H8O6": (2, -0.17),
    "C5H8O7": (1.5, 0.05),
    "C5H8O6": (-3, 0.2),
    "C6H8O5": (1.5, -0.3),
    "C7H8O4": (-0.5, -0.4),
    "C5H6O4": (-7, 0.2),
    "C7H6O3": (1.7, -0.1),
    "C4H6O5": (0.8, 0.15),
    "C5H8O5": (-1.8, -0.35),
    "C7H8O": (-1.8, -0.35),
    "C7H6O2": (1.4, 0),
    "C5H8O3": (1.4, -0.35),
    "C7H8O2": (1.4, -0.35),
    "C6H6O": (-2, -0.4),
    "C4H6O3": (-2.5, -0.32),
    "C6H6O2": (1, -0.3),
    "C7H6O": (-5, 0.2),
    "C5H6O3": (-6, 0.19),
    "C7H8O3": (-2, 0.17),
    "C4H6O4": (0.8, 0.1),
}'''

np.random.seed(0)
labeled_formulas = set()


for src, sub in vk_df_clean.groupby("source"):
    for _, row in sub.iterrows():
        # Determine marker color
        if src == "particle phase":
            color = "grey"  # all particle phase points
        elif src == "Benzyl alcohol":
            color = "k"     # black
        else:
            # gas phase points: use carbon-number colors
            cnum = row["Cnum"]
            color = colors_number_stacked.get(cnum, "lightgray")

        # Scatter point
        ax2.scatter(
            row["Mass"],
            row["OSC"],
            color=color,
            marker=markers_md.get(src, "o"),  # marker by phase
            s=sizes.get(src, 300),
            edgecolors="black",
            linewidths=0.4
        )

        # Add label
        formula_label = extract_formula(row["col_header"])
        if formula_label in labeled_formulas:
            continue

        if formula_label in manual_positions:
            offset_x, offset_y = manual_positions[formula_label]
        else:
            offset_x = np.random.uniform(1.4, 1.5)
            offset_y = np.random.uniform(-0.035, 0.05)

        ax2.text(
            row["Mass"] + offset_x,
            row["OSC"] + offset_y,
            convert_to_latex(formula_label),
            fontsize=14,
            color="k"
        )
        labeled_formulas.add(formula_label)
        
        
        

# create legend elements based on your marker and size dictionaries
phase_legend_elements = [
    Line2D([0], [0], marker=markers_md[phase], color='w', label=phase,
           markerfacecolor='k' if phase == 'Benzyl alcohol' else '#F4F4F4',
           markersize=np.sqrt(sizes[phase]), markeredgecolor='black')
    for phase in markers_md
]
ax1.text(-0.07, 1.05, "(a)", transform = ax1.transAxes, fontsize = 20)#, weight = "bold")
ax2.text(-0.07, -0.23, "(b)", transform = ax1.transAxes, fontsize = 20)#, weight = "bold")

# add the legend to ax2
ax2.legend(handles=phase_legend_elements, loc='lower right', fontsize=18, ncol = 2)
ax2.set_xlabel("Exact Mass (m/z)", fontsize=18)
ax2.set_ylabel("$\overline{\mathrm{OS}}_c$", fontsize=18)
ax2.tick_params(axis='both', labelsize=18)

plt.tight_layout()
plt.show()

#carbon_grouped_df.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/MainText/Figure1/20220425_Cnum_concentration.csv")

#%%
#SI version of the figure without the mass defect bit
fig, ax1 = plt.subplots(figsize=(8, 6), dpi=300)

#top plot: carbon# stackedplot
ordered_labels = []
if "Benzyl alcohol" in carbon_grouped_df.columns:
    ordered_labels.append("Benzyl alcohol")

ordered_labels += [
    c for c in carbon_grouped_df.columns
    if c not in ["Benzyl alcohol", "SOA", "SOA_i"]
]

if "SOA_i" in carbon_grouped_df.columns:
    ordered_labels.append("SOA_i")
elif "SOA" in carbon_grouped_df.columns:
    ordered_labels.append("SOA")

ax1.stackplot(carbon_grouped_df.index, *[carbon_grouped_df[col] for col in ordered_labels],
    labels=ordered_labels, colors=[colors_number_stacked.get(lbl.split('_')[0], 'lightgray') for lbl in ordered_labels])


#add error bars!!
last_row = carbon_grouped_df.iloc[-4]
stack_heights = np.cumsum(last_row.values)
total_height = last_row[ordered_labels].sum()
errors = np.array([group_errors.get(lbl, 0.5) for lbl in ordered_labels])  # fallback 0.5 if missing
total_error = np.sqrt(np.sum(errors**2))

#positions error bars (outside right edge)
xmax = carbon_grouped_df.index[-1]
base_outside_x = 2.04e11
x_offset_step = 0.0065 * xmax

stack_colors = [colors_number_stacked.get(lbl, "lightgray") for lbl in labels]

for i, (height, err, color) in enumerate(zip(stack_heights, errors, stack_colors)):
    lower_err = min(err, height)
    staggered_x = base_outside_x + i * x_offset_step
    ax1.errorbar(staggered_x, height, yerr=[[lower_err], [err]], fmt='none',
                 ecolor=color, capsize=1, elinewidth=3, clip_on=False)
overall_x = base_outside_x + len(stack_heights) * x_offset_step + x_offset_step
ax1.errorbar(overall_x, total_height, yerr=total_error, fmt='o',
             color='#D4D4D4', ecolor='#474747', capsize=1, elinewidth=5, markersize=3.5, clip_on=False)

ax1.set_xlabel("OH Exposure (molecules cm⁻³ s)", fontsize = 16)
ax1.set_ylabel("Concentration (ppb C)", fontsize = 16)
ax1.set_xlim(0, 2e11)
ax1.set_ylim(0, 540)

secax = ax1.secondary_xaxis("top", functions=(OHexp_to_age, age_to_OHexp))
secax.set_xlabel("Photochemical Age (hr)", fontsize = 16)

legend_labels = []
for lbl in ordered_labels:
    if lbl in ["SOA_i", "SOA"]:
        legend_labels.append("SOA")
    else:
        legend_labels.append(convert_to_latex(lbl))
        
#ax1.legend(legend_labels, ncol=2, loc="upper right", fontsize = 16)
from matplotlib.ticker import MultipleLocator

ax1.xaxis.set_major_locator(MultipleLocator(2.5e10))
secax.xaxis.set_major_locator(MultipleLocator(5)) 
formatter_x = CustomScalarFormatter(useMathText=True)
formatter_x.set_scientific(True)   # prevent sci notation
secax.xaxis.set_major_formatter(formatter_x)
ax1.xaxis.set_major_formatter(formatter_x)

plt.tight_layout()
plt.show()
#%%
#this is to calculate carbon closure numbers
'''
BnOH = carbon_grouped_df['Benzyl alcohol']
BnOH0 = BnOH.iloc[0]
BnOH_rxt = BnOH0 - BnOH

prods = carbon_grouped_df['C2'] + carbon_grouped_df['C3'] + carbon_grouped_df['C4'] + carbon_grouped_df['C5'] + carbon_grouped_df['C6'] + carbon_grouped_df['C7'] + carbon_grouped_df['SOA'] 

prods_plus = 0
prods_minus = 0

for col in carbon_grouped_df:
    if col == "Benzyl alcohol":
        continue  # skip BA
    err = group_err_rel.get(col, 0.0)  # fallback to 0 if missing
    prods_plus += carbon_grouped_df[col] * (1 + err)
    prods_minus += carbon_grouped_df[col] * (1 - err)



ratio = prods/BnOH_rxt
ratio_p = prods_plus/BnOH_rxt
ratio_m = prods_minus/BnOH_rxt
#print(ratio)

n_consec = 3   # require 3 consecutive points below 1
tol = 1.0

valid = ratio.replace([np.inf, -np.inf], np.nan)
valid_p = ratio_p.replace([np.inf, -np.inf], np.nan)
valid_m = ratio_m.replace([np.inf, -np.inf], np.nan)
valid = valid.dropna()
valid_p = valid_p.dropna()
valid_m = valid_m.dropna()

def find_failure(valid, label, tol=1.0, n_consec=3):
    if len(valid) == 0:
        print(f"No valid ratio values for {label}.")
        return None

    below = (valid < tol).astype(int)
    consec = below.rolling(window=n_consec, min_periods=n_consec).sum()
    hits = np.where(consec >= n_consec)[0]

    if hits.size > 0:
        first_idx = hits[0]
        exposure_at_failure = valid.index[first_idx]
        age_at_failure = OHexp_to_age(exposure_at_failure)

        print(
            f"{label}: First photochemical age with "
            f"{n_consec} consecutive ratio < {tol}: "
            f"{age_at_failure:.2f} h"
        )

        return exposure_at_failure

    else:
        print(f"{label}: No run of {n_consec} consecutive points with ratio < {tol}.")
        return None

# Apply to all three
find_failure(valid, "central", tol=tol, n_consec=n_consec)
find_failure(valid_p, "plus error", tol=tol, n_consec=n_consec)
find_failure(valid_m, "minus error", tol=tol, n_consec=n_consec)



last_idx = ratio.index[-1]

ratio_central_last = ratio.iloc[-1]
ratio_plus_last = ratio_p.iloc[-1]
ratio_minus_last = ratio_m.iloc[-1]

print(f"Central ratio at final OH: {ratio_central_last:.4f}")
print(f"Plus-error ratio at final OH: {ratio_plus_last:.4f}")
print(f"Minus-error ratio at final OH: {ratio_minus_last:.4f}")
'''
#%%
BnOH = carbon_grouped_df['Benzyl alcohol']
BnOH0 = BnOH.iloc[0]
BnOH_rxt = BnOH0 - BnOH

product_cols = ['C2','C3','C4','C5','C6','C7','SOA']
prods = carbon_grouped_df[product_cols].sum(axis=1)
prods_plus = 0
prods_minus = 0

for col in product_cols:
    err = group_err_rel.get(col, 0.0)
    prods_plus  += carbon_grouped_df[col] * (1 + err)
    prods_minus += carbon_grouped_df[col] * (1 - err)
    
ratio = prods/BnOH_rxt
ratio_p = prods_plus/BnOH_rxt
ratio_m = prods_minus/BnOH_rxt

n_consec = 2   # require 3 consecutive points below 1
tol = 1.0

valid = ratio.replace([np.inf, -np.inf], np.nan)
valid_p = ratio_p.replace([np.inf, -np.inf], np.nan)
valid_m = ratio_m.replace([np.inf, -np.inf], np.nan)
valid = valid.dropna()
valid_p = valid_p.dropna()
valid_m = valid_m.dropna()

def find_failure(valid, label, tol=1.0, n_consec=3):
    if len(valid) == 0:
        print(f"No valid ratio values for {label}.")
        return None, None

    below = (valid < tol).astype(int)
    consec = below.rolling(window=n_consec, min_periods=n_consec).sum()
    hits = np.where(consec >= n_consec)[0]

    if hits.size > 0:
        first_idx = hits[0]

        exposure_at_failure = valid.index[first_idx]
        age_at_failure = OHexp_to_age(exposure_at_failure)
        ratio_at_failure = valid.iloc[first_idx]

        print(
            f"{label}: Failure at {age_at_failure:.2f} h "
            f"(ratio = {ratio_at_failure:.3f})"
        )

        return exposure_at_failure, ratio_at_failure

    else:
        print(f"{label}: No run of {n_consec} consecutive points with ratio < {tol}.")
        return None, None
# Apply to all three
fail_c_exp, fail_c_ratio = find_failure(valid,   "central",     tol=tol, n_consec=n_consec)
fail_p_exp, fail_p_ratio = find_failure(valid_p, "plus error",  tol=tol, n_consec=n_consec)
fail_m_exp, fail_m_ratio = find_failure(valid_m, "minus error", tol=tol, n_consec=n_consec)

#this bit finds the error in the failure point
fail_exp, fail_ratio = find_failure(valid, "central", tol=tol, n_consec=n_consec)
fail_age = OHexp_to_age(fail_exp)
ratio_c = ratio.loc[fail_exp]
ratio_p_at_c = ratio_p.loc[fail_exp]
ratio_m_at_c = ratio_m.loc[fail_exp]
err_plus  = ratio_p_at_c - ratio_c
err_minus = ratio_c - ratio_m_at_c
print("\n"
    f"At photochemical age {fail_age:.2f} h:\n"
    f"Ratio = {ratio_c:.3f} "
    f"+{err_plus:.3f} / -{err_minus:.3f}"
"\n")


last_idx = ratio.index[-1]

ratio_central_last = ratio.iloc[-1]
ratio_plus_last = ratio_p.iloc[-1]
ratio_minus_last = ratio_m.iloc[-1]

print(f"Central ratio at final OH: {ratio_central_last:.4f}")
print(f"Plus-error ratio at final OH: {ratio_plus_last:.4f}")
print(f"Minus-error ratio at final OH: {ratio_minus_last:.4f}")




target_age = 5  # hours
target_exp = age_to_OHexp(target_age)

nearest_idx = np.abs(ratio.index.values - target_exp).argmin()
nearest_exp = ratio.index[nearest_idx]

ratio_c_20 = ratio.iloc[nearest_idx]
ratio_p_20 = ratio_p.iloc[nearest_idx]
ratio_m_20 = ratio_m.iloc[nearest_idx]

err_plus_20  = ratio_p_20 - ratio_c_20
err_minus_20 = ratio_c_20 - ratio_m_20

age_at_point = OHexp_to_age(nearest_exp)

print("\n""\n"
      "BEHOLD: CARBON CLOSURE \n" 
    f"At photochemical age {age_at_point:.2f} h:\n"
    f"Ratio = {ratio_c_20:.3f} "
    f"+{err_plus_20:.3f} / -{err_minus_20:.3f}"
)

target_age = 20.0  # hours
target_exp = age_to_OHexp(target_age)

nearest_idx = np.abs(ratio.index.values - target_exp).argmin()
nearest_exp = ratio.index[nearest_idx]

ratio_c_20 = ratio.iloc[nearest_idx]
ratio_p_20 = ratio_p.iloc[nearest_idx]
ratio_m_20 = ratio_m.iloc[nearest_idx]

err_plus_20  = ratio_p_20 - ratio_c_20
err_minus_20 = ratio_c_20 - ratio_m_20

age_at_point = OHexp_to_age(nearest_exp)

print(
    f"At photochemical age {age_at_point:.2f} h:\n"
    f"Ratio = {ratio_c_20:.3f} "
    f"+{err_plus_20:.3f} / -{err_minus_20:.3f}"
"\n"
"\n")



peak_idx = prods.idxmax()
age_at_peak = OHexp_to_age(peak_idx)

# Amount of C7 at that peak
C7_peak = carbon_grouped_df.loc[peak_idx, 'C7']

# BnOH reacted at that peak
BnOH_rxt_peak = BnOH_rxt.loc[peak_idx]

# Relative error for C7
C7_err = group_err_rel.get('C7', 0.0)

# Compute ratio and error bounds
ratio_C7 = C7_peak / BnOH_rxt_peak
ratio_C7_plus = (C7_peak * (1 + C7_err)) / BnOH_rxt_peak
ratio_C7_minus = (C7_peak * (1 - C7_err)) / BnOH_rxt_peak

print("\n")
print(f"C7 peaks at photochemical age: {age_at_peak:.2f} h")
print(f"C7 yield at peak total products: {ratio_C7:.3f}")
print(f"Yield + error: {ratio_C7_plus:.3f}")
print(f"yield - error: {ratio_C7_minus:.3f}")
print(f"Peak C7: {C7_peak:.3f} +/- {C7_err: 2f} ppbC")

        
#%%
#this colors the mass defect plot by assigned pathway
colors = {
    "Benzyl alcohol": "#000000",
    "BnAld_VWL": "#B83030",
    "HBnOH": "#52A7F1",
    "BPR Frag": "#5817AC",
    "Phen_VWL": "#B3D90B",
    "Unassigned GP Decomp Products": "#D4D4D4",
    "BPR RO2 accretion": "#6A0101",
    "BPR Ring-retaining": "#E28113"}

markers = {
    "gas phase": "o",   # circles
    "particle phase": "^"}   # triangles}

manual_positions_pathway = {
    "C7H10O7": (-3, -0.22),
    "C7H8O7" : (0.5, 0.01),
    "C5H8O4": (-5, 0.13),
    "C5H4O3": (-4.0, 0.35),
    "C7H6O3": (1.4, -0.05),
    "C6H6O2": (-18, 0),
    "C5H8O3": (2, -0.47),
    "NO2C6H5O": (1.6, 0.01),
    "C2H6O2": (-2, 0.15),
    "C7H6O": (1.4, -0.01),
    "C4H6O6": (1.8, 0.01),
    "C7H7NO4": (0, -0.45),
    "C4H6O3": (0.5, 0.05),
    "C7H8O5": (0.5, -0.13),
    "C7H8O6": (0.5, -0.13),
    "C7H10O8": (-3, -0.15),
    "C8H12O7": (-3, -0.15),
    }#,   # (x_offset, y_offset)

# Plot with different markers
fig, ax = plt.subplots(figsize=(16,8), dpi=300)

labeled_formulas = set()

for (src, grp), sub in vk_df_clean.groupby(['source','group']):
    ax.scatter(sub['Mass'], sub['OSC'],
               label=f"{grp}",# ({src})",
               color=colors.get(grp, 'black'),
               marker=markers.get(src, 'x'),
               s=150)
    
    # add labels with small jitter
    for i, row in sub.iterrows():
        formula_label = extract_formula(row["col_header"])
        if formula_label in labeled_formulas:
            continue

        if formula_label in manual_positions_pathway:
            offset_x, offset_y = manual_positions_pathway[formula_label]
        else:
            offset_x = np.random.uniform(2, 2.1)
            offset_y = np.random.uniform(-0.035, 0.05)
        ax.text(
            row["Mass"] + offset_x,
            row["OSC"] + offset_y,
            convert_to_latex(formula_label),
            fontsize=13,
            color="k")
        labeled_formulas.add(formula_label)
  

ax.set_xlabel('Molecular Mass (amu)', fontsize=14)
ax.set_ylabel('Average Carbon Oxidation State', fontsize=14)
ax.legend(fontsize = 13, loc = "upper left")
plt.tight_layout()
plt.show()

#%%

# === Smoothing parameters ===
window = 5  # rolling window size

# === (1) SOA concentration ===
soa_smooth_conc = carbon_grouped_df['SOA'].rolling(window=window, center=True, min_periods=1).mean()

max_idx_conc = soa_smooth_conc.idxmax()
max_soa_conc = soa_smooth_conc[max_idx_conc]
max_oh_exp_conc = OH_exp_for_aerosol[soa_smooth_conc.index.get_loc(max_idx_conc)]
max_age_conc = OHexp_to_age(max_oh_exp_conc)

print(f"Smoothed max SOA conc: {max_soa_conc:.4f} at OH exposure: {max_idx_conc}")
print(f"  → Corresponding photochemical age: {max_age_conc:.2f} hours")


# === (2) SOA mass yield ===
yields = SOA_df['SOAAMSAVGpwl_vocobs_yield']
yields_smooth = yields.rolling(window=window, center=True, min_periods=1).mean()
OH_exp_yield = SOA_df['OHexp']

max_idx_yield = yields_smooth.idxmax()
max_soa_yield = yields_smooth[max_idx_yield]
max_oh_exp_yield = OH_exp_yield[yields_smooth.index.get_loc(max_idx_yield)]
max_age_yield = OHexp_to_age(max_oh_exp_yield)

print(f"Smoothed max SOA yield: {max_soa_yield:.4f} at OH exposure: {max_oh_exp_yield}")
print(f"  → Corresponding photochemical age: {max_age_yield:.2f} hours")

target_age = 20  # hours
target_OH = age_to_OHexp(target_age)
yield_at_20h_smooth = np.interp(target_OH, OH_exp_yield, yields_smooth)

print(f"Smoothed SOA yield at {target_age} hours: {yield_at_20h_smooth:.4f}")


# === Combined subplot figure ===
fig, axes = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

# --- Top: SOA concentration ---
ax1 = axes[0]
ax1.scatter(carbon_grouped_df.index, carbon_grouped_df['SOA'], label='Original SOA conc', alpha=0.3)
ax1.scatter(carbon_grouped_df.index, soa_smooth_conc, color='steelblue', label='Smoothed SOA conc')
ax1.scatter(max_idx_conc, max_soa_conc, color='red', zorder=5, label='Max SOA conc')
ax1.set_ylabel("SOA Concentration (ppb C)")
ax1.legend()


# Add secondary x-axis for photochemical age
ax1_top = ax1.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
ax1_top.set_xlabel("Photochemical Age (hours)")

# --- Bottom: SOA mass yield ---
ax2 = axes[1]
ax2.scatter(OH_exp_yield, yields, color='m', label='SOA yield', alpha=0.3)
ax2.scatter(OH_exp_yield, yields_smooth, color='m', label='Smoothed SOA yield')
ax2.scatter(max_oh_exp_yield, max_soa_yield, color='c', zorder=5, label='Max SOA yield')
ax2.set_xlabel("OH Exposure (molec cm$^{-3}$ s)")
ax2.set_ylabel("SOA Mass Yield")
ax2.set_ylim(0, 0.5)
ax2.legend()

ax2_top = ax2.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
#ax2_top.set_xlabel("Photochemical Age (hours)")

plt.tight_layout()
plt.show()