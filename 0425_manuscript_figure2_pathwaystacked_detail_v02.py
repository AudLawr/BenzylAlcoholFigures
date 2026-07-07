# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 15:43:40 2025

@author: audlaw
"""

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

mpl.rcParams.update({'font.size': 16})

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


#this makes the origins 0 instead of 0.0
sf = ScalarFormatter(useMathText=True)
def format_ticks(value, pos):
    if value == 0:
        return "0"
    return sf.format_data(value)

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
NH4_ppb_df = NH4_ppb_df.resample('min').mean()   #resamples to get the same time index for NH4 and ptr data
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

PTR_ppb_df = PTR_ppb_df.resample('min').mean()
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

#corrected_nh4.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/20220425_benzyl_alcohol_50_200_bg_sub_r3_ppb.csv")

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

merged_df_corr.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/20220425_benzyl_alcohol_50_200_bg_sub_r3_ppb.csv")

#plt.scatter(OH_exp, merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total'], label = 'HBnOH')
plt.scatter(OH_exp, merged_df_corr['m85_C4H4O2H_ppb'], label = 'C4H4O2')
plt.scatter(OH_exp, merged_df_corr['130.049870;(NH4)C5H4O3+_ppb_total'], label = 'C5H4O3')
plt.show()

#%%
#%%
#**************************************************************
#branching fraction corrections baybeeeeee
#u need to come back to this one, its not accurate
merged_df_ohexp = merged_df_corr.copy()
merged_df_ohexp['OH_exp'] = OH_exp
merged_df_ohexp.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS01/20220425_gasphase_ppb.csv")

OH_rates = {"Benzyl alcohol" : 2.8e-11,
            "HBnOH" : 5.59e-11, 
            "BnAld": 1.29e-11,
            "Phen": 2.83e-11, 
            "Frag": 3.45e-11} #come back to this one

bnoh0 = merged_df_corr['C7H8O_ppb'].iloc[0]
bnoh_ratio = merged_df_corr['C7H8O_ppb'] / bnoh0
hbnoh_ratio = merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']/merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total'].iloc[0]
    
HBnOH_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["HBnOH"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["HBnOH"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
HBnOH_CFh = ((OH_rates["Benzyl alcohol"] - OH_rates["HBnOH"])/OH_rates["Benzyl alcohol"]) * ((1-hbnoh_ratio)/(hbnoh_ratio**(OH_rates["HBnOH"]/OH_rates["Benzyl alcohol"]) - hbnoh_ratio))
BnAld_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["BnAld"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["BnAld"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
Phen_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["Phen"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["Phen"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))

HBnOH_OHcorr = merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']*HBnOH_CF
BnAld_OHcorr = merged_df_corr['m107_C7H6OH_ppb_matt']*BnAld_CF
Phen_OHcorr = merged_df_corr['m95_C6H6OH_ppb']*Phen_CF
#THIS IS LOOKING AT CORRECTING THE BRANCHING FRACTIONS, NEED TO COME BACK TO IT

fig, ax = plt.subplots(figsize = (10, 6), dpi = 300)
ax.scatter(OH_exp, merged_df_corr['m107_C7H6OH_ppb_matt'], label = "bnald", color = "#D95A30")
ax.scatter(OH_exp, BnAld_OHcorr, label = "bnald corr for OH ox", color = "#E28413")
ax.scatter(OH_exp, merged_df_corr['m95_C6H6OH_ppb'], label = "phenol", color = "#377216")
ax.scatter(OH_exp, merged_df_corr['m95_C6H6OH_ppb']*Phen_CF, label = "phenol corr for OHox", color = "#6CAC34")
#ax.scatter(OH_exp, merged_df_corr['m61_C2H4O2H_ppb'], label = "c2", color = "grey")
ax.scatter(OH_exp, merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']*HBnOH_CF, label = 'hbnoh corrected for OH oxidation', color = "#4589D7")
ax.scatter(OH_exp, merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total'], label = 'hbnoh', color = "#2835C7")

plt.xlim(0)
plt.ylim(0, 20)
plt.xlabel("OH Exposure (molecules cm⁻³ s)")
plt.ylabel("Concentration (ppb)")
plt.legend(loc='upper right', fontsize = 12)
plt.title(titles)
plt.show()

fig,ax = plt.subplots(figsize = (10, 6), dpi = 300)
#ax.scatter(OH_exp, merged_df_corr['174.076084;(NH4)C7H8O4+_ppb_total'], label = 'epox', color = 'plum')
#ax.scatter(OH_exp, merged_df_corr['m73_C3H4O2H_ppb'], label = 'frag', color = 'c')
ax.scatter(OH_exp, merged_df_corr['140.070605;(NH4)C7H6O2+_ppb_total'], label = 'C7H6O2 nh4', color = 'steelblue')
ax.scatter(OH_exp, merged_df_corr['m123_C7H6O2H_ppb'], label = 'C7H6O2 ptr', color = 'orange')
#ax.scatter(OH_exp, merged_df_corr['157.060769;(NH4)(NO2)C6H5O+_ppb_total'], label = 'C6H5NO3 nh4', color = 'steelblue')
#ax.scatter(OH_exp, merged_df_corr['m140_C6H5NO3H_ppb'], label = 'C6H5NO3 ptr', color = 'orange')

plt.legend()
plt.show


fig,ax = plt.subplots(figsize = (10, 6), dpi = 300)
#ax.scatter(OH_exp, merged_df_corr['m107_C7H6OH_ppb_matt'], label = 'C7H6O ptr', color = 'orange')
ax.scatter(OH_exp, merged_df_corr['m140_C6H5NO3H_ppb'], label = 'C6H5NO3 ptr', color = 'm')
ax.scatter(OH_exp, merged_df_corr['157.060769;(NH4)(NO2)C6H5O+_ppb_total'], label = 'C6H5O(NO2) nh4', color = 'orange')
plt.legend()
plt.show

hbnoh_Y = HBnOH_OHcorr.nlargest(10).mean()/bnoh0
bnald_Y = BnAld_OHcorr.nlargest(10).mean()/bnoh0
phen_Y = Phen_OHcorr.nlargest(10).mean()/bnoh0

bnoh_rxt = bnoh0 - merged_df_corr['C7H8O_ppb']
bnald_y2 = BnAld_OHcorr/bnoh_rxt*100
hbnoh_y2 = HBnOH_OHcorr/bnoh_rxt*100
phen_y2 = Phen_OHcorr/bnoh_rxt*100


print(f"HBnOH Yield: {hbnoh_Y * 100:.1f}%")
print(f"BnAld Yield: {bnald_Y * 100:.1f}%")
print(f"Phenol Yield: {phen_Y * 100:.1f}%")

firstgen_frag = {"BPR Frag": ['80.070605;(NH4)C2H6O2+_ppb_total', '118.049870;(NH4)C4H4O3+_ppb_total', '120.065520;(NH4)C4H6O3+_ppb_total', '132.065520;(NH4)C5H6O3+_ppb_total', '134.044784;(NH4)C4H4O4+_ppb_total', '134.081170;(NH4)C5H8O3+_ppb_total', '148.060434;(NH4)C5H6O4+_ppb_total', '168.050263;(NH4)C4H6O6+_ppb_total', '150.076084;(NH4)C5H8O4+_ppb_total', '174.076084;(NH4)C7H8O4+_ppb_total']}
grouped_frag = {group: merged_df_corr[species].sum(axis=1) for group, species in firstgen_frag.items()}

# Convert to DataFrame for plotting
grouped_frag = pd.DataFrame(grouped_frag)
grouped_frag['OH_exp'] = OH_exp  # Add elapsed time as a column
grouped_frag.set_index('OH_exp', inplace=True)  # Set elapsed time as index

#calcs branching fraction for frag, need to revisit this
frag_y = grouped_frag['BPR Frag'].nlargest(5).mean()/bnoh0
print(f"Frag Yield: {frag_y * 100:.1f}%")



#%%

# Function to extract carbon count from species name
def extract_carbon_count(species_name):
    match = re.search(r'C(\d+)', species_name)  # Look for 'C' followed by a number
    return int(match.group(1)) if match else 1  # Default to 1 if no match found

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


#%%
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
#merging the VWL corrected df
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
'''
groups_corr_SOA = {
    "Benzyl alcohol": ['C7H8O_ppb'],
    "BnAld": ['m107_C7H6OH_ppb_matt', 'm123_C7H6O2H_ppb', '156.065520;(NH4)C7H6O3+_ppb_total', '130.049870;(NH4)C5H4O3+_ppb_total'],
    "HBnOH": ["142.086255;(NH4)C7H8O2+_ppb_total", "128.070605;(NH4)C6H6O2+_ppb_total", '158.081170;(NH4)C7H8O3+_ppb_total',
              'm123_C7H6O2H_ppb', '156.065520;(NH4)C7H6O3+_ppb_total'], #'174.076084;(NH4)C7H8O4+_ppb_total'],
    "Nitroaromatic": ['m95_C6H6OH_ppb', "128.070605;(NH4)C6H6O2+_ppb_total", "157.060769;(NH4)(NO2)C6H5O+_ppb_total", 'm170_C7H7NO4H_ppb'], #'173.055683;(NH4)C6H5NO2O2+_ppb'],
    "BPR Decomp": ['80.070605;(NH4)C2H6O2+_ppb_total', '118.049870;(NH4)C4H4O3+_ppb_total', '120.065520;(NH4)C4H6O3+_ppb_total', 
                 '134.044784;(NH4)C4H4O4+_ppb_total', '136.060434;(NH4)C4H6O4+_ppb_total', '168.050263;(NH4)C4H6O6+_ppb_total', 
                 '132.065520;(NH4)C5H6O3+_ppb_total', '134.081170;(NH4)C5H8O3+_ppb_total', '148.060434;(NH4)C5H6O4+_ppb_total', 
                 '150.076084;(NH4)C5H8O4+_ppb_total', '174.076084;(NH4)C7H8O4+_ppb_total'],
    "Unassigned": ['m45_C2H4OH_ppb', 'm61_C2H4O2H_ppb', 'm59_C3H6OH_ppb', 'm73_C3H4O2H_ppb','m71_C4H6OH_ppb', 'm85_C4H4O2H_ppb'],#, 'm47_C2H6OH_ppb'], # 'm33_CH4OH_ppb'],
    "SOA": ['OC (ppb)']}
'''

groups_corr_SOA = { #im just rearranging the fragment groups
    "Benzyl alcohol": ['C7H8O_ppb'],
    "BnAld": ['m107_C7H6OH_ppb_matt', 'm123_C7H6O2H_ppb', '156.065520;(NH4)C7H6O3+_ppb_total', '130.049870;(NH4)C5H4O3+_ppb_total'],
    "HBnOH": ["142.086255;(NH4)C7H8O2+_ppb_total", "128.070605;(NH4)C6H6O2+_ppb_total", '158.081170;(NH4)C7H8O3+_ppb_total',
              'm123_C7H6O2H_ppb', '156.065520;(NH4)C7H6O3+_ppb_total'], #'174.076084;(NH4)C7H8O4+_ppb_total'],
    "Nitroaromatic": ['m95_C6H6OH_ppb', "128.070605;(NH4)C6H6O2+_ppb_total", "m140_C6H5NO3H_ppb", 'm170_C7H7NO4H_ppb'], #'157.060769;(NH4)(NO2)C6H5O+_ppb_total', '173.055683;(NH4)C6H5NO2O2+_ppb'],m140_C6H5NO3H_ppb
    "BPR Decomp": ['106.049870;(NH4)C3H4O3+_ppb_total', 'm85_C4H4O2H_ppb', '118.049870;(NH4)C4H4O3+_ppb_total', '122.044784;(NH4)C3H4O4+_ppb_total', '134.044784;(NH4)C4H4O4+_ppb_total', 
                   '136.060434;(NH4)C4H6O4+_ppb_total', '168.050263;(NH4)C4H6O6+_ppb_total', '132.065520;(NH4)C5H6O3+_ppb_total',
                   '148.060434;(NH4)C5H6O4+_ppb_total', '150.076084;(NH4)C5H8O4+_ppb_total', '174.076084;(NH4)C7H8O4+_ppb_total'],
    "Unassigned": ['m45_C2H4OH_ppb', 'm61_C2H4O2H_ppb', 'm59_C3H6OH_ppb', 'm73_C3H4O2H_ppb', #'80.070605;(NH4)C2H6O2+_ppb_total', 'm71_C4H6OH_ppb',
                   '120.065520;(NH4)C4H6O3+_ppb_total', '134.081170;(NH4)C5H8O3+_ppb_total'],#, 'm47_C2H6OH_ppb'], # 'm33_CH4OH_ppb'],
    "SOA": ['OC (ppb)']}


colors = {
    "Benzyl alcohol": "k",
    "BnAld": "#B83030", 
    "HBnOH": "#52A7F1", 
    "Nitroaromatic": "#B3D90B", 
    "BPR Decomp": "#5817AC", 
    "Unassigned": "#D4D4D4",
    "SOA": "#ACACAC",
    "SOA_i": "#7F7F7F"}



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

#plt.scatter(OH_exp_for_aerosol, ppb_C_df_with_aerosol['174.076084;(NH4)C7H8O4+_ppb_total'], label = 'C7H8O4')
#plt.scatter(OH_exp_for_aerosol, ppb_C_df_with_aerosol['m47_C2H6OH_ppb'], label = 'C2H6O')
#plt.legend()
#plt.show()

def normalize(df_to_normalize):#this functin normalizes the data
    normalized_df = (df_to_normalize-df_to_normalize.min())/(df_to_normalize.max()-df_to_normalize.min())
    return normalized_df

normal_ppb_C_df_with_aerosol = normalize(ppb_C_df_with_aerosol)

fig, ax = plt.subplots(figsize=(10, 6))
plt.scatter(OH_exp_for_aerosol, normal_ppb_C_df_with_aerosol['134.044784;(NH4)C4H4O4+_ppb_total'], label = 'C4H4O4')
plt.scatter(OH_exp_for_aerosol, normal_ppb_C_df_with_aerosol['136.060434;(NH4)C4H6O4+_ppb_total'], label = 'C4H6O4')
plt.scatter(OH_exp_for_aerosol, normal_ppb_C_df_with_aerosol['168.050263;(NH4)C4H6O6+_ppb_total'], label = 'C4H6O6')
plt.xlabel('OH exp molec cm$^{-3}$ s')
plt.ylabel('normalized concentration')
plt.xlim(0)
plt.ylim(0)
plt.legend()
plt.show()

#%%

# --- Build grouped_aero_df (as before) ---
grouped_aero = {group: ppb_C_df_with_aerosol[species].sum(axis=1)
                for group, species in groups_corr_SOA.items()}
grouped_aero_df = pd.DataFrame(grouped_aero)

# Attach OH exposure and ensure sorted index
grouped_aero_df['OH_exp'] = OH_exp_for_aerosol
grouped_aero_df.set_index('OH_exp', inplace=True)
grouped_aero_df = grouped_aero_df.sort_index()

# --- Align SOA to same OH index ---
SOA_OH_exp = (np.array(elapsed_minutes(pd.Series(ppb_C_df_with_aerosol.index)), dtype=float)
              * OH_exp_slope)
SOA_aligned = pd.Series(ppb_C_df_with_aerosol["OC (ppb)"].values, index=SOA_OH_exp).sort_index()
grouped_aero_df["SOA"] = SOA_aligned.reindex(grouped_aero_df.index, method="nearest")

# Replace zeros with NaN, then interpolate only those NaNs
grouped_aero_df["SOA"] = grouped_aero_df["SOA"].replace(0, np.nan)
grouped_aero_df["SOA_i"] = grouped_aero_df["SOA"].interpolate(method='linear')

# Identify gas-phase columns (exclude SOA and SOA_i)
gas_phase_cols = [c for c in grouped_aero_df.columns if c not in ["SOA", "SOA_i"]]

# Mask both SOA and SOA_i when all gas-phase values are NaN or zero
mask = (grouped_aero_df[gas_phase_cols].isna() | (grouped_aero_df[gas_phase_cols] == 0)).all(axis=1)
grouped_aero_df.loc[mask, ["SOA", "SOA_i"]] = np.nan

# Reorder for clarity
cols = [c for c in grouped_aero_df.columns if c not in ["SOA", "SOA_i"]] + ["SOA", "SOA_i"]
grouped_aero_df = grouped_aero_df[cols]

# Fill NaNs in gas-phase for stackplot use
stackplot_df = grouped_aero_df.copy()
stackplot_df[gas_phase_cols] = stackplot_df[gas_phase_cols].fillna(0)

# Diagnostics
print("Interpolated SOA_i count:", grouped_aero_df['SOA_i'].notna().sum())
print("Masked rows:", mask.sum())



# --- Diagnostics (helpful to print after running) ---
print("Total rows:", len(grouped_aero_df))
print("SOA original: NaNs =", grouped_aero_df["SOA"].isna().sum())
print("SOA_i (interpolated): NaNs =", grouped_aero_df["SOA_i"].isna().sum())
print("Rows where all gas-phase are NaN/zero (mask.sum):", mask.sum())
print("Example rows (masked rows):")
print(grouped_aero_df.loc[mask, ["SOA", "SOA_i"]].head())



stackplot_df = grouped_aero_df.copy()
stackplot_df[gas_phase_cols] = stackplot_df[gas_phase_cols].fillna(0) # Fill gas-phase NaNs with 0 for stacking
stackplot_df["SOA_i"] = stackplot_df["SOA_i"].fillna(method='ffill').fillna(method='bfill')# Create plotting-friendly SOA_i (interpolated) to avoid NaN gaps
plot_columns = gas_phase_cols + ["SOA"] # define columns to plot: gas-phase + SOA_i only

# --- Mask SOA_i_plot where all gas-phase species are zero ---
mask = (stackplot_df[gas_phase_cols] == 0).all(axis=1)  # True where all gas-phase = 0
stackplot_df.loc[mask, "SOA_i"] = 0  # or np.nan, depending on whether you want stackplot to ignore it



#%%
#finding uncorrected peaks for each section
BnOH0 = grouped_aero_df['Benzyl alcohol'].iloc[0]

# Optional: choose smoothing window size
window = 3  # e.g., 3-point rolling average
apply_smoothing = True

peak_info = []

for col in grouped_aero_df.columns:
    if col == "Benzyl alcohol":
        continue  # skip the reactant itself

    y = grouped_aero_df[col] / BnOH0  # yield relative to initial BnOH

    # Smooth if desired
    if apply_smoothing:
        y = y.rolling(window=window, center=True, min_periods=1).mean()

    # Skip empty or all-NaN data
    if y.isna().all():
        print(f"{col}: all NaN, skipping.")
        continue

    # Find OH exposure and value at peak
    max_idx = y.idxmax()
    max_val = y.loc[max_idx]

    print(f"{col}: Peak yield = {max_val:.4f} at OH exposure = {max_idx:.3e}")
    peak_info.append((col, max_idx, max_val))

# Convert to DataFrame for easy viewing or saving
peak_df = pd.DataFrame(peak_info, columns=["Species", "OH_exposure_at_peak", "Peak_yield"])

#%%



#relative error section
rel_err_dict = rel_err_df.set_index("compound")["rel_error"].to_dict()


#add OH exposure as index to df_filtered beacuse apparently its still on datetime
ppb_C_df_with_aerosol_with_OH = ppb_C_df_with_aerosol.copy()
ppb_C_df_with_aerosol_with_OH["OH_exp"] = OH_exp_for_aerosol
ppb_C_df_with_aerosol_with_OH = ppb_C_df_with_aerosol_with_OH.set_index("OH_exp")

fig, ax = plt.subplots(figsize=(10, 6))

rel_err = rel_err_dict.get('128.070605;(NH4)C6H6O2+_ppb_total', 0.2)
yerr = rel_err * ppb_C_df_with_aerosol_with_OH['128.070605;(NH4)C6H6O2+_ppb_total']
y_low = ppb_C_df_with_aerosol_with_OH['128.070605;(NH4)C6H6O2+_ppb_total'] - yerr
y_high = ppb_C_df_with_aerosol_with_OH['128.070605;(NH4)C6H6O2+_ppb_total'] + yerr

plt.scatter(ppb_C_df_with_aerosol_with_OH.index, ppb_C_df_with_aerosol_with_OH['128.070605;(NH4)C6H6O2+_ppb_total'], label = "C$_6$H$_6$O$_2$")
plt.axhline(y= 0.06120, linestyle='--', c = "r", label = "detection limit")  #hoz line at catechol DL
ax.fill_between(
    ppb_C_df_with_aerosol_with_OH.index, y_low, y_high,
    color='gray',
    alpha=0.2,
    label="uncertainty")

plt.legend()
plt.xlabel("OH Exposure (molecules cm⁻³ s)")
plt.ylabel("concentration (ppb)")
plt.show()




fig, ax = plt.subplots(figsize=(10, 6))

col = '128.070605;(NH4)C6H6O2+_ppb_total'
y_raw = ppb_C_df_with_aerosol_with_OH[col]
y = y_raw / bnoh0
x = ppb_C_df_with_aerosol_with_OH.index

rel_err = rel_err_dict.get(col, 0.2)
yerr = rel_err * y
y_low = y - yerr
y_high = y + yerr


top_n = 10


# --- top points ---
top_idx = y.nlargest(top_n).index
y_top = y.loc[top_idx]
x_top = top_idx
top_avg = y_top.mean()

yerr_top = rel_err * y_top
mean_err = np.sum(yerr_top) / len(y_top)

print(f"Top {top_n} average: {top_avg:.4f} ± {mean_err:.4f}")


# uncertainty
ax.fill_between(x, y_low, y_high,
                color='gray', alpha=0.2,
                label="uncertainty")

# main data
ax.scatter(x, y, label="C$_6$H$_6$O$_2$ / BnOH$_0$")

# highlight top points
ax.scatter(x_top, y_top,
           color='gold',
           edgecolor='k',
           s=120,
           zorder=5,
           label=f"Top {top_n}")

# average line
ax.axhline(top_avg,
           color='gold',
           linestyle='--',
           linewidth=2,
           label=f"Top {top_n} avg")

# detection limit
ax.axhline(y=0.06120 / bnoh0,
           linestyle='--', c="r",
           label="detection limit")

ax.legend()
ax.set_ylabel("Fraction of initial carbon")

plt.show()

#%%
#propagate errors at the last OH exposure point
last_oh_exp = ppb_C_df_with_aerosol_with_OH.index[-1]  # <- use the DataFrame with OH index
group_errors = {}

group_contribs = {}  # store individual contributions

for group, species_list in groups_corr_SOA.items():
    variances = []
    contribs = {}  # store each species' abs error
    for sp in species_list:
        if sp in ppb_C_df_with_aerosol.columns and sp in rel_err_dict:
            conc_val = ppb_C_df_with_aerosol_with_OH.loc[last_oh_exp, sp]
            print(conc_val)
            rel_err = rel_err_dict[sp]
            abs_err = conc_val * rel_err
            variances.append(abs_err**2)
            contribs[sp] = abs_err
    if variances:
        group_errors[group] = np.sqrt(sum(variances))
        group_contribs[group] = contribs

# Example: inspect
for group, contrib in group_contribs.items():
    print(f"Group: {group}")
    for sp, err in contrib.items():
        print(f"  {sp}: {err:.4f}")
    print(f"  => Total propagated error: {group_errors[group]:.4f}\n")


#manually adding in SOA error based off of what tucker said shantanu said the error was (40%)
if "SOA" in grouped_aero_df.columns:
    last_soa = grouped_aero_df.loc[grouped_aero_df.index[-1], "SOA"]
    group_errors["SOA"] = last_soa * 0.4
    
#manually adding in SOA error based off of what tucker said shantanu said the error was (40%)
if "SOA_i" in grouped_aero_df.columns:
    last_soa = grouped_aero_df.loc[grouped_aero_df.index[-1], "SOA_i"]
    group_errors["SOA_i"] = last_soa * 0.4
    
#error for the frac carbon stackedplots
frac_carb = ppb_C_df_with_aerosol/260

ppb_C_df_with_aerosol_with_OH_frac = frac_carb.copy()
ppb_C_df_with_aerosol_with_OH_frac["OH_exp"] = OH_exp_for_aerosol
ppb_C_df_with_aerosol_with_OH_frac = ppb_C_df_with_aerosol_with_OH_frac.set_index("OH_exp")

#ppb_C_df_with_aerosol_with_OH_frac.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/MainText/Figure2/20220425_TEST.csv")

plt.scatter(ppb_C_df_with_aerosol_with_OH_frac.index, ppb_C_df_with_aerosol_with_OH_frac['m123_C7H6O2H_ppb'], label = "C7H6O2") #140.070605;(NH4)C7H6O2+_ppb_total
plt.scatter(ppb_C_df_with_aerosol_with_OH_frac.index, ppb_C_df_with_aerosol_with_OH_frac['156.065520;(NH4)C7H6O3+_ppb_total'], label = "C7H6O3")
plt.legend()
plt.show()
#%%
#propagate errors at the last OH exposure point
last_oh_exp_frac = ppb_C_df_with_aerosol_with_OH_frac.index[-1]  # <- use the DataFrame with OH index
group_errors_frac = {}

group_contribs_frac = {}  # store individual contributions

for group, species_list in groups_corr_SOA.items():
    variances_frac = []
    contribs_frac = {}  # store each species' abs error
    for sp in species_list:
        if sp in ppb_C_df_with_aerosol.columns and sp in rel_err_dict:
            conc_val_frac = ppb_C_df_with_aerosol_with_OH_frac.loc[last_oh_exp, sp]
            print(conc_val_frac)
            rel_err_frac = rel_err_dict[sp]
            abs_err_frac = conc_val_frac * rel_err_frac
            variances_frac.append(abs_err_frac**2)
            contribs_frac[sp] = abs_err_frac
    if variances_frac:
        group_errors_frac[group] = np.sqrt(sum(variances_frac))
        group_contribs_frac[group] = contribs_frac

# Example: inspect
for group, contrib in group_contribs_frac.items():
    print(f"Group: {group}")
    for sp, err in contrib.items():
        print(f"  {sp}: {err:.4f}")
    print(f"  => Total propagated error: {group_errors_frac[group]:.4f}\n")


#manually adding in SOA error based off of what tucker said shantanu said the error was (40%)
if "SOA" in grouped_aero_df.columns:
    last_soa = grouped_aero_df.loc[grouped_aero_df.index[-1], "SOA"]
    group_errors_frac["SOA"] = last_soa * 0.4

if "SOA_i" in grouped_aero_df.columns:
    last_soa = grouped_aero_df.loc[grouped_aero_df.index[-1], "SOA_i"]
    group_errors_frac["SOA_i"] = last_soa * 0.4



'''        
#this just shows what the errors are for each group
for group in grouped_aero_df.columns:
    if group in group_errors:
        val = grouped_aero_df.loc[last_idx, group]
        err = group_errors[group]
        print(f"{group}: {val:.3f} ± {err:.3f}")
'''





#%%

# --- Plot stackplot ---
fig, ax = plt.subplots(figsize=(10, 6))
labels = [c for c in stackplot_df.columns if c not in ["SOA", "SOA_i"]] + ["SOA_i"]
plot_columns = gas_phase_cols + ["SOA_i"]
stack_colors = [colors.get(lbl, "lightgray") for lbl in labels]



ax.stackplot(stackplot_df.index, stackplot_df[plot_columns].T, labels=plot_columns, colors=stack_colors)

# --- Compute last point ±10% error ---
last_row = grouped_aero_df.iloc[-1]
stack_heights = np.cumsum(last_row.values)
errors = np.array([group_errors.get(lbl, 0.5) for lbl in labels])
# X positioning for error bars (staggered outside right edge)
xmax = grouped_aero_df.index[-1]
base_outside_x = xmax
x_offset_step = 0.015 * xmax

#stack_colors = [colors.get(lbl, "lightgray") for lbl in labels]

for i, (height, err, color) in enumerate(zip(stack_heights, errors, stack_colors)):
    lower_err = min(err, height)
    staggered_x = base_outside_x + i * x_offset_step
    ax.errorbar(staggered_x, height, yerr=[[lower_err], [err]], fmt='none',
                ecolor=color, capsize=2, elinewidth=3, clip_on=False)


# secondary x-axis: Photochemical Age
def OHexp_to_age(OH_exp):
    return OH_exp/1.5e6/3600

def age_to_OHexp(age):
    return age * 1.5e6 * 3600

secax = ax.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
secax.set_xlabel("Photochemical Age (hr)")
ax.set_xlabel("OH Exposure (molecules cm⁻³ s)")


ax.set_xlim(0, 1.15e11)
ax.set_ylabel("Concentration (ppb C)")
ax.set_ylim(0, 310)
#ax.set_title("0425 Benzyl Alcohol VWL corrected", fontsize=17)
ax.legend(fontsize=12, ncol=2)  # two-column legend
plt.tight_layout()
#plt.title(titles)
plt.show()


#%%

#uhh so this does a plot for the species in each group
group_cmaps = {
    "BnAld": cm.Reds,
    "HBnOH": cm.Blues,
    "Nitroaromatic": cm.Greens,
    "BPR Decomp": cm.Purples,
    "Unassigned": cm.Greys,
    "SOA": cm.Greys,
    "SOA_i": cm.Greys}

# groups to plot (exclude BnOH_VWL)
plot_groups = [g for g in groups_corr_SOA.keys() if g != "BnOH_VWL"]

top_row_groups = ["BnAld", "HBnOH"]
bottom_row_groups = ["Nitroaromatic", "BPR Decomp", "Unassigned"]


def clean_label(col_name):
    # Try first regex: things like (fg)(formula)+
    match = re.search(r"(?:\([^\)]+\))?(\([A-Za-z0-9]+\))?([A-Za-z0-9]+)\+", col_name)
    if match:
        fg = match.group(1) if match.group(1) else ""
        formula = match.group(2)
        # remove trailing H (and optional digit(s))
        formula = re.sub(r'H\d*$', '', formula)
        return fg + formula

    # Try second regex: things like _C7H7NO4H_ppb
    match = re.search(r"_(C[0-9A-Za-z]+)_ppb", col_name)
    if match:
        formula = match.group(1)
        formula = re.sub(r'H\d*$', '', formula)
        return formula

    # Fallback
    formula = col_name.replace("_ppb_total", "").replace("_ppb_matt", "").strip()
    formula = re.sub(r'H\d*$', '', formula)
    return formula


#%%
'''
#**********************************************************
#figure 2: stacked plots for days
#these are in ppbC with errors applied to the final point in the series
#**********************************************************


fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=False)
axes2d = np.array(axes).reshape(2, 3)

# Parameters for error bar positioning
xmax = OH_exp_for_aerosol[-1]  # last OH exposure
base_outside_x = 1.18e11
x_offset_step = 0.015 * xmax  # horizontal spacing for staggered error bars


# --- (1) First stackedplot goes in top-left (0,0) ---
ax_main = axes2d[0, 0]

labels = grouped_aero_df.columns

ax_main.stackplot(stackplot_df.index, stackplot_df[plot_columns].T, labels=plot_columns, colors=stack_colors)

# Error bars at second-to-last OH exposure point
stack_heights = stackplot_df.iloc[-1].values
errors = [group_errors.get(lbl, 0) for lbl in labels]

# Compute cumulative heights to match the stackplot top
cumulative_heights = np.cumsum(stack_heights)

for i, (height, cum_height, err, color) in enumerate(zip(stack_heights, cumulative_heights, errors, stack_colors)):
    if err == 0:
        continue
    err = abs(err)
    if cum_height >= 0:
        lower_err = min(err, abs(cum_height))
    else:
        lower_err = min(err, abs(cum_height))
    staggered_x = base_outside_x + i * x_offset_step
    # y-value for error bar = top of this stack layer
    ax_main.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                     fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

ax_main.set_title("Total", fontsize=16)
ax_main.set_xlim(0, 1.15e11)
ax_main.set_ylim(0, 310)
ax_main.legend(fontsize = 10)
#ax_main.set_ylabel("ppb C")

# --- (2) Group subplots go in remaining 5 slots ---
placement = {
    "BnAld": axes2d[0, 1],
    "HBnOH": axes2d[0, 2],
    "Nitroaromatic": axes2d[1, 0],
    "BPR Decomp": axes2d[1, 1],
    "Unassigned": axes2d[1, 2]}

for group_name, ax in placement.items():
    if group_name not in groups_corr_SOA:
        continue
    species_list = groups_corr_SOA[group_name]
    species_in_df = [s for s in species_list if s in ppb_C_df_with_aerosol.columns]
    if not species_in_df:
        continue

    cmap = group_cmaps.get(group_name, cm.viridis)
    species_colors = [cmap(i) for i in np.linspace(0.3, 0.9, len(species_in_df))]

    ax.stackplot(
        OH_exp_for_aerosol,
        [ppb_C_df_with_aerosol[s] for s in species_in_df],
        labels = [convert_to_latex(clean_label(s)) for s in species_in_df],
        colors=species_colors)
    
    # error bars for individual species
    species_heights  = ppb_C_df_with_aerosol[species_in_df].iloc[-1].values
    errors = [group_contribs[group_name].get(sp, 0) for sp in species_in_df]
   # Compute cumulative heights for the stack
    cumulative_heights = np.cumsum(species_heights)

    for i, (height, cum_height, err, color) in enumerate(zip(species_heights, cumulative_heights, errors, species_colors)):
        if err == 0:
            continue
        err = abs(err)
        if cum_height >= 0:
            lower_err = min(err, abs(cum_height))
        else:
            lower_err = min(err, abs(cum_height))
        #lower_err = min(err, height)
        staggered_x = base_outside_x + i * x_offset_step
        ax.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                    fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

    
    
    ax.set_title(group_name, fontsize=16)
    ax.set_xlim(0, 1.16e11)

    if group_name in ["BnAld", "HBnOH"]:
        ax.set_ylim(0, 60)
    else:
        ax.set_ylim(0, 21)

    if group_name == "Unassigned":
        ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(-0.001, 0.8))
    else:
        ax.legend(fontsize=9)
        
        


# global x/y labels
fig.text(0.5, 0.04, "OH Exposure (molecules cm⁻³ s)", ha="center", fontsize=16)
fig.text(0.05, 0.5, "Concentration (ppb C)", va="center", rotation="vertical", fontsize=16)

plt.xlim(0, 1.16e11)
plt.tight_layout(rect=[0.05, 0.05, 1, 1], w_pad=0, h_pad=-0.5)
plt.show()




#%%
#plots fraction of carbon with errors applied to the final point
fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=False)
axes2d = np.array(axes).reshape(2, 3)

# Parameters for error bar positioning
xmax = OH_exp_for_aerosol[-1]  # last OH exposure
base_outside_x = 1.18e11
x_offset_step = 0.015 * xmax  # horizontal spacing for staggered error bars


# --- (1) First stackedplot goes in top-left (0,0) ---
ax_main = axes2d[0, 0]

labels = grouped_aero_df.columns

ax_main.stackplot(stackplot_df.index, stackplot_df[plot_columns].T, labels=plot_columns, colors=stack_colors)

# Error bars at second-to-last OH exposure point
stack_heights = stackplot_df.iloc[-1].values
errors = [group_errors.get(lbl, 0) for lbl in labels]

# Compute cumulative heights to match the stackplot top
cumulative_heights = np.cumsum(stack_heights)

for i, (height, cum_height, err, color) in enumerate(zip(stack_heights, cumulative_heights, errors, stack_colors)):
    if err == 0:
        continue
    lower_err = min(err, height)
    staggered_x = base_outside_x + i * x_offset_step
    # y-value for error bar = top of this stack layer
    ax_main.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                     fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

ax_main.set_title("Total", fontsize=16)
ax_main.set_xlim(0, 1.15e11)
ax_main.set_ylim(0, 310)
ax_main.legend(fontsize = 10)
#ax_main.set_ylabel("ppb C")

# --- (2) Group subplots go in remaining 5 slots ---
placement = {
    "BnAld": axes2d[0, 1],
    "HBnOH": axes2d[0, 2],
    "Nitroaromatic": axes2d[1, 0],
    "BPR Decomp": axes2d[1, 1],
    "Unassigned": axes2d[1, 2]}



for group_name, ax in placement.items():
    if group_name not in groups_corr_SOA:
        continue
    species_list = groups_corr_SOA[group_name]
    species_in_df = [s for s in species_list if s in frac_carb.columns]
    if not species_in_df:
        continue

    cmap = group_cmaps.get(group_name, cm.viridis)
    species_colors = [cmap(i) for i in np.linspace(0.3, 0.9, len(species_in_df))]

    ax.stackplot(
        OH_exp_for_aerosol,
        [frac_carb[s] for s in species_in_df],
        labels = [convert_to_latex(clean_label(s)) for s in species_in_df],
        colors=species_colors)

    # error bars for individual species
    species_heights  = frac_carb[species_in_df].iloc[-1].values
    errors = [group_contribs_frac[group_name].get(sp, 0) for sp in species_in_df]
   # Compute cumulative heights for the stack
    cumulative_heights = np.cumsum(species_heights)

    for i, (height, cum_height, err, color) in enumerate(zip(species_heights, cumulative_heights, errors, species_colors)):
        if err == 0:
            continue
        err = abs(err)
        if cum_height >= 0:
            lower_err = min(err, abs(cum_height))
        else:
            lower_err = min(err, abs(cum_height))
        staggered_x = base_outside_x + i * x_offset_step
        ax.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                    fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

    
    
    ax.set_title(group_name, fontsize=16)
    ax.set_xlim(0, 1.16e11)

    if group_name in ["BnAld", "HBnOH"]:
        ax.set_ylim(0, 0.25)
    else:
        ax.set_ylim(0, 0.08)

    if group_name == "Unassigned":
        ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(-0.001, 0.785))
    else:
        ax.legend(fontsize=9)
        
        


# global x/y labels
fig.text(0.5, 0.04, "OH Exposure (molecules cm⁻³ s)", ha="center", fontsize=16)
fig.text(0.05, 0.78, "Concentration (ppb C)", va="center", rotation="vertical", fontsize=16)
fig.text(0.05, 0.3, "Fraction Carbon ([product]$_t$/[BnOH]$_0$)", va="center", rotation="vertical", fontsize=16)

plt.xlim(0, 1.16e11)
plt.tight_layout(rect=[0.05, 0.05, 1, 1], w_pad=1.5, h_pad=-0.5)
plt.show()

'''

#%%
#THIS IS THE MANUSCRIPT ONE
average_OH = 1.5e6  

# secondary x-axis: Photochemical Age
def OHexp_to_age(OH_exp):
    return OH_exp/1.5e6/3600

def age_to_OHexp(age):
    return age * 1.5e6 * 3600

def add_photochem_age_axis(ax, label=True):
    secax = ax.secondary_xaxis(
        'top',
        functions=(OHexp_to_age, age_to_OHexp)
    )
    if label:
        secax.set_xlabel(" ", fontsize=14)
    secax.tick_params(labelsize=16)
    return secax

class CustomScalarFormatter(ScalarFormatter): #this sets the xaxis origin to 0 and makes sci notation 10^11 instead of e11
    def __call__(self, x, pos=None):
        if x == 0.0:
            return "0"
        return super().__call__(x, pos)

#makes yaxis origin 0 and all values 2 decimal places
def format_y_ticks(value, pos):
    if value == 0:
        return "0"
    return f"{value:.2f}"  # two decimal places for all other ticks

# dictionary to store results
max_summary = {}

for lbl in labels:
    # find the OH exposure where this species reaches its max
    peak_idx = stackplot_df[lbl].idxmax()
    peak_val = stackplot_df.loc[peak_idx, lbl]
    
    # convert OH exposure → photochemical age (hours)
    photochem_age_hr = peak_idx / average_OH / 3600
    
    # store results
    max_summary[lbl] = {
        "OH_exposure_max": peak_idx,
        "value_max": peak_val,
        "photochemical_age_hr": photochem_age_hr
    }

# convert to DataFrame for convenience
max_summary_df = pd.DataFrame(max_summary).T
print(max_summary_df)


fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharex=False, dpi=300)
axes2d = np.array(axes).reshape(2, 3)

# Parameters for error bar positioning
xmax = OH_exp_for_aerosol[-1]  # last OH exposure
base_outside_x = 1.18e11
x_offset_step_tot = 0.01 * xmax  # horizontal spacing for staggered error bars on the total plot
x_offset_step_detail = 0.015 * xmax #horizontal spacing for the pathway detail plots
# --- (1) First stackedplot goes in top-left (0,0) ---
ax_main = axes2d[0, 0]
labels = []
for col in plot_columns:  # plot_columns = gas_phase_cols + ["SOA_i"]
    if col == "SOA_i":
        labels.append("SOA")   # show as "SOA" in the legend
    else:
        labels.append(col)

ax_main.stackplot(stackplot_df.index, stackplot_df[plot_columns].T, labels=labels, colors=stack_colors)
#stackplot_df.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/MainText/Figure2/20220425_pathway_concentration.csv")

errors = [group_errors.get(lbl, 0) for lbl in labels]

# Compute error bar positions based on each species' own max
for i, (lbl, err, color) in enumerate(zip(labels, errors, stack_colors)):
    if err == 0:
        continue

    # Find where this species hits its maximum
    species_series = stackplot_df[lbl]
    peak_idx = species_series.idxmax()
    peak_val = species_series.loc[peak_idx]

    # Compute cumulative stack height up to this species at the same OH exposure
    cumulative_heights = stackplot_df.loc[peak_idx].cumsum()
    cum_height = cumulative_heights.iloc[i]

    # Determine lower error (so it doesn’t go below the base of its own layer)
    lower_err = min(err, peak_val)

    # Position error bar slightly outside to the right
    staggered_x = base_outside_x + i * x_offset_step_tot

    # Plot error bar at correct y (cumulative) height
    ax_main.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                     fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

    #Overall stack error for total plot
    total_stack = stackplot_df[plot_columns].sum(axis=1)
    max_idx_total = total_stack.idxmax()
    max_total = total_stack.loc[max_idx_total]
    
    #combine errore
    overall_err_total = np.sqrt(np.sum(np.array(errors)**2))
    
    #positioning
    overall_x = base_outside_x + len(labels) * x_offset_step_tot * 1.2
    
    ax_main.errorbar(overall_x, max_total, yerr=[[min(overall_err_total, max_total)], [overall_err_total]],
                     fmt='none', ecolor='#474747', capsize=1, elinewidth=7, clip_on=False)




add_photochem_age_axis(ax_main, label=True)

ax_main.set_xlim(0, 1.15e11)
ax_main.set_ylim(0, 310)
ax_main.legend(fontsize=10)

#this changes the 1e11 to x10^11
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))
ax_main.xaxis.set_major_formatter(formatter)
ax_main.xaxis.get_offset_text().set_visible(True)


# --- (2) Group subplots go in remaining 5 slots ---
placement = {
    "BnAld": axes2d[0, 1],
    "HBnOH": axes2d[0, 2],
    "Nitroaromatic": axes2d[1, 0],
    "BPR Decomp": axes2d[1, 1],
    "Unassigned": axes2d[1, 2]
}

for group_name, ax in placement.items():
    if group_name not in groups_corr_SOA:
        continue
    species_list = groups_corr_SOA[group_name]
    species_in_df = [s for s in species_list if s in frac_carb.columns]
    if not species_in_df:
        continue

    cmap = group_cmaps.get(group_name, cm.viridis)
    species_colors = [cmap(i) for i in np.linspace(0.3, 0.9, len(species_in_df))]

    ax.stackplot(
        OH_exp_for_aerosol,
        [frac_carb[s] for s in species_in_df],
        labels=[convert_to_latex(clean_label(s)) for s in species_in_df],
        colors=species_colors)

    # --- Align error bars with real stack peak ---
    stack_sum = frac_carb[species_in_df].sum(axis=1)
    max_idx = stack_sum.idxmax()  # OH exposure where this group peaks
   

    # Heights at that OH exposure
    species_heights = frac_carb.loc[max_idx, species_in_df].values
    errors = [group_contribs_frac[group_name].get(sp, 0) for sp in species_in_df]

    cumulative_heights = np.cumsum(species_heights)

    for i, (height, cum_height, err, color) in enumerate(zip(species_heights, cumulative_heights, errors, species_colors)):
        if err == 0:
            continue
        err = abs(err)
        if cum_height >= 0:
            lower_err = min(err, abs(cum_height))
        else:
            lower_err = min(err, abs(cum_height))
        staggered_x = base_outside_x + i * x_offset_step_detail
        ax.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                    fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)
        
        
        # ---- Overall stack error for this group ----
        stack_sum = frac_carb[species_in_df].sum(axis=1)
        max_idx = stack_sum.idxmax()
        max_val = stack_sum.loc[max_idx]
        
        species_errors = [group_contribs_frac[group_name].get(sp, 0)
                          for sp in species_in_df]
        
        overall_err = np.sqrt(np.sum(np.array(species_errors)**2))
        
        overall_x = base_outside_x + len(species_in_df) * x_offset_step_detail * 1.2
        
        ax.errorbar(overall_x, max_val,
                    yerr=[[min(overall_err, max_val)], [overall_err]], fmt='none',
                    ecolor='#474747', capsize=1, elinewidth=7, clip_on=False)


#    ax.set_title(group_name, fontsize=16)
    ax.set_xlim(0, 1.16e11)
    
    
    #FORMATTING THE ORIGIN AND THE SCI NOTATION
    #this changes the 1e11 to x10^11
    '''formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((0, 0))
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.get_offset_text().set_visible(True)
    '''
    formatter_x = CustomScalarFormatter(useMathText=True) # makes x 0 just 0
    formatter_x.set_scientific(True)   # prevent sci notation
    ax_main.xaxis.set_major_formatter(formatter_x) #overall plot
    ax.xaxis.set_major_formatter(formatter_x) #detail plots
    ax.yaxis.set_major_formatter(FuncFormatter(format_y_ticks)) #detail plot yaxes
    
    
    add_photochem_age_axis(ax, label=False)
    
    if group_name in ["BnAld", "HBnOH"]:
        ax.set_ylim(0, 0.25)
    else:
        ax.set_ylim(0, 0.08)

    if group_name == "Unassigned":
        ax.legend(title = group_name, title_fontsize=11, fontsize=10, loc='center left', bbox_to_anchor=(-0.001, 0.755))
    elif group_name == "BnAld":
        ax.legend(title=group_name, title_fontsize=11, fontsize=10, loc='center left', bbox_to_anchor=(-0.001, 0.823))
    elif group_name == "BPR Decomp":   #this one gets 2 columns
        ax.legend(title=group_name,
                  title_fontsize=11,
                  fontsize=10,
                  ncol=2)
    else:
        ax.legend(title=group_name, title_fontsize=11, fontsize=10)


    
#labels because we have to add these i guess
fig.text(0.5, 0.055, "OH Exposure (molecules cm⁻³ s)", ha="center", fontsize=16)
fig.text(0.5, 0.97, "Photochemical age (hr)", ha="center", fontsize=16)
fig.text(0.05, 0.78, "Concentration (ppb C)", va="center", rotation="vertical", fontsize=16)
fig.text(0.05, 0.3, "Fraction Carbon", va="center", rotation="vertical", fontsize=16) # ([product]$_t$/[BnOH]$_0$)
fig.text(0.365, 0.78, "Fraction Carbon", va="center", rotation="vertical", fontsize=16)
# fig.text(0.1, 0.365, "365 PARTY GIRL", va="center", rotation="horizontal", fontsize=160, weight = "bold", color = "#69FF23") # ([product]$_t$/[BnOH]$_0$)



plt.xlim(0, 1.16e11)
plt.tight_layout(rect=[0.05, 0.05, 1, 1], w_pad=0, h_pad=0)
plt.show()

#%%



#%%
#loops through and plots them separately
# =====================================================
# (1) Main stacked plot (Total)
# =====================================================


fig, ax_main = plt.subplots(figsize=(7, 7), dpi=300)



labels = ["SOA" if col == "SOA_i" else col for col in plot_columns]
ax_main.stackplot(stackplot_df.index, stackplot_df[plot_columns].T, labels=labels, colors=stack_colors)

errors = [group_errors.get(lbl, 0) for lbl in labels]
xmax = OH_exp_for_aerosol[-1]
base_outside_x = 1.18e11
x_offset_step_tot = 0.01 * xmax

# Error bars for main plot
for i, (lbl, err, color) in enumerate(zip(labels, errors, stack_colors)):
    if err == 0:
        continue
    err = abs(err)
    if cum_height >= 0:
        lower_err = min(err, abs(cum_height))
    else:
        lower_err = min(err, abs(cum_height))
    series = stackplot_df[lbl]
    peak_idx = series.idxmax()
    peak_val = series.loc[peak_idx]
    cum_height = stackplot_df.loc[peak_idx].cumsum().iloc[i]
    #lower_err = min(err, peak_val)
    staggered_x = base_outside_x + i * x_offset_step_tot
    ax_main.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                     fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

# apply the wrapped formatter

# Keep normal formatter
class CustomScalarFormatter(ScalarFormatter):
    def __call__(self, x, pos=None):
        if x == 0.0:
            return "0"
        return super().__call__(x, pos)
 


formatter_x = CustomScalarFormatter(useMathText=True)
formatter_x.set_scientific(True)   # prevent sci notation
ax_main.xaxis.set_major_formatter(formatter_x)



ax_main.set_xlim(0, 1.15e11)
ax_main.set_ylim(0, 310)
ax_main.legend(fontsize=10)
ax_main.set_xlabel("OH Exposure (molecules cm⁻³ s)", fontsize=14)
ax_main.set_ylabel("Concentration (ppb C)", fontsize=14)
ax_main.set_title("Total Carbon Distribution", fontsize=16)
plt.tight_layout()
plt.show()

#%%

plt.scatter(OH_exp_for_aerosol, frac_carb['128.070605;(NH4)C6H6O2+_ppb_total'])

#%%
# =====================================================
# (2) Loop over 5 group-specific plots
# =====================================================
remove_species = [
    '''
    "134.044784;(NH4)C4H4O4+_ppb_total",
    "168.050263;(NH4)C4H6O6+_ppb_total",
    '150.076084;(NH4)C5H8O4+_ppb_total',
    '174.076084;(NH4)C7H8O4+_ppb_total', 
    # add more here…'''
]


# Groups to iterate over
#group_order = ["BnAld", "HBnOH", "Nitroaromatic", "BPR Decomp", "Unassigned"]
group_order = ["Nitroaromatic"]

for group_name in group_order:
    if group_name not in groups_corr_SOA:
        continue

    species_list = groups_corr_SOA[group_name]
    species_in_df = [s for s in species_list if s in frac_carb.columns and s not in remove_species]
    n_species = len(species_in_df)
    if n_species == 0:
        continue

    cmap = group_cmaps.get(group_name, cm.viridis)
    color_list = [cmap(i) for i in np.linspace(0.3, 0.9, n_species)]

    # --- Loop through 1 product → all products ---
    for n in range(1, n_species + 1):
        subset = species_in_df[:n]
        subset_colors = color_list[:n]

        fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
        ax.stackplot(
            OH_exp_for_aerosol,
            [frac_carb[s] for s in subset],
            labels=[convert_to_latex(clean_label(s)) for s in subset],
            colors=subset_colors
        )

        # --- Add error bars for this subset ---
        stack_sum = frac_carb[subset].sum(axis=1)
        max_idx = stack_sum.idxmax()
        species_heights = frac_carb.loc[max_idx, subset].values
        errors = [group_contribs_frac[group_name].get(sp, 0) for sp in subset]
        cumulative_heights = np.cumsum(species_heights)

        for i, (height, cum_height, err, color) in enumerate(zip(species_heights, cumulative_heights, errors, subset_colors)):
            if err == 0:
                continue
            err = abs(err)
            if cum_height >= 0:
                lower_err = min(err, abs(cum_height))
            else:
                lower_err = min(err, abs(cum_height))
            staggered_x = base_outside_x + i * x_offset_step_detail
            ax.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                        fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

        # --- Axis formatting ---
        ax.set_xlim(0, 1.16e11)
        ax.set_ylim(0)#, 0.25 if group_name in ["BnAld", "HBnOH"] else 0.08)
        ax.legend(fontsize=18)
        ax.set_xlabel("OH Exposure (molecules cm⁻³ s)", fontsize=16)
        ax.set_ylabel("Fraction Carbon", fontsize=20)
        
        class CustomScalarFormatter(ScalarFormatter):
            def __call__(self, x, pos=None):
                if x == 0.0:
                    return "0"
                return super().__call__(x, pos)
         
        def format_y_ticks(value, pos):
            if value == 0:
                return "0"
            return f"{value:.3f}"  # two decimal places for all other ticks

        formatter_x = CustomScalarFormatter(useMathText=True)
        formatter_x.set_scientific(True)   # prevent sci notation
        ax.xaxis.set_major_formatter(formatter_x)
        ax.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))
        

        
        ax.tick_params(axis='both', which='both', labelsize=20)
        
        #ax.set_title(f"{group_name} Pathway: First {n} Products", fontsize=15)

        # --- Add top x-axis for photochemical age ---
        secax = ax.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
        secax.set_xlabel("Photochemical Age (hours)", fontsize=20)

        plt.tight_layout()
        plt.show()
        
#%%

group_order = ["BPR Decomp"]

for group_name in group_order:
    if group_name not in groups_corr_SOA:
        continue

    species_list = groups_corr_SOA[group_name]
    species_in_df = [s for s in species_list if s in frac_carb.columns and s not in remove_species]
    selected_species = [
     'm85_C4H4O2H_ppb',
     '134.044784;(NH4)C4H4O4+_ppb_total',
     '136.060434;(NH4)C4H6O4+_ppb_total',
     '132.065520;(NH4)C5H6O3+_ppb_total',
    ]
    
    # keep only selected species that exist
    species_in_df = [s for s in selected_species if s in species_in_df]
    
    n_species = len(species_in_df)
    if n_species == 0:
        continue

    cmap = group_cmaps.get(group_name, cm.viridis)
    color_list = [cmap(i) for i in np.linspace(0.3, 0.9, n_species)]

    # --- Loop through 1 product → all products ---
    for n in range(1, n_species + 1):
        subset = species_in_df[:n]
        subset_colors = color_list[:n]

        fig, ax = plt.subplots(figsize=(10, 7), dpi=300)
        ax.stackplot(
            OH_exp_for_aerosol,
            [frac_carb[s] for s in subset],
            labels=[convert_to_latex(clean_label(s)) for s in subset],
            colors=subset_colors
        )

        # --- Add error bars for this subset ---
        stack_sum = frac_carb[subset].sum(axis=1)
        max_idx = stack_sum.idxmax()
        species_heights = frac_carb.loc[max_idx, subset].values
        errors = [group_contribs_frac[group_name].get(sp, 0) for sp in subset]
        cumulative_heights = np.cumsum(species_heights)

        for i, (height, cum_height, err, color) in enumerate(zip(species_heights, cumulative_heights, errors, subset_colors)):
            if err == 0:
                continue
            err = abs(err)
            if cum_height >= 0:
                lower_err = min(err, abs(cum_height))
            else:
                lower_err = min(err, abs(cum_height))
            staggered_x = base_outside_x + i * x_offset_step_detail
            ax.errorbar(staggered_x, cum_height, yerr=[[lower_err], [err]],
                        fmt='none', ecolor=color, capsize=1, elinewidth=3.5, clip_on=False)

        # --- Axis formatting ---
        ax.set_xlim(0, 1.16e11)
        ax.set_ylim(0)#, 0.25 if group_name in ["BnAld", "HBnOH"] else 0.08)
        ax.legend(fontsize=18)
        ax.set_xlabel("OH Exposure (molecules cm⁻³ s)", fontsize=16)
        ax.set_ylabel("Fraction Carbon", fontsize=20)
        
        class CustomScalarFormatter(ScalarFormatter):
            def __call__(self, x, pos=None):
                if x == 0.0:
                    return "0"
                return super().__call__(x, pos)
         
        def format_y_ticks(value, pos):
            if value == 0:
                return "0"
            return f"{value:.4f}"  # two decimal places for all other ticks

        formatter_x = CustomScalarFormatter(useMathText=True)
        formatter_x.set_scientific(True)   # prevent sci notation
        ax.xaxis.set_major_formatter(formatter_x)
        ax.yaxis.set_major_formatter(FuncFormatter(format_y_ticks))
        

        
        ax.tick_params(axis='both', which='both', labelsize=20)
        
        #ax.set_title(f"{group_name} Pathway: First {n} Products", fontsize=15)

        # --- Add top x-axis for photochemical age ---
        secax = ax.secondary_xaxis('top', functions=(OHexp_to_age, age_to_OHexp))
        secax.set_xlabel("Photochemical Age (hours)", fontsize=20)

        plt.tight_layout()
        plt.show()