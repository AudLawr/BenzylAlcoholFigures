# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:05:20 2026

@author: audlaw
"""

#plots figure S3

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import re
from matplotlib.ticker import ScalarFormatter


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

#**************************************************************
#branching fraction corrections baybeeeeee
#u need to come back to this one, its not accurate

OH_rates = {"Benzyl alcohol" : 2.8e-11,
            "HBnOH" : 5.59e-11, 
            "BnAld": 1.29e-11,
            "Phen": 2.83e-11, 
            "Frag": 3.45e-11} #come back to this one

bnoh0 = merged_df_corr['C7H8O_ppb'].iloc[0]
bnoh_ratio = merged_df_corr['C7H8O_ppb'] / bnoh0
hbnoh_ratio = merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']/merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total'].iloc[0]
    
#HBnOH_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["HBnOH"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["HBnOH"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
HBnOH_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["HBnOH"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["HBnOH"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))

HBnOH_CFh = ((OH_rates["Benzyl alcohol"] - OH_rates["HBnOH"])/OH_rates["Benzyl alcohol"]) * ((1-hbnoh_ratio)/(hbnoh_ratio**(OH_rates["HBnOH"]/OH_rates["Benzyl alcohol"]) - hbnoh_ratio))
BnAld_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["BnAld"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["BnAld"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
Phen_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["Phen"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["Phen"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))

HBnOH_OHcorr = merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']*HBnOH_CF
BnAld_OHcorr = merged_df_corr['m107_C7H6OH_ppb_matt']*BnAld_CF
Phen_OHcorr = merged_df_corr['m95_C6H6OH_ppb']*Phen_CF
#THIS IS LOOKING AT CORRECTING THE BRANCHING FRACTIONS, NEED TO COME BACK TO IT

fig, ax = plt.subplots(figsize = (10, 6), dpi = 300)
ax.scatter(OH_exp, merged_df_corr['m107_C7H6OH_ppb_matt'], label = "Benzaldehyde", color = "#D95A30")
ax.scatter(OH_exp, BnAld_OHcorr, label = "Benzaldehyde corr for OH ox", color = "#E28413")
ax.scatter(OH_exp, merged_df_corr['m95_C6H6OH_ppb'], label = "Phenol", color = "#377216")
ax.scatter(OH_exp, merged_df_corr['m95_C6H6OH_ppb']*Phen_CF, label = "Phenol corr for OHox", color = "#6CAC34")
#ax.scatter(OH_exp, merged_df_corr['m61_C2H4O2H_ppb'], label = "c2", color = "grey")

#plt.scatter(OH_exp, NH4_ppb_df_deltaC['142.086255;(NH4)C7H8O2+_ppb'], label = "C7H8O2 no VWL")


'''
fig,ax = plt.subplots(figsize = (10, 6), dpi = 300)
#ax.scatter(OH_exp, merged_df_corr['174.076084;(NH4)C7H8O4+_ppb_total'], label = 'epox', color = 'plum')
#ax.scatter(OH_exp, merged_df_corr['m73_C3H4O2H_ppb'], label = 'frag', color = 'c')
ax.scatter(OH_exp, merged_df_corr['140.070605;(NH4)C7H6O2+_ppb_total'], label = 'C7H6O2 nh4', color = 'steelblue')
ax.scatter(OH_exp, merged_df_corr['m123_C7H6O2H_ppb'], label = 'C7H6O2 ptr', color = 'orange')
plt.legend()
plt.show()

bnoh_pred = bnoh0 * np.exp(-OH_rates["Benzyl alcohol"] * OH_exp)
bnoh_react = bnoh0 - bnoh_pred
bnoh_tot = merged_df_corr['C7H8O_ppb'] + bnoh_react
plt.scatter(OH_exp, bnoh_tot)
plt.ylim(0,50)
plt.show()
'''
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
grouped_frag = grouped_frag.where(grouped_frag >= 0.1).dropna()

Frag_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["Frag"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["Frag"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
Frag_OHcorr = grouped_frag['BPR Frag']*Frag_CF

OH_exp_S = pd.Series(OH_exp)
grouped_frag['OH_exp'] = OH_exp_S.loc[grouped_frag.index]
grouped_frag.set_index('OH_exp', inplace=True)


plt.scatter(grouped_frag.index, grouped_frag, label = "Decomp prods", color = "#4B2194")
plt.scatter(OH_exp, Frag_OHcorr, label = "Decomp prods corr for OH ox", color = "#8657DE")
ax.scatter(OH_exp, merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total'], label = 'Hydroxybenzyl alcohol', color = "#2835C7")
ax.scatter(OH_exp, merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']*HBnOH_CF, label = 'Hydroxybenzyl alcohol corrected for OH oxidation', color = "#4589D7")
#calcs branching fraction for frag, need to revisit this
#frag_y = grouped_frag['BPR Frag'].nlargest(10).mean()/bnoh0
frag_y = Frag_OHcorr.nlargest(10).mean()/bnoh0
print(f"Frag Yield: {frag_y * 100:.1f}%")

#thop ten points
'''
Frag_OHcorr.index = OH_exp
top10_frag = Frag_OHcorr.nlargest(10)
plt.scatter(top10_frag.index, top10_frag, s=120,
    facecolors='none', edgecolors='black', linewidths=2)

HBnOH_OHcorr.index = OH_exp
top10_hbnoh = HBnOH_OHcorr.nlargest(10)
plt.scatter(top10_hbnoh.index, top10_hbnoh, s=120,
    facecolors='none', edgecolors='black', linewidths=2)
'''

plt.xlim(0)
plt.ylim(0, 15)
plt.xlabel("OH Exposure (molecules cm⁻³ s)")
plt.ylabel("Concentration (ppb)")
plt.legend(loc='upper right', fontsize = 11, ncol = 2)
#plt.title(titles)
plt.show()

#%%

#yield the way they do in beuconsejo

bnoh0 = merged_df_corr['C7H8O_ppb'].iloc[0]
bnoh_ratio = merged_df_corr['C7H8O_ppb'] / bnoh0
bnoh_rxt = bnoh0 - merged_df_corr['C7H8O_ppb']


bnald_yield_B = merged_df_corr['m107_C7H6OH_ppb_matt']/bnoh_rxt
BnAld_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["BnAld"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["BnAld"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
bnald_yield_B_corr = bnald_yield_B * BnAld_CF

hbnoh_yield_B = merged_df_corr['142.086255;(NH4)C7H8O2+_ppb_total']/bnoh_rxt
HBnOH_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["HBnOH"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["HBnOH"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
hbnoh_yield_B_corr = hbnoh_yield_B * HBnOH_CF

phen_yield_B = merged_df_corr['m95_C6H6OH_ppb']/bnoh_rxt
Phen_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["Phen"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["Phen"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
phen_yield_B_corr = phen_yield_B * Phen_CF

firstgen_frag = {"BPR Frag": ['80.070605;(NH4)C2H6O2+_ppb_total', '118.049870;(NH4)C4H4O3+_ppb_total', '120.065520;(NH4)C4H6O3+_ppb_total', '132.065520;(NH4)C5H6O3+_ppb_total', '134.044784;(NH4)C4H4O4+_ppb_total', '134.081170;(NH4)C5H8O3+_ppb_total', '148.060434;(NH4)C5H6O4+_ppb_total', '168.050263;(NH4)C4H6O6+_ppb_total', '150.076084;(NH4)C5H8O4+_ppb_total', '174.076084;(NH4)C7H8O4+_ppb_total']}
grouped_frag_B = {group: merged_df_corr[species].sum(axis=1) for group, species in firstgen_frag.items()}
grouped_frag_B = pd.DataFrame(grouped_frag_B)
grouped_frag_B = grouped_frag_B.where(grouped_frag_B >= 0.01).dropna()
Frag_yield_B = grouped_frag_B["BPR Frag"] / bnoh_rxt
Frag_CF = ((OH_rates["Benzyl alcohol"] - OH_rates["Frag"])/OH_rates["Benzyl alcohol"]) * ((1-bnoh_ratio)/(bnoh_ratio**(OH_rates["Frag"]/OH_rates["Benzyl alcohol"]) - bnoh_ratio))
frag_yield_B_corr = Frag_yield_B*Frag_CF


corr_yields = {
    "BnAld": bnald_yield_B_corr,
    "HBnOH": hbnoh_yield_B_corr,
    "Phen": phen_yield_B_corr, 
    "Decomp": frag_yield_B_corr, }

target = 0.5 * bnoh0

idx_65 = (bnoh_rxt - target).abs().idxmin()

bnald_val = bnald_yield_B_corr.loc[idx_65]
hbnoh_val = hbnoh_yield_B_corr.loc[idx_65]
phen_val = phen_yield_B_corr.loc[idx_65]

# print results
print("Yields at 50% BnOH remaining:")
print(f"BnAld corrected yield: {bnald_val:.4f}")
print(f"HBnOH corrected yield: {hbnoh_val:.4f}")
print(f"Phen corrected yield: {phen_val:.4f}")
print("\n")

for name, series in corr_yields.items():
    avg = series.tail(25).mean()
    
    print(f"{name} corrected yield (last 25 avg): {avg:.4f}")
    
print(type(OH_exp), OH_exp.shape)
print(type(Frag_yield_B), Frag_yield_B.shape)

fig, ax = plt.subplots(figsize = (10, 6), dpi = 300)
x_50 = OH_exp[idx_65]
plt.scatter(OH_exp, bnald_yield_B, label="Benzaldehyde yield", color="#D95A30")
plt.scatter(OH_exp, bnald_yield_B_corr, label="Benzaldehyde yield OH corr", color="#E28413")
plt.scatter(OH_exp, phen_yield_B, label="Phenol yield", color="#377216")
plt.scatter(OH_exp, phen_yield_B_corr, label="Phenol yield OH corr", color="#6CAC34")
plt.scatter(OH_exp, Frag_yield_B, label = "Decomp prods", color = "#4B2194")
plt.scatter(OH_exp, frag_yield_B_corr, label = "Decomp prods corr for OH ox", color = "#8657DE")
plt.scatter(OH_exp, hbnoh_yield_B, label="Hydroxybenzyl alcohol yield", color="#2835C7")
plt.scatter(OH_exp, hbnoh_yield_B_corr, label="Hydroxybenzyl alcohol yield OH corr", color="#4589D7")

# vertical line for 50% reacted
plt.axvline(x=x_50, color="black", linestyle="--", label="50% BnOH reacted")

plt.xlim(0)
plt.ylim(0, 1)
plt.ylabel("Yield ([Prod]t/[BnOH reacted])", fontsize = 14)
plt.xlabel("OH Exposure (molecules cm⁻³ s)", fontsize = 14)
plt.legend(fontsize=12)

plt.show()
