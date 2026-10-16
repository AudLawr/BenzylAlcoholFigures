# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:54:16 2024

@author: Willis Lab
"""

#script adapted from Tucker Melles with contributions from Huiying Luo

import pandas as pd
import glob
from datetime import datetime, timedelta
import datetime as dt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import re
df_pwl = pd.read_excel("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/wall_loss/PWL_results20231103.xlsx")
def Gauss(x, a1, b1, c1):
    y = a1 * np.exp(-((x - b1) / c1) ** 2)
    return y

pwlx1 = np.nanmean(df_pwl.iloc[2, 1:3]) # mean of AMS and SMPS estimates
pwlx2 = np.nanmean(df_pwl.iloc[2, 3:5])
if np.isnan(pwlx1) and np.isnan(pwlx2):
    pass
else:
    if np.isnan(pwlx1):
        pwlx1 = pwlx2
    elif np.isnan(pwlx2):
        pwlx2 = pwlx1


EESI_data = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/20220425_BnOH_EESI_CES.csv")#, skiprows=1)
EESI_time = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/20220425_BnOH_EESI_time.csv")
#EESI_zeros = pd.read_csv("E:/CSU_research/VCP/VCP_data/NOAA_EESI/20220321_carbitol/igor_exp/2022.03.21_Workup_TM_2_zeros.csv")
A = pd.concat([EESI_data, EESI_time], axis=1)
#A = A[A["'USER_DO0 monitor []'"] != 1]
A.rename(columns={'t_start_MDT': 't_series'}, inplace=True)
print(A.columns)
## estimate SOA6 (the AMS+SMPS mean PWL corrected SOA)

SOA6 = np.zeros(len(A.t_series))
SOA6[:] = np.nan

C = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/wall_loss/20220425_SMPS_dNdlogDp_raw.csv", skiprows=0)

on = '2022-04-25 10:45:00'
date_Lon = datetime.strptime(on, '%Y-%m-%d %H:%M:%S')
off = '2022-04-25 15:53:00'
date_Loff = datetime.strptime(off, '%Y-%m-%d %H:%M:%S')

t3 = date_Lon
t4 = date_Loff
A['t_series'] = pd.to_datetime(A['t_series'])
i3 = np.min(np.where(A.t_series >= t3))


#%%
columns = C.filter(like='.').columns
dp = np.zeros(len(columns))
dp[:] = np.nan
for i in range(len(columns)):
    temp = re.findall(r"[-+]?(?:\d*\.*\d+)", columns[i])
    dp[i] = float(temp[0])


ntime_smps = C.shape[0]  # every 8-9 min
nsize_smps_full = len(dp)
smps_data = C.iloc[:, 1: nsize_smps_full + 1]
bin_cutoff = 4
nsize_smps = nsize_smps_full - bin_cutoff

# for raw calculation in Step 1
number = np.zeros(ntime_smps)
surfarea = np.zeros(ntime_smps)
volume = np.zeros(ntime_smps)
numbermat = np.zeros([ntime_smps, nsize_smps])
volumemat = np.zeros([ntime_smps, nsize_smps])

x2 = 10 ** np.arange(1, 3.301, 0.1)
nsize_smps2 = len(x2)
dp2 = x2

# for fitted data in Steps 3 and 5
volume2 = np.zeros(ntime_smps)
volumemat2 = np.zeros([ntime_smps, nsize_smps2])
number2 = np.zeros(ntime_smps)
numbermat2 = np.zeros([ntime_smps, nsize_smps2])

# for excess data in Step 6
number_exc = np.zeros(ntime_smps)
volume_exc = np.zeros(ntime_smps)
stime = np.zeros(ntime_smps)

#%%

SOA6susp = np.zeros(len(A.t_series))
MassAS=np.zeros(len(A.t_series))
kwallfinal_t = np.zeros(len(A.t_series))
i5 = np.min(np.where(A.t_series >= t3))

# Step 1 AS mass conc. at t0
i = i5 
dens_ammsulf = 1.77 
seed_at_LO = np.nanmean(volume2[i3 - 1:i3 + 2])
MassAS[i]=seed_at_LO*dens_ammsulf
test=np.zeros(len(A.t_series))
test[i]=1
SOA6susp[i]=0
SOAonwall6=0
print(pwlx1)
print(pwlx2)
indices = np.where(A.t_series >= t4)[0]
if len(indices) > 0:
    i6 = np.min(np.where(A.t_series >= t4))
else: #AMS off before end
    t4 = A.t_series.iloc[-1]#TM the last value if NaT for some reason so used
    print("t4 set to:", t4)  # Print the value of t4 for debugging purposes
    i6 = np.min(np.where(A.t_series >= t4))
tsa = np.array([t.hour + t.minute / 60 + t.second / 3600 for t in A.t_series])
sus_amt = A['m213_C7H10NaO6_cps_bg_corr'].values
wall_amt2 = 0
total_amt = np.zeros(len(sus_amt))
total_amt[i] = sus_amt[i]
for i in range(i5 + 1, i6 + 1):
    kwallfinal_t[i] = pwlx1 + (pwlx2 - pwlx1) / (tsa[i6] - tsa[i5]) * (tsa[i] - tsa[i5]) # PWL interpolation. i5 is the index for lights on time; tsa = np.array([t.hour + t.minute / 60 for t in A.t_series]) is the time array in hours for AMS data
    test[i]=test[i - 1] - kwallfinal_t[i] * test[i - 1] * (tsa[i] - tsa[i - 1])
    wall_amt = sus_amt[i - 1] * kwallfinal_t[i] * (tsa[i] - tsa[i - 1])
    total_amt[i] = (sus_amt[i] + wall_amt)
    
test_2 = 1/test



A.set_index('t_series', inplace = True)
#A['m213_C7H10NaO6_cps_bg_corr_pwlcorr'] = A['m213_C7H10NaO6_cps_bg_corr']*test_2
for column in A.columns:
    A[column + '_pwlcorr'] = A[column] * test_2
    
# Replace inf with NaN
A.replace([np.inf, -np.inf], np.nan, inplace=True)
    
A.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/20220425_BnOH_EESI_CES_PWLcorr.csv")
    
    
A['m213_C7H10NaO6_cps_bg_corr_2'] = total_amt
A_1min = A.resample('min').mean()
A_1min = A_1min[(A_1min.index >= date_Lon) & (A_1min.index <= date_Loff)]


elapsed_time = A_1min.index - A_1min.index.min()
elapsed_time_seconds = elapsed_time.total_seconds()
A_1min.index = elapsed_time_seconds
A_1min['OH exp'] = A_1min.index*284516505/60

fig1, ax = plt.subplots(figsize=(8,6), dpi = 800)
ax.scatter(A_1min['OH exp'], A_1min['m213_C7H10NaO6_cps_bg_corr'], label = 'C$_7$H$_{10}$O$_6$ uncorrected', s=12, c = "steelblue")
ax.scatter(A_1min['OH exp'], A_1min['m213_C7H10NaO6_cps_bg_corr_pwlcorr'], label = 'C$_7$H$_{10}$O$_6$ corrected', s=12, c = "darkorange")
#ax.scatter(A_1min.index, A_1min['221.063173;C6H14O7Na+_corr_2'], label = 'corrected 2', s=5)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.legend(loc='upper right', fontsize=14, frameon=False)
ax.set_xlim(0, 1.0E11)  # Set x-axis range from 0 to 1.0E11
ax.set_ylim(0, None)  


ax.set_ylabel('EESI intensity (c/s)', fontsize=14)
ax.set_xlabel('OH exposure molec cm$^{-3}$ s' , fontsize=14)




#A_1min.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/EESI/20220425_BnOH_EESI_CES_PWLcorr.csv")
