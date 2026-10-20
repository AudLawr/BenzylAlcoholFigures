# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:17:57 2026

@author: audlaw
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

#this plots figure S4 in the most inefficient way possible

# Load csvs
SMPS_0328 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220328_BenzylAlcohol/Python Analylsis/0328SMPSTS_v0.csv")
AMSSUM_0328 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220328_BenzylAlcohol/Python Analylsis/0328SOAamssum.csv")
SMPS_0425 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SMPSTS_v0.csv")
AMSSUM_0425 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SOAamssum.csv")
SMPS_0706 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220706_benzyl_alcohol/Python_analysis/0706SMPSTS_v0.csv")
AMSSUM_0706 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220706_benzyl_alcohol/Python_analysis/0706SOAamssum.csv")
SMPS_0708 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220708_benzyl_alcohol/Python_analysis/0708SMPSTS_v0.csv")
AMSSUM_0708 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220708_benzyl_alcohol/Python_analysis/0708SOAamssum.csv")
SMPS_0722 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220722_benzyl_alcohol/Python_analysis/0722SMPSTS_v0.csv")
AMSSUM_0722 = pd.read_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220722_benzyl_alcohol/Python_analysis/0722SOAamssum.csv")
charan_2020 = pd.read_csv("C:/Users/audlaw/git_push_pull/BnOH_manuscript/charan_2020_fig2.csv")

# Convert time columns to datetime
SMPS_0328["Time"] = pd.to_datetime(SMPS_0328["Time"])
AMSSUM_0328["Time"] = pd.to_datetime(AMSSUM_0328["Time"])
SMPS_0425["Time"] = pd.to_datetime(SMPS_0425["Time"])
AMSSUM_0425["Time"] = pd.to_datetime(AMSSUM_0425["Time"])
SMPS_0706["Time"] = pd.to_datetime(SMPS_0706["Time"])
AMSSUM_0706["Time"] = pd.to_datetime(AMSSUM_0706["Time"])
SMPS_0708["Time"] = pd.to_datetime(SMPS_0708["Time"])
AMSSUM_0708["Time"] = pd.to_datetime(AMSSUM_0708["Time"])
SMPS_0722["Time"] = pd.to_datetime(SMPS_0722["Time"])
AMSSUM_0722["Time"] = pd.to_datetime(AMSSUM_0722["Time"])

# Select only the following columns from the soa csv
AMSSUM_0328 = AMSSUM_0328[["Time", "SOA_AVGpwl", "OHexp", "photohravg", "SOAAMSAVGpwl_vocobs_yield"]]
AMSSUM_0425 = AMSSUM_0425[["Time", "SOA_AVGpwl", "OHexp", "photohravg", "SOAAMSAVGpwl_vocobs_yield"]]
AMSSUM_0706 = AMSSUM_0706[["Time", "SOA_AVGpwl", "OHexp", "photohravg", "SOAAMSAVGpwl_vocobs_yield"]]
AMSSUM_0708 = AMSSUM_0708[["Time", "SOA_AVGpwl", "OHexp", "photohravg", "SOAAMSAVGpwl_vocobs_yield"]]
AMSSUM_0722 = AMSSUM_0722[["Time", "SOA_AVGpwl", "OHexp", "photohravg", "SOAAMSAVGpwl_vocobs_yield"]]

# Merge based on time column, keeping only matching times
merged_0328 = pd.merge(SMPS_0328, AMSSUM_0328, on="Time", how="inner")
merged_0425 = pd.merge(SMPS_0425, AMSSUM_0425, on="Time", how="inner")
merged_0708 = pd.merge(SMPS_0708, AMSSUM_0708, on="Time", how="inner")
merged_0722 = pd.merge(SMPS_0722, AMSSUM_0722, on="Time", how="inner")

#the 0706 file is merged on the nearest minute timestamp because the surface area data come from the TSI and not the SMPS so the SA data are off by 10s
merged_0706 = pd.merge(SMPS_0706, AMSSUM_0706, left_on = SMPS_0706['Time'].dt.floor('min'), right_on=AMSSUM_0706['Time'].dt.floor('min'), how = "inner")
merged_0706.drop(columns=['key_0'], errors='ignore', inplace=True) #drop extra merging columns
merged_0706.drop(columns=['Time_y'], errors='ignore', inplace=True)
merged_0706.rename(columns = {'Time_x': 'Time'}, inplace = True)
merged_0706.dropna(inplace=True)#drop nans

# Save the merged DataFrame
merged_0328.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220328_BenzylAlcohol/Python Analylsis/0328SMPSTS_v1.csv", index=False)
merged_0425.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/0425SMPSTS_v1.csv", index=False)
merged_0706.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220706_benzyl_alcohol/Python_analysis/0706SMPSTS_v1.csv", index=False)
merged_0708.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220708_benzyl_alcohol/Python_analysis/0708SMPSTS_v1.csv", index=False)
merged_0722.to_csv("C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220722_benzyl_alcohol/Python_analysis/0722SMPSTS_v1.csv", index=False)

merged_csvs = {'0328':merged_0328, '0425':merged_0425, '0706':merged_0706, '0708': merged_0708, '0722': merged_0722}


#import variables
OHexp_0328 = merged_0328['OHexp']
SOA_AVGpwl_0328 = merged_0328["SOA_AVGpwl"]
SA_0328 = merged_0328["area_fit"]
phohr_0328 = merged_0328['photohravg']
SOA_yield_0328 = merged_0328['SOAAMSAVGpwl_vocobs_yield']

OHexp_0425 = merged_0425['OHexp']
SOA_AVGpwl_0425 = merged_0425["SOA_AVGpwl"]
SA_0425 = merged_0425["area_fit"]
phohr_0425 = merged_0425['photohravg']
SOA_yield_0425 = merged_0425['SOAAMSAVGpwl_vocobs_yield']

OHexp_0706 = merged_0706['OHexp']
SOA_AVGpwl_0706 = merged_0706["SOA_AVGpwl"]
SA_0706 = merged_0706["area_fit"]
phohr_0706 = merged_0706['photohravg']
SOA_yield_0706 = merged_0706['SOAAMSAVGpwl_vocobs_yield']

OHexp_0708 = merged_0708['OHexp']
SOA_AVGpwl_0708 = merged_0708["SOA_AVGpwl"]
SA_0708 = merged_0708["area_fit"]
phohr_0708 = merged_0708['photohravg']
SOA_yield_0708 = merged_0708['SOAAMSAVGpwl_vocobs_yield']

OHexp_0722 = merged_0722['OHexp']
SOA_AVGpwl_0722 = merged_0722["SOA_AVGpwl"]
SA_0722 = merged_0722["area_fit"]
phohr_0722 = merged_0722['photohravg']
SOA_yield_0722 = merged_0722['SOAAMSAVGpwl_vocobs_yield']

dfs = {
       "SA_0328": SA_0328,
       "SA_0425": SA_0425,
       "SA_0706": SA_0706,
       "SA_0708": SA_0708,
       "SA_0722": SA_0722}

for name, df in dfs.items():
    mean_val = df.mean()
    initial_val = df.iloc[0]
    
    print(f"{name}:initial = {initial_val},  mean = {mean_val}")


#%%

SOA_yield_0328 = SOA_yield_0328.interpolate(method='linear')
SOA_yield_0425 = SOA_yield_0425.interpolate(method='linear')
SOA_yield_0706 = SOA_yield_0706.interpolate(method='linear')
SOA_yield_0708 = SOA_yield_0708.interpolate(method='linear')
SOA_yield_0722 = SOA_yield_0722.interpolate(method='linear')

# SOA yield and OHexp colored by sus SA
fig, ax12 = plt.subplots(figsize=(14,6))

# Normalize colorbar to apply the same scale to all scatter plots
norm = mcolors.Normalize(vmin=min(SA_0328.min(), SA_0425.min(), SA_0706.min(), SA_0708.min(), SA_0722.min()), 
                          vmax=max(SA_0328.max(), SA_0425.max(), SA_0706.max(), SA_0708.max(), SA_0722.max()))


def add_uncertainty_band(ax, x, y, color, alpha=0.2):
    # sort so fill_between doesn't zig-zag
    idx = np.argsort(x)
    x_sorted = np.array(x)[idx]
    y_sorted = np.array(y)[idx]

    lower = 0.6 * y_sorted   # -40%
    upper = 1.4 * y_sorted   # +40%

    ax.fill_between(
        x_sorted,
        lower,
        upper,
        color="gray",
        alpha=alpha,
        linewidth=0,
        zorder=0,
        label = "uncertainty")

# pick colors consistent with viridis (or approximate by sampling cmap)
cmap = plt.cm.viridis
colors = [cmap(0.15), cmap(0.35), cmap(0.55), cmap(0.75), cmap(0.9)]

#plot the stupid uncertainty bands FIRST
add_uncertainty_band(ax12, OHexp_0328, SOA_yield_0328, colors[0])
add_uncertainty_band(ax12, OHexp_0425, SOA_yield_0425, colors[1])
add_uncertainty_band(ax12, OHexp_0706, SOA_yield_0706, colors[2])
add_uncertainty_band(ax12, OHexp_0708, SOA_yield_0708, colors[3])
add_uncertainty_band(ax12, OHexp_0722, SOA_yield_0722, colors[4])


scatter_0328 = ax12.scatter(
    OHexp_0328, SOA_yield_0328, 
    c=SA_0328, cmap="viridis", edgecolor="k", norm=norm, label="Expt. 1", marker = "o", s = 62)
scatter_0425 = ax12.scatter(
    OHexp_0425, SOA_yield_0425, 
    c=SA_0425, cmap="viridis", edgecolor="k", norm=norm, label="Expt. 2", marker = "^", s = 75)
scatter_0706 = ax12.scatter(
    OHexp_0706, SOA_yield_0706, 
    c=SA_0706, cmap="viridis", edgecolor="k", norm=norm, label="Expt. 3", marker = "P", s = 95)
scatter_0708 = ax12.scatter(
    OHexp_0708, SOA_yield_0708, 
    c=SA_0708, cmap="viridis", edgecolor="k", norm=norm, label="Expt. 4", marker = ",", s = 60)
scatter_0722 = ax12.scatter(
    OHexp_0722, SOA_yield_0722, 
    c=SA_0722, cmap="viridis", edgecolor="k", norm=norm, label="Expt. 5", marker = "D", s = 60)

# Label axes
ax12.set_xlabel(r"OH Exposure (molec $cm^{-3} s$)", fontsize=17)
ax12.set_ylabel(r"SOA Mass Yield", fontsize=17)

# Set X-axis limit
ax12.set_xlim(0)

# Adjust tick font size for axes
ax12.tick_params(axis='both', labelsize=17)

# Create secondary X-axis at the top
ax13_ = ax12.twiny()
ax13_.set_xlim(0, 30)
ax13_.set_xlabel(r"Photochemical Age (hr)", fontsize=17)

# Adjust tick font size for secondary axis
ax13_.tick_params(axis='x', labelsize=16)

# Color bar (shared for both scatter plots)
cbar = plt.colorbar(scatter_0328, ax=ax12)
cbar.set_label(r"Suspended Aerosol Surface Area ($\mu m^2 cm^{-3}$)", fontsize=15)

handles, labels = ax12.get_legend_handles_labels()

unique = dict(zip(labels, handles))  # removes duplicates

ax12.legend(unique.values(), unique.keys(), markerscale=1.5, fontsize=16)

plt.show()

#%%
'''
# Add experiment date columns
SOA_yield_0328["exp_date"] = "20220328"
SOA_yield_0425["exp_date"] = "20220425"
SOA_yield_0706["exp_date"] = "20220706"
SOA_yield_0708["exp_date"] = "20220708"
SOA_yield_0722["exp_date"] = "20220722"
'''
#combine those suckers
all_SOA_yield_df = pd.concat([
    pd.DataFrame({
        "exp_date": "20220328",
        "OHexp": OHexp_0328.to_numpy(),
        "SOA_yield": SOA_yield_0328.to_numpy(),
        "SOA_AVGpwl": SOA_AVGpwl_0328.to_numpy(),
        "SA": SA_0328.to_numpy(),
        "photohr": phohr_0328.to_numpy(),
    }),

    pd.DataFrame({
        "exp_date": "20220425",
        "OHexp": OHexp_0425.to_numpy(),
        "SOA_yield": SOA_yield_0425.to_numpy(),
        "SOA_AVGpwl": SOA_AVGpwl_0425.to_numpy(),
        "SA": SA_0425.to_numpy(),
        "photohr": phohr_0425.to_numpy(),
    }),

    pd.DataFrame({
        "exp_date": "20220706",
        "OHexp": OHexp_0706.to_numpy(),
        "SOA_yield": SOA_yield_0706.to_numpy(),
        "SOA_AVGpwl": SOA_AVGpwl_0706.to_numpy(),
        "SA": SA_0706.to_numpy(),
        "photohr": phohr_0706.to_numpy(),
    }),

    pd.DataFrame({
        "exp_date": "20220708",
        "OHexp": OHexp_0708.to_numpy(),
        "SOA_yield": SOA_yield_0708.to_numpy(),
        "SOA_AVGpwl": SOA_AVGpwl_0708.to_numpy(),
        "SA": SA_0708.to_numpy(),
        "photohr": phohr_0708.to_numpy(),
    }),

    pd.DataFrame({
        "exp_date": "20220722",
        "OHexp": OHexp_0722.to_numpy(),
        "SOA_yield": SOA_yield_0722.to_numpy(),
        "SOA_AVGpwl": SOA_AVGpwl_0722.to_numpy(),
        "SA": SA_0722.to_numpy(),
        "photohr": phohr_0722.to_numpy(),
    }),

], ignore_index=True)

# Save to CSV
all_SOA_yield_df.to_csv(
    r"C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/Dataverse_for_upload/SI/FigureS04/SOA_yield_all_experiments.csv",
    index=False)

#%%

tmp = merged_0706[['OHexp',
                   'SOAAMSAVGpwl_vocobs_yield',
                   'area_fit',
                   'photohravg']]

print(tmp.shape)

print(tmp[tmp.isna().any(axis=1)])
print(tmp.tail(10))

