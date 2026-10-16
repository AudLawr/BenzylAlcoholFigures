# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 10:11:36 2024

@author: trmel
"""

#this is adapted from tucker melles' script

import math
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['font.size']= 22
mpl.rcParams['lines.linewidth']= 2
mpl.rcParams['xtick.major.size']= 10
mpl.rcParams['xtick.minor.size']= 5
mpl.rcParams['xtick.major.width']= 1
mpl.rcParams['xtick.minor.width']= 1
mpl.rcParams['ytick.major.size']= 10
mpl.rcParams['ytick.minor.size']= 5
mpl.rcParams['ytick.major.width']= 1
mpl.rcParams['ytick.minor.width']= 1
mpl.rcParams['axes.labelsize'] = 22


#%%
#alright simplifying the heck out of this 
#with a mere three reactions for three structures, so nine i guess
#just starting with the bprs

T = 298.15 #temp in K
avo = 6.022e23 #avogadros number
R = 8.314 # R in m3 * Pa K-1 mol-1
P = 85000  # FC pressure in Pa (850 mbar)
#P = 104000

OH_conc = 1.5E6

#HO2_ppb = 0.35 #this is from MCM
HO2_ppb = 0.41 #0425
#HO2_ppb = 0.46
NO_ppb = 1.1 #0425
#NO_ppb = 1.2
BnOH_ppb = 38

NO_conc = ((NO_ppb*10**-9*avo*P)/(R*T))/100**3 #in molec/cm3
HO2_conc = ((HO2_ppb*10**-9*avo*P)/(R*T))/100**3 #in molec/cm3

BnOH_conc = ((NO_ppb*10**-9*avo*P)/(R*T))/100**3


#%%
parent_conc = 1.35E12

ortho_conc = parent_conc*0.15
ipso_conc = parent_conc*0.60

#ka is the ipso BPR
ka1 = 0.95
ka2 = (2.44E-12*math.exp(360/T)*NO_conc)
ka3 = (2.58E-13*math.exp(1300/T)*HO2_conc)

ka_tot = ka1 + ka2 + ka3

ka1_frac = ka1/ka_tot
ka2_frac = ka2/ka_tot
ka3_frac = ka3/ka_tot

#kb is the ortho BPR CHOH
kb1 = 20
kb2 = 2.59E-12*math.exp(360/T)*NO_conc
kb3 = 2.58E-13*math.exp(1300/T)*HO2_conc

kb_tot = kb1 + kb2 + kb3

kb1_frac = kb1/kb_tot
kb2_frac = kb2/kb_tot
kb3_frac = kb3/kb_tot

#kc is the ortho BPR CH2
kc1 = 114
kc2 = 2.59E-12*math.exp(360/T)*NO_conc
kc3 = 2.58E-13*math.exp(1300/T)*HO2_conc

kc_tot = kc1 + kc2 + kc3

kc1_frac = kc1/kc_tot
kc2_frac = kc2/kc_tot
kc3_frac = kc3/kc_tot

#kd is the para BPR 

kd2 = 2.59E-12*math.exp(360/T)*NO_conc
kd3 = 2.58E-13*math.exp(1300/T)*HO2_conc

kd_tot = kd2 + kd3

kd2_frac = kd2/kd_tot
kd3_frac = kd3/kd_tot


print(f"Pseudo first order rate for Ipso OH BPR with NO: {ka2:.4f} s-1")
print(f"Pseudo first order rate for Ipso OH BPR with HO2: {ka3:.4f} s-1")
print(f"Pseudo first order rate for Ortho OH BPRs with NO: {kb2:.4f} s-1")
print(f"Pseudo first order rate for Ortho OH BPRs with HO2: {kb3:.4f} s-1")

ipso_BPR_kbi = ka2+ka3 #in /seconds
ortho_BPR_kbi = kb2+kb3

ipso_BPR_tbi = 1/ipso_BPR_kbi
ortho_BPR_tbi = 1/ortho_BPR_kbi 

print("\n")
print(f"Bimolecular rate of ipso BPR with NO and HO2 = {ipso_BPR_kbi:.4f} s-1")
print(f"bimolecular lifetime = {ipso_BPR_tbi:.4f} s")

print("\n")
print(f"Bimolecular rate or ortho BPRs with NO and HO2 = {ortho_BPR_kbi:.4f} s-1")
print(f"bimolecular lifetime = {ortho_BPR_tbi:.4f} s")
print("\n")

scenarios = {
    "ipso-OH meta-OO CH$_2$": { #ka ipso CH2 (blue); colors go from lightest to dark, with H-shift matching fig 5 color
        "rates": {"H-shift": ka1_frac, "NO": ka2_frac, "HO$_2$": ka3_frac},
        "colors": {"H-shift": "#4AC6FF", "NO": "#25A9E6", "HO$_2$": "#008BCC"}},
    "ortho-OH ortho-OO CHOH": { #ortho CHOH (teal)
        "rates": {"H-shift": kb1_frac, "NO": kb2_frac, "HO$_2$": kb3_frac},
        "colors": {"H-shift": "#65DDC9", "NO": "#41B29F", "HO$_2$": "#1D8674"}},
    "ortho-OH ortho-OO CH$_2$": {#ortho CH2 (green)
        "rates": {"H-shift": kc1_frac, "NO": kc2_frac, "HO$_2$": kc3_frac},
        "colors": {"H-shift": "#B3D90B", "NO": "#A2C40A", "HO$_2$": "#90AE09"}},
    "ortho-OH para-OO BPR": { #ka ipso CH2 (blue); colors go from lightest to dark, with H-shift matching fig 5 color
        "rates": {"NO": kd2_frac, "HO$_2$": kd3_frac},
        "colors": {"NO": "#1AAE87", "HO$_2$": "#158E6E"}}}


#%%
manual_labels = {
    ("ortho-OH ortho-OO CHOH", "NO"): (-1, 0.5),
    ("ortho-OH ortho-OO CHOH", "HO$_2$"): (-1, 1),
    ("ortho-OH ortho-OO CH$_2$", " "): (-1, 0.5), #blank labels for these because they dont show up
    ("ortho-OH ortho-OO CH$_2$", "  "): (-1, 1)}

fig, ax = plt.subplots(figsize=(10, 8))

for i, (exp_name, data) in enumerate(scenarios.items()):
    rates = data["rates"]
    colors = data["colors"]
    left = 0
    for path, val in rates.items():
        ax.barh(exp_name, val, left=left, color=colors[path], edgecolor="white")
        # Only annotate manually if this segment is in the manual_labels list
        if (exp_name, path) in manual_labels:
            x_offset, y_offset = manual_labels[(exp_name, path)]
            ax.annotate(path,
                        xy=(left + val, i),       # point at end of segment
                        xytext=(x_offset*50, y_offset*50),  # scaled for visibility
                        textcoords='offset points',
                        ha='left', va='center',
                        fontsize=13, weight='bold', color='w',
                        arrowprops=dict(arrowstyle="->", color='w', lw=1.2))

        else:
            # Automatic label inside bar if segment is large enough
            if val > 0.1:  # threshold for visibility
                ax.text(left + val/2, i, path,
                        ha='center', va='center',
                        fontsize=13, weight='bold', color='white')
        
        left += val

ax.set_xlim(0, 1)
ax.set_xlabel("Fraction of BPR Fate")
plt.tight_layout()
plt.show()

#%%

# define structure formation fractions (sum to 1)
formation_fracs = {
    "ipso-OH meta-OO CH$_2$": 0.40/0.4375,
    "ortho-OH ortho-OO CHOH": 0.001/0.4375,
    "ortho-OH ortho-OO CH$_2$": 0.006/0.4375,  
    "ortho-OH para-OO BPR": 0.03/0.4375       
}

manual_labels = {
    ("ortho-OH ortho-OO CHOH", " "): (-1, 0.5),
    ("ortho-OH ortho-OO CHOH", "  "): (-1, 1),
    ("ortho-OH ortho-OO CH$_2$", "H-shift"): (-1.1, 1.1), #blank labels for these because they dont show up
    ("ortho-OH ortho-OO CH$_2$", "  "): (-1, 1), 
    ("ortho-OH para-OO BPR", "HO$_2$"): (-0.4, -0.5),
    ("ortho-OH para-OO BPR", "NO"): (0.01, 0.5),
    }


#%%
# define structure formation fractions (sum to 1)
formation_fracs = {
    "ipso-OH meta-OO CH$_2$": 0.40/0.4375,
    "ortho-OH ortho-OO CHOH": 0.001/0.4375,
    "ortho-OH ortho-OO CH$_2$": 0.006/0.4375,  
    "ortho-OH para-OO BPR": 0.03/0.4375       
}


manual_labels2 = {
    #("ortho-OH ortho-OO CHOH", " "): (-1, 4),
    #("ortho-OH ortho-OO CHOH", "  "): (-1, 1),
    ("ortho-OH ortho-OO CH$_2$", "1,5-H-shift"): (1, 1), #blank labels for these because they dont show up
    #("ortho-OH ortho-OO CH$_2$", "  "): (-1, 1), 
    #("ortho-OH para-OO BPR", "HO$_2$"): (-0.4, -0.5),
    #("ortho-OH para-OO BPR", "NO"): (0.01, 0.5),
    }



#build a list of segments for the overall stacked bar
segments = []

for struct_name, data in scenarios.items():
    rates = data["rates"]
    colors = data["colors"]
    form_frac = formation_fracs.get(struct_name, 0)
    
    for path, val in rates.items():
        weighted_val = val * form_frac
        segments.append({
            "structure": struct_name,
            "path": path,
            "value": weighted_val,
            "color": colors[path]
        })

# sort segments by structure order (optional, for consistent layering)
structure_order = ["ipso-OH meta-OO CH$_2$", "ortho-OH ortho-OO CHOH", "ortho-OH ortho-OO CH$_2$", "ortho-OH para-OO BPR"]
segments.sort(key=lambda x: structure_order.index(x["structure"]))


# plot single stacked bar (vertical)
fig, ax = plt.subplots(figsize=(4, 10))  # taller for vertical

bottom = 0  # start stacking from bottom
for i, seg in enumerate(segments):
    ax.bar(" ", seg["value"], bottom=bottom, color=seg["color"], width=0.5)

    # check if this structure-pathway pair is manually labeled
    key = (seg["structure"], seg["path"])
    if key in manual_labels2:
        x_offset, y_offset = manual_labels2[key]
        ax.annotate(
            seg["path"],
            xy=(0, bottom + seg["value"] / 2),        # center of bar segment
            xytext=(x_offset, y_offset),   # tweak offset for vertical orientation
            textcoords="offset points",
            ha="center", va="bottom",
            fontsize=18, weight="bold", color="k",
            arrowprops=dict(arrowstyle="->", color="w", lw=2),
            annotation_clip=False, clip_on = False, zorder = 20
        )
    else:
        # automatic label for large segments
        if seg["value"] > 0.03:
            ax.text(
                0, bottom + seg["value"] / 2,
                f"{seg['path']}",
                ha="center", va="center",
                fontsize=18, weight="bold", color="white"
            )

    bottom += seg["value"]
'''
ax.text(
    0.5, -0.05, "Estimated Fractional BPR Fate",  # move below the bar
    transform=ax.transAxes,
    fontsize=16,
    ha="center",
    bbox=dict(facecolor="white", edgecolor="white", boxstyle="round,pad=0.3")
)
'''
ax.set_ylim(0, 1)

ax.set_ylabel("Estimated Fractional BPR Fate")
ax.set_xticks([])  # hide x-axis tick
plt.tight_layout()
plt.show()
#%%

#this is to calculate RO2 specifications 
# just shows how high the RO2 needs to be if the H-shift rate is correct (blue) or if the HO2 conc is correct (green)
RO2RO2_rate = 1.2e-13 #this is the rate of cyclohexane RO2 + CH3O2


RO2_ppb_HO2 = 103
RO2_conc_HO2 = ((RO2_ppb_HO2*10**-9*avo*P)/(R*T))/100**3 #in molec/cm3
ka4_HO2 = RO2RO2_rate * RO2_conc_HO2

RO2_ppb_Hs = 383
RO2_conc_Hs = ((RO2_ppb_Hs*10**-9*avo*P)/(R*T))/100**3 #in molec/cm3
ka4_Hs = RO2RO2_rate * RO2_conc_Hs

#this is what [RO2] needs to be for RO2+RO2 to equal Hshift for ipso-OH meta RO2
Hshift_RO2_ppb = ((ka1/RO2RO2_rate)*R*T*10**9*100**3)/(avo*P)
Hshift_RO2_mcm = ka1/RO2RO2_rate

HO2_RO2_ppb = ((ka3/RO2RO2_rate)*R*T*10**9*100**3)/(avo*P)
HO2_RO2_mcm = ka3/RO2RO2_rate

species_HO2 = ["HO$_2$ + RO$_2$", "RO$_2$ + RO$_2$"]
species_Hs = ["H-shift", "RO$_2$ + RO$_2$"]
values_HO2 = [ka3, ka4_HO2]   # first RO2 value
values_Hs = [ka1, ka4_Hs]



RO2_conc_amel = 7e8 #molecule/cm3
RO2RO2_rate_pred = ka3/RO2_conc_amel
RO2RO2_rate_pred_2 = ka1/RO2_conc_amel


# Plot
plt.figure(figsize=(10, 6))
plt.barh(species_Hs, values_Hs, alpha=1, color = "#257995")
plt.barh(species_HO2, values_HO2, color = "#9EC925")

y = 1
plt.text(0.5, y, rf"[RO$_2$] = {Hshift_RO2_mcm:.2e}" "\n" rf"molecules/cm$^3$", va = "center", color = "#B6EDFF", fontsize = 10)
plt.text(0.05, y, rf"[RO$_2$] = {HO2_RO2_mcm:.2e}" "\n" rf"molecules/cm$^3$", va = "center", color = "#276C09", fontsize = 10)

plt.xlabel("(Pseudo) first order rate constant s$^{-1}$", fontsize = 14)
plt.xlim(0,1)
plt.ylabel("Fate", fontsize = 14)
plt.tight_layout()

plt.show()

print("  ")
print(f"If H-shift rate is correct, [RO₂] =  {Hshift_RO2_ppb: .2f} ppb or {Hshift_RO2_mcm: .2e} molecules/cm³" "\n" f"and RO2 rate = {RO2RO2_rate_pred_2: .2e} cm³/molecules/s")
print("  ")
print(f"If HO₂ ammount is correct, [RO₂] =  {HO2_RO2_ppb: .2f} ppb or {HO2_RO2_mcm: .2f} molecules/cm³" "\n" f"and RO2 rate = {RO2RO2_rate_pred: .2e} cm³/molecules/s")
