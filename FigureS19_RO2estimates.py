# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:00:43 2026

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

mpl.rcParams.update({'font.size': 20})

titles = "0425"

OH_values = {
    "20220328": 6.9e6,
    "20220425": 6.4e6,
    "20220706": 9.1e6,
    "20220708": 1.1e7,
    "20220722": 7.7e6
}

RO2_0328file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220328_BenzylAlcohol/Python Analylsis/GECKO_03282022_RO2_SOA.csv"
RO2_0425file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/GECKO_04252022_RO2_SOA.csv"
RO2_0706file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220706_benzyl_alcohol/Python_analysis/GECKO_07062022_RO2_SOA.csv"
RO2_0708file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220708_benzyl_alcohol/Python_analysis/GECKO_07082022_RO2_SOA.csv"
RO2_0722file = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220722_benzyl_alcohol/Python_analysis/GECKO_07222022_RO2_SOA.csv"

RO2_0328 = pd.read_csv(RO2_0328file, header = 0)
RO2_0425 = pd.read_csv(RO2_0425file, header = 0)
RO2_0706 = pd.read_csv(RO2_0706file, header = 0)
RO2_0708 = pd.read_csv(RO2_0708file, header = 0)
RO2_0722 = pd.read_csv(RO2_0722file, header = 0)

RO2_0328_OH = RO2_0328['Time'] * OH_values['20220328']
RO2_0425_OH = RO2_0425['Time'] * OH_values['20220425']
RO2_0706_OH = RO2_0706['Time'] * OH_values['20220706']
RO2_0708_OH = RO2_0708['Time'] * OH_values['20220708']
RO2_0722_OH = RO2_0722['Time'] * OH_values['20220722']


fig, ax = plt.subplots(figsize=(10, 6), sharex=False)

ax.scatter(RO2_0328_OH, RO2_0328['RO2'], c='#DBE5FA', label = "expt 1", s = 50)
ax.scatter(RO2_0425_OH, RO2_0425['RO2'], c='#94AEF0', label = "expt 2", s = 50)
ax.scatter(RO2_0706_OH, RO2_0706['RO2'], c='#5F83E7', label = "expt 3", s = 50)
ax.scatter(RO2_0708_OH, RO2_0708['RO2'], c='#1D49C3', label = "expt 4", s = 50)
ax.scatter(RO2_0722_OH, RO2_0722['RO2'], c='#0F2563', label = "expt 5", s = 50)

ax.set_xlabel("OH Exposure ((molec cm$^{-3}$) s)")
ax.set_ylabel("RO$_2$ (molec cm$^{-3}$)")
ax.legend(ncol = 2)