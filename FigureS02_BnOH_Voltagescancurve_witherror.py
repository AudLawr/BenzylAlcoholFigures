# -*- coding: utf-8 -*-
"""
Created on Thu May 16 08:39:08 2024

@author: audlaw
"""
#this script is intended to take KE50s and sensitivites that were calculated from a series of tank calibrations
#and use them to generate a cal curve. These data have been fit using both a hill fit equation and an ODR linear fit
#if you want to plot your own data on the curves, change the data_file_exp to your own file path 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.odr
from scipy.odr import Model, Data, ODR, RealData

#importing and extracting data
data_file = "voltagescancurve_input.csv"
data = pd.read_csv(data_file, header=0)

#extract variables
xerr=data['KE_50_cm_err']
yerr=data['Serr_t']
x = data['KE_50_cm']
y = data['S']

#putting together ODR linear fit
def f(B, x):
    return B[0]*x + B[1]

linear = Model(f)
mydata = RealData(x, y, sx=xerr, sy=yerr)
myodr = ODR(mydata, linear, beta0=[1., 2.])
myoutput = myodr.run()
#myoutput.pprint()   #can uncomment this line if you want it to print all the stats

#extract info from odr fit
slope, intercept = myoutput.beta
slope_err, intercept_err = myoutput.sd_beta

#calculating r-squared like a plebian
y_pred = slope * x + intercept # predicted y-values using odr fit
ss_res = np.sum((y - y_pred) ** 2) #residual sum of squares
ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
r_squared_odr = 1 - (ss_res / ss_tot) #r-squared

#plotting the ODR line w the scatter data and piecewise for the lowest and highest x values
x_intercept = -intercept / slope if slope != 0 else 0 #calc x-intercept
x_min, x_max = min(x), max(x)  #x-range for the regression bit
y_min = slope * x_min + intercept  #calculate y at the smallest x

piecewise_x = np.linspace(0, 2 * x_max, 300) # piecewise x values
piecewise_y = np.piecewise(
    piecewise_x,
    [piecewise_x < x_intercept, (piecewise_x >= x_intercept) & (piecewise_x <= x_max), piecewise_x > x_max],
    [0, lambda x: slope * x + intercept, lambda x: slope * x_max + intercept]) #piecewise y values

plt.scatter(x, y, color='#437C90', edgecolor='k', alpha=0.9, label='tank cal data') #plot the scatter plot
plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color='#437C90', alpha=0.5, label='error')

plt.plot(piecewise_x, piecewise_y, '--', color='grey', label='ODR', linewidth = 2.5) #plot the piecewise regression line

print("Optimized ODR parameters:")
print(f"y = ({slope:.3f} ± {slope_err:.3f})x + ({intercept:.3f} ± {intercept_err:.3f})")
print(f"R-squared: {r_squared_odr:.3f}")
print(" ")


#hill fit plotting
def func(x, a, b, c):
    return (a * (x ** c)) / ((b ** c) + (x ** c))

# Define the model for x
x_model = np.linspace(0.0, 0.6, 100)
popt, pcov = curve_fit(func, x, y, bounds=([0, 0, 0], [2000, 0.6, 10]))
plt.plot(x_model, (popt[0] * (x_model ** popt[2])) / ((popt[1] ** popt[2]) + (x_model ** popt[2])), 
         color='k', alpha = 0.7, linewidth =2.5, label = "Hill equation")
perr = np.sqrt(np.diag(pcov)) # Extract standard deviations from the covariance matrix to propagate error

# add compounds names
offset = 50
for i, txt in enumerate(data['compound']):
    plt.annotate(txt, (data['KE_50_cm'][i], data['S'][i] + offset), fontsize = 13)

    
plt.xlabel('KE50, cm', fontsize=15)
plt.ylabel('Sensitivity', fontsize=15)
plt.xlim([0, 0.6])
plt.ylim([0, 2500])
plt.tick_params(axis='both', which='major', labelsize=15)

print("Optimized Hill fit parameters:")
print(f"a: {popt[0]:.2f} ± {perr[0]:.2f}")
print(f"b: {popt[1]:.2f} ± {perr[1]:.2f}")
print(f"c: {popt[2]:.2f} ± {perr[2]:.2f}")

plt.legend(fontsize=13)
plt.show()


#%%
#this section allows you to see where your experimental points land on the sensitivity curves

data_file_exp = "C:/Users/audlaw/OneDrive - Colostate/Desktop/CSU_Research/DataAnalysis/SCENTS/benzyl_alcohol_AL/20220425_benzyl_alcohol/Python analysis/20220425_voltagescan_sensitivity_20260215.csv"
data_exp = pd.read_csv(data_file_exp, header=0)
fig, ax = plt.subplots(figsize = (10,7))
plt.plot(x_model, (popt[0] * (x_model ** popt[2])) / ((popt[1] ** popt[2]) + (x_model ** popt[2])), '-', color='k', alpha = 0.8, label = "Hill fit", linewidth = 2.5) #plotting hill fit curve
plt.scatter(data['KE_50_cm'], data['S'], color='#27327F', s = 60, edgecolor = "k", label = "Gas cylinder cal.") # plotting tank cal scatter data
plt.errorbar(data['KE_50_cm'], data['S'], xerr=xerr, yerr=yerr, fmt='o', color='#437C90', alpha=0.5, label='error')

# Annotate compounds
offset = 50

plt.plot(piecewise_x, piecewise_y, '--', color='k', alpha = 0.7, label='ODR-based linear fit', linewidth = 2.5) #plot the piecewise regression line

for i, x in enumerate(data_exp['KE_50_cm']):
    label = 'Exp. 2 KE$_{50,cm}$ values' if i == 0 else ""  # only label the first one
    plt.axvline(x=x, alpha=0.5, color='#E07808', linewidth=1.5, label=label)

# Annotate compounds
offset = 50
for i, txt in enumerate(data['compound']):
    x = data['KE_50_cm'][i]
    y = data['S'][i]

    if txt == "acrylonitrile":  # the one you want moved
        plt.annotate(txt, (x + 0.01, y - offset), fontsize=13)  # shift left
    else:
        plt.annotate(txt, (x, y + offset), fontsize=13)

eq_text = (rf"Optimized ODR Parameters" + "\n"
    rf"$Sens. = ({slope:.3f} \pm {slope_err:.3f})KE_{{50,cm}} + ({intercept:.3f} \pm {intercept_err:.3f})$" + "\n"
    rf"$R^2 = {r_squared_odr:.3f}$")
plt.text(
    0.02, 0.97,  # position in axes coords
    eq_text,
    transform=plt.gca().transAxes,
    fontsize=12,
    verticalalignment='top',
    bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'))


plt.xlabel('KE$_{50,cm}$ (eV)', fontsize=15)
plt.ylabel('Sensitivity (cps ppbv$^{-1}$)', fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.xlim([0, 0.6])
plt.ylim([0, 2500])
plt.legend(fontsize = 12)
plt.show()
