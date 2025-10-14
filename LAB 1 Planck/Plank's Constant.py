# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 10:00:19 2025

@author: firew
"""
import numpy as np
from scipy import odr 
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sci


#Column Headers
V_r = "V1"
I_r = "I1"
V_y = "V2"
I_y = "I2"
V_g = "V3"
I_g = "I3"
V_b = "V4"
I_b = "I4"
V_o = "V5"
I_o = "I5"
WL  = 'Wavelengths'
FWHM = 'FWHM'

#Error for I and V, as well as list of errors and intercepts
V_err = 0.0001
I_err = 0.0010
inter = []
errors=[]
data = pd.read_excel('planck.xlsx').dropna()
 # Requires that xlsx file is in same dir as code.
                                             # Dropna removes NaN rows which cause the ODR module to malfunction. this again took a while to figure out.

def solve(m,c): #solving for the x intercept
    return -c/m

def propagate(x,y,m,c): #calculating the absolute uncertainty in the x intercept 
    x_rel = x/m
    y_rel= y/c
    z= (np.sqrt((x_rel**2) + (y_rel**2)))*(-c/m)
    return z 
    
    
def f(B,x):
    return B[0]*x + B[1]


for i in range(1, 6): #probably a better way to generate each graph but im bad at coding so here we are.
    plt.clf()
    match i:
        case 1:
            x = data[V_r]
            y = data[I_r]
            text = 'Red'
        case 2:
            x = data[V_y]
            y = data[I_y]
            text = 'Yellow'
        case 3:
            x = data[V_g]
            y = data[I_g]
            text = 'Green'
        case 4:
            x = data[V_b]
            y = data[I_b]
            text = 'Blue'
        case 5:
            x = data[V_o]
            y = data[I_o]
            text = 'Orange'  

    model = odr.Model(f)
    data2 = odr.RealData(x,y,V_err,I_err)
    fit = odr.ODR(data2, model, beta0 = [1,0]).run()
    a,b = fit.beta #slope = a, b = intercept
    a_err, b_err = fit.sd_beta #error of each a and b
    inter.append(solve(a,b))
    errors.append(propagate(a_err, b_err, a, b)) #appending threshold voltages and their uncertainties to a list for use in next graph
    print(f"Slope = {a:3f} +- {a_err:.3f}")
    print(f"Intercept = {b:.3f} +- {b_err:.3f}")
    plt.errorbar(x, y,xerr=100*V_err, yerr=1000*I_err, fmt ='o', label = 'data') #error bars enlarged for visual purposes
    plt.plot(x,f([a,b],x), label=f'y = {a:.2f}x + {b:.2f}' , color = text.lower())
    plt.legend
    plt.xlabel('V')
    plt.ylabel('I')
    plt.title("Linear I-V for " + text + " L.E.D (100x horizontal 1000x vertical error bars)")
    plt.grid(True)
    plt.show
    plt.savefig(text)
    
plt.clf()
planckData = { #creating data frame for plotting the Vt vs 1/lambda graph
    "Threshold Voltages": inter,
    "Voltage Errors": errors,
    'Wavelengths': data[WL],
    'FWHM': data[FWHM]
    }
PlanckDF = pd.DataFrame(planckData)
wave_recip = 'Reciprocal of the Wavelength'
lambda_sd = 'Standard deviation of 1/lambda'
PlanckDF[lambda_sd] = (PlanckDF[FWHM]*10**-9)/2.355 * ((PlanckDF['Wavelengths']*(10**-9))**-2) # calculating uncertainty in 1/lambda according to formulas in back of lab manual
PlanckDF[wave_recip] = 1/(PlanckDF['Wavelengths']*(10**-9))
x = PlanckDF[wave_recip]
y = PlanckDF["Threshold Voltages"]
x_err = PlanckDF['Standard deviation of 1/lambda'].tolist() # converting to list from panda series is required for ODR to work, this took me too long to trouble shoot
y_err = PlanckDF["Voltage Errors"].tolist()
print(type(y_err))
model2 = odr.Model(f)
data3 = odr.RealData(x,y)
fit = odr.ODR(data3, model2, beta0 = [1,0]).run()
c,d = fit.beta #slope = c, d = intercept
c_err, d_err = fit.sd_beta #error of each c and d 
print(f"Slope = {c:20f} +- {c_err:.20f}")
print(f"Intercept = {d:.3f} +- {d_err:.3f}")
plt.errorbar(x, y,xerr=x_err, yerr=y_err, fmt ='o', label = 'data')
plt.plot(x,f([c,d],x), label=f'y = {c:.2f}x + {d:.2f}' , color = 'red')
plt.legend
plt.xlabel('1/λ  (1/m)')
plt.ylabel('Threshold Voltage (Volts)')
plt.title("Threshold Voltage vs 1/λ")
plt.grid(True)
plt.show
name = 'plankgraph'
plt.savefig('plankgraph')

print(f"h ={c*(sci.constants.e)/(sci.constants.c)} +- {c_err}")
