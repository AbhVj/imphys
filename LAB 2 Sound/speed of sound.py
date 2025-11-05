import numpy as np 
import pandas as pd
from scipy import odr
import matplotlib.pyplot as plt
import math 


#Column Headers
dist = 'Distance'
phase = 'PD'
s_dev = 'SD'

data = pd.read_excel('C:/Users/firew/Documents/GitHub/imphys/LAB 2 Sound/Sound.xlsx').dropna() #ur directory, dropna to prevent NaN errors

print(data)
freq = 20000 #frequency is 20khz

d_error = 0.0005 #=- 0.0005 in the callipers, errors in the phase are given by the s.dev

def f(B,x): #Linear fit function
    return B[0]*x + B[1]

phase_error = data[s_dev].tolist() #ODR cant parse errors stored in a panda series, convert to list
x = data[dist]
y = data[phase]
print(x)
print(y)
model = odr.Model(f)
data2 = odr.RealData(x,y, sx = d_error, sy = phase_error)
fit = odr.ODR(data2, model, beta0 = [1,0]).run()
a,b = fit.beta #slope = a, b = intercept
a_err, b_err = fit.sd_beta #error of each a and b
print(f"Slope = {a:3f} +- {a_err:.3f}")
print(f"Intercept = {b:.3f} +- {b_err:.3f}")
plt.errorbar(x, y,xerr=d_error, yerr=phase_error , fmt ='o', label = 'data') 
plt.plot(x,f([a,b],x), label=f'y = {a:.2f}x + {b:.2f}' , color = 'blue')
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Phase Difference (rad)')
plt.title("Phase Difference vs Detector Distance")
plt.grid(True)
plt.savefig('C:/Users/firew/Documents/GitHub/imphys/LAB 2 Sound/Figures/Sound Graph')
plt.show
v = (2*math.pi/a) *freq #a = 2pi/lambda, v = f(lambda)
dv = v *(a_err/a)
print(f"The value of the speed of sound is {v} +- {dv}")
print((343.3-v)/dv)
