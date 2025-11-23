import numpy as np 
import pandas as pd
from scipy import odr
import matplotlib.pyplot as plt
import math 


#Column Headers
P_p= 'P_p'
P_s = 'P_s'



data = pd.read_excel('C:/Users/firew/Documents/GitHub/imphys/Faraday/Data/Faraday.xlsx').dropna() #ur directory, dropna to prevent NaN errors




def f(B,x): #Linear fit function
    return B[0]*x + B[1]
sP_p = data['sP_p'].tolist()
sP_s = data['sP_s'].tolist()

x = data[P_s]
y = data[P_p]
print(x)
print(y)
model = odr.Model(f)
data2 = odr.RealData(x,y, sx = sP_s, sy = sP_p)
fit = odr.ODR(data2, model, beta0 = [1,0]).run()
a,b = fit.beta #slope = a, b = intercept
a_err, b_err = fit.sd_beta #error of each a and b
print(f"Slope = {a:3f} +- {a_err:.3f}")
print(f"Intercept = {b:.3f} +- {b_err:.3f}")
plt.errorbar(x, y,xerr=sP_s, yerr=sP_p , fmt ='o', label = 'data') 
plt.plot(x,f([a,b],x), label=f'y = {a:.2f}x + {b:.2f}' , color = 'blue')
plt.legend()
plt.xlabel('Power in the Secondary Coil (W)')
plt.ylabel('Power in the Primary Coil (W)')
plt.title("Power in the Primary Coil  vs Power in the Secondary Coil (of a transformer) ")
plt.grid(True)
plt.savefig('C:/Users/firew/Documents/GitHub/imphys/LAB 2 Sound/Figures/Sound Graph')
eff = 1/a
d_eff = b_err/b *eff
print(f"The value of the speed of sound is {eff} +- {d_eff}")
plt.show()


