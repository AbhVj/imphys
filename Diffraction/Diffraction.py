import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
#everythings in mm
lam = 670e-6          # wavelength in metres (670 nm)
dlam = 1e-6           # wavelength uncertainty (±1 nm)

f = 500             # focal length in metres (500 mm)
sf = 0.5          # focal length uncertainty (example ±5 mm)

def linear(x,m,b):                            #linear model
    return m*x +b

fig, ax = plt.subplots()
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

df = pd.read_csv('C:/Users/firew/Documents/GitHub/imphys/Diffraction/Data/double 500 polarised 2.csv') # your directory
x = df['Distance_(mm)'] -3.5# the minus value depends on centre, best to analyze the particular image again to get the exact peak intensity
# purpose of minus: shift x to align m=0 at x=0
y = df['Gray_Value']
df2 = pd.read_excel('C:/Users/firew/Documents/GitHub/imphys/Diffraction/Data/Double.xlsx')
x2 =df2['Min']
y2 =df2['X']


s= np.ones_like(y2)/0.0502
fit2, lin_cov2 = np.polyfit(x2,y2 ,deg =1, w=s, cov = True) 
g = fit2[0]
dg= np.sqrt(np.diag(lin_cov2)[0])
d1 = lam * f / g
frac_unc = np.sqrt((dlam/lam)**2 + (sf/f)**2 + (dg/g)**2)
sd = d1* frac_unc
print(g, dg)
print(f"Slit separation d = ({d1:.4f} ± {sd:.4f}) mm")
print(f"Fractional uncertainty = {frac_unc*100:.2f}%")


plt.plot(x2,g*x2 + fit2[1], label ='fit')
plt.ylabel('Distance (mm)')
plt.xlabel('Minima Number')
plt.title('Distance of Minima (mm) vs Minima Number')
plt.errorbar(x2, y2, yerr=0.0502 , fmt ='o', label = 'data') 
ax.text(0.05, 0.95, f"d = ({d1:.3e} ± {sd:.3e})mm", transform=ax.transAxes, fontsize=13,verticalalignment='top', bbox=props)
plt.legend()
plt.grid(True)
plt.savefig('C:/Users/firew/Documents/GitHub/imphys/Diffraction/Figures/DSPD', dpi=1000)
plt.show()
plt.clf()
# fit function
# from csv file, x=5.1562 y=44.2716 is the point with max I --> i.e. I_0 = I_max / 4 = 11.0679
l = 670*10**(-6) # lambda in mm
f = 500 # focal length used 160mm
I_0 = 15
def function(x,a,d,c):
    return 4*I_0* (np.sin(np.pi*a*x/(l*f)) / (np.pi*a*x/(l*f)))**2 * (np.cos(np.pi*d*x/(l*f)))**2 + c

#initial guesses (in mm)
guess_a = 0.05
guess_d = 0.24
guess_c = 1
guess = [guess_a, guess_d, guess_c]

# Set reasonable bounds for the parameters (a, d) in mm
#bounds = ([0.001, 0.001], [1.0, 1.0]) # a and d must be positive, and reasonable maximum is 1mm

#fit
fit, fit_cov = curve_fit(function, x, y, p0=guess, maxfev=5000, sigma=0.005)

sf = np.sqrt(np.diag(fit_cov)[1])   # uncertainty in d
sa = np.sqrt(np.diag(fit_cov)[0])   # uncertainty in a

d = fit[1]
a = fit[0]

print('The slit sep is: ', d, '+/-', sf)
print('The width is :', a, '+/-', sa)

data_fit = function(x, *fit)

plt.plot(x, data_fit, label='fit', linewidth=2, zorder=3)

plt.errorbar(
    x, y,
    xerr=0.005,
    fmt='o',
    color='gray',
    alpha=0.6,
    markersize=4,
    label='data',
    zorder=2
)

plt.xlabel('Distance (mm)')
plt.ylabel('Intensity (arbitrary units)')
plt.title('Intensity profile for Double Slit System')

# ---- Round uncertainties to 1 significant figure ----
sf_1sf = float(f"{sf:.1g}")
sa_1sf = float(f"{sa:.1g}")

# Determine decimal places needed
d_dp = max(0, -int(np.floor(np.log10(sf_1sf))))
a_dp = max(0, -int(np.floor(np.log10(sa_1sf))))

# ---- Add text box with results ----
textstr = (
    rf"$d = ({d:.{d_dp}f} \pm {sf_1sf:.{d_dp}f})\ \mathrm{{mm}}$" "\n"
    rf"$a = ({a:.{a_dp}f} \pm {sa_1sf:.{a_dp}f})\ \mathrm{{mm}}$"
)

plt.text(
    0.02, 0.95, textstr,
    transform=plt.gca().transAxes,
    fontsize=11,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.85)
)
plt.grid(True)
plt.legend()
plt.savefig('C:/Users/firew/Documents/GitHub/imphys/Diffraction/Figures/DS500P')
plt.show()