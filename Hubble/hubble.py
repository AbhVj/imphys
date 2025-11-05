import pandas as pd 
import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

h_beta = 486.1 *(10**-9) # wavelength of H beta line 
c_l = 2.9979 *(10**8)    # speed of light

spec_data = pd.read_csv("C:/Users/firew/Documents/GitHub/imphys/Hubble/Data/SpectralData_Hbeta.csv", skiprows=4)        #skipping comments at the start of the csv and txt, filestring needs adjustment based on location of data
dist_data = pd.read_csv('C:/Users/firew/Documents/GitHub/imphys/Hubble/Data/dist_data.txt', delimiter='	', skiprows=5)
instrumentErrors = dist_data[dist_data.Instrument_response == "E"]['Observation_number'].values.tolist()                # List of all observation numbers that are instrumental errors

dist_data = dist_data[dist_data.Instrument_response != "E"]   #clearing the distance data of rows that contain instrumental errors 
for col in spec_data.columns:                                 #clearing the spectral data of instrumental errors by looping through all spectral data headers
    if int(float(col)) in instrumentErrors:                   # pandas stores duplicate headers as incremented by 0.1 : eg, the second instance of 47700 is denoted as 47700.1 , in order to check whether both columns are errors its necessary to convert -
        #print(col)                                           the observation number to an integer (to remove the increment). If the observation number of a column is in the errors list, it is removed. 
        spec_data.drop(col, axis = 1, inplace = True)         #dropping column if observation number corresponds to an error
        

def gauss(x,a, mu, sig, m, c ):             #gaussian function with linear background 
    gauss = a * np.exp(-(x - mu)**2 / (2 * sig**2))
    line = m * x + c
    return gauss + line

def guess(x,y):                             #returns a guess for the parameters of the noisy gaussian: estimates the linear background, subtracts it from the y data so that parameters can be estimnated.
    grad_guess = (y[-1]-y[0])/(x[-1]-x[0])  #gradient between end and beginning data points
    c_guess = y[0]-grad_guess*(x[0])         #solving for y intercept
    line = grad_guess*x + c_guess
    y_corr = y - line                       #y_corr is the data minus the estimated linear background
    a_guess = y_corr.max()                  # the maximum of the corrected data is guessed to be the amplitude, the x position of the largest value is guessed to be the mean
    mu_guess = x[np.where(y_corr == a_guess)[0][0]] 
    sig_guess = np.std(x)                    # guess for sigma is the std of the x data
    return [a_guess, mu_guess, sig_guess , grad_guess, c_guess]

def vel(f,sf):                              #function that takes in the frequency and its uncertainty and returns velocity and uncertainty in the velocity. (uses relativistic doppler equation)
    lam = (c_l)/f
    s_lam = (sf/f)  * lam                   # relative uncertainty is the same in wavelength and frequency
    v = ((lam/h_beta)**2-1)/(1+(lam/h_beta)**2) * c_l
    r = lam/h_beta                   
    sv = ((c_l*4*r*s_lam)/(1+r**2)**2)/h_beta        #error propagation formula determined through the partial derivative of v w.r.t lambda, r is used to simplify the expression
    return [v,sv]

def linear(x,m,):                            #linear model
    return m*x 

observ_nums = []
velocities = []
vel_err = []
distances=[]
print(len(spec_data.columns))
for i  in range(0,len(spec_data.columns), 2):                           #looping through the spectral data indexes, skipping two columns each iteration as each pair of columns yields the frequency and intensity data
    x = spec_data.iloc[:, i].to_numpy()             #grabbing x data from one column
    y = spec_data.iloc[:, i+1].to_numpy()           #grabbing y data from the adjacent column
    fit, fit_cov = curve_fit(gauss, x , y, guess(x,y))
    f = fit[1]
    sf = np.sqrt(np.diag(fit_cov)[1]) #scipy.optimise.curve_fit gives uncertainties as the sqrt of the diagonals of the covariance matrix
    observ_nums.append(spec_data.columns[i])
    v, sv = vel(f,sf)
    velocities.append(v)
    vel_err.append(sv)

hubbleData = {                                            #Setting up data frame for production of later plots and storing of velocities, uncertainties and the corresponding observation numbers
    "Observation Number": observ_nums,
    "Velocities": velocities, 
    "St.Deviation of Velocity": vel_err
    }
hubbleDF = pd.DataFrame(hubbleData)                       #converting above dictionary to a dataframe
for i in hubbleDF["Observation Number"]:                  #looping through all the observation numbers in new dataframe
    distances.append(dist_data[dist_data.Observation_number == int(i)]["Distance"].values[0])  #finding corresponding distance  (in dist_data) to ith observation numbeer, appending it to a list

hubbleDF['Distance'] = distances
x = hubbleDF['Distance']
y = hubbleDF['Velocities']/1000 #converting to km/s  (same with errors below)
y_error = hubbleDF['St.Deviation of Velocity']/1000 
print(y_error)

#two fitting aproaches can be taken, np.polyfit or curve_fit.
#Polyfit with degree 1 is currently used.

fit, lin_cov = curve_fit(linear, x,y, sigma=y_error) 
fit2, lin_cov2 = np.polyfit(x,y ,deg =1, w = 1/y_error, cov = True) 
m = fit2[0]
m_err = np.sqrt(np.diag(lin_cov2)[0])
print(f"Slope = {fit2[0]:20f} +- {m_err:.20f}")

#Plotting data and fits 
fig, ax = plt.subplots()
plt.errorbar(x, y, yerr=y_error, fmt ='o', label = 'Raw Data')
plt.plot(x,linear(x, m), label=f'V = {m:.2f}D ' , color = 'red')
plt.legend()
plt.xlabel('Distance  (Mpc)')
plt.ylabel('Redshift Inferred Velocity (km/s)')
plt.title("Redshifted Velocity vs Distance from Earth")
plt.grid(True)

#Code in the next 2 lines adapted from: https://matplotlib.org/stable/gallery/text_labels_and_annotations/placing_text_boxes.html
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, f"H_0 = ({m:.3f} ± {m_err:.3f}) km/s/Mpc", transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
plt.savefig('C:/Users/firew/Documents/GitHub/imphys/Hubble/Figures/hubble')
plt.show()

