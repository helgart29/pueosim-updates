import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy import integrate
import scipy.integrate as integrate
from scipy.integrate import simpson, dblquad
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
import pandas as pandasForSortingCSV 
from csv import writer
import astropy.units as u
from astropy.time import Time
from datetime import datetime, timezone
from astropy.coordinates import EarthLocation, Angle, SkyCoord, AltAz
from datetime import timedelta

"Antenna Gain Data"
#reading antenna gain file
data = pd.read_csv("Copy_RG_65in_to600MHz_Standardard_Hport_Zport.csv")

#filtering so only get values for desired frequency
frequency = 0.50#GHz
filtered_df = data[ data ['Freq (GHz)'] == frequency]

#adding gain_phi and gain_theta
summed_gain = filtered_df[' Gain_phi (in linear unit) (this is copol)'] + filtered_df[' Gain_theta']
summed_gain = summed_gain.to_numpy()

#theta, phi, and summed gain to numpy array
zen = filtered_df[' Theta (zenith)'].to_numpy()
azi = filtered_df [' Phi (azimuth)'].to_numpy()
gain = np.column_stack((zen, azi, summed_gain))

#gain interpolate
phi_grid, theta_grid = np.meshgrid(np.arange(0,365,5), np.arange(0,185,5))
gain_grid = griddata((azi, zen), summed_gain, (phi_grid, theta_grid))

# df = pd.DataFrame(gain_grid)
# df.to_csv('gain pt 2.csv', index=False)

#interpolate gain plot 
# g_dBi = 10*np.log10(summed_gain) #convert from linear to dBi
# g_dBi_grid = griddata((azi,zen), g_dBi, (phi_grid, theta_grid))
# plt.figure(figsize=(10, 8))
# plt.contourf(phi_grid, theta_grid, g_dBi_grid, levels=100, cmap='viridis')
# plt.colorbar(label='Gain (dBi)')
# plt.xlabel('Azimuth (degrees)')
# plt.ylabel('Zenith (degrees)')
# plt.title('Interpolated Gain at 100MHz')
# plt.savefig('interp gain.png')

"AltAz data"
#import ZenAz and temp data; sort in ascending zenith values
ZenAz = pandasForSortingCSV.read_csv("Avg_Elev_12hr_500MHz.csv")
# ZenAz.sort_values(ZenAz.columns[0], axis=0, inplace=True) ##only sorts zenith column
ZenAz_sorted = ZenAz.sort_values(by=['Zen', 'Az']) #sorts zenith and then sorts azimuth accordingly

#put data into arrays
zen_T = ZenAz_sorted['Zen'].to_numpy()
az_T = ZenAz_sorted['Az'].to_numpy()
T_b = ZenAz_sorted['Sky'].to_numpy()

#temp interpolate
az_grid, zen_grid = np.meshgrid(np.arange(0,365,5), np.arange(0,185,5))
# T_b_grid = griddata((az_T, zen_T), T_b, (az_grid, zen_grid), method='nearest')

T_b_grid = griddata((az_T, zen_T), T_b, (az_grid, zen_grid), method='linear')
grid_fill = griddata((az_T, zen_T), T_b, (az_grid, zen_grid), method = 'nearest')
nan_indices = np.isnan(T_b_grid) 
T_b_grid[nan_indices] = grid_fill[nan_indices] #filling in any Nan spaces with nearest interpolation method

# df = pd.DataFrame(T_b_grid)
# df.to_csv('ice temp nearest.csv', index=False)


"Antenna Temp Integral"
#finding effective antenna aperture
c = 2.998 * pow(10,8)
freq = 500#MHz
Ae_num = c*c*gain_grid
Ae_den = freq*freq*4*math.pi
A_e = Ae_num / Ae_den

#calulate antenna temp
Effec_A = A_e
Bright_temp = T_b_grid

bounds_spacing = np.deg2rad(5)

#compute integral
def coef_find(theta, phi):
    coef_grid = gain_grid*Bright_temp
    coef = coef_grid[int(np.rad2deg(theta)/5), int(np.rad2deg(phi)/5)]
    return coef
T_A = dblquad(lambda theta, phi: coef_find(theta, phi)*np.sin(theta), 0, 2*np.pi, lambda theta:0, lambda theta: np.pi)

mult = 1/(4*np.pi)
ant_temp = T_A[0]*mult
print(ant_temp)

# #create csv file to store antenna temp output
# # only need to do for the first frequency and then just append future outputs
# data = {
#     'Frequency (MHz)': [freq],
#     'Antenna Temperature (K)': [ant_temp],
# }

# df = pd.DataFrame(data)
# df.to_csv('Ant_Temp_13hr.csv', index=False)

# append antenna temp output to csv file
freq_temp = [freq, ant_temp]
#hr_temp = [hour, ant_temp]
with open('Ant_Temp_12hr.csv', 'a', newline='') as f_object:
    writer_object = writer(f_object)
    writer_object.writerow(freq_temp)
    f_object.close