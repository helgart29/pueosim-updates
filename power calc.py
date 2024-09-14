import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import scipy
from scipy.interpolate import make_interp_spline, BSpline
import pandas as pandasForSortingCSV 
from csv import writer

"Plot antenna temp data"
#reading antenna gain file and converting to numpy arrays
data = pd.read_csv("Ant_Temp_pt1.csv")
freq = data['Frequency (MHz)'].to_numpy()
ant_temp = data['Antenna Temperature (K)'].to_numpy()

data = pd.read_csv("Ant_Temp_pt2.csv")
freq2 = data['Frequency (MHz)'].to_numpy()
ant_temp2 = data['Antenna Temperature (K)'].to_numpy()

data = pd.read_csv("Ant_Temp_pt3.csv")
freq3 = data['Frequency (MHz)'].to_numpy()
ant_temp3 = data['Antenna Temperature (K)'].to_numpy()

data = pd.read_csv("Ant_Temp_pt4.csv")
freq4 = data['Frequency (MHz)'].to_numpy()
ant_temp4 = data['Antenna Temperature (K)'].to_numpy()

data = pd.read_csv("Ant_Temp_pt5.csv")
freq5 = data['Frequency (MHz)'].to_numpy()
ant_temp5 = data['Antenna Temperature (K)'].to_numpy()

data = pd.read_csv("Ant_Temp_Avg.csv")
freq_avg = data['Frequency (MHz)'].to_numpy()
ant_temp_avg = data['Antenna Temperature (K)'].to_numpy()


#call the plot
# fig, ax = plt.subplots()
# ax.plot(freq, ant_temp, 'tab:purple')
# ax.set_xticks([50,100,200,300,400,500])
# ax.set_xlabel('Frequency [MHz]')
# ax.set_ylabel('Antenna Temperature [K])')
# ax.set_title('Antenna Temperature for the Fourth Point in the Ascent')
# plt.savefig('Antenna Temperature for the Fourth Point in the Ascent.png')

#combined plot
fig, ax = plt.subplots()
ax.plot(freq, ant_temp, 'cornflowerblue', label='Elevation = 11.5 km')
ax.plot(freq2, ant_temp2, 'mediumpurple', label='Elevation = 18.5 km')
ax.plot(freq3, ant_temp3, 'darkorchid', label='Elevation = 27.6 km')
ax.plot(freq4, ant_temp4, 'indigo', label='Elevation = 34.7 km')
ax.plot(freq5, ant_temp5, 'plum', label='Elevation = 39.3 km')
plt.axhline(y=135, color='midnightblue', linestyle='--')
#ax.text(350, 300, r'$T_{ice}$' ' / 2 = 135 K', fontsize=14, bbox=dict(facecolor='thistle', alpha=0.5))
plt.legend(fontsize=14)
ax.set_xticks([50,100,200,300,400,500])
plt.xticks(fontsize=13)
plt.yticks(fontsize=12.6)
ax.set_xlabel('Frequency [MHz]', fontsize=14.5)
ax.set_ylabel('Antenna Temperature [K]', fontsize=14.5)
ax.set_title('Expected ' r'$T_{A}$' ' for Flight Ascent', fontsize=18)
plt.savefig('Poster Ascent Antenna Temp.png', dpi=1000)

data = pd.read_csv("Ant_Temp_LST_Avg.csv")
time = data['Time (Hr)'].to_numpy()
ant_temp = data['Antenna Temperature (K)'].to_numpy()

fig, ax = plt.subplots()
ax.plot(time, ant_temp, 'tab:purple', label = 'Average Elevation = 38 km')
ax.set_xticks([0,5,10,15,20,23])
plt.legend(fontsize=14)
plt.xticks(fontsize=13)
plt.yticks(fontsize=12.7)
ax.set_xlabel('Local Sidereal Time [Hr]', fontsize=14.5)
ax.set_ylabel('Antenna Temperature [K]', fontsize=14.5)
#ax.set_title('Fig. 7. T r\N{SUBSCRIPT A} Generated at 100MHz for a Sidereal Day', fontsize=13)
ax.set_title(r'$T_{A}$' ' Generated at 100 MHz for a Sidereal Day', fontsize = 16)
plt.savefig('LST Antenna Temp.png', dpi=1000)

"Power and Voltage Calculation"
#power
k_B = 1.389649 * pow(10,-23) #Boltz constant [J/K]
delt_freq = 10 * pow(10,6) #10 MHz
power = k_B*ant_temp_avg*delt_freq

# coord = np.column_stack((freq, power))
# data = {
#     'Frequency (MHz)': coord[:,0],
#     'power (W)': coord[:,1],
# }

# df = pd.DataFrame(data)
# df.to_csv('pt_1_power.csv', index=False)

# #voltage as function of frequency
v_freq = np.sqrt(50 * k_B * ant_temp_avg * delt_freq)

coord2 = np.column_stack((freq_avg, v_freq))
data = {
    'Frequency (MHz)': coord2[:,0],
    'Voltage as func of frequency (V)': coord2[:,1],
}

df = pd.DataFrame(data)
df.to_csv('avg_voltage.csv', index=False)

data = pd.read_csv("pt_1_voltage.csv")
freq = data['Frequency (MHz)'].to_numpy()
volt1 = data['Voltage as func of frequency (V)'].to_numpy()

data = pd.read_csv("pt_2_voltage.csv")
freq2 = data['Frequency (MHz)'].to_numpy()
volt2 = data['Voltage as func of frequency (V)'].to_numpy()

data = pd.read_csv("pt_3_voltage.csv")
freq3 = data['Frequency (MHz)'].to_numpy()
volt3 = data['Voltage as func of frequency (V)'].to_numpy()

data = pd.read_csv("pt_4_voltage.csv")
freq4 = data['Frequency (MHz)'].to_numpy()
volt4 = data['Voltage as func of frequency (V)'].to_numpy()

data = pd.read_csv("pt_5_voltage.csv")
freq5 = data['Frequency (MHz)'].to_numpy()
volt5 = data['Voltage as func of frequency (V)'].to_numpy()

data = pd.read_csv("avg_voltage.csv")
freq6 = data['Frequency (MHz)'].to_numpy()
volt6 = data['Voltage as func of frequency (V)'].to_numpy()

# #plot voltage vs frequency
# fig, ax = plt.subplots()
# ax.plot(freq6, volt6*pow(10,6), 'tab:purple')
# ax.set_xticks([50,100,200,300,400,500])
# plt.xticks(fontsize=13)
# plt.yticks(fontsize=13)
# ax.set_xlabel('Frequency [MHz]', fontsize=14.5)
# ax.set_ylabel('Voltage [μV]', fontsize=14.5)
# ax.set_title('Fig. 8. Voltage as a Function of Frequency', fontsize=16)
# plt.savefig('Voltage for Avg Elev.png', dpi=1000)

# #overlaying voltage vs frequency
# fig, ax = plt.subplots()
# ax.plot(freq, volt1, 'teal', label='Ascent Pt 1')
# ax.plot(freq2, volt2, 'deepskyblue', label='Ascent Pt 2')
# ax.plot(freq3, volt3, 'lightseagreen', label='Ascent Pt 3')
# ax.plot(freq4, volt4, 'steelblue', label='Ascent Pt 4')
# ax.plot(freq5, volt5, 'dodgerblue', label='Ascent Pt 5')
# ax.plot(freq6, volt6, 'powderblue', label='Ascent Pt 6')
# plt.legend(fontsize=12)
# ax.set_xticks([50,100,200,300,400,500])
# plt.xticks(fontsize=13)
# plt.yticks(fontsize=13)
# ax.set_xlabel('Frequency [MHz]', fontsize=14.5)
# ax.set_ylabel('Voltage [V]', fontsize=14.5)
# ax.set_title('Voltage as Function of Frequency for a Nominal Flight', fontsize=16)
# plt.savefig('Ascent Volt.png')

data = pd.read_csv("Vrms for avg elev.csv")
time = data['Time (Hr)'].to_numpy()
v_rms = data['Vrms (μV)'].to_numpy()

spline = make_interp_spline(time, v_rms, k=3)
time_smooth = np.linspace(time.min(), time.max(), 300)
v_rms_smooth = spline(time_smooth)

fig, ax = plt.subplots()
#ax.plot(time, v_rms, label='original')
ax.plot(time_smooth, v_rms_smooth/(np.sqrt(2)), 'tab:purple')
plt.axhline(y=7.8, color='midnightblue', label = 'Before PyGDSM ' r'$V_{RMS}$' '=7.8μV')
#ax.annotate('Before PyGDSM ' r'$V_{RMS}$' '=7.8μV', xy=(0, 0.5), xytext=(10, 0.52), textcoords='offset points',
#            ha='right', fontsize=12, color='r')
# plt.legend(fontsize=12)
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.tight_layout(pad=2.1, w_pad=0.5, h_pad=1.0)
ax.set_xlabel('Local Sidereal Time [Hr]', fontsize=16)
ax.set_ylabel(r'$V_{RMS}$' ' [μV]', fontsize=16)
ax.set_title(r'$V_{RMS}$' ' for Half a Sidereal Day', fontsize=18)
plt.savefig('smoothed Vrms for Avg Elev.png', dpi = 1000)