import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, cheby1
import math
import scipy.stats as stats
from scipy.interpolate import interp1d
import pandas as pd
from csv import writer

Re = 6399.0  # polar radius of curvature in km.
# Re = 6356.755 # polar radius in km
c_km = 299792.458  # speed of light in km/s
k_b = 1.38064852e-23  # boltzmann's constant Watts / Hz / K
Z_0 = 376.730  # the impedance of free-space in Ohms

Z_A = 50.0

# and that they see a 50 ohm load
Z_L = 50.0

# we currently assume that the antenna sees 50% sky and 50% ice
sky_frac = 0.5

# the reflection of the sky OFF the ice and up into the antenna beam
# reduces the brightness temperature of the ice (see ITU-R noise reference)
# this is an average scaling factor across the band to match the reflected
# sky temperature to a more detailed simulation here:
# https://www.phys.hawaii.edu/~prechelt/presentations/thermal_noise.pdf
CFr = 0.65

def butter_bandpass(lowcut, highcut, fs, order=2):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band') #butter
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=2):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def read_S11(freq, filename):

    freq=0.001*freq
    freq[freq>0.79]=0.79
    data=np.loadtxt(filename,skiprows=0,delimiter=',')
    r_s11 = data[:,1]
    i_s11 = data[:,2]
    s11_c=r_s11+i_s11*1j

    s11=20*np.log10(np.abs(s11_c))
    # plt.plot(data[:,0],s11)
    # plt.savefig('freq.png')
    trans=1-10**(0.1*s11)

    f=interp1d(data[:,0],trans,'linear')

    return f(freq)

def add_Gaussian_noise(h_phi,sigma):
    noise=stats.norm.rvs(scale=sigma,loc=0,size=h_phi.size)
    h_phi_noise=noise+h_phi
    return h_phi_noise

def add_Rayleigh_noise(h_phi, sigma):
    amplitude = stats.rayleigh.rvs(loc=0, scale=sigma, size=h_phi.size)
    phase = np.random.uniform(0, 2 * np.pi, size=h_phi.size)    
    noise = np.real(amplitude * np.exp(1j * phase))
    h_phi_noise=noise+h_phi
    return h_phi_noise

def add_Rician_noise(h_phi, sigma):
    noise = stats.rice.rvs(loc=h_phi, scale=sigma, size=h_phi.size)
    ###in work
    h_phi_noise=noise+h_phi
    return h_phi_noise
    
def noise_voltage(freqs: np.ndarray):
    """
    Return the noise voltage (in V) at given frequencies.

    Parameters
    ----------
    freqs: np.ndarray
        The frequencies to calculate at (in GHz).
    -------
    noise: np.ndarray
        The noise voltages at each frequency (in V).
    """
    #convert the freq to MHz
    freqs=freqs*1000 


    # and the system noise (from RF chain) temperature in Kelvin
    Tsys = get_Tsys(freqs)

    # combine to make the total noise temperature
    Tcomb = Tsys + noise_T(freqs)
 
    # get the bandwidth at each frequency step (Hz) The len(freqs) is for fft normalization factor 
    bw = 1e6 * (freqs[2] - freqs[1])* len(freqs)
    #print("bw=",bw)
    # and now convert this into a combined noise voltage
    Vn = np.sqrt(Tcomb * bw * k_b * Z_L) * pow(10,6)

    # and we are done
    return Vn

def noise_T(freqs: np.ndarray):
    """
    Return the noise Temperature from environment (sky+ice) (in V) at given frequencies.

    Parameters
    ----------
    freqs: np.ndarray
        The frequencies to calculate at (in MHz).
    -------
    noise_T: np.ndarray
        The Antenna noise Temperaure at each frequency (in V).
    """
    data = pd.read_csv("Ant_Temp_Avg.csv")
    frequencies = data['Frequency (MHz)'].to_numpy()
    ant_temp = data['Antenna Temperature (K)'].to_numpy()
    interp_temp = interp1d(frequencies, ant_temp, kind='linear', fill_value="extrapolate")
    Tcomb = interp_temp(freqs)
    Tcomb[Tcomb<0] = 0

    # filtered_df = data[data['Frequency (MHz)'].isin(freqs)]
    # Tcomb = filtered_df['Antenna Temperature (K)'].to_numpy()
    # fig, ax = plt.subplots()
    # ax.plot(freqs, Tcomb)
    # plt.savefig(f'interp1d.png')

    return Tcomb


def get_Tice(freqs: np.ndarray) -> np.ndarray:

    return 270.0 * np.ones_like(freqs)

def get_Tsys(freqs: np.ndarray, flight: int = 3) -> np.ndarray:
    """
    return RF chain noise 
    """
    return 100.0 * np.ones_like(freqs)


def get_freqs(samplingrate,sampling_points) -> np.ndarray: #sampling rate in GHz
    df = samplingrate/sampling_points # 3, 1.5 GHz/ #points
    return np.arange(0+df, samplingrate/2.0+df, df)

# def get_freqs() -> np.ndarray:
#     return np.arange(50, 510,10) #frequencies from 50-500MHz, step size of 10


def get_noise_spectra(freqs):
    
    # The sqrt(float(len(freqs))) is for fft normalization purpose
    Vn_spectrum = noise_voltage(freqs)*np.sqrt(float(len(freqs))) 
    
    # and sample the Rayleigh amplitudes
    amplitude = stats.rayleigh.rvs(loc=0, scale=Vn_spectrum, size=freqs.size)

    # and sample random phases
    phase = np.random.uniform(0, 2 * np.pi, size=freqs.size)
         
    return amplitude * np.exp(1j * phase)


def get_noise_t(samplingrate, time):
    freq = get_freqs(samplingrate, len(time))
    #freq = get_freqs()
    noise = get_noise_spectra(freq)
    noise = np.insert(noise,0,0)
    s = np.fft.irfft(noise)
    return s
    

def get_scaling_factor(h_phi, noise_RMS, SNR):
    peak2peak = np.max(h_phi) - np.min(h_phi)
    return SNR*noise_RMS/(peak2peak/2)

def add_noise(time, samplingrate, h_phi, scaling, filterRangelow, filterRangehi): # filterRange in GHz
    
    step_size = time[1] - time[0]
    
    extended_time = np.concatenate((time, time[-1] + np.arange(1, 200 + 1) * step_size))

    noise = get_noise_t(samplingrate, extended_time)
    noise_t_filter = butter_bandpass_filter(noise, filterRangelow, filterRangehi, 1.0 / (time[1] - time[0]), order=4)
    
    # skip first 200 events and adjust length to the len(h_phi)
    noise_t_filter = noise_t_filter[-1*len(h_phi):]
    #plt.plot(extended_time[200:],noise_t_filter)
    #plt.show()
    #plt.close()

    #print("add noise")    
    return noise_t_filter+scaling*h_phi, noise_t_filter

def interpolate_timesampling(ts,h_phi,rate):#sampling rate at GHz
        ts_sample=np.arange(0,100,1.0/rate) 
        h_phi_sample=np.interp(ts_sample,ts,h_phi)
        return ts_sample,h_phi_sample

def GetNoiseRMS(samplingrate, ts_sample, N, filterlow, filterhi): # filter in GHz
    assemble=[]
    for i in range(N):
        noise_t = get_noise_t(samplingrate,ts_sample)
        noise_t_filter = butter_bandpass_filter(noise_t, filterlow, filterhi, 1.0 / (ts_sample[1] - ts_sample[0]), order=4)

        assemble.extend(noise_t_filter)

    return np.std(assemble)

def GetNoiseRMS_sum(samplingrate, ts_sample, N, filterlow, filterhi, Nant, Digitize, num_bits, ceiling): # filter in GHz
    Sum=[]
    for n_ant in range(Nant):
        assemble=[]
        for i in range(N):
            noise_t = get_noise_t(samplingrate, ts_sample)
            noise_t_filter = butter_bandpass_filter(noise_t, filterlow, filterhi, 1.0 / (ts_sample[1] - ts_sample[0]), order=4)
            
            if Digitize == True:
                noise_t_filter = digitize_waveform(noise_t_filter, num_bits, ceiling) 

            assemble.extend(noise_t_filter)
        if n_ant==0:
            Sum = np.array(assemble)
        else:
            Sum = Sum + np.array(assemble)
    return np.std(Sum)


def digitize_waveform(voltage_values, num_bits, ceiling):

    max_voltage = ceiling
    min_voltage = -1*ceiling

    voltage_range = max_voltage - min_voltage

    # Define the number of bits for quantization
    max_levels = 2 ** (num_bits-1)

    # Calculate the quantization step size
    quantization_step = voltage_range / (max_levels*2)

    # Digitize the waveform into 5 bits
    digitized_waveform = np.zeros(len(voltage_values))
    for i, voltage in enumerate(voltage_values):
        # Quantize the voltage value
        quantized_value = int(voltage / quantization_step)

        # Clip the value to ensure it fits within the 5-bit range
        quantized_value = min(quantized_value, max_levels - 1)
        quantized_value = max(quantized_value, -max_levels)

        digitized_waveform[i] = quantized_value

    return digitized_waveform



if __name__=="__main__":
    '''
    An example to generate a time domain noise waveform, and Vrms
    '''
    # Sampling rate in GHz
    samplingrate=3
    # Time sample in ns
    ts_sample = np.arange(0,300,1/samplingrate) # time sample in ns
    
    # Generate noise time domain waveform
    
    noise_time = get_noise_t(samplingrate, ts_sample)
	
    # apply a filter
    
    noise_time = butter_bandpass_filter(noise_time, 0.05, 0.5, 1.0 / (ts_sample[1] - ts_sample[0]), order=4)

    plt.plot(noise_time, 'midnightblue')
    plt.xlabel('Time [ns]', fontsize=16.0)
    plt.ylabel('Voltage [μV]', fontsize=16.0)
    # plt.show()
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.tight_layout(pad=2.1, w_pad=0.5, h_pad=1.0)
    plt.title('Voltage as a Function of Time', fontsize = 18)
    plt.savefig('v(t) for avg elev.png', dpi=1000)

   
    # Calculate Vrms, repeat noise waveform N times to make it accurate
    N=100
    print("Vrms=",GetNoiseRMS(samplingrate, ts_sample, N, 0.05, 0.5)) 
    hour = 12

    # data = {
    # 'Time (Hr)': [hour],
    # 'Vrms (μV)': [GetNoiseRMS(samplingrate, ts_sample, N, 0.05, 0.5)],
    # }

    # df = pd.DataFrame(data)
    # df.to_csv('Vrms for avg elev.csv', index=False)

    # hr_volt = [hour, GetNoiseRMS(samplingrate, ts_sample, N, 0.05, 0.5)]
    # with open('Vrms for avg elev.csv', 'a', newline='') as f_object:
    #     writer_object = writer(f_object)
    #     writer_object.writerow(hr_volt)
    #     f_object.close
