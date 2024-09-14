import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy.interpolate import interp1d
import os


def read_gain_Interpolate(freq, theta, phi, filename):
    data = np.loadtxt(filename,delimiter=',',skiprows=1)
    mask = (data[:,1] < theta+1) & (data[:,1] > theta-1) & (data[:,2] < phi+1) & (data[:,2] > phi-1)
    data_theta_phi = data[mask]
    Freq = data_theta_phi[:,0]
    G_phi = data_theta_phi[:,3]
    G_theta = data_theta_phi[:,4]

    f_phi = interp1d(Freq, G_phi,'linear')
    f_theta = interp1d(Freq, G_theta,'linear')

    return f_phi(freq), f_theta(freq)

def read_Hpol_Hgain(freq, theta, phi, filenameZ):
    Zport = read_gain_Interpolate(freq, theta, phi, filenameZ)
    return 10*np.log10(Zport[0])

def read_Hpol_Vgain(freq, theta, phi, filenameZ):
    Zport = read_gain_Interpolate(freq, theta, phi, filenameZ)
    return 10*np.log10(Zport[1])



if __name__=="__main__":
    freq=0.105 #GHz
    theta=90 # deg
    phi=0  #deg
    Colpol_gain = read_Hpol_Hgain(freq,theta,phi, "Copy_RG_65in_to600MHz_Standardard_Hport_Zport.csv") # in dBi
    Xpol_gain = read_Hpol_Vgain(freq,theta,phi, "Copy_RG_65in_to600MHz_Standardard_Hport_Zport.csv") # in dBi
    print("Colpol G=", Colpol_gain,'dBi', "Xpol gain=", Xpol_gain,"dBi")