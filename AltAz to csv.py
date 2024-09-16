import os
import numpy as np
import healpy as hp
import h5py
import pygdsm
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.table import Table
from datetime import datetime
from scipy.interpolate import griddata
from scipy.interpolate import interp2d, Rbf
import pandas as pandasForSortingCSV 

#ordering RaDec data in ascending Ra
RaDec = pandasForSortingCSV.read_csv("RaDec_pt1_50MHz.csv")
RaDec.sort_values(RaDec.columns[0], axis=0, inplace=True)

#converting data into arrays
ra = RaDec['RA'].to_numpy()
dec = RaDec['Dec'].to_numpy()
temp = RaDec['Sky'].to_numpy()

ra = ra[::75]
dec = dec[::75]
temp = temp[::75]

#finding the lowest horizon angle
horz_angle = np.degrees(-np.arccos((6378)/(6378+38.93628049)))
lowest_zen = 90 + abs(horz_angle)

# #ascent pt 2
# location = EarthLocation(lon = 169.572 * u.deg, lat = -77.924 * u.deg, height = 18482.9268 * u.m)

# #ascent pt 3
# location = EarthLocation(lon = 169.314 * u.deg, lat = -78.031 * u.deg, height = 27539.6342 * u.m)

# #ascent pt 4
# location = EarthLocation(lon = 168.044 * u.deg, lat = -78.052 * u.deg, height = 34710.6707 * u.m)

# #ascent pt 5
# location = EarthLocation(lon = 165.726 * u.deg, lat = -78.168 * u.deg, height = 39271.9512 * u.m)

#pt 6
location = EarthLocation(lon = 131.395 * u.deg, lat = -80.15 * u.deg, height = 38936.28049 * u.m)

obstime = datetime(2016,12,3,1,30,0) #utc 

#defining the frame that coord will be converted to
frame = AltAz(obstime=obstime, location=location)

#convert to AltAz
sc = SkyCoord(ra, dec, frame='icrs', unit="deg")
altaz_sc = sc.transform_to(frame)
az = altaz_sc.az.deg
alt = altaz_sc.alt.deg

#convert Alt to zen
zen = 90 - alt

#converting any below horizon points to temp of ice
coord = np.column_stack((zen, az, temp))
coord[np.where(zen>=lowest_zen),2]=270

#naming columns of coord
data = {
    'Zen': np.round(coord[:,0]),
    'Az': coord[:,1],
    'Sky': coord[:,2]
}
#putting AltAz data into csv file
df = pd.DataFrame(data)
df.to_csv('Avg_Elev_15hr_50MHz.csv', index=False)


# "sanity check polar plot"
# alt_grid, azi_grid = np.meshgrid(np.linspace(0,180,180), np.linspace(0,360,360))

# temp_grid = griddata((coord[:,0], coord[:,1]), coord[:,2], (alt_grid, azi_grid))
# # temp_grid = griddata((valid_alt, valid_az), valid_sky, (alt_grid, azi_grid)) #using AltAz coord and corresponding temps to interpolate at grid points

# #convert to radians because polar plot uses radians as input for theta value
# #don't convert altitude bc radius value of polar plot is viewed as a distance not angle
# azi_rad = np.deg2rad(azi_grid) 

# fig, ax = plt.subplots(subplot_kw = {'projection': 'polar'})

# mesh = plt.pcolormesh(azi_rad, alt_grid, temp_grid, cmap = 'viridis')

# ax.set_rlim(bottom=0, top=180)
# ax.set_rticks([90, 60, 30, 0])
# ax.tick_params(axis = 'y', colors='white')
# ax.set_theta_offset(np.pi / 2)
# cbar = plt.colorbar(mesh, ax=ax, pad = 0.2)
# cbar.set_label('Temperature (K)')
# ax.set_title('Observer Sky Map at 50 MHz')

# #markings for top=lowest_horz
# ax.text(0.59,-30, 'Alt (Degrees)')
# ax.text(3.11, -35, 'S')
# ax.text(1.58, -35, 'E')
# ax.text(2.25,-60,'Az (Degrees)')

# plt.savefig('Obs_pt1_0_MHz.png.png')

