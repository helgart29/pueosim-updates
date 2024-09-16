##need PyGDSM in order to run this file

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

from pygdsm import gsm16, GSMObserver16, base_skymodel, GlobalSkyModel16, GlobalSkyModel, GSMObserver
from base_observer import BaseObserver
from gsm08 import GlobalSkyModel, GSMObserver

# #first pt in ascent
# (latitude, longitude, elevation) = ('-77.762','168.888', 11516.1585)

#second pt in ascent
(latitude, longitude, elevation) = ('-77.924','169.572', 18482.9268)

#third pt in ascent
(latitude, longitude, elevation) = ('-78.031','169.314', 27539.6342)

#fourth pt in ascent
(latitude, longitude, elevation) = ('-78.052','168.044', 34710.6707)

#fifth pt in ascent
(latitude, longitude, elevation) = ('-78.168','165.726', 39271.9512)

#sixth pt
(latitude, longitude, elevation) = ('-80.15','131.395', 38936.28049)

ov = GSMObserver()
ov.lon = longitude
ov.lat = latitude
ov.elev = elevation

#start time
#ov.date = datetime(2016,12,2,13,11,0) #utc 

# #first pt in ascent
# ov.date = datetime(2016,12,2,14,6,42)

#second pt in ascent
ov.date = datetime(2016,12,2,15,11,0)

#third pt in ascent
ov.date = datetime(2016,12,2,16,11,0)

#fourth pt in ascent
ov.date = datetime(2016,12,2,17,11,0)

#fifth pt in ascent
ov.date = datetime(2016,12,2,18,11,0)

#sixth pt
ov.data = datetime(2016,12,3, 9,11,0)

ov.generate(50)
sky = ov.view_observed_gsm(logged=False)

#Sanity check of what map output actually is
hp.projview(sky, coord = ['G','C'], unit="T [K]",fontsize=dict(title='large',
                  xlabel='medium', ylabel='medium', 
                  cbar_label='small'),
                  graticule=True, graticule_labels=True, projection_type='mollweide',
                  title = f'Mollweide View from McMurdo Station at 50 MHz', xlabel='Right Ascension', ylabel='Declination')
plt.savefig(f'RaDec_pt1_at_60MHz.png')

#collecting all the pixels and their corresponding coordinates
npix = len(sky)
nside = hp.npix2nside(npix) #returning the resolution of the map
gl, gb = hp.pix2ang(nside, np.arange(npix), lonlat=True) 

#convert galactic coor output to equatorial coord in degrees
gal_sc = SkyCoord(gl, gb, frame = 'galactic', unit = ('deg', 'deg'))
eq_sc = gal_sc.transform_to('icrs')
ra = eq_sc.ra.deg
dec = eq_sc.dec.deg

valid_ra = np.array(ra)
valid_dec = np.array(dec)
valid_sky = np.array(sky)

data = {
    'RA': valid_ra,
    'Dec': valid_dec,
    'Sky': valid_sky
}
df = pd.DataFrame(data)
df.to_csv('RaDec_pt1_60MHz.csv', index=False)
