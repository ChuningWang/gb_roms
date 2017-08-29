'''
Calculate residual flow, interpolate and compare with station data.
'''

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp2d
import xarray as xr
import pyroms

def find_nearest(lat, lon, lat0, lon0):

in_dir = '/glade/p/work/chuning/gb_roms/tides/'

# load velocity data after processed with tide
mod = loadmat(in_dir + 'Tide_model.mat')
stn = loadmat(in_dir + 'Tide_stn.mat')

lat = stn['lat']
lon = stn['lon']
latm = mod['lat']
lonm = mod['lon']

ut = np.zeros((lat.shape))
vt = np.zeros((lat.shape))

i = 0
for j in range(len(lat)):

    u_mod = np.real(mod['Utide'][i, :, :])
    v_mod = np.imag(mod['Utide'][i, :, :])
    v_mod[np.isnan(u_mod)] = np.nan

    # f = interp2d(mod['lon'], mod['lat'], u)
