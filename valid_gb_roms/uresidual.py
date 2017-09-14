'''
Calculate residual flow, interpolate and compare with station data.
'''

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp2d
import pyroms
import xarray as xr

def find_nearest(lat, lon, lat0, lon0):
    dis = (lat-lat0)**2+(lon-lon0)**2
    idx = np.where(dis==dis.min())
    xx = idx[0][0]
    yy = idx[1][0]
    return xx, yy

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

    lat0, lon0 = lat[i, 0], lon[i, 0]
    xx, yy = find_nearest(latm, lonm, lat0, lon0)

    u_mod = np.real(mod['Utide'][i, :, :])
    v_mod = np.imag(mod['Utide'][i, :, :])
    v_mod[np.isnan(u_mod)] = np.nan

    f = interp2d(lonm, latm, u_mod)
