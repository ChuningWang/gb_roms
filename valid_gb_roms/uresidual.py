'''
Calculate residual flow, interpolate and compare with station data.
'''

from glob import glob
import numpy as np
from scipy.io import loadmat
# import hdf5storage
from scipy.interpolate import interp2d
import pyroms
import netCDF4 as nc

def find_nearest(lat, lon, lat0, lon0):
    dis = (lat-lat0)**2+(lon-lon0)**2
    idx = np.where(dis==dis.min())
    xx = idx[0][0]
    yy = idx[1][0]
    return xx, yy

data_dir = '/glade/scratch/chuning/tmpdir_GB-TIDE/outputs/2008/'
in_dir = '/glade/p/work/chuning/gb_roms/tides/'
s0 = 'SEA1009'
grd0 = 'GB_lr'

# load grid info
grd = pyroms.grid.get_ROMS_grid(grd0)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho
h = grd.vgrid.h

# load velocity data after processed with tide
# mod = hdf5storage.loadmat(in_dir + 'Tide_model_lr.mat')
fin = nc.Dataset(in_dir + 'Tide_model_lr.nc', 'r')
time = fin.variables['time'][:]
utide = fin.variables['utide'][:]
vtide = fin.variables['vtide'][:]
ures = fin.variables['ures'][:]
vres = fin.variables['vres'][:]
fin.close()

stn = {}
# stn = loadmat(in_dir + 'Tide_stn.mat')
flist = glob(in_dir + 'SEA*.nc')
for fname in flist:
    stn_name = fname.split('/')[-1].split('.')[0]
    stn[stn_name] = {}
    fin = nc.Dataset(fname, 'r')
    stn[stn_name]['time'] = fin.variables['time'][:]
    stn[stn_name]['lat'] = fin.variables['lat'][:]
    stn[stn_name]['lon'] = fin.variables['lon'][:]
    stn[stn_name]['utide'] = fin.variables['utide'][:]
    stn[stn_name]['vtide'] = fin.variables['vtide'][:]
    stn[stn_name]['ures'] = fin.variables['ures'][:]
    stn[stn_name]['vres'] = fin.variables['vres'][:]
    fin.close()

lat0 = stn[s0]['lat']
lon0 = stn[s0]['lon']
xx, yy = find_nearest(lat, lon, lat0, lon0)

umod = utide[:, xx, yy]
vmod = vtide[:, xx, yy]
tstn = stn[s0]['time']
ustn = stn[s0]['utide']
vstn = stn[s0]['vtide']

# lat = stn['lat']
# lon = stn['lon']
# latm = mod['lat']
# lonm = mod['lon']
# 
# ut = np.zeros((lat.shape))
# vt = np.zeros((lat.shape))
# 
# i = 0
# for j in range(len(lat)):
# 
#     lat0, lon0 = lat[i, 0], lon[i, 0]
#     xx, yy = find_nearest(latm, lonm, lat0, lon0)
# 
#     u_mod = np.real(mod['Utide'][i, :, :])
#     v_mod = np.imag(mod['Utide'][i, :, :])
#     v_mod[np.isnan(u_mod)] = np.nan
# 
#     f = interp2d(lonm, latm, u_mod)
