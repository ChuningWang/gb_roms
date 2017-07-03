import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms
import glob
from gb_toolbox import gb_ctd
from matplotlib.mlab import griddata
from geopy.distance import vincenty

# load data
outputs_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/outputs/spinup/'
fig_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/figs/trans/spinup/'

flist = glob.glob(outputs_dir+'*.nc')
flist = flist[-1:]

depth = 200
dd = 10
zlev = 40
tindex = 0
var = 'salt'
uvar = 'u'
vvar = 'v'
clim = [28, 32]

# load geological coordinates
ctd = gb_ctd.rd_ctd('/Volumes/R1/scratch/chuning/gb_roms/data/ctd/ctd.nc')
lat_ctd = ctd['lat_stn']
lon_ctd = ctd['lon_stn']

station = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 21]
lat_ctd = lat_ctd[station]
lon_ctd = lon_ctd[station]

ct_tr = (len(station)-1)*dd
lat_t = np.zeros(ct_tr)
lon_t = np.zeros(ct_tr)

for i in range(len(station)-1):
    lat_t[i*dd:(i+1)*dd] = np.linspace(lat_ctd[i], lat_ctd[i+1], dd+1)[:-1]
    lon_t[i*dd:(i+1)*dd] = np.linspace(lon_ctd[i], lon_ctd[i+1], dd+1)[:-1]

grd = pyroms.grid.get_ROMS_grid('GB')
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

# interpolate topography
h = grd.vgrid.h
h_tr = griddata(lon.flatten(), lat.flatten(), h.flatten(), lon_t, lat_t, interp='linear').diagonal()

# get rho point depth
s_rho = grd.vgrid.s_rho
z_tr = np.dot(np.matrix(h_tr).T, np.matrix(s_rho))
z_tr = z_tr.T

# calculate distance
dis = np.zeros(lat_t.size)
for i in range(1, lat_t.size):
    dis[i] = vincenty(
                      (lat_t[i-1], lon_t[i-1]),
                      (lat_t[i], lon_t[i])
                     ).meters
dis = np.cumsum(dis)
dis = dis/1000  # [km]
dis = np.tile(dis, (zlev, 1))

# interpolate salinity
s_tr = np.zeros((zlev, ct_tr))

for fn in flist:
    tag = fn.split('/')[-1].split('.')[0]
    print 'processing ' + tag + ' ...'
    fh = nc.Dataset(fn)
    s = fh.variables['salt'][:].squeeze()
    fh.close()

    for i in range(40):
        s_tr[i, :] = griddata(lon.flatten(), lat.flatten(), s[i, :, :].squeeze().flatten(), lon_t, lat_t, interp='linear').diagonal()

plt.pcolor(dis, z_tr, s_tr)
