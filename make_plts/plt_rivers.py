import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

# grd1 = 'GB_USGS'
# # Read grid
# grd = pyroms.grid.get_ROMS_grid(grd1)
# lon = grd.hgrid.lon_rho
# lat = grd.hgrid.lat_rho
# z = grd.vgrid.h
# msk = grd.hgrid.mask_rho

grd1 = 'GB_USGS'
tag = 'GlacierBay_usgs'

grd = pyroms.grid.get_ROMS_grid(grd1)
msk = grd.hgrid.mask_rho

fh = nc.Dataset(data_dir+tag+'_runoff_2008_Hill.nc', 'r')
time = fh.variables['time'][:]
lon = fh.variables['lon_rho'][:]
lat = fh.variables['lat_rho'][:]
r = fh.variables['Runoff'][:]
fh.close

r = np.ma.masked_where(r<0, r)

def plt_rivers(t, r, out_dir):
    plt.pcolormesh(r, cmap='Greens')
    plt.xlim(0, 502)
    plt.ylim(0, 1002)
    plt.clim(0, 2)
    plt.colorbar()
    plt.contour(msk, np.array([0.5, 0.5]), linewidths=0.05, colors='k')
    plt.savefig(out_dir+'figs/rivers/runoff_'+str(int(t))+'.png')
    plt.close()
    return None

# for i in range(len(time)):
# for i in range(3):

from joblib import Parallel, delayed
# Parallel(n_jobs=16)(delayed(plt_rivers)(time[i], r[i, :, :], out_dir) for i in range(len(time)))
Parallel(n_jobs=16)(delayed(plt_rivers)(time[i], r[i, :, :], out_dir) for i in range(150, 165))
