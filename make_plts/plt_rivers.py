import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

if grd1=='GB_USGS':
    tag = 'GlacierBay_usgs'
elif grd1=='GB_lr':
    tag = 'GlacierBay_lr'

dtype = 1

grd = pyroms.grid.get_ROMS_grid(grd1)
lon_grd = grd.hgrid.lon_rho
lat_grd = grd.hgrid.lat_rho
msk = grd.hgrid.mask_rho

if dtype==1:
    # load remapped data
    fh = nc.Dataset(data_dir+tag+'_runoff_2008_Hill.nc', 'r')
    time = fh.variables['time'][:]
    lon = fh.variables['lon_rho'][:]
    lat = fh.variables['lat_rho'][:]
    r = fh.variables['Runoff'][:]
    fh.close
elif dtype==2:
    # load original data
    fh = nc.Dataset(data_dir+'gb_discharge.nc', 'r')
    time = fh.variables['t'][:]
    lon = fh.variables['lon'][:]
    lat = fh.variables['lat'][:]
    r = fh.variables['discharge'][:]
    coast = fh.variables['coast'][:]
    tmsk = (time>=39446) & (time <=39812)
    time = time[tmsk]
    r = r[tmsk, :, :]
    for i in range(len(time)):
        r[i, :, :] = np.ma.masked_where(coast.mask, r[i, :, :])
    fh.close()

r = np.ma.masked_where(r<0, r)
tp, xp, yp = r.shape

def plt_rivers(t, r, out_dir):
    plt.pcolormesh(lon, lat, r, cmap='Greens')
    plt.xlim(-137.5, -135)
    plt.ylim(58., 59.25)
    plt.clim(0, 10)
    plt.colorbar()
    plt.contour(lon_grd, lat_grd, msk, np.array([0.5, 0.5]), linewidths=0.05, colors='k')
    plt.savefig(out_dir+'figs/rivers/runoff_'+str(int(t))+'.png')
    plt.close()
    return None

from joblib import Parallel, delayed
Parallel(n_jobs=16)(delayed(plt_rivers)(time[i], r[i, :, :], out_dir) for i in range(len(time)))
# Parallel(n_jobs=16)(delayed(plt_rivers)(time[i], r[i, :, :], out_dir) for i in range(150, 165))
