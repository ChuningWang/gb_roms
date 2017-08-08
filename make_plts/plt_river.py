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

tag = 'GlacierBay_usgs'

fh = nc.Dataset(data_dir+tag+'_runoff_2000_Hill.nc', 'r')
time = fh.variables['time'][:]
lon = fh.variables['lon_rho'][:]
lat = fh.variables['lat_rho'][:]
r = fh.variables['Runoff'][:]
fh.close

r = np.ma.masked_where(r<0, r)

for i in range(len(time)):
# for i in range(3):
    plt.pcolormesh(r[i, :, :], cmap='Greens')
    plt.xlim(0, 500)
    plt.ylim(0, 1000)
    plt.colorbar()
    plt.savefig(out_dir+'figs/rivers/runoff_'+str(int(time[i]))+'.png')
    plt.close()
