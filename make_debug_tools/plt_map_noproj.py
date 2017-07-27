import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_USGS'
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
h = grd.vgrid.h
msk = grd.hgrid.mask_rho

plt.pcolormesh(h, cmap='Greens')
plt.clim(0, 40)
plt.colorbar()
plt.contour(msk, np.array([0.5, 0.5]), linewidths=0.05, colors='k')
plt.show()
