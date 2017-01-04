import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import pyroms
import pyroms_toolbox

fh = nc.Dataset('../data/GB_tides_otps.nc')
Eamp = fh.variables['tide_Eamp'][:]
fh.close()

grd = pyroms.grid.get_ROMS_grid('GB')

lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

m = Basemap(projection='merc', llcrnrlon=lon.min(), llcrnrlat=lat.min(),
            urcrnrlon=lon.max(), urcrnrlat=lat.max(),
            resolution='f')

m.drawcoastlines()

xh, yh = m(lon, lat)
m.contourf(xh, yh, Eamp[0, :, :].squeeze())
m.colorbar()

plt.savefig('../figs/tide.eps',format='eps')
plt.close()

