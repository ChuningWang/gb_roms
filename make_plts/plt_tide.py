import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import pyroms
import pyroms_toolbox

fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/GB_tides_otps.nc')
name = fh.variables['tide_name'][:]
prd = fh.variables['tide_period'][:]
Eamp = fh.variables['tide_Eamp'][:]
Epha = fh.variables['tide_Ephase'][:]
Cmax = fh.variables['tide_Cmax'][:]
Cmin = fh.variables['tide_Cmin'][:]
Cang = fh.variables['tide_Cangle'][:]
Cpha = fh.variables['tide_Cphase'][:]
fh.close()

consts_num = prd.shape[0]

grd = pyroms.grid.get_ROMS_grid('GB')

lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

for i in range(consts_num):
    m = Basemap(projection='merc', llcrnrlon=lon.min(), llcrnrlat=lat.min(),
                urcrnrlon=lon.max(), urcrnrlat=lat.max(),
                resolution='f')

    m.drawcoastlines()

    xh, yh = m(lon, lat)

    m.contourf(xh, yh, Eamp[i, :, :].squeeze())
    m.contour(xh, yh, Epha[i, :, :].squeeze(), color='k')
    m.colorbar()
    plt.title(name[i])

    plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/tide' + str(i) + '.tiff',format='tiff')
    plt.close()

