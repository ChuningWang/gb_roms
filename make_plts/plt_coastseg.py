# Plot Glacier Bay map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc
import time

# Read topography
fh = nc.Dataset('../data/ARDEMv2.0.nc', mode='r')
lon = fh.variables['lon'][:]
lon = lon-360
lat = fh.variables['lat'][:]
z = fh.variables['z'][:]
fh.close()

# Read grid
fh = nc.Dataset('../data/GB_grd.nc', mode='r')
lon_psi = fh.variables['lon_psi'][:]
lat_psi = fh.variables['lat_psi'][:]
h = fh.variables['h'][:]
msk = fh.variables['mask_psi'][:]
fh.close()

plt.close()
fig = plt.figure()

lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)
# lat_0 = 58.7

lon_min = -137.5
lon_max = -134.5
lon_0 = 0.5 * (lon_min + lon_max)
# lon_0 = -136.0

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='i')

m.drawcoastlines()
plt.show(block=False)

# for verts in m.coastsegs:
#     vertss = np.asarray(verts)
#     m.plot(vertss[:, 0], vertss[:, 1], 'b')
#     plt.draw()
#     time.sleep(0.5)


