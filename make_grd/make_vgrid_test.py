import pyroms
import pyroms_toolbox
import bathy_smoother

import netCDF4
import numpy as np

from scipy.interpolate import griddata
import matplotlib.pyplot as plt

lat_min = 57.
lat_max = 60.
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -138.
lon_max = -134.
lon_0 = 0.5 * (lon_min + lon_max)

# generate the bathy
bathydir = '../../data/ARDEMv2.0.nc'
fh = netCDF4.Dataset(bathydir, mode='r')
topo = fh.variables['z'][:]
lons = fh.variables['lon'][:] 
lons[lons>180] = lons[lons>180]-360  # lons from -180 to 180
lats = fh.variables['lat'][:] 
fh.close()

msk1 = (lons>lon_min) & (lons<lon_max)
msk2 = (lats>lat_min) & (lats<lat_max)

lons = lons[msk1]
lats = lats[msk2]
topo = topo[msk2,:][:,msk1]

# depth positive
topo = -topo

# fix minimum depth
hmin = 1.
topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

hgrd = pyroms.grid.get_ROMS_hgrid('GB')

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
h = griddata((lon.flatten(), lat.flatten()), topo.flatten(), (hgrd.lon_rho, hgrd.lat_rho), method='cubic')
hraw = h.copy()

# insure that depth is always deeper than hmin
h = pyroms_toolbox.change(h, '<', hmin, hmin)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

# ----------------------------------------------------------------------------------------------------------
# vertical grd
theta_s = 8.0
theta_b = 2.0
Tcline = 10
N = 40

vgrd = pyroms.vgrid.s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=hraw)


