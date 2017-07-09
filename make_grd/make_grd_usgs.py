import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4 as nc

import pyroms
import pyroms_toolbox
import bathy_smoother

from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import path

grd1 = 'GB3'
out_dir = '/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/'
grd_name = 'GlacierBay_usgs'
hgrd = pyroms.grid.get_ROMS_hgrid(grd1)
vgrd = pyroms.grid.get_ROMS_vgrid(grd1)

# generate the bathymetry
h = pyroms.grid.get_ROMS_vgrid(grd1).h
water = hgrd.mask_rho

# fix bathymetry with USGS and NOAA data
bathydir = '/Volumes/R1/scratch/chuning/data/bathy/'

fh = nc.Dataset(bathydir + 'bathy_noaa.nc', 'r')
lon1 = fh.variables['lon'][:]
lat1 = fh.variables['lat'][:]
z1 = fh.variables['z'][:]
fh.close()

fh = nc.Dataset(bathydir + 'bathy_usgs.nc', 'r')
lon2 = fh.variables['lon'][:][1::3, 1::3]
lat2 = fh.variables['lat'][:][1::3, 1::3]
z2 = fh.variables['z'][:][1::3, 1::3]
fh.close()

msk = ~z2.mask
lon2 = lon2[msk]
lat2 = lat2[msk]
z2 = z2[msk]

lon = np.concatenate((lon1, lon2))
lat = np.concatenate((lat1, lat2))
z = np.concatenate((z1, z2))

# load grid boundary
fh = nc.Dataset('bdry.nc')
x = fh.variables['lon'][:]
y = fh.variables['lat'][:]
fh.close()

p0 = [(x[i], y[i]) for i in range(len(x))]

p = path.Path(p0)
pc = p.contains_points(np.array([hgrd.lon_rho.flatten(), hgrd.lat_rho.flatten()]).T).reshape(h.shape)
 
# interpolate new bathymetry
hc = h.copy()
hc[(water==1) & pc] = griddata((lon.flatten(), lat.flatten()), z.flatten(),
                               (hgrd.lon_rho[(water==1) & pc], hgrd.lat_rho[(water==1) & pc]), method='linear')
# hgrd.h = hc
# 
# # mask out small channels
# msk0 = grd.hgrid.mask_rho
# msk0[:230, :75] = 0
# msk0[:150, :150] = 0
# grd.hgrid.mask_rho = msk0
# 
# # shapiro filter
# h0 = pyroms_toolbox.shapiro_filter.shapiro2(h0, 32)
# 
# # ----------------------------------------------------------------------------------------------------------
# grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
# # write grid file
# pyroms.grid.write_ROMS_grid(grd, filename=out_dir + grd_name + '_grd.nc')
# 
