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

import read_host_info
sv = read_host_info.read_host_info()
bathy_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_USGS'
grd_name = 'GlacierBay_usgs'
hgrd = pyroms.grid.get_ROMS_hgrid(grd1)
vgrd = pyroms.grid.get_ROMS_vgrid(grd1)

# generate the bathymetry
h = vgrd.h
water = hgrd.mask_rho

# fix bathymetry with USGS and NOAA data
fh = nc.Dataset(bathy_dir + 'bathy_noaa.nc', 'r')
lon1 = fh.variables['lon'][:]
lat1 = fh.variables['lat'][:]
z1 = fh.variables['z'][:]
fh.close()

fh = nc.Dataset(bathy_dir + 'bathy_usgs.nc', 'r')
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
 
# copy variables
hc = h.copy()
mskc = hgrd.mask_rho.copy()

# interpolate new bathymetry
hc[(water==1) & pc] = griddata((lon.flatten(), lat.flatten()), z.flatten(),
                               (hgrd.lon_rho[(water==1) & pc], hgrd.lat_rho[(water==1) & pc]), method='linear')

# ------------------------------------------------
# mask out small channels
mskc[:230, :75] = 0
mskc[:150, :150] = 0

# ------------------------------------------------
# use a 2D filter to smooth locally
from scipy.ndimage import uniform_filter
h1 = hc[:500, :220]
hc[:500, :220] = uniform_filter(h1, size=5)
h1 = hc[260:360, 440:470]
hc[260:360, 440:470] = uniform_filter(h1, size=5)
h1 = hc[0:140, 220:380]
hc[0:140, 220:380] = uniform_filter(h1, size=5)
h1 = hc[790:820, 0:25]
hc[790:820, 0:25] = uniform_filter(h1, size=5)
h1 = hc[575:585, 225:235]
hc[575:585, 225:235] = uniform_filter(h1, size=3)
h1 = hc[910:930, 85:95]
hc[910:930, 85:95] = uniform_filter(h1, size=5)

# ------------------------------------------------
# # locally smooth first blow-up point
# x0 = 339
# y0 = 112
# dy = 28
# dx = 28
# 
# h = hc[x0-dx:x0+dx, y0-dy:y0+dy]
# msk1 = mskc[x0-dx:x0+dx, y0-dy:y0+dy]
# 
# # smooth bathymetry
# rx0_max = 0.3
# RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, msk1)
# print 'Max Roughness value is: ', RoughMat.max()
# hs = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(msk1, h, rx0_max)
# RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hs, msk1)
# print 'Max Roughness value is: ', RoughMat.max()
# 
# hc[x0-dx:x0+dx, y0-dy:y0+dy] = hs

# ------------------------------------------------
# final smooth
# smooth bathymetry
rx0_max = 0.35
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hc, mskc)
print 'Max Roughness value is: ', RoughMat.max()
hc = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(mskc, hc, rx0_max)
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hc, mskc)
print 'Max Roughness value is: ', RoughMat.max()

# shapiro filter
hc = pyroms_toolbox.shapiro_filter.shapiro2(hc, 32)

# insure that depth is always deeper than hmin
hc = pyroms_toolbox.change(hc, '<', vgrd.hmin, vgrd.hmin)

vgrd.h = hc
hgrd.mask_rho = mskc

# ----------------------------------------------------------------------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
# write grid file
pyroms.grid.write_ROMS_grid(grd, filename = out_dir + 'grd/' + grd_name + '_grd.nc')

