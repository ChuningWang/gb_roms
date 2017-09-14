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

# ------------------------------------------------------------------------
grd1 = 'GB_lr'
grd2 = 'GB_300m_orig'
grd_name = 'GlacierBay_lr'

hgrd = pyroms.grid.get_ROMS_hgrid(grd2)
water = hgrd.mask_rho

# ------------------------------------------------------------------------
# mask out small channels
water[190:260, :10] = 0
water[:100, :55] = 0
water[90:125, 10:25] = 0
water[:45, 95:120] = 0
water[:15, 115:135] = 0
water[330:360, 210:] = 0

# ------------------------------------------------------------------------
# msk_c = np.loadtxt('mask_change_GB_lr.txt')
# for i in range(len(msk_c)):
#     hgrd.mask_rho[int(msk_c[i, 1]), int(msk_c[i, 0])] = msk_c[i, 2]
# print 'mask done...'
 
# ------------------------------------------------------------------------
# defind the boundary of mapping domain
lat_min = 57.
lat_max = 60.
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -138.
lon_max = -134.
lon_0 = 0.5 * (lon_min + lon_max)

# ------------------------------------------------------------------------
# generate the base bathymetry
fh = nc.Dataset(bathy_dir + 'ARDEMv2.0.nc', mode='r')
h0 = fh.variables['z'][:]
lon0 = fh.variables['lon'][:] 
lon0[lon0>180] = lon0[lon0>180]-360  # lons from -180 to 180
lat0 = fh.variables['lat'][:] 
fh.close()

msk1 = (lon0>lon_min) & (lon0<lon_max)
msk2 = (lat0>lat_min) & (lat0<lat_max)

lon0 = lon0[msk1]
lat0 = lat0[msk2]
h0 = h0[msk2,:][:,msk1]

# depth positive
h0 = -h0

# fix minimum depth
hmin = -1  # allow dry_wet
hmax = 425
# h0 = pyroms_toolbox.change(h0, '<', hmin, hmin)

# interpolate new bathymetry
lon0, lat0 = np.meshgrid(lon0, lat0)
h = griddata((lon0.flatten(), lat0.flatten()), h0.flatten(), (hgrd.lon_rho, hgrd.lat_rho), method='linear')
print 'griddata from ARDEM done...'

# ------------------------------------------------------------------------
# fix bathymetry with USGS and NOAA data
fh = nc.Dataset(bathy_dir + 'bathy_noaa.nc', 'r')
lon1 = fh.variables['lon'][:]
lat1 = fh.variables['lat'][:]
h1 = fh.variables['z'][:]
fh.close()

bdry = np.loadtxt('bdry_usgs.txt')
p0 = [(bdry[i, 0], bdry[i, 1]) for i in range(len(bdry[:, 0]))]
p = path.Path(p0)
pc = ~p.contains_points(np.array([lon1, lat1]).T) 

lon1 = lon1[pc]
lat1 = lat1[pc]
h1 = h1[pc]

fh = nc.Dataset(bathy_dir + 'bathy_usgs.nc', 'r')
lon2 = fh.variables['lon'][:][1::3, 1::3]
lat2 = fh.variables['lat'][:][1::3, 1::3]
h2 = fh.variables['z'][:][1::3, 1::3]
fh.close()

msk = ~h2.mask
lon2 = lon2[msk]
lat2 = lat2[msk]
h2 = h2[msk]

lon0 = np.concatenate((lon1, lon2))
lat0 = np.concatenate((lat1, lat2))
h0 = np.concatenate((h1, h2))

# load grid boundary
bdry = np.loadtxt('bdry.txt')
x = bdry[:, 0]
y = bdry[:, 1]

p0 = [(x[i], y[i]) for i in range(len(x))]

p = path.Path(p0)
pc = p.contains_points(np.array([hgrd.lon_rho.flatten(), hgrd.lat_rho.flatten()]).T).reshape(h.shape)
 
# interpolate new bathymetry
h[(water==1) & pc] = griddata((lon0.flatten(), lat0.flatten()), h0.flatten(),
                               (hgrd.lon_rho[(water==1) & pc], hgrd.lat_rho[(water==1) & pc]), method='linear')

# copy raw bathymetry
hraw = h.copy()

# ------------------------------------------------------------------------
# locally constrain hmin at some location
hmin0 = 5  # m
h1 = h[:, :145]
h1[h1<hmin0] = hmin0
h[:, :145] = h1

h1 = h[:180, :]
h1[h1<hmin0] = hmin0
h[:180, :] = h1

h1 = h[182:188, 150:158]
h1[h1<hmin0] = hmin0
h[182:188, 150:158] = h1

h1 = h[327:, 162:200]
h1[h1<hmin0] = hmin0
h[327:, 162:200] = h1

h1 = h[250:, :190]
h1[h1<hmin0] = hmin0
h[250:, :190] = h1

# # constrain hmin at river discharge points
# hmin0 = 20
# h1 = h[780:805, 230:255]
# h1[h1<hmin0] = hmin0
# h[780:805, 230:255] = h1
# 
# h1 = h[920:940, 90:110]
# h1[h1<hmin0] = hmin0
# h[920:940, 90:110] = h1
# 
# h1 = h[920:940, 90:110]
# h1[h1<hmin0] = hmin0
# h[920:940, 90:110] = h1

# ------------------------------------------------------------------------
# use a 2D filter to smooth locally
from scipy.ndimage import uniform_filter
h1 = h[:250, :110]
h[:250, :110] = uniform_filter(h1, size=3)
# h1 = h[300:320, 300:380]
# h[300:320, 300:380] = uniform_filter(h1, size=5)
# h1 = h[260:360, 440:470]
# h[260:360, 440:470] = uniform_filter(h1, size=5)
# h1 = h[0:140, 220:380]
# h[0:140, 220:380] = uniform_filter(h1, size=5)
# h1 = h[790:820, 0:25]
# h[790:820, 0:25] = uniform_filter(h1, size=5)
# h1 = h[575:585, 225:235]
# h[575:585, 225:235] = uniform_filter(h1, size=3)
# h1 = h[910:930, 85:95]
# h[910:930, 85:95] = uniform_filter(h1, size=5)
# 
# # ------------------------------------------------------------------------
# # deal with shallow water regions
# def local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.3):
# 
#     h1 = h[xmin:xmax, ymin:ymax]
#     msk1 = water[xmin:xmax, ymin:ymax]
# 
#     # smooth bathymetry
#     # rx0_max = 0.3
#     RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h1, msk1)
#     print 'Max Roughness value is: ', RoughMat.max()
#     hs = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(msk1, h1, rx0_max)
#     RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hs, msk1)
#     print 'Max Roughness value is: ', RoughMat.max()
# 
#     h[xmin:xmax, ymin:ymax] = hs
# 
# 
# xmin = 630
# xmax = 685
# ymin = 380
# ymax = 500
# local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.25)
# 
# # xmin = 653
# # xmax = 663
# # ymin = 411
# # ymax = 425
# # local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.05)
# 
# xmin = 370
# xmax = 420
# ymin = 285
# ymax = 355
# local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.25)
# 
# xmin = 420
# xmax = 490
# ymin = 310
# ymax = 380
# local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.25)
# 
# xmin = 0
# xmax = 350
# ymin = 0
# ymax = 502
# local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.20)
# 
# xmin = 0
# xmax = 500
# ymin = 0
# ymax = 220
# local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.20)

# ------------------------------------------------------------------------
# shapiro filter
h = pyroms_toolbox.shapiro_filter.shapiro2(h, 32)

# final smooth
# smooth bathymetry
rx0_max = 0.35
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, water)
print 'Max Roughness value is: ', RoughMat.max()
h = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(water, h, rx0_max)
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, water)
print 'Max Roughness value is: ', RoughMat.max()

# insure that depth is always deeper than hmin
h = pyroms_toolbox.change(h, '<', hmin, hmin)
# constrain maximum depth
h = pyroms_toolbox.change(h, '>', hmax, hmax)
# fix depth of land points
h[water==0] = hmin

# # ------------------------------------------------------------------------
# # change bathymetry by hand using h_change.txt
# msk_c = np.loadtxt('h_change.txt')
# for i in range(len(msk_c)):
#     h[msk_c[i, 1], msk_c[i, 0]] = msk_c[i, 2]

# ------------------------------------------------------------------------
# redesign the vertical coordinate
theta_b = 2.0
theta_s = 8.0
Tcline = 10
N = 40
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

# ------------------------------------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
# write grid file
pyroms.grid.write_ROMS_grid(grd, filename = out_dir + 'grd/' + grd_name + '_grd.nc')