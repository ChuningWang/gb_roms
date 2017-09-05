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
grd1 = 'GB_hr'
grd2 = 'GB_150m_orig'
grd_name = 'GlacierBay_hr'

hgrd = pyroms.grid.get_ROMS_hgrid(grd2)
# laod from mask_change.txt
msk_c = np.loadtxt('mask_change_GB_hr.txt')
for i in range(len(msk_c)):
    hgrd.mask_rho[int(msk_c[i, 1]), int(msk_c[i, 0])] = msk_c[i, 2]
print 'mask done...'

water = hgrd.mask_rho

# ------------------------------------------------------------------------
# mask out small channels
# Lisianski Inlet
water[:245, :75] = 0
water[:150, :150] = 0

# Charpentier Inlet
# water[585:623, 170:200] = 0

# Adams Inlet (part)
water[718:725, 420:430] = 0
water[660:685, 480:490] = 0

# Tenakee Inlet
water[:100, 190:225] = 0
water[:40, 225:250] = 0
water[:20, 250:260] = 0

# other places
water[380:420, :20] = 0
water[480:500, 370:390] = 0
water[420:440, 150:170] = 0
water[840:860, 250:260] = 0

# mask out adams inlet
# water[620:720, 400:] = 0

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
lon3 = fh.variables['lon'][:]
lat3 = fh.variables['lat'][:]
h3 = fh.variables['z'][:]
fh.close()

bdry = np.loadtxt('bdry_usgs.txt')
p0 = [(bdry[i, 0], bdry[i, 1]) for i in range(len(bdry[:, 0]))]
p = path.Path(p0)
pc = ~p.contains_points(np.array([lon3, lat3]).T) 

lon3 = lon3[pc]
lat3 = lat3[pc]
h3 = h3[pc]

fh = nc.Dataset(bathy_dir + 'bathy_usgs.nc', 'r')
lon2 = fh.variables['lon'][:][1::3, 1::3]
lat2 = fh.variables['lat'][:][1::3, 1::3]
h2 = fh.variables['z'][:][1::3, 1::3]
fh.close()

msk = ~h2.mask
lon2 = lon2[msk]
lat2 = lat2[msk]
h2 = h2[msk]

lon0 = np.concatenate((lon3, lon2))
lat0 = np.concatenate((lat3, lat2))
h0 = np.concatenate((h3, h2))

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
h1 = h[:, :300]
h1[h1<hmin0] = hmin0
h[:, :300] = h1

h1 = h[:380, :]
h1[h1<hmin0] = hmin0
h[:380, :] = h1

h1 = h[630:, :]
h1[h1<hmin0] = hmin0
h[630:, :] = h1

hmax0 = 5  # m
h1 = h[340:350, 290:310]
h1[h1<hmax0] = hmax0
h[340:350, 290:310] = h1

h1 = h[340:350, 330:340]
h1[h1<hmax0] = hmax0
h[340:350, 330:340] = h1

c0 = [(-135.81307137660747, 58.385557938193749),
      (-135.8021646671254, 58.38488902601744),
      (-135.79529747967371, 58.385557938193749),
      (-135.78964214883118, 58.387007247909068),
      (-135.78600657900381, 58.38578090891918),
      (-135.77348406070959, 58.38767616008537),
      (-135.75974968580624, 58.390574779516015),
      (-135.75571016377583, 58.3904632941533),
      (-135.7452074064968, 58.392358545319489),
      (-135.73874417124819, 58.392247059956773),
      (-135.73308884040563, 58.391578147780471),
      (-135.72662560515698, 58.391132206329601),
      (-135.71450703906581, 58.391466662417756),
      (-135.70763985161412, 58.391912603868619),
      (-135.69632918992903, 58.393361913583945),
      (-135.69188571569558, 58.393919340397531),
      (-135.64987468657947, 58.39860172563165),
      (-135.65108654318857, 58.400608462160555),
      (-135.66522487029496, 58.407409069286302),
      (-135.68097900621351, 58.408412437550751),
      (-135.70400428178678, 58.404287479130218),
      (-135.76298130343056, 58.396929445190892),
      (-135.81145556779529, 58.38767616008537),
      (-135.81185951999834, 58.385000511380163)]

p = path.Path(c0)
pc = p.contains_points(np.array([hgrd.lon_rho.flatten(), hgrd.lat_rho.flatten()]).T).reshape(h.shape)
pc = pc & (h>hmax0)
h[pc] = hmax0

# Secret Bay
hmin0 = 10
h1 = h[400:440, 300:320]
h1[h1<hmin0] = hmin0
h[400:440, 300:320] = h1

h1 = h[400:420, 320:325]
h1[h1<hmin0] = hmin0
h[400:420, 320:325] = h1

h1 = h[380:415, 320:340]
h1[h1<hmin0] = hmin0
h[380:415, 320:340] = h1

h1 = h[420:460, 350:360]
h1[h1<hmin0] = hmin0
h[420:460, 350:360] = h1

# adams inlet
hmin0 = 10  # m
h1 = h[660:720, 390:]
h1[h1<hmin0] = hmin0
h[660:720, 390:] = h1

h1 = h[680:700, 380:390]
h1[h1<hmin0] = hmin0
h[680:700, 380:390] = h1

# ------------------------------------------------------------------------
# deal with shallow water regions
def local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.3):

    h1 = h[xmin:xmax, ymin:ymax]
    msk1 = water[xmin:xmax, ymin:ymax]

    # smooth bathymetry
    RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h1, msk1)
    print 'Max Roughness value is: ', RoughMat.max()
    hs = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(msk1, h1, rx0_max)
    RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hs, msk1)
    print 'Max Roughness value is: ', RoughMat.max()

    h[xmin:xmax, ymin:ymax] = hs

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
