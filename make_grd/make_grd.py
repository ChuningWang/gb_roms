# Use gridgen to generate a coarse, unfixed grid. Since this step is particularly slow, I split it into a
# separate script for efficiency purpose.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

import pyroms
import pyroms_toolbox
import bathy_smoother

from scipy.interpolate import griddata
import matplotlib.pyplot as plt

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_300m_orig'
grd_name = 'GlacierBay_300m_orig'
tag = ''
bathy_dir = in_dir + 'ARDEMv2.0.nc'
out_file = out_dir + 'grd/' + grd_name + '_grd' + tag + '.nc'

# ----------------------------------------------------------------------------------------------------------
# grid dimension
Lg = 250   # horizontal
Mg = 500  # vertical

# Lg = 250   # horizontal
# Mg = 500   # vertical

# ----------------------------------------------------------------------------------------------------------
# defind the boundary of mapping domain
lat_min = 57.
lat_max = 60.
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -138.
lon_max = -134.
lon_0 = 0.5 * (lon_min + lon_max)

# ----------------------------------------------------------------------------------------------------------
# These coords are handpicked using the Boundary Interactor.

lon_bry = np.array([    -137.35,    -136.30,    -135.00,    -136.05])
lat_bry = np.array([    59.05,      57.80,      58.05,      59.25  ])
beta = np.array([ 1., 1., 1., 1.])

# # other test cases
# ------------------------------------------------------------------------
# lon_bry = np.array([-137.347, -136.974, -137.136, -137.239,
#                     -137.275, -137.245, -137.112, -136.932,
#                     -136.739, -136.522, -136.133, -134.987,
#                     -134.963, -134.945, -134.963, -135.017,
#                     -135.102, -135.974])
# 
# lat_bry = np.array([ 59.051,  58.401,  58.313,  58.202,
#                      58.091,  57.992,  57.887,  57.816,
#                      57.778,  57.752,  57.752,  58.033,
#                      58.091,  58.151,  58.202,  58.240,
#                      58.272,  59.221])
# 
# beta = np.array([ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,
#                   1.,  0.,  0.,  0.,  0.,  0.,  1.])

# ------------------------------------------------------------------------
# lon_bry = np.array([-137.347, -136.974, -137.136, -137.239,
#                     -137.275, -137.245, -137.112, -136.932,
#                     -136.739, -136.522, -136.133, -136.000,
#                     -135.900, -135.500, -135.300, -134.987,
#                     -134.963, -134.945, -134.963, -135.017,
#                     -135.102, -135.974])

# lat_bry = np.array([ 59.051,  58.401,  58.313,  58.202,
#                      58.091,  57.992,  57.887,  57.816,
#                      57.778,  57.752,  57.752,  57.800,
#                      58.050,  58.000,  58.025,  58.033,
#                      58.091,  58.151,  58.202,  58.240,
#                      58.272,  59.221])

# beta = np.array([ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,
#                   0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  1.])

# ----------------------------------------------------------------------------------------------------------
# generate hgrid
bdryInteractor = 2
if bdryInteractor == 1:
    # use boundary interactor
    m = Basemap(projection='lcc', width = 12000000, height = 9000000,
                lat_1 = 30, lat_2 = 70, lat_0=lat_0, lon_0=lon_0,
                resolution='f')
    xp, yp = m(lon_bry, lat_bry)
    plt.figure()
    m.drawcoastlines()
    # plt.show()
    bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mg+3,Lg+3), proj=m)
    hgrd =bry.grd
elif bdryInteractor == 2:
    # use gridgen directly
    m = Basemap(projection='lcc', width = 8000000, height = 6000000,
                lat_1 = 30, lat_2 = 70, lat_0=lat_0, lon_0=lon_0,
                resolution='f')
    xp, yp = m(lon_bry, lat_bry)
    hgrd = pyroms.grid.Gridgen(lon_bry, lat_bry, beta, (Mg+3, Lg+3), proj=m)
else:
    # load grid that has been generated before
    hgrd = pyroms.grid.get_ROMS_hgrid(grd1)

print 'hgrid generated'

if (bdryInteractor == 1) | (bdryInteractor == 2):

    lonv, latv = m(hgrd.x_vert, hgrd.y_vert, inverse=True)
    hgrd = pyroms.grid.CGrid_geo(lonv, latv, m)

    # generate the mask
    for verts in m.coastsegs:
        # hgrd.mask_polygon(verts)
        if np.shape(verts)[0] == 2:
            verts.append(verts[0])

        hgrd.mask_polygon(verts)

print 'mask done'

# generate the bathy
fh = netCDF4.Dataset(bathy_dir, mode='r')
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
hmin = 10  # allow dry_wet
topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
h = griddata((lon.flatten(), lat.flatten()), topo.flatten(), (hgrd.lon_rho, hgrd.lat_rho), method='linear')

print 'griddata done...'

# save raw bathymetry
hraw = h.copy()

# smooth bathymetry
rx0_max = 0.35
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()
h = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()

# ensure that depth is always deeper than hmin
h = pyroms_toolbox.change(h, '<', hmin, hmin)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

hgrd.h = h

# ----------------------------------------------------------------------------------------------------------

# vertical coordinate
theta_b = 2.0
theta_s = 8.0
Tcline = 10
N = 40
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid file
pyroms.grid.write_ROMS_grid(grd, filename=out_file)

