import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

import pyroms
import pyroms_toolbox
import bathy_smoother

from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import path

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB'
grd_name = 'GlacierBay'
tag = ''
bathy_dir = in_dir + 'Glacier_Bay_grd_kate.nc'
out_file = out_dir + 'grd/' + grd_name + '_grd' + tag + '.nc'

# ----------------------------------------------------------------------------------------------------------
# grid dimension
Lg = 500   # horizontal
Mg = 1000  # vertical

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

lon_bry = np.array([    -137.40,    -136.30,    -135.00,    -136.10])
lat_bry = np.array([    59.10,      57.80,      58.05,      59.35])
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
bdryInteractor = -1
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
    m = Basemap(projection='lcc', width = 12000000, height = 9000000,
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

# using masking GUI to change land mask
GUImsk = 2
if GUImsk == 1:
    m = Basemap(projection='lcc', llcrnrlon=lon_min, llcrnrlat=lat_min,
                urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
                resolution='f')
    plt.ion()
    plt.figure()
    m.drawcoastlines()
    plt.show()
    pyroms.grid.edit_mask_mesh(hgrd, proj=m)
    # or (this is faster 'I guess')
    # coast = pyroms.utility.get_coast_from_map(m)
    # pyroms.grid.edit_mask_mesh(hgrd, coast=coast)
elif GUImsk == 2:
    # laod from mask_change.txt
    msk_c = np.loadtxt('mask_change.txt')
    for i in range(len(msk_c)):
        hgrd.mask_rho[msk_c[i, 1], msk_c[i, 0]] = msk_c[i, 2]

print 'mask done...'

# generate the bathy
h = pyroms.grid.get_ROMS_vgrid(grd1).h

fh = netCDF4.Dataset(bathydir, mode='r')
topo = fh.variables['h'][:]
lon = fh.variables['lon_rho'][:]
lat = fh.variables['lat_rho'][:]
lon[lon>180] = lon[lon>180]-360  # lon from -180 to 180
fh.close()

# fix minimum depth
hmin = 10  # allow dry_wet
topo = pyroms_toolbox.change(topo, '<', hmin, hmin)

# interpolate Kate's grid onto this one
# find the grid boundary
xi, eta = topo.shape
p0 = []
idx0 = []
for i in range(eta-1):
    p0.append((lon[0, i], lat[0, i]))
    idx0.append((0, i))
for i in range(xi-1):
    p0.append((lon[i, -1], lat[i, -1]))
    idx0.append((i, eta-1))
for i in range(eta-1):
    p0.append((lon[-1, eta-i-1], lat[-1, eta-i-1]))
    idx0.append((xi-1, eta-i-1))
for i in range(xi-1):
    p0.append((lon[xi-i-1, 0], lat[xi-i-1, 0]))
    idx0.append((xi-i-1, 0))

p = path.Path(p0)

pc = p.contains_points(np.array([hgrd.lon_rho.flatten(), hgrd.lat_rho.flatten()]).T).reshape(h.shape)

# interpolate new bathymetry
h[pc] = griddata((lon.flatten(), lat.flatten()), topo.flatten(), (hgrd.lon_rho[pc], hgrd.lat_rho[pc]), method='linear')

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

# insure that depth is always deeper than hmin
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

# # longitude between 0 and 360
# grd.hgrid.lon_vert = grd.hgrid.lon_vert + 360
# grd.hgrid.lon_rho = grd.hgrid.lon_rho + 360
# grd.hgrid.lon_u = grd.hgrid.lon_u + 360
# grd.hgrid.lon_v = grd.hgrid.lon_v + 360
# grd.hgrid.lon_psi = grd.hgrid.lon_psi + 360

# smooth the grid, mask out small channels
h0 = grd.vgrid.h
msk0 = grd.hgrid.mask_rho

# ------------------------------------------------
# mask out small channels
msk0[:230, :75] = 0
msk0[:150, :150] = 0

# ------------------------------------------------
# shapiro filter
h0 = pyroms_toolbox.shapiro_filter.shapiro2(h0, 32)

# write grid file
pyroms.grid.write_ROMS_grid(grd, filename=out_file)

