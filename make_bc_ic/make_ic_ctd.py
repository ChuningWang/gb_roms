'''
The goal of this script is to reconstruct monthly T and S at all depth as initial condition
for Glacier Bay ROMS model. The purpose of using CTD instead of SODA data is that SODA data is
Too coarse for detailed mapping of T and S inside the bay.
'''

import sys
import numpy as np
import pyroms
import pyroms_toolbox
from scipy.interpolate import interp1d
from  scipy import ndimage

import netCDF4 as nc

from ocean_toolbox import ctd

def floodFill(c, r, mask):
    '''
    produce a mask array that identifies the flood of ctd stations.
    '''

    # number of seeds
    seednum = c.size
    # cells already filled
    filled = set()
    # cells to fill
    fill = [set() for _ in range(seednum)] 
    for i in range(seednum):
        fill[i].add((c[i], r[i]))

    eta = mask.shape[0]
    xi = mask.shape[1]

    # Our output inundation array
    flood = np.ones_like(mask, dtype=np.int8)*-1
    # Loop through and modify the cells which
    # need to be checked.

    while any(fill):
        fill_p1 = [set() for _ in range(seednum)]

        # loop through the stations
        for i in range(seednum):

            while fill[i]:
                # Grab a cell
                x, y = fill[i].pop()
                if x == eta or y == xi or x < 0 or y < 0:
                    # Boundary, don't fill
                    continue
                if mask[int(x)][int(y)] == 1:
                    # Do fill
                    flood[int(x)][int(y)] = i
                    filled.add((x, y))
                    # Check neighbors for 1 values
                    west = (x-1, y)
                    east = (x+1, y)
                    north = (x, y-1)
                    south = (x, y+1)
                    if west not in filled:
                        fill_p1[i].add(west)
                    if east not in filled:
                        fill_p1[i].add(east)
                    if north not in filled:
                        fill_p1[i].add(north)
                    if south not in filled:
                        fill_p1[i].add(south)

        fill = fill_p1[:]

    return flood

def nan_gaussian_filter(U, sigma):
    '''
    http://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    '''

    V = U.copy()
    msk = U!=U
    V[msk] = 0
    VV = ndimage.gaussian_filter(V, sigma = sigma)

    W = 0*U.copy()+1
    W[msk] = 0
    WW = ndimage.gaussian_filter(W, sigma = sigma)

    Z=VV/WW
    Z[msk] = np.NaN
    return Z

# ----------------------------------------------------------------
# preparation
# output dir
import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
dst_dir = sv['out_dir']
soda_dir = sv['soda_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

month = 4
src_grd_name = '2008_04_15_CTD_floodFill'

# Load target grid and land mask
grd = pyroms.grid.get_ROMS_grid(grd1)
msk = grd.hgrid.mask
h = grd.vgrid.h
Cs_r = grd.vgrid.Cs_r
zlev = grd.vgrid.N
lat_grd = grd.hgrid.lat_rho
lon_grd = grd.hgrid.lon_rho
lat_grd[msk==0] = np.NaN
lon_grd[msk==0] = np.NaN

eta, xi = msk.shape

# ----------------------------------------------------------------
# Load CTD climatology
info = {'data_dir': in_dir + 'ctd_raw/',
        'file_dir': in_dir,
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
       }
gb_ctd = ctd.ctd(info)
gb_ctd()

# Extract data
stn = gb_ctd.climatology['station']
lat = gb_ctd.climatology['lat']
lon = gb_ctd.climatology['lon']
z = gb_ctd.climatology['z']
temp = gb_ctd.climatology['temp']
salt = gb_ctd.climatology['salt']

# replicate station 00 to station 24
stn = np.append(stn, 24)
lat = np.append(lat, 58.35)
lon = np.append(lon, -136.2)
temp = np.dstack((temp, temp[:, :, 0:1]))
salt = np.dstack((salt, salt[:, :, 0:1]))

# ----------------------------------------------------------------
# flood CTD stations to the whole grid
seednum = stn.size
c = np.zeros(seednum)
r = np.zeros(seednum)

for i in range(seednum):
    dis = (lat_grd-lat[i])**2+(lon_grd-lon[i])**2
    c[i], r[i] = np.where(dis == np.nanmin(dis))

temp_ini = np.zeros((1, zlev, eta, xi))*np.NaN
salt_ini = np.zeros((1, zlev, eta, xi))*np.NaN

# actual flood
fl = floodFill(c, r, msk)

for i in range(eta):
    for j in range(xi):
        if fl[i, j] != -1:
            salt_ij = salt[:, month, fl[i, j]]
            temp_ij = temp[:, month, fl[i, j]]

            z_grd = -h[i, j]*Cs_r
            salt_ini[0, :, i, j] = interp1d(z, salt_ij)(z_grd)
            temp_ini[0, :, i, j] = interp1d(z, temp_ij)(z_grd)

# ----------------------------------------------------------------
# post-process
# filter the results
for i in range(zlev):
    temp_ini[:, i, :, :] = nan_gaussian_filter(temp_ini[:, i, :, :], 5.0)
    salt_ini[:, i, :, :] = nan_gaussian_filter(salt_ini[:, i, :, :], 5.0)

# mask data below water depth
temp_ini = np.ma.masked_invalid(temp_ini)
salt_ini = np.ma.masked_invalid(salt_ini)

# set zeta, u, v to zeros
msk_rho = np.zeros((1, zlev, grd.hgrid.mask_rho.shape[0], grd.hgrid.mask_rho.shape[1]))
msk_u = np.zeros((1, zlev, grd.hgrid.mask_u.shape[0], grd.hgrid.mask_u.shape[1]))
msk_v = np.zeros((1, zlev, grd.hgrid.mask_v.shape[0], grd.hgrid.mask_v.shape[1]))

for i in range(zlev):
    msk_rho[0, i, :, :] = grd.hgrid.mask_rho
    msk_u[0, i, :, :] = grd.hgrid.mask_u
    msk_v[0, i, :, :] = grd.hgrid.mask_v

zeta = np.zeros(msk_rho.shape)
u = np.zeros(msk_u.shape)
v = np.zeros(msk_v.shape)

zeta = np.ma.masked_where(msk_rho==0, zeta)
u = np.ma.masked_where(msk_u==0, u)
v = np.ma.masked_where(msk_v==0, v)

ubar = u.mean(axis=1)
vbar = v.mean(axis=1)

# ----------------------------------------------------------------
# write into nc file
ic_file = dst_dir + 'bc_ic/' + grd.name + '_ic_' + src_grd_name + '.nc'
spval = -1.0e20
class nctime(object):
    pass

nctime.long_name = 'time'
nctime.units = 'days since 1900-01-01 00:00:00'
pyroms_toolbox.nc_create_roms_file(ic_file, grd, nctime)

fh = nc.Dataset(ic_file, 'r+')
fh.createVariable('zeta', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['zeta'].long_name = 'free-surface'
fh.variables['zeta'].units = 'meter'
fh.variables['zeta'].field = 'free-surface, scalar, series'
fh.variables['zeta'][:] = zeta

fh.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['temp'].long_name = 'potential temperature'
fh.variables['temp'].units = 'Celsius'
fh.variables['temp'].field = 'temperature, scalar, series'
fh.variables['temp'][:] = temp_ini

fh.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['salt'].long_name = 'salinity'
fh.variables['salt'].units = 'PSU'
fh.variables['salt'].field = 'salinity, scalar, series'
fh.variables['salt'][:] = salt_ini

fh.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
fh.variables['u'].long_name = '3D u-momentum component'
fh.variables['u'].units = 'meter second-1'
fh.variables['u'].field = 'u-velocity, scalar, series'
fh.variables['u'][:] = u

fh.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
fh.variables['v'].long_name = '3D v-momentum component'
fh.variables['v'].units = 'meter second-1'
fh.variables['v'].field = 'v-velocity, scalar, series'
fh.variables['v'][:] = v

fh.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
fh.variables['ubar'].long_name = '2D u-momentum component'
fh.variables['ubar'].units = 'meter second-1'
fh.variables['ubar'].field = 'ubar-velocity, scalar, series'
fh.variables['ubar'][:] = ubar

fh.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
fh.variables['vbar'].long_name = '2D v-momentum component'
fh.variables['vbar'].units = 'meter second-1'
fh.variables['vbar'].field = 'vbar-velocity, scalar, series'
fh.variables['vbar'][:] = vbar

fh.close()
