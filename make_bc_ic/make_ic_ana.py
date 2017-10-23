'''
generate uniform IC for a test run.
'''

import sys
import numpy as np
import pyroms
import pyroms_toolbox
from datetime import datetime

import netCDF4 as nc

# ----------------------------------------------------------------
# preparation
import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
dst_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

month = 6
my_year = 2008
s0 = 30
t0 = 10
zeta0 = 0
u0 = 0
v0 = 0
src_grd_name = str(my_year) + '_' + "%02d" % month + '_uniform'

# Load target grid and land mask
grd = pyroms.grid.get_ROMS_grid(grd1)
msk = grd.hgrid.mask
eta, xi = msk.shape

# ----------------------------------------------------------------
# write into nc file
ic_file = dst_dir + 'bc_ic/' + grd.name + '_ic_' + src_grd_name + '.nc'
spval = -1.0e20
class nctime(object):
    pass

nctime.long_name = 'time'
nctime.units = 'days since 1900-01-01 00:00:00'
pyroms_toolbox.nc_create_roms_file(ic_file, grd, nctime)

ocean_time = nc.date2num(datetime(my_year, month, 15), nctime.units)

fh = nc.Dataset(ic_file, 'r+')

fh.variables['ocean_time'][:] = ocean_time

fh.createVariable('zeta', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['zeta'].long_name = 'free-surface'
fh.variables['zeta'].units = 'meter'
fh.variables['zeta'].field = 'free-surface, scalar, series'
fh.variables['zeta'][:] = zeta0

fh.createVariable('temp', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['temp'].long_name = 'potential temperature'
fh.variables['temp'].units = 'Celsius'
fh.variables['temp'].field = 'temperature, scalar, series'
fh.variables['temp'][:] = t0

fh.createVariable('salt', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
fh.variables['salt'].long_name = 'salinity'
fh.variables['salt'].units = 'PSU'
fh.variables['salt'].field = 'salinity, scalar, series'
fh.variables['salt'][:] = s0

fh.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
fh.variables['u'].long_name = '3D u-momentum component'
fh.variables['u'].units = 'meter second-1'
fh.variables['u'].field = 'u-velocity, scalar, series'
fh.variables['u'][:] = u0

fh.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
fh.variables['v'].long_name = '3D v-momentum component'
fh.variables['v'].units = 'meter second-1'
fh.variables['v'].field = 'v-velocity, scalar, series'
fh.variables['v'][:] = v0

fh.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
fh.variables['ubar'].long_name = '2D u-momentum component'
fh.variables['ubar'].units = 'meter second-1'
fh.variables['ubar'].field = 'ubar-velocity, scalar, series'
fh.variables['ubar'][:] = u0

fh.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
fh.variables['vbar'].long_name = '2D v-momentum component'
fh.variables['vbar'].units = 'meter second-1'
fh.variables['vbar'].field = 'vbar-velocity, scalar, series'
fh.variables['vbar'][:] = v0

fh.close()
