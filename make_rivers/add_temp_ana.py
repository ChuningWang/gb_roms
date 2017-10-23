import numpy as np
from scipy import signal
import netCDF4 as nc
from datetime import datetime, timedelta
import sys

import pyroms
from ocean_toolbox import ctd

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

my_year = 2008
tag = 'Hill_ana'
t0 = 10.
s0 = 0.

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'
out = nc.Dataset(out_file, 'a', format='NETCDF3_64BIT')
river_time = out.variables['river_time'][:]
out.close()

river_temp = t0*np.ones(river_time.shape)
river_salt = s0*np.ones(river_time.shape)

# create file with all the objects
fh = nc.Dataset(out_file, 'a', format='NETCDF3_64BIT')

# temp = fh.createVariable('river_temp', 'f8', ('river_time'))
temp = fh.variables['river_temp']
temp.long_name = 'river runoff temperature'
temp.units = 'Celcius'
temp.time = 'river_time'

# salt = fh.createVariable('river_salt', 'f8', ('river_time'))
salt = fh.variables['river_salt']
salt.long_name = 'river runoff salinity'
salt.units = 'PSU'
salt.time = 'river_time'

fh.variables['river_temp'][:] = river_temp
fh.variables['river_salt'][:] = river_salt

print('add temp, salt done')

fh.close()
