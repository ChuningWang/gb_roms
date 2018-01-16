import numpy as np
from scipy import signal
import netCDF4
from datetime import datetime, timedelta
import sys

import pyroms

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
tag = 'Hill'

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

out = netCDF4.Dataset(out_file, 'a', format='NETCDF3_64BIT')
river_temp = out.variables['river_temp'][:]
out.close()

# create file with all the objects
out = netCDF4.Dataset(out_file, 'a', format='NETCDF3_64BIT')

dye1 = out.createVariable('river_dye_01', 'f8', ('river_time'))
dye1.long_name = 'river dye 01'
dye1.units = 'kilogram meter-3'
dye1.time = 'river_time'

dye2 = out.createVariable('river_dye_02', 'f8', ('river_time'))
dye2.long_name = 'river dye 02'
dye2.units = 'kilogram meter-3'
dye2.time = 'river_time'

dye3 = out.createVariable('river_dye_03', 'f8', ('river_time'))
dye3.long_name = 'river dye 03'
dye3.units = 'kilogram meter-3'
dye3.time = 'river_time'

# out.variables['river_tracer_time'][:] = river_time
out.variables['river_dye_01'][:] = np.zeros(river_temp.shape)
out.variables['river_dye_02'][:] = np.zeros(river_temp.shape)
out.variables['river_dye_03'][:] = np.ones(river_temp.shape)

print('add dye done')

out.close()
