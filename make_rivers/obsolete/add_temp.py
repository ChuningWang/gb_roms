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
tag = 'Hill'

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'
out = nc.Dataset(out_file, 'a', format='NETCDF3_64BIT')
river_time = out.variables['river_time'][:]
out.close()
river_datetime = nc.num2date(river_time, 'days since 1900-01-01')
river_yearday = np.array([i.timetuple().tm_yday for i in river_datetime])

# Read river temperature
info = {'data_dir': '/glade/p/work/chuning/data/ctd_raw/',
        'file_dir': '/glade/p/work/chuning/data/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['temp'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
       }

c = ctd.ctd(info)
c()

temp = c.climatology['temp'][:, :, 12]
temp = temp[:5, :].mean(axis=0)
time = c.climatology['time']
temp = np.concatenate((np.array([temp[-1]]), temp, np.array([temp[0]])))
time = np.concatenate((np.array([time[-1]-365.]), time, np.array([time[0]+365.])))

river_temp = np.interp(river_yearday, time, temp)
river_salt = np.zeros(river_temp.shape)

# create file with all the objects
fh = nc.Dataset(out_file, 'a', format='NETCDF3_64BIT')

temp = fh.createVariable('river_temp', 'f8', ('river_time'))
temp.long_name = 'river runoff temperature'
temp.units = 'Celcius'
temp.time = 'river_time'

salt = fh.createVariable('river_salt', 'f8', ('river_time'))
salt.long_name = 'river runoff salinity'
salt.units = 'PSU'
salt.time = 'river_time'

fh.variables['river_temp'][:] = river_temp
fh.variables['river_salt'][:] = river_salt

print('add temp, salt done')

fh.close()
