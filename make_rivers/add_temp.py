import numpy as np
from scipy import signal
import netCDF4
from datetime import datetime, timedelta
import sys

import pyroms
from gb_toolbox import gb_ctd

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

if len(sys.argv)>0:
    grd1 = sys.argv[1]
else:
    grd1 = 'GB_lr'

my_year = 2008
tag = 'Hill'

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

temp_file = in_dir + 'ctd.nc'

out = netCDF4.Dataset(out_file, 'a', format='NETCDF3_64BIT')
river_time = out.variables['river_time'][:]
out.close()

# Read the river temperatures
ctd = gb_ctd.rd_ctd(temp_file)

# calculate climatology
clm_ctd = gb_ctd.cal_ctd_climatology(ctd)
clm_t = np.mean(clm_ctd['t'][:, 12, 1:5], axis=1)
clm_time = np.array([15, 45, 74, 105, 136, 166, 197, 227, 258, 288, 319, 349])

clm_t = np.concatenate((clm_t, clm_t, clm_t))
clm_time = np.concatenate((clm_time-365, clm_time, clm_time+365))

river_source_time_raw = np.arange(-300, 366+300)
river_temp_raw = np.interp(river_source_time_raw, clm_time, clm_t)

b, a = signal.butter(8, 0.015)
river_temp_raw2 = signal.filtfilt(b, a, river_temp_raw)

river_source_time = np.arange(366)+1
river_temp = np.interp(river_source_time, river_source_time_raw, river_temp_raw2)

river_yearday = np.array([(datetime(1900, 1, 1)+timedelta(i)).timetuple().tm_yday for i in river_time])
river_temp = np.interp(river_yearday, river_source_time, river_temp)
river_salt = np.zeros(river_temp.shape)

# create file with all the objects
out = netCDF4.Dataset(out_file, 'a', format='NETCDF3_64BIT')

# out.createDimension('river_tracer_time', len(river_time))
# 
# times = out.createVariable('river_tracer_time', 'f8', ('river_tracer_time'))
# times.units = 'day'
# times.cycle_length = 365.25
# times.long_name = 'river tracer time'

temp = out.createVariable('river_temp', 'f8', ('river_time'))
temp.long_name = 'river runoff temperature'
temp.units = 'Celcius'
temp.time = 'river_time'

salt = out.createVariable('river_salt', 'f8', ('river_time'))
salt.long_name = 'river runoff salinity'
salt.units = 'PSU'
salt.time = 'river_time'

# out.variables['river_tracer_time'][:] = river_time
out.variables['river_temp'][:] = river_temp
out.variables['river_salt'][:] = river_salt

print('add temp, salt done')

out.close()
