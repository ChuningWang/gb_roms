import netCDF4 as nc
import numpy as np

rstfile = '/glade/p/work/chuning/gb_roms/bc_ic/GB-SPINUP_rst.nc'
fh = nc.Dataset(rstfile, 'a')
t0 = fh.variables['ocean_time'][0]
t1 = fh.variables['ocean_time'][1]
if t0>t1:
    fh.variables['ocean_time'][0] = 36527*24*60*60
    fh.variables['ocean_time'][1] = fh.variables['ocean_time'][0]-24*60*60
else:
    fh.variables['ocean_time'][1] = 36527*24*60*60
    fh.variables['ocean_time'][0] = fh.variables['ocean_time'][0]-24*60*60
fh.close()
