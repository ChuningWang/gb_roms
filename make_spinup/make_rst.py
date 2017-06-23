import netCDF4 as nc
import numpy as np

rstfile = '/glade/scratch/chuning/tmpdir_GB-SPINUP/GB-SPINUP_rst.nc'
fh = nc.Dataset(rstfile)
t0 = fh.variables['ocean_time'][0]
t1 = fh.variables['ocean_time'][1]
if t0>t1:
    fh.variables['ocean_time'][0] = 36526*24*60*60
    fh.variables['ocean_time'][1] = fh.variables['ocean_time'][0]-24*60*60
else:
    fh.variables['ocean_time'][1] = 36526*24*60*60
    fh.variables['ocean_time'][0] = fh.variables['ocean_time'][0]-24*60*60

fh.close()
