import netCDF4 as nc
import numpy as np

rstfile = '/glade/scratch/chuning/tmpdir_GB-SPINUP/GB-SPINUP_rst.nc'
fh = nc.Dataset(rstfile)
fh.variables['ocean_time'][0] = fh.variables['ocean_time'][1]+24*60*60
fh.close()
