import netCDF4 as nc
import numpy as np

icfile = '/glade/p/work/chuning/gb_roms/bc_ic/GlacierBay_usgs_ic_2000_01_03_SODA3.3.1.nc'
fh = nc.Dataset(icfile, 'a')
fh.variables['ocean_time'][0] = 36527.
fh.close()
