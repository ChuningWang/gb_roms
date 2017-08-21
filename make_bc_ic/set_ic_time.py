import netCDF4 as nc
import numpy as np

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']

tag='2008_06_04'
icfile = out_dir + 'bc_ic/GlacierBay_usgs_ic_'+tag+'_SODA3.3.1_0.25.nc'
fh = nc.Dataset(icfile, 'a')
# fh.variables['ocean_time'][0] = 36527.
fh.variables['ocean_time'][0] = np.floor(fh.variables['ocean_time'][0])
fh.close()
