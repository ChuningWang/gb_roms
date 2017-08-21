import numpy as np
import netCDF4
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_USGS'

my_year = 2008
grd = pyroms.grid.get_ROMS_grid(grd1)
tag = 'Hill'
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

# Set the vertical distribution of the river transport.

out = netCDF4.Dataset(out_file, 'a', format='NETCDF3_64BIT')
N = len(out.dimensions['s_rho'])
Nr = len(out.dimensions['river'])

vshape = np.zeros([N, Nr])
for k in range(N):
    vshape[k,:] = k

area = sum(vshape[:,0])
vshape = (1.0/area)*vshape
print vshape[:,0]

out.variables['river_Vshape'][:] = vshape
out.close()
