import re
import numpy as np
import netCDF4
import sys
import pdb

my_year = 2000
grd = pyroms.grid.get_ROMS_grid('GB3')
tag = 'Hill'
out_dir = '/glade/p/work/chuning/gb_roms/frc/'
out_file = out_dir + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

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
