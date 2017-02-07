import re
import numpy as np
import netCDF4
import sys
import pdb

outfile = sys.argv[1]

# Set the vertical distribution of the river transport.

out = netCDF4.Dataset(outfile, 'a', format='NETCDF3_64BIT')
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
