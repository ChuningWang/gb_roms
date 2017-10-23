import numpy as np
import netCDF4
import pyroms
import sys

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

my_year = 2008
discharge_depth = 3
grd = pyroms.grid.get_ROMS_grid(grd1)
h = grd.vgrid.h
Cs_r = grd.vgrid.Cs_r

tag = 'Hill_ana'
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

# Set the vertical distribution of the river transport.

fh = netCDF4.Dataset(out_file, 'a', format='NETCDF3_64BIT')
N = len(fh.dimensions['s_rho'])
Nr = len(fh.dimensions['river'])
eta = fh.variables['river_Eposition'][:]
xi = fh.variables['river_Xposition'][:]
sign = fh.variables['river_sign'][:]
direction = fh.variables['river_direction'][:]

hh = np.zeros(Nr)
for i in range(Nr):
    if sign[i] == 1:
        hh[i] = h[eta[i], xi[i]]
    else:
        if direction[i] == 0:
            hh[i] = h[eta[i], xi[i]-1]
        else:
            hh[i] = h[eta[i]-1, xi[i]]

vshape = np.zeros((N, Nr))
for i in range(Nr):
    z = -Cs_r*hh[i]
    msk = z <= discharge_depth
    weight = np.arange(sum(msk))
    vshape[msk, i] = weight
    vshape[:, i] = vshape[:, i]/sum(vshape[:, i])

fh.variables['river_Vshape'][:] = vshape
fh.close()
