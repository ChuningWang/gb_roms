from shutil import copy2

import numpy as np
import netCDF4 as nc
import sys

import pyroms

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

# load Glacier Bay grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

tag = 'Hill'
my_year = 2008
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

sens = [0.5, 0.1]
for i in sens:
    out_file_i = out_dir + 'frc/' + grd.name + '_rivers_' + \
            str(my_year) + '_' + tag + '_' + "%4.2f" % i + '.nc'
    copy2(out_file, out_file_i)

    # scale the new river file
    fin = nc.Dataset(out_file_i, 'a')
    fin.variables['river_transport'][:] = fin.variables['river_transport'][:]*i
    fin.close()
