import subprocess
import os
import sys
import numpy as np
import pyroms
import pyroms_toolbox
import netCDF4 as nc

from remap import remap
from remap_uv import remap_uv

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
data_dir = sv['soda_dir']

dst_dir = out_dir + 'bc_ic/'

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

my_year = 2008
tag='2008_05_30'

# filein=data_dir+str(my_year)+'/'+'soda3.3.1_5dy_ocean_reg_'+tag+'.nc'
filein=data_dir+str(my_year)+'/'+'soda3.3.1_5dy_ocean_or_'+tag+'.nc'

# load grids
# src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', xrange=(400, 500), yrange=(180, 280))
src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir + 'grid/SODA3_0.25deg_grid.nc', 
                                                      name='SODA3.3.1_0.25', xrange=(550, 600), yrange=(750, 800))
dst_grd = pyroms.grid.get_ROMS_grid(grd1)

print '\nBuild IC file from %s' %filein

zeta_dst_file = dst_dir + dst_grd.name + '_ic_zeta_' + tag + '_' + src_grd.name + '.nc'
temp_dst_file = dst_dir + dst_grd.name + '_ic_temp_' + tag + '_' + src_grd.name + '.nc'
salt_dst_file = dst_dir + dst_grd.name + '_ic_salt_' + tag + '_' + src_grd.name + '.nc'
u_dst_file    = dst_dir + dst_grd.name + '_ic_u_'    + tag + '_' + src_grd.name + '.nc'
v_dst_file    = dst_dir + dst_grd.name + '_ic_v_'    + tag + '_' + src_grd.name + '.nc'

# remap ssh
zeta = remap('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

# reload grid with zeta (more accurate)
dst_grd = pyroms.grid.get_ROMS_grid(grd1, zeta=zeta)

# regrid temp, salt and velocities
remap('temp',filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
remap('salt',filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
remap_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

# merge file
ic_file = dst_dir + dst_grd.name + '_ic_' + tag + '_' + src_grd.name + '.nc'

command1 = 'mv '      + zeta_dst_file + ' '    + ic_file
command2 = 'ncks -A ' + temp_dst_file + ' -o ' + ic_file
command3 = 'ncks -A ' + salt_dst_file + ' -o ' + ic_file
command4 = 'ncks -A ' + u_dst_file    + ' -o ' + ic_file
command5 = 'ncks -A ' + v_dst_file    + ' -o ' + ic_file

subprocess.call(command1,shell=True)
subprocess.call(command2,shell=True)
subprocess.call(command3,shell=True)
subprocess.call(command4,shell=True)
subprocess.call(command5,shell=True)

# clean up
os.remove(temp_dst_file)
os.remove(salt_dst_file)
os.remove(u_dst_file)
os.remove(v_dst_file)

# set ic file time
icfile = out_dir + 'bc_ic/' + dst_grd.name + '_ic_'+tag+'_SODA3.3.1_0.25.nc'
fh = nc.Dataset(icfile, 'a')
fh.variables['ocean_time'][0] = np.floor(fh.variables['ocean_time'][0])
fh.close()
