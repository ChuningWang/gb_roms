import subprocess
import os
import numpy as np

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

import read_host_info
sv = read_host_info.read_host_info()
dst_dir = sv['out_dir']
data_dir = sv['soda_dir']

grd1 = 'GB_USGS'

my_year = 2000

data_dir_year = data_dir + str(my_year) + '/'

filelst = subprocess.check_output(['ls', data_dir_year]).replace('/n',' ').split()

src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(data_dir+'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', xrange=(400, 500), yrange=(180, 280) )
dst_grd = pyroms.grid.get_ROMS_grid(grd1)

for filein in filelst:
    tag = filein.replace('soda3.3.1_5dy_ocean_reg_','').replace('.nc','')
    print '\nBuild OBC file for time %s' %filein
    zeta_dst_file = dst_dir + 'temp/' + dst_grd.name + '_bdry_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + 'temp/' +  dst_grd.name + '_bdry_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + 'temp/' +  dst_grd.name + '_bdry_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + 'temp/' +  dst_grd.name + '_bdry_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + 'temp/' +  dst_grd.name + '_bdry_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap ssh
    zeta = remap_bdry('ssh', data_dir_year + filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

    # reload grid with zeta (more accurate)
    dst_grd = pyroms.grid.get_ROMS_grid(grd1, zeta=zeta)

    # regrid temp, salt and velocities
    remap_bdry('temp',data_dir_year + filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
    remap_bdry('salt',data_dir_year + filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
    remap_bdry_uv(data_dir_year + filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + 'temp/' + dst_grd.name + '_bdry_' + tag + '_' + src_grd.name + '.nc'

    command1 = 'mv '      + zeta_dst_file + ' '    + bdry_file
    command2 = 'ncks -A ' + temp_dst_file + ' -o ' + bdry_file
    command3 = 'ncks -A ' + salt_dst_file + ' -o ' + bdry_file
    command4 = 'ncks -A ' + u_dst_file    + ' -o ' + bdry_file
    command5 = 'ncks -A ' + v_dst_file    + ' -o ' + bdry_file

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

# concatenate
cmd1 = 'ncrcat -O ' + dst_dir + 'temp/' + dst_grd.name + '_bdry_' + str(my_year) + '_*_' + src_grd.name + '.nc ' + \
        dst_dir + dst_grd.name + '_bdry_' + str(my_year) + '_' + src_grd.name + '.nc',
subprocess.call(cmd1, shell=True)
subprocess.call('rm ' + dst_dir + 'temp/*', shell=True)
