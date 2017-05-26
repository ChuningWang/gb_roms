import numpy as np
import netCDF4 as nc
from datetime import datetime

import subprocess
import os

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv


def fix_ic():

    t0 = (datetime(2000, 1, 3, 0, 0, 0)-datetime(1900, 1, 1, 0, 0, 0)).days
    fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_spinup_roms/data/GlacierBay_ic_2000.nc', 'a')
    fh.variables['ocean_time'][:] = t0
    fh.close()


def make_bc():

    filein = '/Volumes/P1/Data/SODA/SODA_3.3.1/2000/soda3.3.1_5dy_ocean_reg_2000_01_03.nc'
    dst_dir= '/Volumes/R1/scratch/chuning/gb_spinup_roms/data/'
    src_grd_dir = '/Volumes/P1/Data/SODA/SODA_3.3.1/grid/SODA3_0.5deg_grid.nc'

    tag = '2000'

    src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(src_grd_dir, name='SODA3.3.1', xrange=(400, 500), yrange=(180, 280) )
    dst_grd = pyroms.grid.get_ROMS_grid('GB')

    zeta_dst_file = dst_dir + 'temp/' + dst_grd.name + '_bdry_zeta_' + tag + '_' + src_grd.name + '.nc'
    temp_dst_file = dst_dir + 'temp/' +  dst_grd.name + '_bdry_temp_' + tag + '_' + src_grd.name + '.nc'
    salt_dst_file = dst_dir + 'temp/' +  dst_grd.name + '_bdry_salt_' + tag + '_' + src_grd.name + '.nc'
    u_dst_file    = dst_dir + 'temp/' +  dst_grd.name + '_bdry_u_'    + tag + '_' + src_grd.name + '.nc'
    v_dst_file    = dst_dir + 'temp/' +  dst_grd.name + '_bdry_v_'    + tag + '_' + src_grd.name + '.nc'

    # remap ssh
    zeta = remap_bdry('ssh', filein, src_grd, dst_grd, zeta_dst_file, dst_dir=dst_dir)

    # reload grid with zeta (more accurate)
    dst_grd = pyroms.grid.get_ROMS_grid('GB', zeta=zeta)

    # regrid temp, salt and velocities
    remap_bdry('temp',filein, src_grd, dst_grd, temp_dst_file, dst_dir=dst_dir)
    remap_bdry('salt',filein, src_grd, dst_grd, salt_dst_file, dst_dir=dst_dir)
    remap_bdry_uv(filein, src_grd, dst_grd, u_dst_file, v_dst_file, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + dst_grd.name + '_bdry_' + tag + '.nc'

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


def fix_bc():

    k = 5
    dt = 7

    # make a few copies of the boundary condition file
    dst_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/data/'
    tag = 'GlacierBay_bdry_2000'

    t0 = (datetime(2000, 1, 3, 0, 0, 0)-datetime(1900, 1, 1, 0, 0, 0)).days

    for i in range(k):
        subprocess.call('cp ' + dst_dir + tag + '.nc ' + dst_dir + 'temp/' + tag + '.' + str(i) + '.nc', shell=True)
        # change time stamp
        fh = nc.Dataset(dst_dir + 'temp/' + tag + '.' + str(i) + '.nc', 'a')
        fh.variables['ocean_time'][:] = t0 + dt*i
        fh.close()

    # concatenate
    cmd1 = 'ncrcat -O ' + dst_dir + 'temp/' + tag + '.*.nc ' + \
            dst_dir + tag + '_repeat.nc',
    subprocess.call(cmd1, shell=True)
    subprocess.call('rm ' + dst_dir + 'temp/*', shell=True)


def fix_frc():

    dir_in = '/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/frc/'
    dir_out = '/Volumes/R1/scratch/chuning/gb_spinup_roms/data/frc/'
    tag = '2000'

    k = 28
    dt = 1

    vlist = ['Qair', 'rain', 'snow', 'lwrad_down', 'swrad', 'Pair', 'Tair', 'Uwind', 'Vwind']
    vlist2 = ['qair', 'rain', 'snow', 'lrf', 'srf', 'pair', 'tair', 'wind', 'wind']

    for i in range(len(vlist)):
        fileout = dir_out + vlist[i] + '_' + tag + '.nc'
        command1 = 'ncks -d ' + vlist2[i] + '_time,16,23 ' + dir_in + vlist[i] + '_2000_JRA55v1.1.nc' + ' ' + \
                   fileout
        subprocess.call(command1, shell=True)

        for j in range(k):
            subprocess.call('cp ' + fileout + ' ' + dir_out + vlist[i] + '_' + tag + '.' + "%02d" % j + '.nc', shell=True)
            # change time stamp
            fh = nc.Dataset(dir_out + vlist[i] + '_' + tag + '.' + "%02d" % j + '.nc', 'a')
            fh.variables[vlist2[i]+'_time'][:] = fh.variables[vlist2[i]+'_time'][:] + dt*j
            fh.close()

        # concatenate
        cmd1 = 'ncrcat -O ' + dir_out + vlist[i] + '_' + tag + '.*.nc -o ' + \
                dir_out + vlist[i] + '_' + tag + '.nc',
        subprocess.call(cmd1, shell=True)
        subprocess.call('rm ' + dir_out + vlist[i] + '_' + tag + '.*.nc', shell=True)


def fix_river():

    filein = '/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/GlacierBay_rivers.nc'
    dir_out = '/Volumes/R1/scratch/chuning/gb_spinup_roms/data/'
    fileout = dir_out + 'GlacierBay_rivers_2000.nc'
    cmd1 = 'ncks -d river_time,3,3 ' + filein + ' ' + fileout

    subprocess.call(cmd1, shell=True)

    k = 5
    dt = 7

    for j in range(k):
        subprocess.call('cp ' + fileout + ' ' + dir_out + 'temp/GlacierBay_rivers_2000' + '.' + "%02d" % j + '.nc', shell=True)
        fh = nc.Dataset(dir_out + 'temp/GlacierBay_rivers_2000' + '.' + "%02d" % j + '.nc', 'a')
        fh.variables['river_time'][:] = fh.variables['river_time'][:] + dt*j
        fh.close()

    # concatenate
    cmd1 = 'ncrcat -O ' + dir_out + 'temp/GlacierBay_rivers_2000.*.nc -o ' + \
            fileout,
    subprocess.call(cmd1, shell=True)
    subprocess.call('rm ' + dir_out + 'temp/*.nc', shell=True)


fix_river()
