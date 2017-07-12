import numpy as np
import netCDF4 as nc
from datetime import datetime

import subprocess
import os

import pyroms
import pyroms_toolbox

from remap_bdry import remap_bdry
from remap_bdry_uv import remap_bdry_uv

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

bc_ic_dir = out_dir + 'bc_ic/'
frc_dir = out_dir + 'frc/'
tag = 'SODA3.3.1'
tag1 = 'JRA55v1.1'
tag2 = '2000_01_03'
tag3 = 'Hill'
grd1 = 'GB_USGS'
my_year = 2000

t0 = (datetime(2000, 1, 3, 0, 0, 0)-datetime(1900, 1, 1, 0, 0, 0)).days
grd = pyroms.grid.get_ROMS_grid(grd1)

fix_ic = 0
fix_bc = 0
fix_frc = 1
fix_river = 0


if fix_ic == 1:

    fh = nc.Dataset(bc_ic_dir + grd.name + '_ic_' + tag2 + '_' + tag + '.nc', 'a')
    fh.variables['ocean_time'][:] = t0
    fh.close()

if fix_bc == 1:

    k = 10
    dt = 7

    in_file = bc_ic_dir + grd.name + '_bdry_' + str(my_year) + '_' + tag + '.nc'
    # in_file = bc_ic_dir + 'GlacierBay' + '_bdry_' + str(my_year) + '_' + tag + '.nc'
    out_file = bc_ic_dir + grd.name + '_bdry_spinup_' + tag + '.nc'

    for i in range(k):

        out_file_temp = bc_ic_dir + 'temp/' + "%02d" % i + '.nc'
        command1 = 'ncks -d ocean_time,0,0 ' + in_file + ' ' + out_file_temp
        subprocess.call(command1, shell=True)
        # change time stamp
        fh = nc.Dataset(out_file_temp, 'a')
        fh.variables['ocean_time'][:] = t0 + dt*i
        fh.close()

    # concatenate
    cmd1 = 'ncrcat -O ' + bc_ic_dir + 'temp/*.nc -o ' + out_file,
    subprocess.call(cmd1, shell=True)
    subprocess.call('rm ' + bc_ic_dir + 'temp/*.nc', shell=True)

if fix_frc == 1:

    k = 70
    dt = 1

    vlist = ['Qair', 'rain', 'snow', 'lwrad_down', 'swrad', 'Pair', 'Tair', 'Uwind', 'Vwind']
    vlist2 = ['qair', 'rain', 'snow', 'lrf', 'srf', 'pair', 'tair', 'wind', 'wind']

    for i in range(len(vlist)):
        in_file = frc_dir + vlist[i] + '_' + str(my_year) + '_' + tag1 + '.nc'
        out_file = frc_dir + vlist[i] + '_spinup_' + tag1 + '.nc'
        command1 = 'ncks -d ' + vlist2[i] + '_time,15,22 ' + in_file + ' ' + out_file
        subprocess.call(command1, shell=True)

        for j in range(k):
            out_file_temp = frc_dir + 'temp/' + vlist[i] + '_' + tag1 + '_' + "%02d" % j + '.nc'
            subprocess.call('cp ' + out_file + ' ' + out_file_temp, shell=True)
            # change time stamp
            fh = nc.Dataset(out_file_temp, 'a')
            fh.variables[vlist2[i]+'_time'][:] = fh.variables[vlist2[i]+'_time'][:] + dt*j
            fh.close()

        # concatenate
        cmd1 = 'ncrcat -O ' + frc_dir + 'temp/' + vlist[i] + '_' + tag1 + '_*.nc -o ' + \
                out_file,
        subprocess.call(cmd1, shell=True)
        subprocess.call('rm ' + frc_dir + 'temp/' + vlist[i] + '_' + tag1 + '_*.nc', shell=True)

if fix_river == 1:

    in_file = frc_dir + grd.name + '_rivers_' + str(my_year) + '_' + tag3 + '.nc'
    out_file = frc_dir + grd.name + '_rivers_spinup_' + tag3 + '.nc'
    cmd1 = 'ncks -d river_time,2,2 ' + in_file + ' ' + out_file

    subprocess.call(cmd1, shell=True)

    k = 10
    dt = 7

    for j in range(k):
        out_file_temp = frc_dir + 'temp/rivers_' + "%02d" % j + '.nc'
        subprocess.call('cp ' + out_file + ' ' + out_file_temp, shell=True)
        fh = nc.Dataset(out_file_temp, 'a')
        fh.variables['river_time'][:] = fh.variables['river_time'][:] + dt*j
        fh.close()

    # concatenate
    cmd1 = 'ncrcat -O ' + frc_dir + 'temp/rivers_*.nc -o ' + \
            out_file,
    subprocess.call(cmd1, shell=True)
    subprocess.call('rm ' + frc_dir + 'temp/*.nc', shell=True)

