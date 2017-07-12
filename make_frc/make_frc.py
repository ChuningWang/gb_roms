import subprocess
import os
import numpy as np

import pyroms
import pyroms_toolbox

import netCDF4 as nc

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
dst_dir = sv['out_dir']

grd1 = 'GB_USGS'

my_year = 2000
tag = 'JRA55v1.1'
vlist = ['Qair', 'rain', 'snow', 'lwrad_down', 'swrad', 'Pair', 'Tair', 'Uwind', 'Vwind']
vlist2 = ['huss_10m', 'prrn', 'prsn', 'rlds', 'rsds', 'psl', 'tas_10m', 'uas_10m', 'vas_10m']
vlist3 = ['qair', 'rain', 'snow', 'lrf', 'srf', 'pair', 'tair', 'wind', 'wind']

dst_grd = pyroms.grid.get_ROMS_grid(grd1)

data_dir =  in_dir + 'jra55do/'
data_dir_year = data_dir + str(my_year) + '/'

filelst = subprocess.check_output(['ls', data_dir_year]).replace('/n',' ').split()

# for filein in filelst:
for var in vlist:
    filein = var + '*.nc'
    fileout = var + '.' + str(my_year) + '.' + tag + '.nc'
    command1 = 'ncks -d latitude,255,270 -d longitude,390,410 ' \
                + data_dir_year + filein + ' -O ' + dst_dir + 'temp/' + fileout
    subprocess.call(command1, shell=True)

kk = 0
for var in vlist:
    filein = dst_dir + 'temp/' + var + '*.nc'
    fileout = dst_dir + 'frc/' + var + '_' + str(my_year) + '_' + tag + '.nc' 
    subprocess.call('mv ' + filein + ' ' + fileout, shell=True)

    subprocess.call('ncrename -d longitude,lon ' + fileout, shell=True)
    subprocess.call('ncrename -d latitude,lat ' + fileout, shell=True)
    subprocess.call('ncrename -d time,' + vlist3[kk] + '_time ' + fileout, shell=True)

    subprocess.call('ncrename -v longitude,lon ' + fileout, shell=True)
    subprocess.call('ncrename -v latitude,lat ' + fileout, shell=True)
    subprocess.call('ncrename -v time,' + vlist3[kk] + '_time ' + fileout, shell=True)
    subprocess.call('ncrename -v ' + vlist2[kk] + ',' + var + ' ' + fileout, shell=True)

    fh = nc.Dataset(fileout, 'a')
    fh.variables[var].coordinates = 'lon lat'
    fh.variables['lon'][:] = fh.variables['lon'][:]-360
    fh.close()

    kk = kk+1

# convert Tair unit from Kelvin to Celcius
fileout = dst_dir + 'frc/' + 'Tair' + '_' + str(my_year) + '_' + tag + '.nc' 
fh = nc.Dataset(fileout, 'a')
fh.variables['Tair'][:] = fh.variables['Tair'][:]-273.15
fh.close()

