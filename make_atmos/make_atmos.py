import subprocess
import os
import numpy as np

import pyroms
import pyroms_toolbox

import netCDF4 as nc

my_year = 2000
tag = 'JRA55v0.8'
vlist = ['q_10', 'rain', 'rlds', 'rsds', 'slp', 'snow', 't_10', 'u_10', 'v_10']

dst_grd = pyroms.grid.get_ROMS_grid('GB')

data_dir = '/Volumes/R1/scratch/chuning/gb_roms/data/jra55do/'
data_dir_year = data_dir + str(my_year) + '/'
dst_dir='/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/'

filelst = subprocess.check_output(['ls', data_dir_year]).replace('/n',' ').split()

# for filein in filelst:
for var in vlist:
    filein = var + '*.nc'
    fileout = var + '.' + str(my_year) + '.' + tag + '.nc'
    command1 = 'ncks -d latitude,255,270 -d longitude,390,410 ' \
                + data_dir_year + filein + ' -O ' + dst_dir + 'temp/' + fileout
    subprocess.call(command1, shell=True)

# merge file
atmos_file = dst_dir + dst_grd.name + '_atmos_' + str(my_year) + '_' + tag + '.nc'

kk = 0
for var in vlist:
    filein = dst_dir + 'temp/' + var + '*.nc'
    if kk==0:
        subprocess.call('mv ' + filein + ' ' + atmos_file, shell=True)
    else:
        subprocess.call('ncks -A ' + filein + ' -o ' + atmos_file ,shell=True)
    kk = kk+1

# clean up
subprocess.call('rm ' + dst_dir + 'temp/*', shell=True)

# change variable names, unit conversion
subprocess.call('ncrename -d longitude,lon ' + atmos_file, shell=True)
subprocess.call('ncrename -d latitude,lat ' + atmos_file, shell=True)
subprocess.call('ncrename -d time,ocean_time ' + atmos_file, shell=True)

subprocess.call('ncrename -v latitude,lat ' + atmos_file, shell=True)
subprocess.call('ncrename -v longitude,lon ' + atmos_file, shell=True)

subprocess.call('ncrename -v huss_10m,Qair ' + atmos_file, shell=True)
subprocess.call('ncrename -v prrn,rain ' + atmos_file, shell=True)
subprocess.call('ncrename -v prsn,snow ' + atmos_file, shell=True)
subprocess.call('ncrename -v rlds,lwrad_down ' + atmos_file, shell=True)
subprocess.call('ncrename -v rsds,swrad ' + atmos_file, shell=True)
subprocess.call('ncrename -v psl,Pair ' + atmos_file, shell=True)
subprocess.call('ncrename -v tas_10m,Tair ' + atmos_file, shell=True)
subprocess.call('ncrename -v uas_10m,Uwind ' + atmos_file, shell=True)
subprocess.call('ncrename -v vas_10m,Vwind ' + atmos_file, shell=True)

fh = nc.Dataset(atmos_file, 'a') 
time = fh.variables['time'][:]
nt = len(time)
fh.createDimension('qair_time', nt)
fh.createDimension('rain_time', nt)
fh.createDimension('snow_time', nt)
fh.createDimension('lrf_time', nt)
fh.createDimension('srf_time', nt)
fh.createDimension('pair_time', nt)
fh.createDimension('tair_time', nt)
fh.createDimension('wind_time', nt)

fh.createVariable('qair_time', 'f8', ('qair_time'))
fh.createVariable('rain_time', 'f8', ('rain_time'))
fh.createVariable('snow_time', 'f8', ('snow_time'))
fh.createVariable('lrf_time', 'f8', ('lrf_time'))
fh.createVariable('srf_time', 'f8', ('srf_time'))
fh.createVariable('pair_time', 'f8', ('pair_time'))
fh.createVariable('tair_time', 'f8', ('tair_time'))
fh.createVariable('wind_time', 'f8', ('wind_time'))

fh.variables['qair_time'][:] = time
fh.variables['rain_time'][:] = time
fh.variables['snow_time'][:] = time
fh.variables['lrf_time'][:] = time
fh.variables['srf_time'][:] = time
fh.variables['pair_time'][:] = time
fh.variables['tair_time'][:] = time
fh.variables['wind_time'][:] = time
fh.close()
