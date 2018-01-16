import os
import subprocess
import fnmatch
from datetime import datetime

import numpy as np
import netCDF4 as nc

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['out_dir'] + 'frc/'
out_dir = sv['out_dir'] + 'frc_clim/'

# -------------------- my inputs ---------------------------------------
my_year = 2008
tag = 'JRA55v1.1'

lon_min = 222
lon_max = 226
lat_min = 57
lat_max = 60

# lon_min = 0
# lon_max = 360
# lat_min = -90
# lat_max = 90

vlist = ['Qair', 'rain', 'snow', 'lwrad_down', 'swrad', 'Pair', 'Tair', 'Uwind', 'Vwind']
vlist2 = ['huss_10m', 'prrn', 'prsn', 'rlds', 'rsds', 'psl', 'tas_10m', 'uas_10m', 'vas_10m']
vlist3 = ['qair', 'rain', 'snow', 'lrf', 'srf', 'pair', 'tair', 'wind', 'wind']

# -------------------- funcitonals -------------------------------------
def recursive_glob(in_dir, pattern='*'):
    ''' Search recursively for files matching a specified pattern. '''

    matches = []
    for root, dirnames, filenames in os.walk(in_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

# -------------------- load data and calculate clim --------------------

for i2 in range(len(vlist3)):
# for i2 in range(1):
    flist = recursive_glob(data_dir, vlist[i2]+'*.nc')
    file_out = out_dir+vlist[i2]+'_clim.nc'

    # create output file
    subprocess.call('cp ' + flist[0] + ' ' + file_out, shell=True)
    fout = nc.Dataset(out_dir+vlist[i2]+'_clim.nc', 'a')

    for i, fn in enumerate(flist):
        fin = nc.Dataset(fn, 'r')
        if i == 0:
            time = fin.variables[vlist3[i2]+'_time'][:]
            data = fin.variables[vlist[i2]][:]
        else:
            time = np.concatenate((time, fin.variables[vlist3[i2]+'_time'][:]))
            data = np.concatenate((data, fin.variables[vlist[i2]][:]))
        fin.close()

    # # create dimensions
    # fout.createDimension('time', None)
    # fout.createDimension('latitude', eta)
    # fout.createDimension('longitude', xi)

    # # create variables
    # fout.createVariable('time', 'f8', ('time'))
    # fout.createVariable('latitude', 'f8', ('latitude'))
    # fout.createVariable('longitude', 'f8', ('longitude'))
    # fout.createVariable(vlist3[i2], 'f8', ('time', 'latitude', 'longitude'))

    pytime = nc.num2date(time, 'days since 1900-01-01')
    tbase = np.array([datetime(pyt.year, 1, 1) for pyt in pytime])
    tbase = nc.date2num(tbase, 'days since 1900-01-01')
    tbase2 = nc.date2num(datetime(my_year, 1, 1), 'days since 1900-01-01')
    yearday = time - tbase
    # time_clim = np.arange(0, 366, 0.125)
    time_clim = np.unique(yearday)
    data_clim = []

    for tc in time_clim:
        msk = yearday == tc
        data_clim.append(data[msk, :, :].mean(axis=0))

    data_clim = np.array(data_clim)

    fout.variables[vlist3[i2]+'_time'][:] = time_clim + tbase2
    fout.variables[vlist[i2]][:] = data_clim

    fout.close()
    print(out_dir+vlist[i2]+'_clim.nc'+' done...')

