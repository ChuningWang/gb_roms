import numpy as np
import netCDF4 as nc
import csv

import os
import fnmatch

def recursive_glob(rootdir='.', pattern='*'):
    """Search recursively for files matching a specified pattern."""

    matches = []
    for root, dirnames, filenames in os.walk(rootdir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

import read_host_info
sv = read_host_info.read_host_info()
bathy_dir = sv['in_dir']

# first, glob all data path
bathy_raw_dir = bathy_dir + 'bathy/raw/'
flist = recursive_glob(rootdir=bathy_raw_dir, pattern='*.xyz')

# read csv files
lat = []
lon = []
z = []
for fn in flist:
    f = open(fn, 'rb')
    rd = csv.reader(f, delimiter='\t')
    next(rd, None)
    for row in rd:
        lon.append(float(row[1]))
        lat.append(float(row[2]))
        z.append(float(row[3]))

lon = np.array(lon)
lat = np.array(lat)
z = np.array(z)

fh = nc.Dataset(bathy_dir+'bathy_noaa.nc', 'w')
fh.createDimension('ct')
fh.createVariable('lon', 'd', ('ct'))
fh.createVariable('lat', 'd', ('ct'))
fh.createVariable('z', 'd', ('ct'))
fh.variables['lon'][:] = lon
fh.variables['lat'][:] = lat
fh.variables['z'][:] = z
fh.close()
