import re
import numpy as np
import netCDF4
import sys
import pyroms

def add_to_lists(pairs, i, j, sign, dir):
    x1, y1 = pairs[0]

    for it in range(1,len(pairs)):
        x2, y2 = pairs[it]

	if x2 > x1:
	# negative v-velocity
	    i.append(x1)
	    j.append(y1)
	    sign.append(-1)
	    dir.append(1)
	elif x1 > x2:
	# positive v-velocity
	    i.append(x2)
	    j.append(y1)
	    sign.append(1)
	    dir.append(1)
	elif y2 > y1:
	# positive u-velocity
	    i.append(x1)
	    j.append(y1)
	    sign.append(1)
	    dir.append(0)
	elif y1 > y2:
	# negative u-velocity
	    i.append(x1)
	    j.append(y2)
	    sign.append(-1)
	    dir.append(0)
	x1 = x2
	y1 = y2

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

my_year = 2008
tag = 'Hill'
grd1 = 'GB_USGS'

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

msk_dir = home_dir + 'git/gb_roms/make_rivers/maskedge.out'

# We need to parse the output of the maskedge program for two
# different purposes:
#  1. Create the rivers file for ROMS, at least the locations part.
#  2. Create a scrip grid file for the river locations.
# This routine will only do #1 (so far).

# Read the landmask boundaries
f = open(msk_dir, 'r')
pairs = []
# Eat first line so we don't trigger the add_to_lists routine
f.readline()

# These are for the ROMS sources file
i = []
j = []
sign = []
dir = []

for line in f:
    a, b, c = re.split('\s+', line)
    if a=='-10':
	# wrap up object
	add_to_lists(pairs, i, j, sign, dir)
    elif (a=='-1' or a=='-3'):
	# wrap up object
	add_to_lists(pairs, i, j, sign, dir)
	# start new object
        pairs = []
    else:
        pairs.append([int(a),int(b)]) 

# create file with all the objects
out = netCDF4.Dataset(out_file, 'w', format='NETCDF3_64BIT')
out.type = 'ROMS RIVERS file'
out.title = 'Glacier Bay'
out.source = 'David Hill and Jordan Beamer'

out.createDimension('river_time', None)
out.createDimension('river', len(i))
out.createDimension('s_rho', 30)

times = out.createVariable('river_time', 'f8', ('river_time'))
times.units = 'days since 1900-01-01 00:00:00'
times.long_name = 'river runoff time'

river = out.createVariable('river', 'i4', ('river'))
river.long_name = 'river runoff identification number'
out.variables['river'][:] = range(1,len(i)+1)

flag = out.createVariable('river_sign', 'f8', ('river'))
flag.long_name = 'river directional sign'
out.variables['river_sign'][:] = sign

xi = out.createVariable('river_Xposition', 'i4', ('river'))
xi.long_name = 'river XI-position at RHO-points'
xi.valid_min = 1
xi.valid_max = 793    # WARNING - hardcoded Lm+1
out.variables['river_Xposition'][:] = i

eta = out.createVariable('river_Eposition', 'i4', ('river'))
eta.long_name = 'river ETA-position at RHO-points'
eta.valid_min = 1
eta.valid_max = 361    # WARNING - hardcoded Mm+1
out.variables['river_Eposition'][:] = j

dirs = out.createVariable('river_direction', 'i4', ('river'))
dirs.long_name = 'river runoff direction'
out.variables['river_direction'][:] = dir

vshape = out.createVariable('river_Vshape', 'f8', ('s_rho', 'river'))
vshape.long_name = 'river runoff mass transport vertical profile'

trans = out.createVariable('river_transport', 'f8', ('river_time', 'river'))
trans.long_name = 'river runoff vertically integrated mass transport'
trans.units = 'meter3 second-1'
trans.time = 'river_time'

#temp = out.createVariable('river_temp', 'f8', ('river_time', 's_rho', 'river'))
#temp.long_name = 'river runoff potential temperature'
#temp.units = 'Celsius'
#temp.time = 'river_time'
#
#salt = out.createVariable('river_salt', 'f8', ('river_time', 's_rho', 'river'))
#salt.long_name = 'river runoff salinity'
#salt.time = 'river_time'

out.close()
