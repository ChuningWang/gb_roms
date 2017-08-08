import numpy as np
import netCDF4 as netCDF
from datetime import datetime

import pyroms
import pyroms_toolbox

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

grd1 = 'GB_USGS'

my_year = 2008

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

# load 2-dimentional discharge data 
print 'Load discharge data'
tag = 'Hill'
in_file = in_dir + grd.name + '_runoff_' + str(my_year) + '_' + tag + '.nc'
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

nc_data = netCDF.Dataset(in_file, 'r')
nc_rivers = netCDF.Dataset(out_file, 'a')
data = nc_data.variables['Runoff'][:]
time = nc_data.variables['time'][:]
nc_data.close()
sign = nc_rivers.variables['river_sign'][:]
xi = nc_rivers.variables['river_Xposition'][:]
eta = nc_rivers.variables['river_Eposition'][:]
dir = nc_rivers.variables['river_direction'][:]

# define some variables
nt = data.shape[0]
Nr = sign.shape[0]
Mp, Lp = grd.hgrid.mask_rho.shape
runoff = np.zeros((Nr))
count = np.zeros(grd.hgrid.mask_rho.shape, dtype=np.int32)

# from a Python forum - create an array of lists
filler = np.frompyfunc(lambda x: list(), 1, 1)
rivers = np.empty((Mp, Lp), dtype=np.object)
filler(rivers, rivers)

# pdb.set_trace()

for k in range(Nr):
    if (sign[k]==1):
        count[eta[k],xi[k]] += 1
	rivers[eta[k],xi[k]].append(k)
    elif (sign[k]==-1 and dir[k]==0):
        count[eta[k],xi[k]-1] += 1
	rivers[eta[k],xi[k]-1].append(k)
    elif (sign[k]==-1 and dir[k]==1):
        count[eta[k]-1,xi[k]] += 1
	rivers[eta[k]-1,xi[k]].append(k)

nct=0
for t in range(nt):
    print 'Remapping runoff for time %f' %time[t]
    for j in range(Mp):
        for i in range(Lp):
	    for n in range(count[j,i]):
	        frac = 1.0/count[j,i]
		k = rivers[j,i][n]
		runoff[k] = frac*data[t,j,i]
	        
    # write data in destination file
    nc_rivers.variables['river_transport'][nct] = runoff
    nc_rivers.variables['river_time'][nct] = time[nct]
    nct = nct + 1

# close netcdf file
nc_rivers.close()

