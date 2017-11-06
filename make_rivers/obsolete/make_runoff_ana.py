import numpy as np
import netCDF4 as nc
from datetime import datetime
import sys

import pyroms
import pyroms_toolbox

# You need to edit and run compute_daitren_remap_weights.py first
# then edit and run make_runoff.py to generate the runoff file.
# In make_runoff.py you can tune the number of cell defining the
# littoral band where you will have a runoff value (look for "width"
# variable) and the area over which the runoff is spread in order
# homogenize and to avoid large runoff value (variable "rspread").

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

script_dir = home_dir + 'git/gb_roms/make_rivers/'
dis_file = in_dir + 'gb_discharge.nc'

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

# my inputs
my_year = 2008
rspread = 2
spval = -1e30

# Select time range (days since 1900)
t_base = datetime(1900, 01, 01)
t_ini = datetime(my_year, 01, 01)
t_end = datetime(my_year+1, 01, 01)
# t_ini = datetime(my_year, 06, 01)
# t_end = datetime(my_year, 06, 05)

# load Glacier Bay grid object
grd = pyroms.grid.get_ROMS_grid(grd1)
lat_grd = grd.hgrid.lat_rho
lon_grd = grd.hgrid.lon_rho
msk = grd.hgrid.mask_rho
Mp, Lp = msk.shape

# load 2-dimentional interannual discharge data 
print 'Load interannual discharge data'
fh = nc.Dataset(dis_file, 'r')
time = fh.variables['t'][:]
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]
t1 = (t_ini-t_base).days
t2 = (t_end-t_base).days
mskt = (time>=t1) & (time<=t2)
time = time[mskt]
data = fh.variables['discharge'][mskt, :, :]
fh.close()

# the coordinate of Hill's grid is slightly off - correct it here
lon = lon-0.010
lat = lat-0.010

# define some variables
wts_file = script_dir + 'remap_weights_runoff_to_' + grd.name + '_conservative_nomask.nc'
nt = len(time)

# specify output file
out_dir = in_dir
tag = 'Hill'
out_file = out_dir + grd.name + '_runoff_' + str(my_year) + '_' + tag + '.nc'

# get littoral (here 1 cells wide, with diagonals)
print 'get littoral points'
width = 1
idx = []
idy = []
mskc = msk.copy()

for w in range(width):
    lit = pyroms_toolbox.get_littoral2(mskc)
    idx.extend(lit[0])
    idy.extend(lit[1])

littoral_idx = (np.array(idx), np.array(idy))
latc = lat_grd[littoral_idx]
lonc = lon_grd[littoral_idx]
mskc = np.zeros(msk.shape)
mskc[littoral_idx] = 1

# find out which discharge point each grid point belongs to
cc = np.zeros(msk.shape)
for i in range(Mp):
    for j in range(Lp):
        dis = (latc - lat_grd[i, j])**2 + (lonc - lon_grd[i, j])**2
        cc[i, j] = np.argmin(dis)

# initiate
runoff_raw_nc = np.zeros((nt, Mp, Lp))
runoff_nc = np.zeros((nt, Mp, Lp))

# mask invalid values
# runoff_nc = np.ma.masked_where(np.tile(mskc == 0, (nt, 1, 1)), runoff_nc)

for t in range(nt):
    print 'Remapping runoff for time %f' %time[t]
    # conservative horizontal interpolation using scrip
    runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
                                           spval=spval)
    runoff = np.zeros(runoff_raw.shape)
    runoff = np.ma.masked_where(mskc == 0, runoff)

    for i in range(len(littoral_idx[0])):
        runoff[littoral_idx[0][i], littoral_idx[1][i]] = np.sum(runoff_raw[cc == i])

    # spread the runoff
    for i in range(len(littoral_idx[0])):
        xx, yy = littoral_idx[0][i], littoral_idx[1][i]
        xmin = np.max((xx-rspread, 0))
        xmax = np.min((xx+rspread+1, Mp))
        ymin = np.max((yy-rspread, 0))
        ymax = np.min((yy+rspread+1, Lp))
        runoff_nc[t, xx, yy] = np.mean(runoff[xmin:xmax, ymin:ymax])

    # write data in destination file
    runoff_raw_nc[t, :, :] = runoff_raw

# create runoff file
print 'create runoff file'
fh = nc.Dataset(out_file, 'w', format='NETCDF3_64BIT')
fh.Description = 'Hill & Beamer Glacier Bay freshwater runoff, 1980-2007'
fh.Author = 'make_runoff.py'
fh.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
fh.title = 'Hill Glacier Bay runoff'

# creat dimensions and variables
fh.createDimension('eta_rho', grd.hgrid.mask_rho.shape[0])
fh.createDimension('xi_rho', grd.hgrid.mask_rho.shape[1])
fh.createDimension('runoff_time', None)

print 'create lon-lat variables'
fh.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
fh.variables['lon_rho'].long_name = 'longitude of RHO-points'
fh.variables['lon_rho'].units = 'degree_east'
fh.variables['lon_rho'].field = 'lon_rho, scalar'
fh.variables['lon_rho'][:] = grd.hgrid.lon_rho

fh.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
fh.variables['lat_rho'].long_name = 'latitude of RHO-points'
fh.variables['lat_rho'].units = 'degree_north'
fh.variables['lat_rho'].field = 'lat_rho, scalar'
fh.variables['lat_rho'][:] = grd.hgrid.lat_rho

print 'create time variables'
fh.createVariable('time', 'f8', ('runoff_time'))
fh.variables['time'].long_name = 'time'
fh.variables['time'].units = 'days since 1900-01-01 00:00:00'
fh.variables['time'][:] = time

print 'create runoff variables'
fh.createVariable('Runoff_raw', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
fh.variables['Runoff_raw'].long_name = 'Hill River Runoff raw'
fh.variables['Runoff_raw'].missing_value = str(spval)
fh.variables['Runoff_raw'].units = 'kg/s/m^2'
fh.variables['Runoff_raw'][:] = runoff_raw_nc

fh.createVariable('Runoff', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
fh.variables['Runoff'].long_name = 'Hill River Runoff'
fh.variables['Runoff'].missing_value = str(spval)
fh.variables['Runoff'].units = 'm^3/s'
fh.variables['Runoff'][:] = runoff_nc

fh.close()
