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

my_year = 2008

# Select time range (days since 1900)
t_base = datetime(1900, 01, 01)
t_ini = datetime(my_year, 01, 01)
t_end = datetime(my_year+1, 01, 01)

# load 2-dimentional interannual discharge data 
print 'Load interannual discharge data'
fh = nc.Dataset(dis_file, 'r')
time = fh.variables['t'][:]
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]
t1 = (t_ini-t_base).days
t2 = (t_end-t_base).days
msk = (time>=t1) & (time<=t2)
time = time[msk]
data = fh.variables['discharge'][msk, :, :]
fh.close()

# the coordinate of Hill's grid is slightly off - correct it here
lon = lon-0.010
lat = lat-0.005

# load Glacier Bay grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

# define some variables
wts_file = script_dir + 'remap_weights_runoff_to_' + grd.name + '_conservative_nomask.nc'
nt = time.shape[0]
Mp, Lp = grd.hgrid.mask_rho.shape
spval = -1e30
runoff_raw = np.zeros((Mp,Lp))
runoff = np.zeros((Mp,Lp))
rspread = 3

# specify output file
out_dir = in_dir
tag = 'Hill'
out_file = out_dir + grd.name + '_runoff_' + str(rspread) + '_' + str(my_year) + '_' + tag + '.nc'

# get littoral (here 1 cells wide, with diagonals)
print 'get littoral points'
width = 1
idx = []
idy = []
maskl = grd.hgrid.mask_rho.copy()

for w in range(width):
    lit = pyroms_toolbox.get_littoral2(maskl)
    idx.extend(lit[0])
    idy.extend(lit[1])
    maskl[lit] = 0

littoral_idx = (np.array(idx), np.array(idy))
maskl = np.zeros(grd.hgrid.mask_rho.shape)
maskl[littoral_idx] = 1
mask_idx = np.where(grd.hgrid.mask_rho == 0)

# initiate
runoff_raw_nc = np.zeros((nt, grd.hgrid.mask_rho.shape[0], grd.hgrid.mask_rho.shape[1]))
runoff_spread_nc = np.zeros((nt, grd.hgrid.mask_rho.shape[0], grd.hgrid.mask_rho.shape[1]))

# # use parallel programming here
# from joblib import Parallel, delayed
# 
# def remap_runoff(time, data, wts_file, spval, mask_idx):
#     print 'Remapping runoff for time %f' %time
#     # conservative horizontal interpolation using scrip
#     runoff_raw = pyroms.remapping.remap(data, wts_file, spval=spval)
#     runoff_raw[mask_idx] = spval
#     return runoff_raw
# 
# def spread_runoff(time, runoff_raw, spval, littoral_idx, maskl, grd, rspread, mask_idx):
#     print 'Spreading runoff for time %f' %time
#     idx = np.where(runoff_raw != 0)
#     runoff = pyroms_toolbox.move_runoff(runoff_raw, \
#                   np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
#                   grd.hgrid.x_rho, grd.hgrid.y_rho, grd.hgrid.dx, grd.hgrid.dy)
# 
#     # spread the runoff within the littoral band
#     runoff_spread = np.zeros((runoff.shape))
#     idx = np.where(runoff != 0)
#     for p in range(np.size(idx,1)):
#         j = range(max(0,idx[0][p]-rspread), min(Mp-1,idx[0][p]+rspread+1))
#         i = range(max(0,idx[1][p]-rspread), min(Lp-1,idx[1][p]+rspread+1))
#         ji = np.meshgrid(j,i)
#         sidx = np.where(maskl[ji] == 1)
#         nbpt = np.size(sidx) / 2
#         rpt = runoff[idx[0][p],idx[1][p]] * grd.hgrid.dx[idx[0][p],idx[1][p]] * grd.hgrid.dy[idx[0][p],idx[1][p]]
#         rpt = rpt / nbpt
#         for pa in range(nbpt):
#             pai = sidx[0][pa]+ji[1].min()
#             paj = sidx[1][pa]+ji[0].min()
#             runoff_spread[paj, pai] = runoff_spread[paj, pai]+rpt/(grd.hgrid.dx[paj, pai]*grd.hgrid.dy[paj, pai])
#     runoff_spread[mask_idx] = spval
#     return runoff_spread
# 
# runoff_raw_list = Parallel(n_jobs=8)(delayed(remap_runoff)(time[i], data[i, :, :], wts_file, spval, mask_idx) for i in range(nt))
# for i in range(nt):
#     runoff_raw_nc[i, :, :] = runoff_raw_list[i]
# 
# runoff_spread_list = Parallel(n_jobs=8)(delayed(spread_runoff)(time[i], runoff_raw_nc[i, :, :], spval, littoral_idx, maskl, grd, rspread, mask_idx) for i in range(nt))
# for i in range(nt):
#     runoff_spread_nc[i, :, :] = runoff_spread_list[i]

for t in range(nt):
    print 'Remapping runoff for time %f' %time[t]
    # conservative horizontal interpolation using scrip
    runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
                                           spval=spval)
    idx = np.where(runoff_raw != 0)
    runoff = pyroms_toolbox.move_runoff(runoff_raw, \
                  np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
                  grd.hgrid.x_rho, grd.hgrid.y_rho, grd.hgrid.dx, grd.hgrid.dy)

    # spread the runoff within the littoral band
    runoff_spread = np.zeros((Mp,Lp))
    idx = np.where(runoff != 0)
    for p in range(np.size(idx,1)):
        j = range(max(0,idx[0][p]-rspread), min(Mp-1,idx[0][p]+rspread+1))
        i = range(max(0,idx[1][p]-rspread), min(Lp-1,idx[1][p]+rspread+1))
        ji = np.meshgrid(j,i)
        sidx = np.where(maskl[ji] == 1)
        nbpt = np.size(sidx) / 2
        rpt = runoff[idx[0][p],idx[1][p]] * grd.hgrid.dx[idx[0][p],idx[1][p]] * grd.hgrid.dy[idx[0][p],idx[1][p]]
        rpt = rpt / nbpt
        for pa in range(nbpt):
            pai = sidx[0][pa] + ji[1].min()
            paj = sidx[1][pa] + ji[0].min()
            runoff_spread[paj, pai] = runoff_spread[paj, pai] + \
                           rpt / (grd.hgrid.dx[paj, pai] * grd.hgrid.dy[paj, pai])


    # check spval
    runoff_raw[mask_idx] = spval
    runoff_spread[mask_idx] = spval

    # write data in destination file
    runoff_raw_nc[t, :, :] = runoff_raw
    runoff_spread_nc[t, :, :] = runoff_spread

# mask invalid values
runoff_raw_nc = np.ma.masked_where(runoff_raw_nc==spval, runoff_raw_nc)
runoff_spread_nc = np.ma.masked_where(runoff_spread_nc==spval, runoff_spread_nc)

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
fh.variables['Runoff'][:] = runoff_spread_nc

fh.close()
