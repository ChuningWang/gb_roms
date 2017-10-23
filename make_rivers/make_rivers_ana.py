import numpy as np
import netCDF4 as nc
from datetime import datetime
import sys

import pyroms
import pyroms_toolbox

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

dis_file = in_dir + 'gb_discharge.nc'
script_dir = home_dir + 'git/gb_roms/make_rivers/'
my_year = 2008
rspread = 2
spval = -1e30

# Select time range (days since 1900)
t_base = datetime(1900, 01, 01)
# t_ini = datetime(my_year, 01, 01)
# t_end = datetime(my_year+1, 01, 01)
t_ini = datetime(my_year, 06, 01)
t_end = datetime(my_year, 06, 03)

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

# remap weights file
wts_file = script_dir + 'remap_weights_runoff_to_' + grd.name + '_conservative_nomask.nc'

# number of total time entry
nt = len(time)

# specify output file
tag = 'Hill_ana'
out_file = out_dir + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

# get littoral (here 1 cells wide, with diagonals)
print 'get littoral points'

lit = pyroms_toolbox.get_littoral2(msk)
idx = lit[0]
idy = lit[1]
nrivers = len(idx)

latc = lat_grd[lit]
lonc = lon_grd[lit]
mskc = np.zeros(msk.shape)
mskc[lit] = 1

# loop through littoral points to find Epos, Xpos
epos = []
xpos = []
rdir = []
sign = []
count = np.zeros((Mp, Lp))
for i in range(nrivers):
    ix, iy = idx[i], idy[i]
    # for interior, look through 4 directions
    # for boundaries, only look at possible directions
    # first, identify if neighber cells are valid
    snwe = [ix != 0, ix != Mp-1, iy != 0, iy != Lp-1]

    # then identigy xpos, epos, dir and sign
    if snwe[0]:
        if msk[idx[i]-1, idy[i]] == 0:
            epos.append(idx[i])
            xpos.append(idy[i])
            sign.append(1)
            rdir.append(1)
    if snwe[1]:
        if msk[idx[i]+1, idy[i]] == 0:
            epos.append(idx[i]+1)
            xpos.append(idy[i])
            sign.append(-1)
            rdir.append(1)
    if snwe[2]:
        if msk[idx[i], idy[i]-1] == 0:
            epos.append(idx[i])
            xpos.append(idy[i])
            sign.append(1)
            rdir.append(0)
    if snwe[3]:
        if msk[idx[i], idy[i]+1] == 0:
            epos.append(idx[i])
            xpos.append(idy[i]+1)
            sign.append(-1)
            rdir.append(0)

    # ms = msk[idx[i]-1, idy[i]]
    # mn = msk[idx[i]+1, idy[i]]
    # mw = msk[idx[i], idy[i]-1]
    # me = msk[idx[i], idy[i]+1]
    # if ms == 0:
    #     epos.append(idx[i])
    #     xpos.append(idy[i])
    #     sign.append(1)
    #     rdir.append(1)
    # if mn == 0:
    #     epos.append(idx[i])
    #     xpos.append(idy[i]+1)
    #     sign.append(-1)
    #     rdir.append(1)
    # if mw == 0:
    #     epos.append(idx[i])
    #     xpos.append(idy[i])
    #     sign.append(1)
    #     rdir.append(0)
    # if me == 0:
    #     epos.append(idx[i]+1)
    #     xpos.append(idy[i])
    #     sign.append(-1)
    #     rdir.append(0)

# # find out which discharge point each grid point belongs to
# cc = np.zeros(msk.shape)
# for i in range(Mp):
#     for j in range(Lp):
#         dis = (latc - lat_grd[i, j])**2 + (lonc - lon_grd[i, j])**2
#         cc[i, j] = np.argmin(dis)
# 
# # initiate
# runoff_raw_nc = np.zeros((nt, Mp, Lp))
# runoff_nc = np.zeros((nt, Mp, Lp))
# 
# # mask invalid values
# # runoff_nc = np.ma.masked_where(np.tile(mskc == 0, (nt, 1, 1)), runoff_nc)
# 
# for t in range(nt):
#     print 'Remapping runoff for time %f' %time[t]
#     # conservative horizontal interpolation using scrip
#     runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
#                                            spval=spval)
#     runoff = np.zeros(runoff_raw.shape)
#     runoff = np.ma.masked_where(mskc == 0, runoff)
# 
#     for i in range(len(littoral_idx[0])):
#         runoff[littoral_idx[0][i], littoral_idx[1][i]] = np.sum(runoff_raw[cc == i])
# 
#     # spread the runoff
#     for i in range(len(littoral_idx[0])):
#         xx, yy = littoral_idx[0][i], littoral_idx[1][i]
#         xmin = np.max((xx-rspread, 0))
#         xmax = np.min((xx+rspread+1, Mp))
#         ymin = np.max((yy-rspread, 0))
#         ymax = np.min((yy+rspread+1, Lp))
#         runoff_nc[t, xx, yy] = np.mean(runoff[xmin:xmax, ymin:ymax])
# 
#     # write data in destination file
#     runoff_raw_nc[t, :, :] = runoff_raw
