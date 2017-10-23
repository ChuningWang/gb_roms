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
t0 = 10.
s0 = 0.

# Select time range (days since 1900)
t_base = datetime(1900, 01, 01)
t_ini = datetime(my_year, 01, 01)
t_end = datetime(my_year+1, 01, 01)
# t_ini = datetime(my_year, 06, 01)
# t_end = datetime(my_year, 06, 03)

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
coast = fh.variables['coast'][:]
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
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

# get littoral (here 1 cells wide, without diagonals)
print 'get littoral points'

lit = pyroms_toolbox.get_littoral2(msk)
coord = [[lit[0][i], lit[1][i]] for i in range(len(lit[0]))]
# check for repeated elements
coord2 = []
for i in coord:
    if i not in coord2:
        coord2.append(i)

coord = np.array(coord2)
idx = coord[:, 0]
idy = coord[:, 1]
nrivers = len(coord)

latc = lat_grd[idx, idy]
lonc = lon_grd[idx, idy]
mskc = np.zeros(msk.shape)
mskc[idx, idy] = 1

# get remap river runoff
# find out which discharge point each grid point belongs to
cc = np.zeros((Mp, Lp))
for i in range(Mp):
    for j in range(Lp):
        dis = (latc - lat_grd[i, j])**2 + (lonc - lon_grd[i, j])**2
        cc[i, j] = np.argmin(dis)

# initiate
runoff_raw_nc = np.zeros((nt, Mp, Lp))
runoff_nc = np.zeros((nt, Mp, Lp))

for t in range(nt):
    # remap river runoff
    print 'Remapping runoff for time %f' %time[t]
    # conservative horizontal interpolation using scrip
    runoff_raw = pyroms.remapping.remap(data[t, :, :], wts_file, spval=spval)

    runoff = np.zeros((Mp, Lp))
    runoff2 = np.zeros((Mp, Lp))
    runoff = np.ma.masked_where(mskc == 0, runoff)

    # converge runoff to coastal cells
    for i in range(nrivers):
        runoff[idx[i], idy[i]] = np.sum(runoff_raw[cc == i])

    # spread the runoff
    for i in range(nrivers):
        xmin = np.max((idx[i]-rspread, 0))
        xmax = np.min((idx[i]+rspread+1, Mp))
        ymin = np.max((idy[i]-rspread, 0))
        ymax = np.min((idy[i]+rspread+1, Lp))
        runoff_nc[t, idx[i], idy[i]] = np.mean(runoff[xmin:xmax, ymin:ymax])

# ---------------------------------------------------------------------
# convert runoff to transport
# loop through littoral points to find Epos, Xpos
idx2 = []
idy2 = []
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
            idx2.append(idx[i])
            idy2.append(idy[i])
            epos.append(idx[i])
            xpos.append(idy[i])
            sign.append(1)
            rdir.append(1)
            count[ix, iy] = count[ix, iy] + 1
    if snwe[1]:
        if msk[idx[i]+1, idy[i]] == 0:
            idx2.append(idx[i])
            idy2.append(idy[i])
            epos.append(idx[i]+1)
            xpos.append(idy[i])
            sign.append(-1)
            rdir.append(1)
            count[ix, iy] = count[ix, iy] + 1
    if snwe[2]:
        if msk[idx[i], idy[i]-1] == 0:
            idx2.append(idx[i])
            idy2.append(idy[i])
            epos.append(idx[i])
            xpos.append(idy[i])
            sign.append(1)
            rdir.append(0)
            count[ix, iy] = count[ix, iy] + 1
    if snwe[3]:
        if msk[idx[i], idy[i]+1] == 0:
            idx2.append(idx[i])
            idy2.append(idy[i])
            epos.append(idx[i])
            xpos.append(idy[i]+1)
            sign.append(-1)
            rdir.append(0)
            count[ix, iy] = count[ix, iy] + 1

# summarize river info
nrivers2 = len(epos)
epos = np.array(epos)
xpos = np.array(xpos)
rdir = np.array(rdir)
sign = np.array(sign)
trs = np.zeros((nt, nrivers2))

for t in range(nt):
    print 'Remapping rivers for time %f' %time[t]
    for i in range(nrivers2):
        trs[t, i] = runoff_nc[t, idx2[i], idy2[i]]/count[idx2[i], idy2[i]]

# ---------------------------------------------------------------------
# scale the total discharge to original Hill product
def get_discharge_avgbox(t, lat, lon, trs, coast, box):
    ''' sum up discharge in a region. '''

    discharge = trs.copy()
    from matplotlib import path
    # Get points in boxes
    hydro_box = np.ones(lon.shape)*(-1)
    p0 = path.Path(box)

    shp = coast.shape
    if len(shp) == 2:
        for i in range(lon.shape[0]):
            for j in range(lon.shape[1]):
                if p0.contains_points([(lon[i, j], lat[i, j])]):
                    hydro_box[i, j] = 0
    elif len(shp) == 1:
        for i in range(len(lon)):
            if p0.contains_points([(lon[i], lat[i])]):
                hydro_box[i] = 0

    # Find coastal cells
    hydro_box[coast.mask] = -1

    print 'Sum up data in box...'
    d = np.empty(t.size)
    d[:] = np.NaN
    for i in range(t.size):
        if len(shp) == 2:
            d0 = discharge[i, :, :]
        elif len(shp) == 1:
            d0 = discharge[i, :]
        d0[d0 <= 0] = np.NaN
        d[i] = np.nansum(d0[hydro_box == 0])
    return d

box = np.array([[-137.40, 59.10],
                [-137.00, 58.50],
                [-136.55, 58.30],
                [-136.40, 58.15],
                [-136.00, 57.95],
                [-135.00, 58.05],
                [-136.10, 59.35]])

# calculate total discharge
coast_grd = np.ones(nrivers2)
coast_grd = np.ma.masked_invalid(coast_grd)
d = get_discharge_avgbox(time, lat_grd[idx2, idy2], lon_grd[idx2, idy2], trs, coast_grd, box)
# Hill product total runoff
d2 = get_discharge_avgbox(time, lat, lon, data, coast, box)
rr = np.nanmean(d2/d)
trs = trs*rr

# set discharge direction (sign)
trs = trs*sign

# set temp and salt
river_temp = t0*np.ones(time.shape)
river_salt = s0*np.ones(time.shape)

# save file
# create file with all the objects
fout = nc.Dataset(out_file, 'w', format='NETCDF3_64BIT')
fout.type = 'ROMS RIVERS file'
fout.title = 'Glacier Bay'
fout.source = 'David Hill and Jordan Beamer'

fout.createDimension('river_time', None)
fout.createDimension('river', nrivers2)
fout.createDimension('s_rho', grd.vgrid.N)

times = fout.createVariable('river_time', 'f8', ('river_time'))
times.units = 'days since 1900-01-01 00:00:00'
times.long_name = 'river runoff time'
fout.variables['river_time'][:] = time

river = fout.createVariable('river', 'i4', ('river'))
river.long_name = 'river runoff identification number'
fout.variables['river'][:] = range(1, nrivers2 + 1)

xi = fout.createVariable('river_Xposition', 'i4', ('river'))
xi.long_name = 'river XI-position at RHO-points'
fout.variables['river_Xposition'][:] = xpos

eta = fout.createVariable('river_Eposition', 'i4', ('river'))
eta.long_name = 'river ETA-position at RHO-points'
fout.variables['river_Eposition'][:] = epos

dirs = fout.createVariable('river_direction', 'i4', ('river'))
dirs.long_name = 'river runoff direction'
fout.variables['river_direction'][:] = rdir

flag = fout.createVariable('river_sign', 'f8', ('river'))
flag.long_name = 'river directional sign'
fout.variables['river_sign'][:] = sign

trans = fout.createVariable('river_transport', 'f8', ('river_time', 'river'))
trans.long_name = 'river runoff vertically integrated mass transport'
trans.units = 'meter3 second-1'
trans.time = 'river_time'
fout.variables['river_transport'][:] = trs

vshape = fout.createVariable('river_Vshape', 'f8', ('s_rho', 'river'))
vshape.long_name = 'river runoff mass transport vertical profile'

temp = fout.createVariable('river_temp', 'f8', ('river_time'))
temp.long_name = 'river runoff potential temperature'
temp.units = 'Celsius'
temp.time = 'river_time'
fout.variables['river_temp'][:] = river_temp

salt = fout.createVariable('river_salt', 'f8', ('river_time'))
salt.long_name = 'river runoff salinity'
salt.time = 'river_time'
fout.variables['river_salt'][:] = river_salt

fout.close()
