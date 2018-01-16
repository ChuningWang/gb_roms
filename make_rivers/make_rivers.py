import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from matplotlib import path
import sys

import pyroms
import pyroms_toolbox
from ocean_toolbox import ctd

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
tag = 'Hill'
my_year = 2009
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

# ---------------------------------------------------------------------
# load Glacier Bay grid object
grd = pyroms.grid.get_ROMS_grid(grd1)
lat_grd = grd.hgrid.lat_rho
lon_grd = grd.hgrid.lon_rho
msk = grd.hgrid.mask_rho

out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag + '.nc'

# ---------------------------------------------------------------------
# read in river location, info.
idx = []
idy = []
sign = []
rdir = []
river = []
# read in coast cells
river_num = 0
fin = open('river_cells.txt', 'r')
for line in fin:
    llist = line.rstrip('\n').split(',')
    if int(llist[0]) == -1:
        river_num += 1
        irdir = int(llist[1])
        isign = int(llist[2])
    elif int(llist[0]) == -10:
        break
    else:
        river.append(river_num)
        if isign == 1:
            idx.append(int(llist[0]))
            idy.append(int(llist[1]))
        else:
            if irdir == 1:
                idx.append(int(llist[0])+1)
                idy.append(int(llist[1]))
            else:
                idx.append(int(llist[0]))
                idy.append(int(llist[1])+1)

        sign.append(isign)
        rdir.append(irdir)

river = np.array(river)
idx = np.array(idx)
idy = np.array(idy)
sign = np.array(sign)
rdir = np.array(rdir)
Nriver_grid = len(river)
Nriver = river.max()
lat_river = np.zeros(Nriver)
lon_river = np.zeros(Nriver)
for i in range(Nriver):
    ir = river == i+1
    lat_river[i] = lat_grd[idx[ir], idy[ir]].mean()
    lon_river[i] = lon_grd[idx[ir], idy[ir]].mean()

# ---------------------------------------------------------------------
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

nt = len(time)

# only use data in GB region
box = np.array([[-137.40, 59.10],
                [-137.00, 58.50],
                [-136.55, 58.30],
                [-136.40, 58.15],
                [-136.00, 57.95],
                [-135.00, 58.05],
                [-136.10, 59.35]])

p0 = path.Path(box)

xp, yp = coast.shape
for i in range(xp):
    for j in range(yp):
        if not p0.contains_point((lon[i, j], lat[i, j])):
            coast[i, j] = 0

lat = lat[coast == 1]
lon = lon[coast == 1]
data = data[:, coast == 1]
Nriver_hill = len(lat)

# the coordinate of Hill's grid is slightly off - correct it here
lon = lon-0.010
lat = lat-0.010

# ---------------------------------------------------------------------
# project Hill product to my grid
river_idx = np.zeros(Nriver_hill)
for i in range(Nriver_hill):
    dis = (lon[i] - lon_river)**2 + (lat[i] - lat_river)**2
    river_idx[i] = np.argmin(dis) + 1

# sum up data for each river
trs_river = np.zeros((nt, Nriver))
for i in range(Nriver):
    trs_river[:, i] = data[:, river_idx == i+1].sum(axis=1)

# spread discharge onto grid cells
trs = np.zeros((nt, Nriver_grid))
for i in range(Nriver_grid):
    ngrid = np.sum(river == river[i])
    trs[:, i] = trs_river[:, river[i]-1]/float(ngrid)

trs = trs*sign

# ---------------------------------------------------------------------
# generate analytical temp and salt
# Read in river temperature
info = {'data_dir': '/glade/p/work/chuning/data/ctd_raw/',
        'file_dir': '/glade/p/work/chuning/data/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['temp'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
       }

c = ctd.ctd(info)
c()

# process time
river_datetime = nc.num2date(time, 'days since 1900-01-01')
river_yearday = np.array([i.timetuple().tm_yday for i in river_datetime])

temp = c.climatology['temp'][:, :, 12]
temp = temp[:5, :].mean(axis=0)
ctd_time = c.climatology['time']
temp = np.concatenate((np.array([temp[-1]]), temp, np.array([temp[0]])))
ctd_time = np.concatenate((np.array([ctd_time[-1]-365.]), ctd_time, np.array([ctd_time[0]+365.])))

river_temp = np.interp(river_yearday, ctd_time, temp)
river_salt = np.zeros(river_temp.shape)

# ---------------------------------------------------------------------
# generate v_shape
v_shape = np.ones((grd.vgrid.N, Nriver_grid))/float(grd.vgrid.N)

# ---------------------------------------------------------------------
# save file
# create file with all the objects
fout = nc.Dataset(out_file, 'w', format='NETCDF3_64BIT')
fout.type = 'ROMS RIVERS file'
fout.title = 'Glacier Bay'
fout.source = 'David Hill and Jordan Beamer'

fout.createDimension('river_time', None)
fout.createDimension('river', Nriver_grid)
fout.createDimension('s_rho', grd.vgrid.N)

times = fout.createVariable('river_time', 'f8', ('river_time'))
times.units = 'days since 1900-01-01 00:00:00'
times.long_name = 'river runoff time'
fout.variables['river_time'][:] = time

rivers = fout.createVariable('river', 'i4', ('river'))
rivers.long_name = 'river runoff identification number'
fout.variables['river'][:] = river

eta = fout.createVariable('river_Eposition', 'i4', ('river'))
eta.long_name = 'river ETA-position at RHO-points'
fout.variables['river_Eposition'][:] = idx

xi = fout.createVariable('river_Xposition', 'i4', ('river'))
xi.long_name = 'river XI-position at RHO-points'
fout.variables['river_Xposition'][:] = idy

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
fout.variables['river_Vshape'][:] = v_shape

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
