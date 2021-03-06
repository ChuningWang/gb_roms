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

my_year = 2008
rspread = 3

# load GB grid object
grd = pyroms.grid.get_ROMS_grid(grd1)

# load 2-dimentional discharge data 
print 'Load discharge data'
tag1 = 'Hill'
tag2 = 'Hill'
in_file = in_dir + grd.name + '_runoff_' + str(rspread) + '_' + str(my_year) + '_' + tag1 + '.nc'
out_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_' + tag2 + '.nc'

nc_data = nc.Dataset(in_file, 'r')
nc_rivers = nc.Dataset(out_file, 'a')
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

for k in range(Nr):
    if (sign[k]==1):
        count[eta[k],xi[k]] += 1
        rivers[eta[k],xi[k]].append(k)
    elif (sign[k] == -1 and dir[k] == 0):
        count[eta[k],xi[k]-1] += 1
        rivers[eta[k],xi[k]-1].append(k)
    elif (sign[k] == -1 and dir[k] == 1):
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

nc_rivers.close()

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

# forcing file
fh = nc.Dataset(out_dir + 'frc/GlacierBay_lr_rivers_2008_Hill.nc', 'r')
t = fh.variables['river_time'][:]
epos = fh.variables['river_Eposition'][:]
xpos = fh.variables['river_Xposition'][:]
trs = fh.variables['river_transport'][:]
fh.close

# Hill product
fh = nc.Dataset(in_dir + 'gb_discharge.nc', 'r')
t_h = fh.variables['t'][:]
lat_h = fh.variables['lat'][:]
lon_h = fh.variables['lon'][:]
coast_h = fh.variables['coast'][:]
trs_h = fh.variables['discharge'][:]
fh.close()

mskt = (t_h >= t[0]) & (t_h <= t[-1])
t_h = t_h[mskt]
trs_h = trs_h[mskt, :, :]

lat = np.zeros(xpos.shape)
lon = np.zeros(xpos.shape)

for i in range(len(lat)):
    lat[i] = grd.hgrid.lat[epos[i], xpos[i]]
    lon[i] = grd.hgrid.lon[epos[i], xpos[i]]

coast = np.zeros(lon.shape)
coast = np.ma.masked_invalid(coast)

# calculate total discharge
d = get_discharge_avgbox(t, lat, lon, trs, coast, box)
d_h = get_discharge_avgbox(t_h, lat_h, lon_h, trs_h, coast_h, box)

rr = np.nanmean(d_h/d)
trs = trs*rr

# set discharge direction (sign)
trs = trs*sign

# overwrite river_transport
fh = nc.Dataset(out_dir + 'frc/GlacierBay_lr_rivers_2008_Hill.nc', 'r+')
fh.variables['river_transport'][:] = trs
fh.close()
