import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

fh = nc.Dataset(out_dir + 'frc/GlacierBay_lr_rivers_clim_Hill.nc', 'r')
t = fh.variables['river_time'][:]
epos = fh.variables['river_Eposition'][:]
xpos = fh.variables['river_Xposition'][:]
trs = fh.variables['river_transport'][:]
temp = fh.variables['river_temp'][:]
fh.close

latt = np.zeros(xpos.shape)
lont = np.zeros(xpos.shape)

for i in range(len(xpos)):
    latt[i] = lat[epos[i], xpos[i]]
    lont[i] = lon[epos[i], xpos[i]]

coast = np.zeros(lont.shape)
coast = np.ma.masked_invalid(coast)

# fh = nc.Dataset(data_dir + 'gb_discharge.nc', 'r')
# t_h = fh.variables['t'][:]
# lat_h = fh.variables['lat'][:]
# lon_h = fh.variables['lon'][:]
# coast_h = fh.variables['coast'][:]
# trs_h = fh.variables['discharge'][:]
# fh.close()
# 
# mskt = (t_h >= t[0]) & (t_h <= t[-1])
# t_h = t_h[mskt]
# trs_h = trs_h[mskt, :, :]


def get_discharge_avgbox(t, lat, lon, discharge, coast, box):
    ''' sum up discharge in a region. '''

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

box0 = np.array([[-136.65, 58.65],
                 [-137.30, 58.65],
                 [-137.30, 59.15],
                 [-136.40, 59.15],
                 [-136.20, 58.75]])

box1 = np.array([[-136.20, 58.75],
                 [-136.40, 59.15],
                 [-135.70, 59.15],
                 [-135.60, 58.75]])

box2 = np.array([[-136.65, 58.55],
                 [-136.65, 58.65],
                 [-136.20, 58.75],
                 [-135.60, 58.75],
                 [-135.60, 58.55]])

box3 = np.array([[-136.65, 58.55],
                 [-135.60, 58.55],
                 [-135.60, 58.50],
                 [-136.00, 58.35]])

box4 = np.array([[-136.65, 58.55],
                 [-136.75, 58.25],
                 [-136.50, 57.95],
                 [-135.20, 57.95],
                 [-135.40, 58.55],
                 [-135.60, 58.55],
                 [-135.60, 58.50],
                 [-136.00, 58.35]])

d0 = get_discharge_avgbox(t, latt, lont, trs, coast, box0)
d1 = get_discharge_avgbox(t, latt, lont, trs, coast, box1)
d2 = get_discharge_avgbox(t, latt, lont, trs, coast, box2)
d3 = get_discharge_avgbox(t, latt, lont, trs, coast, box3)
d4 = get_discharge_avgbox(t, latt, lont, trs, coast, box4)

plt.plot(d0, 'k')
plt.plot(d1, 'r')
plt.plot(d2, 'b')
plt.plot(d3, 'g')
plt.plot(d4, 'y')
# plt.plot(d5, 'r')
# plt.plot(d6, 'k')

# plt.plot(d)
# plt.plot(d_h)
