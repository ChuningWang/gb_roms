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

fh = nc.Dataset(out_dir + 'frc/GlacierBay_lr_rivers_2008_Hill_ana2.nc', 'r')
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

box1 = np.array([[  -137.7,     59.1  ],
                 [  -137.7,     59.25 ],
                 [  -136.6,     59.25 ],
                 [  -136.15,    58.725],
                 [  -136.65,    58.725],
                 [  -137.15,    58.65 ]])

box2 = np.array([[  -136.6,     59.25 ],
                 [  -136.3,     59.25 ],
                 [  -135.7,     59.0  ],
                 [  -135.7,     58.725],
                 [  -136.,      58.825],
                 [  -136.15,    58.725],
                 [  -136.3,     58.9  ]])

box3 = np.array([[  -136.15,    58.725],
                 [  -136.,      58.825],
                 [  -135.7,     58.725],
                 [  -135.65,    58.55 ],
                 [  -136.65,    58.55 ],
                 [  -136.65,    58.725]])

box4 = np.array([[  -136.65,    58.55 ],
                 [  -135.65,    58.55 ],
                 [  -135.65,    58.45 ],
                 [  -136.0,     58.375]])

box5 = np.array([[  -136.9,     58.95],
                 [  -137.1,     58.95],
                 [  -137.1,     59.15],
                 [  -136.9,     59.15]])

box6 = np.array([[  -136.4,     58.90],
                 [  -136.8,     58.90],
                 [  -136.8,     59.15],
                 [  -136.4,     59.15]])

box = np.array([[-137.40, 59.10],
                [-137.00, 58.50],
                [-136.55, 58.30],
                [-136.40, 58.15],
                [-136.00, 57.95],
                [-135.00, 58.05],
                [-136.10, 59.35]])

# d1 = get_discharge_avgbox(t, latt, lont, trs, coast, box1)
# d2 = get_discharge_avgbox(t, latt, lont, trs, coast, box2)
# d3 = get_discharge_avgbox(t, latt, lont, trs, coast, box3)
# d4 = get_discharge_avgbox(t, latt, lont, trs, coast, box4)
# d = get_discharge_avgbox(t, latt, lont, trs, coast, box)
d5 = get_discharge_avgbox(t, latt, lont, trs, coast, box5)
d6 = get_discharge_avgbox(t, latt, lont, trs, coast, box6)

# d_h = get_discharge_avgbox(t_h, lat_h, lon_h, trs_h, coast_h, box)

# plt.plot(d1)
# plt.plot(d2)
# plt.plot(d3)
# plt.plot(d4)
plt.plot(d5, 'r')
plt.plot(d6, 'k')

# plt.plot(d)
# plt.plot(d_h)
