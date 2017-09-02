import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_USGS'
grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

fh = nc.Dataset(out_dir+'frc/GlacierBay_usgs_rivers_2008_Hill.nc', 'r')
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

def get_discharge_avgbox(t, lat, lon, discharge, box):
    from matplotlib import path

    # ------------------------------------------------------------------------------------------------------
    # Get points in boxes
    hydro_box = np.ones(lon.shape)*(-1)
    p0 = path.Path(box)

    for i in range(len(lon)):
        if p0.contains_points([(lon[i], lat[i])]):
            hydro_box[i] = 0

    # ------------------------------------------------------------------------------------------------------
    print 'Sum up data in box...'
    d = np.zeros(t.size)
    for i in range(t.size):
        d0 = discharge[i, :]
        d[i] = np.nansum(d0[hydro_box==0])

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

d1 = get_discharge_avgbox(t, latt, lont, trs, box1)
d2 = get_discharge_avgbox(t, latt, lont, trs, box2)
d3 = get_discharge_avgbox(t, latt, lont, trs, box3)
d4 = get_discharge_avgbox(t, latt, lont, trs, box4)

plt.plot(d1)
plt.plot(d2)
plt.plot(d3)
plt.plot(d4)
