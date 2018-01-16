import netCDF4 as nc
import numpy as np
import pyroms
from matplotlib import path

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

box = np.array([[  -136.9,     58.95],
                 [  -137.1,     58.95],
                 [  -137.1,     59.15],
                 [  -136.9,     59.15]])

fh = nc.Dataset(out_dir + 'frc/GlacierBay_lr_rivers_2008_Hill_ana2.nc', 'a')
t = fh.variables['river_time'][:]
epos = fh.variables['river_Eposition'][:]
xpos = fh.variables['river_Xposition'][:]
trs = fh.variables['river_transport'][:]
temp = fh.variables['river_temp'][:]
# fh.close

p0 = path.Path(box)
for i in range(len(epos)):
    eta = epos[i]
    xi = xpos[i]
    if p0.contains_points([(lon[eta, xi], lat[eta, xi])]):
        print eta, xi
        trs[:, i] = trs[:, i]*10

fh.variables['river_transport'][:] = trs
fh.close()
