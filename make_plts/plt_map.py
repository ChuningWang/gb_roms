""" Plot Glacier Bay map """

import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from ocean_toolbox import ctd
from matplotlib.mlab import griddata
import cmocean
import numpy as np
import netCDF4 as nc
import pyroms
import sys

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

# load CTD data
info = {'data_dir': in_dir + 'ctd_raw/',
        'file_dir': in_dir,
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
       }

c = ctd.ctd(info)
c()
lat_ctd = c.data_info['lat']
lon_ctd = c.data_info['lon']

# Read grid
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
z = grd.vgrid.h
msk = grd.hgrid.mask_rho
Mp, Np = msk.shape

plt.close()
fig = plt.figure()

lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

# transect
a0 = []
fh = open('../data/a0.txt')
csvr = csv.reader(fh, delimiter=',')
for line in csvr:
    a0.append(line)
fh.close()
a0 = np.array(a0)
a0 = a0.astype(float)

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f')

# m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='lightgrey')
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.04), labels=[0,0,1,1], fontsize=3, linewidth=.1)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.02), labels=[1,1,0,0], fontsize=3, linewidth=.1)
# mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.2), labels=[0,0,0,1], fontsize=10, linewidth=.1)
# pr = m.drawparallels(np.arange(lat_min, lat_max, 0.1), labels=[1,0,0,0], fontsize=10, linewidth=.1)
setlabelrot(mr,-90)
setlabelrot(pr,0)

x, y = m(lon, lat)
m.contourf(x, y, z, cmap=cmocean.cm.deep)
plt.clim(0, 500)
plt.colorbar()

# plot CTD stations
x2, y2 = m(lon_ctd, lat_ctd)
m.plot(x2, y2, '^k', markersize=.5)

# plot transect
xa0, ya0 = m(a0[:, 0], a0[:, 1])
m.plot(xa0, ya0, '--ok', linewidth=.3, markersize=.3)

# plot boxes
box1 = np.array([[-136.10, 58.40],
                 [-135.975, 58.3125],
                 [-135.925, 58.3125],
                 [-135.875, 58.40],
                 [-135.50, 58.90],
                 [-136.00, 59.20],
                 [-137.25, 59.20],
                 [-137.25, 58.75]])

box2 = np.array([[-136.10, 58.40],
                 [-135.975, 58.3125],
                 [-135.925, 58.3125],
                 [-135.875, 58.40],
                 [-135.50, 58.55],
                 [-135.00, 58.25],
                 [-135.00, 57.90],
                 [-136.75, 57.90],
                 [-136.75, 58.60]])

plt.savefig(out_dir + 'figs/'+grd1+'_map.png', format='png', dpi=600)
plt.close()
