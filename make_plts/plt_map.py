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

# ------------------- data prep ------------------------------------------
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

box0 = np.array([[-136.65, 58.55],
                 [-137.30, 58.65],
                 [-137.30, 59.15],
                 [-135.70, 59.15],
                 [-135.60, 58.75],
                 [-135.60, 58.55]])

box1 = np.array([[-136.65, 58.55],
                 [-135.60, 58.55],
                 [-135.60, 58.50],
                 [-136.00, 58.35]])

box2 = np.array([[-136.65, 58.55],
                 [-136.75, 58.25],
                 [-136.50, 58.05],
                 [-135.20, 58.05],
                 [-135.40, 58.55],
                 [-135.60, 58.55],
                 [-135.60, 58.50],
                 [-136.00, 58.35]])

# load CTD data
info = {'data_dir': in_dir + 'ctd_raw/',
        'file_dir': in_dir,
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes'}

gb_ctd = ctd.ctd(info)
gb_ctd()
lat_ctd = gb_ctd.data_info['lat']
lon_ctd = gb_ctd.data_info['lon']

# Read grid
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
z = grd.vgrid.h
msk = grd.hgrid.mask_rho
Mp, Np = msk.shape

# along channel transect
a0 = []
fh = open('../data/a0.txt')
csvr = csv.reader(fh, delimiter=',')
for line in csvr:
    a0.append(line)
fh.close()
a0 = np.array(a0)
a0 = a0.astype(float)

# cross channel transects
lon_c0 = np.array([lon[378, 100], lon[378, 112]])
lat_c0 = np.array([lat[378, 100], lat[378, 112]])
lon_c1 = np.array([lon[280, 128], lon[280, 182]])
lat_c1 = np.array([lat[280, 128], lat[280, 182]])
lon_c2 = np.array([lon[210, 127], lon[210, 145]])
lat_c2 = np.array([lat[210, 127], lat[210, 145]])
lon_c3 = np.array([lon[129, 100], lon[180, 100]])
lat_c3 = np.array([lat[129, 100], lat[180, 100]])

# ------------------- make plot ------------------------------------------
fig = plt.figure()

lat_min = 58.00
lat_max = 59.20
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f')

m.fillcontinents(color='lightgrey')
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.04), 
                     labels=[0, 0, 1, 1], fontsize=3, linewidth=.1)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.02), 
                     labels=[1, 1, 0, 0], fontsize=3, linewidth=.1)
# mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.2),
#                      labels=[0, 0, 0, 1], fontsize=10, linewidth=.1)
# pr = m.drawparallels(np.arange(lat_min, lat_max, 0.1),
#                      labels=[1, 0, 0, 0], fontsize=10, linewidth=.1)
setlabelrot(mr,-90)
setlabelrot(pr,0)

xlon, ylat = m(lon, lat)
# m.contourf(xlon, ylat, z, 50, cmap=cmocean.cm.deep)
m.pcolormesh(xlon, ylat, z, cmap=cmocean.cm.deep)
plt.clim(0, 500)
plt.colorbar()

# plot CTD stations
x2, y2 = m(lon_ctd, lat_ctd)
m.plot(x2, y2, '^k', markersize=.5)

# plot transect
xtr, ytr = m(a0[:, 0], a0[:, 1])
m.plot(xtr, ytr, '--ok', linewidth=.3, markersize=.3)
xtr, ytr = m(lon_c0, lat_c0)
m.plot(xtr, ytr, '--ok', linewidth=.3, markersize=.3)
xtr, ytr = m(lon_c1, lat_c1)
m.plot(xtr, ytr, '--ok', linewidth=.3, markersize=.3)
xtr, ytr = m(lon_c2, lat_c2)
m.plot(xtr, ytr, '--ok', linewidth=.3, markersize=.3)
xtr, ytr = m(lon_c3, lat_c3)
m.plot(xtr, ytr, '--ok', linewidth=.3, markersize=.3)

# plot boxes
xbox, ybox = m(box0[:, 0], box0[:, 1])
plt.plot(xbox, ybox, '--k', linewidth=0.5)
plt.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)

xbox, ybox = m(box1[:, 0], box1[:, 1])
plt.plot(xbox, ybox, '--k', linewidth=0.5)
plt.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)

xbox, ybox = m(box2[:, 0], box2[:, 1])
plt.plot(xbox, ybox, '--k', linewidth=0.5)
plt.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)

# save figure
plt.savefig(out_dir + 'figs/'+grd1+'_map.png', format='png', dpi=600)
plt.close()
