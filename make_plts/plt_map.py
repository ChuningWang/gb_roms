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

box0 = np.array([[-136.90, 58.55],
                 [-137.30, 58.65],
                 [-137.30, 59.15],
                 [-135.70, 59.15],
                 [-135.60, 58.75],
                 [-135.60, 58.55]])

box1 = np.array([[-136.90, 58.55],
                 [-136.90, 58.45],
                 [-136.00, 58.35],
                 [-135.60, 58.50],
                 [-135.60, 58.55]])

box2 = np.array([[-136.90, 58.45],
                 [-136.65, 58.25],
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

core_stn = [1, 4, 7, 12, 13, 16, 20, 24]

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
lon_c0 = np.array([lon[327, 112], lon[337, 138]])
lat_c0 = np.array([lat[327, 112], lat[337, 138]])
lon_c1 = np.array([lon[211, 127], lon[211, 144]])
lat_c1 = np.array([lat[211, 127], lat[211, 144]])
lon_c2 = np.array([lon[184, 126], lon[175, 145]])
lat_c2 = np.array([lat[184, 126], lat[175, 145]])

# ------------------- make plot ------------------------------------------
fig, ax = plt.subplots()

lat_min = 58.00
lat_max = 59.20
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f', ax=ax)

m.fillcontinents(color='lightgrey')
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.2),
                     labels=[0, 0, 0, 1], fontsize=10, linewidth=.1)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.1),
                     labels=[1, 0, 0, 0], fontsize=10, linewidth=.1)
# mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.04), 
#                      labels=[0, 0, 1, 1], fontsize=3, linewidth=.1)
# pr = m.drawparallels(np.arange(lat_min, lat_max, 0.02), 
#                      labels=[1, 1, 0, 0], fontsize=3, linewidth=.1)

setlabelrot(mr, -30)
setlabelrot(pr, 0)

m.drawmapscale(-137.1, 58.25, lon_0, lat_0, 50,
               barstyle='fancy', labelstyle='simple')

xlon, ylat = m(lon, lat)
# m.contourf(xlon, ylat, z, 50, cmap=cmocean.cm.deep)
pcm = ax.pcolormesh(xlon, ylat, z,
                    vmin=0, vmax=450,
                    alpha=0.8, cmap=cmocean.cm.deep)

# plot colorbar
cbar_ax = fig.add_axes([0.70, 0.50, 0.02, 0.35])
cb = fig.colorbar(pcm, cax=cbar_ax, ticks=np.linspace(0, 450, 10))
cbar_ax.set_ylabel(r'Depth [m]')

# plot CTD stations
x2, y2 = m(lon_ctd, lat_ctd)
ax.plot(x2, y2, '^r', markersize=3, label='CTD Station')
ax.plot(x2[core_stn], y2[core_stn], '^r',
        markeredgecolor='k', markeredgewidth=1,
        markerfacecolor='r', markersize=4,
        label='Core Station')

# # plot transect
xtr, ytr = m(a0[:, 0], a0[:, 1])
ax.plot(xtr, ytr, '--or', linewidth=.2, markersize=.2)
xtr, ytr = m(lon_c0, lat_c0)
ax.plot(xtr, ytr, '--or', linewidth=.2, markersize=.2)
xtr, ytr = m(lon_c1, lat_c1)
ax.plot(xtr, ytr, '--or', linewidth=.2, markersize=.2)
xtr, ytr = m(lon_c2, lat_c2)
ax.plot(xtr, ytr, '--or', linewidth=.2, markersize=.2)

# plot boxes
xbox, ybox = m(box0[:, 0], box0[:, 1])
ax.plot(xbox, ybox, '--k', linewidth=0.5)
ax.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)

xbox, ybox = m(box1[:, 0], box1[:, 1])
ax.plot(xbox, ybox, '--k', linewidth=0.5)
ax.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)

xbox, ybox = m(box2[:, 0], box2[:, 1])
ax.plot(xbox, ybox, '--k', linewidth=0.5)
ax.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)

# plot texts
x1, y1 = m(-135.7, 58.25)
ax.text(x1, y1, 'Icy Strait', fontsize=10, rotation=-30)
x1, y1 = m(-136.5, 58.25)
ax.text(x1, y1, 'Cross Sound', fontsize=10, rotation=30)
x1, y1 = m(-136.8, 58.475)
ax.text(x1, y1, 'Lower Bay', fontsize=10)
x1, y1 = m(-137.1, 58.7)
ax.text(x1, y1, 'Upper Bay', fontsize=10)

# plot legend
ax.legend(loc=3)

# save figure
plt.savefig(out_dir + 'figs/' + grd1 + '_map.png', format='png', dpi=600)
plt.close()
