""" Plot Glacier Bay map """

import csv
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects
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

boxes = {}
boxes['box0'] = np.array([
    [-136.65, 58.65],
    [-137.30, 58.65],
    [-137.30, 59.15],
    [-136.60, 59.15],
    [-136.30, 58.95],
    [-136.20, 58.75]])

boxes['box1'] = np.array([
    [-136.20, 58.75],
    [-136.30, 58.95],
    [-136.60, 59.15],
    [-135.70, 59.15],
    [-135.60, 58.75]])

boxes['box2'] = np.array([
    [-136.65, 58.55],
    [-136.65, 58.65],
    [-136.20, 58.75],
    [-135.60, 58.75],
    [-135.60, 58.55]])

boxes['box3'] = np.array([
    [-136.65, 58.55],
    [-135.60, 58.55],
    [-135.60, 58.50],
    [-136.00, 58.35]])

boxes['box4'] = np.array([
    [-136.65, 58.55],
    [-136.75, 58.25],
    [-136.50, 57.95],
    [-135.20, 57.95],
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
fh = open('../../data/a0.txt')
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
lon_max = -135.3
lon_0 = 0.5 * (lon_min + lon_max)

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f', ax=ax)

m.fillcontinents(color='lightgrey')
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.4), 
                     labels=[0, 0, 0, 1], fontsize=7, linewidth=.1)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.2),
                     labels=[1, 0, 0, 0], fontsize=7, linewidth=.1)

# setlabelrot(mr, -90)
# setlabelrot(pr, 0)

# map scale
m.drawmapscale(-137.1, 58.6, lon_0, lat_0, 50,
               barstyle='fancy', labelstyle='simple',
               fontsize=5)

xlon, ylat = m(lon, lat)
# m.contourf(xlon, ylat, z, 50, cmap=cmocean.cm.deep)
pcm = ax.pcolormesh(xlon, ylat, z,
                    vmin=0, vmax=300,
                    alpha=0.8, cmap=cmocean.cm.deep)

# plot colorbar
cbar_ax = fig.add_axes([0.69, 0.50, 0.015, 0.35])
cbar_ax.tick_params(labelsize=7)
cb = fig.colorbar(pcm, cax=cbar_ax, ticks=np.linspace(0, 450, 10))
cbar_ax.set_ylabel(r'Depth [m]', fontsize=7)

# plot texts
x1, y1 = m(-135.7, 58.30)
txt = ax.text(x1, y1, 'Icy Strait',
              fontsize=6, rotation=-30, color='k',
              ha='center', va='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-136.4, 58.22)
txt = ax.text(x1, y1, 'Cross Sound',
              fontsize=6, rotation=35, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-135.95, 58.475)
txt = ax.text(x1, y1, 'Lower Bay',
              fontsize=6, rotation=-95, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-136.05, 58.65)
txt = ax.text(x1, y1, 'Central Bay',
              fontsize=6, rotation=-60, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-136.53, 58.84)
txt = ax.text(x1, y1, 'West',
              fontsize=5, rotation=-45, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-136.4, 58.79)
txt = ax.text(x1, y1, 'Arm',
              fontsize=5, rotation=-35, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-136.12, 58.94)
txt = ax.text(x1, y1, 'East',
              fontsize=5, rotation=-65, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])
x1, y1 = m(-136.09, 58.85)
txt = ax.text(x1, y1, 'Arm',
              fontsize=5, rotation=-90, color='k',
              horizontalalignment='center', verticalalignment='center')
txt.set_path_effects([path_effects.Stroke(linewidth=0.7, foreground='w'),
                      path_effects.Normal()])

# plot CTD stations
x2, y2 = m(lon_ctd, lat_ctd)
ax.plot(x2, y2, '^r', markersize=3, label='CTD Station')
ax.plot(x2[core_stn], y2[core_stn], '^r',
        markeredgecolor='k', markeredgewidth=1,
        markerfacecolor='r', markersize=4,
        label='Core Station')
for i in range(len(core_stn)):
    ax.text(x2[core_stn[i]]-5000, y2[core_stn[i]]-5000, str(core_stn[i]),
            ha='center', va='center', fontsize=5)

# # plot ADCP stations
# lon_adcp = [-136.03453]
# lat_adcp = [58.467072]
# x3, y3 = m(lon_adcp, lat_adcp)
# ax.plot(x3, y3, 'Hb', markersize=3, label='ADCP Station')

# plot transect
xtr, ytr = m(a0[:, 0], a0[:, 1])
ax.plot(xtr, ytr, '--r', linewidth=.5, label='Main Channel')

# # plot boxes
# for boxi in boxes.keys():
#     xbox, ybox = m(boxes[boxi][:, 0], boxes[boxi][:, 1])
#     ax.plot(xbox, ybox, '--k', linewidth=0.5)
#     ax.plot([xbox[-1], xbox[0]], [ybox[-1], ybox[0]], '--k', linewidth=0.5)
# 
# # plot box numbers
# x1, y1 = m(-137.25, 59.125)
# ax.add_patch(patches.Circle((x1, y1), 4000, fill=False))
# ax.text(x1, y1-500, r'0', fontsize=8, horizontalalignment='center',
#         verticalalignment='center')
# x1, y1 = m(-136.45, 59.125)
# ax.add_patch(patches.Circle((x1, y1), 4000, fill=False))
# ax.text(x1, y1-500, r'1', fontsize=8, horizontalalignment='center',
#         verticalalignment='center')
# x1, y1 = m(-136.6, 58.625)
# ax.add_patch(patches.Circle((x1, y1), 4000, fill=False))
# ax.text(x1, y1-500, r'2', fontsize=8, horizontalalignment='center',
#         verticalalignment='center')
# x1, y1 = m(-136.45, 58.525)
# ax.add_patch(patches.Circle((x1, y1), 4000, fill=False))
# ax.text(x1, y1-500, r'3', fontsize=8, horizontalalignment='center',
#         verticalalignment='center')
# x1, y1 = m(-136.6, 58.475)
# ax.add_patch(patches.Circle((x1, y1), 4000, fill=False))
# ax.text(x1, y1-500, r'4', fontsize=8, horizontalalignment='center',
#         verticalalignment='center')

# plot legend
ax.legend(loc=(0.05, 0.55), fontsize=6)

# plot small map
axz = plt.axes([0.25, 0.05, 0.2, 0.4])
m2 = Basemap(width=6000000, height=6000000, projection='lcc',
             resolution='l', lat_1=45., lat_2=55, lat_0=50, lon_0=-127.,
             ax=axz)
m2.drawmapboundary(fill_color='aqua')
m2.fillcontinents(color='lightgrey')
# m2.etopo()
x0, y0 = m2(lon_0, lat_0)
m2.plot(x0, y0, 'k', marker=(5, 1))

# save figure
plt.savefig(out_dir + 'figs/osm2018/osm_' + grd1 + '_map.png', format='png', dpi=600)
plt.close()
