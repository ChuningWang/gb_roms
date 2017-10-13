# Plot Glacier Bay map
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
c0 = np.array([[-137.05000708,   59.05576767],
               [-137.02431432,   59.03650282],
               [-136.99075294,   59.02081148],
               [-136.97549777,   58.99798771],
               [-136.95566605,   58.99085529],
               [-136.9251557 ,   58.98086989],
               [-136.91295157,   58.96375207],
               [-136.89922191,   58.94806073],
               [-136.89617088,   58.92666345],
               [-136.88396674,   58.90954562],
               [-136.85955846,   58.90526617],
               [-136.83667571,   58.90098671],
               [-136.81379295,   58.89528077],
               [-136.79091019,   58.89528077],
               [-136.7741295 ,   58.89813374],
               [-136.74972123,   58.89813374],
               [-136.72226192,   58.89228077],
               [-136.69937916,   58.88386889],
               [-136.67802192,   58.87816295],
               [-136.66429226,   58.87530998],
               [-136.63835847,   58.86817755],
               [-136.60937365,   58.86389809],
               [-136.58954192,   58.86389809],
               [-136.5605571 ,   58.8539127 ],
               [-136.53614882,   58.8439273 ],
               [-136.50716399,   58.82823596],
               [-136.49190882,   58.82110354],
               [-136.47970468,   58.81610354],
               [-136.45987296,   58.81282408],
               [-136.44156675,   58.80997111],
               [-136.42326055,   58.7974122 ],
               [-136.40953089,   58.7954268 ],
               [-136.38969917,   58.77366195],
               [-136.37291848,   58.77250898],
               [-136.34851021,   58.77260303],
               [-136.330204  ,   58.76560303],
               [-136.29816814,   58.76289709],
               [-136.285964  ,   58.76019115],
               [-136.2768109 ,   58.74835278],
               [-136.26003021,   58.73408793],
               [-136.24630055,   58.71839659],
               [-136.23714745,   58.70127877],
               [-136.20968814,   58.68986689],
               [-136.17002469,   58.67702852],
               [-136.15171849,   58.67417555],
               [-136.11358056,   58.66276367],
               [-136.0906978 ,   58.6399399 ],
               [-136.07696815,   58.61711614],
               [-136.06781504,   58.59001292],
               [-136.05866194,   58.55149781],
               [-136.05561091,   58.51726216],
               [-136.0373047 ,   58.47018815],
               [-136.01899849,   58.44736438],
               [-135.99916677,   58.42168765],
               [-135.97323298,   58.38174606],
               [-135.94272263,   58.35036338],
               [-135.91373781,   58.33467204],
               [-135.86949781,   58.3289661 ],
               [-135.82068126,   58.31898071],
               [-135.77949229,   58.31898071],
               [-135.72609919,   58.31042179],
               [-135.67880816,   58.290451  ],
               [-135.62999161,   58.27618615],
               [-135.56286885,   58.26620075],
               [-135.53083299,   58.26192129],
               [-135.48049093,   58.25478887],
               [-135.44082748,   58.25478887],
               [-135.39658748,   58.24480347],
               [-135.34777093,   58.21912673],
               [-135.30810748,   58.21342079],
               [-135.26844404,   58.21484728],
               [-135.21505094,   58.19059703],
               [-135.17386197,   58.16206732],
               [-135.13419852,   58.13781707],
               [-135.08995853,   58.12069925],
               [-135.05792266,   58.10358142]])


m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f')

# m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='lightgrey')
# mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.05), labels=[0,0,1,1], fontsize=4, linewidth=.1)
# pr = m.drawparallels(np.arange(lat_min, lat_max, 0.025), labels=[1,1,0,0], fontsize=4, linewidth=.1)
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.2), labels=[0,0,0,1], fontsize=10, linewidth=.1)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.1), labels=[1,0,0,0], fontsize=10, linewidth=.1)
setlabelrot(mr,-30)
setlabelrot(pr,-30)

x, y = m(lon, lat)
m.pcolor(x, y, z, cmap=cmocean.cm.deep)
plt.clim(0, 400)
plt.colorbar()

# plot CTD stations
x2, y2 = m(lon_ctd, lat_ctd)
m.plot(x2, y2, '^k', ms=2)

# plot transect
xc0, yc0 = m(c0[:, 0], c0[:, 1])
m.plot(xc0, yc0, '--k', linewidth=.3)

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

# xbox1, ybox1 = m(box1[:, 0], box1[:, 1])
# m.plot(xbox1, ybox1, '--k')
# m.plot([xbox1[0], xbox1[-1]], [ybox1[0], ybox1[-1]], '--k')
# 
# xbox2, ybox2 = m(box2[:, 0], box2[:, 1])
# m.plot(xbox2, ybox2, '--k')
# m.plot([xbox2[0], xbox2[-1]], [ybox2[0], ybox2[-1]], '--k')

plt.savefig(out_dir + 'figs/'+grd1+'_map.tiff', format='tiff', dpi=600)
plt.close()
