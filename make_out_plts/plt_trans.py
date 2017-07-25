import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms
import glob
from gb_toolbox import gb_ctd
from matplotlib.mlab import griddata
from geopy.distance import vincenty

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

grd1 = 'GB_USGS'
model = 'tmpdir_GB-TIDE/outputs/2000/'
clim = [25, 35]

# load data
outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/trans/2000/'

flist = glob.glob(outputs_dir+'*.nc')
# flist = flist[-48:]

depth = 200
dd = 3
zlev = 40
tindex = 0
var = 'salt'
uvar = 'u'
vvar = 'v'
clim = [28, 32]

# load geological coordinates
ctd = gb_ctd.rd_ctd(in_dir + 'ctd.nc')
lat_ctd = ctd['lat_stn']
lon_ctd = ctd['lon_stn']

station = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 21]
lat_ctd = lat_ctd[station]
lon_ctd = lon_ctd[station]

c0 = np.array([[-137.04719708,   59.05076767],
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

lon_ctd = c0[:, 0]
lat_ctd = c0[:, 1]


ct_tr = (len(lon_ctd)-1)*dd
lat_t = np.zeros(ct_tr)
lon_t = np.zeros(ct_tr)

for i in range(len(lon_ctd)-1):
    lat_t[i*dd:(i+1)*dd] = np.linspace(lat_ctd[i], lat_ctd[i+1], dd+1)[:-1]
    lon_t[i*dd:(i+1)*dd] = np.linspace(lon_ctd[i], lon_ctd[i+1], dd+1)[:-1]

grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

# interpolate topography
h = grd.vgrid.h
msk = grd.hgrid.mask_rho
h_tr = griddata(lon.flatten(), lat.flatten(), h.flatten(), lon_t, lat_t, interp='linear').diagonal()

# get rho point depth
s_rho = grd.vgrid.s_rho
z_tr = np.dot(np.matrix(h_tr).T, np.matrix(s_rho))
z_tr = z_tr.T

# calculate distance
dis = np.zeros(lat_t.size)
for i in range(1, lat_t.size):
    dis[i] = vincenty(
                      (lat_t[i-1], lon_t[i-1]),
                      (lat_t[i], lon_t[i])
                     ).meters
dis = np.cumsum(dis)
dis = dis/1000  # [km]
dis = np.tile(dis, (zlev, 1))

# interpolate salinity
s_tr = np.zeros((zlev, ct_tr))

plt.switch_backend('Agg')

# plot transaction on map
fig = plt.figure()
plt.pcolor(lon, lat, h, cmap='Greens')
plt.clim(0, 400)
plt.contour(lon, lat, msk, np.array([0.5, 0.5]), colors='k')
plt.plot(c0[:, 0], c0[:, 1], '--.k', ms=3)
plt.savefig(fig_dir + 'map_transac.png')
plt.close()

for fn in flist:
    tag = fn.split('/')[-1].split('.')[0]
    print 'processing ' + tag + ' ...'
    fh = nc.Dataset(fn)
    s = fh.variables['salt'][:].squeeze()
    fh.close()

    for i in range(40):
        s_tr[i, :] = griddata(lon.flatten(), lat.flatten(), s[i, :, :].squeeze().flatten(), lon_t, lat_t, interp='linear').diagonal()

    plt.pcolormesh(dis, z_tr, s_tr)
    plt.clim(clim[0], clim[1])
    plt.colorbar()
    plt.savefig(fig_dir + var + '_' + tag + '.png')
    plt.close()
