'''
Contour transect along Glacier Bay for initial condition.

2017/10/12
Use nearest neighbor instead griddata to get the transect data.
'''

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms
import glob
from matplotlib.colors import LogNorm
import cmocean
from geopy.distance import vincenty
from matplotlib import path
import sys

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

file_in = out_dir + 'bc_ic/GlacierBay_lr_ic_2008_06_15_CTD_floodFill.nc'

# my inputs
plt_contourf = 1
varlist = ['salt', 'temp', 'dye_01', 'dye_03']
varlist = ['salt']
dd = 3
depth1 = 450
depth0 = 50

# dicts for variable clims, colormaps, and other properties
clim = {'temp': [2, 10],
        'salt': [15, 35],
        'zeta': [-3, 3],
        'wetdry_mask_rho': [0, 1],
        'dye_01': [0, 1],
        'dye_02': [0, 1],
        'dye_03': [0, 1],
        'tke': [1e-5, 1e0],
        'gls': [1e-6, 1e-4],
       }

cmap = {'temp': cmocean.cm.thermal,
        'salt': cmocean.cm.haline,
        'zeta': cmocean.cm.balance,
        'wetdry_mask_rho': cmocean.cm.balance,
        'dye_01': cmocean.cm.matter,
        'dye_02': cmocean.cm.matter,
        'dye_03': cmocean.cm.matter,
        'tke': cmocean.cm.matter,
        'gls': cmocean.cm.matter,
       }

var_omega = ['tke', 'gls']
var_log = ['tke', 'gls']

grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)

zlev = grd.vgrid.N

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
               [-136.38969917,   58.78166195],
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

lon_ct = c0[:, 0]
lat_ct = c0[:, 1]

ct_tr = (len(lon_ct)-1)*dd
lon_tr = np.zeros(ct_tr)
lat_tr = np.zeros(ct_tr)

for i in range(len(lon_ct)-1):
    lon_tr[i*dd:(i+1)*dd] = np.linspace(lon_ct[i], lon_ct[i+1], dd+1)[:-1]
    lat_tr[i*dd:(i+1)*dd] = np.linspace(lat_ct[i], lat_ct[i+1], dd+1)[:-1]

# read grid information
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

# interpolate topography
h = grd.vgrid.h
ang = grd.hgrid.angle_rho
msk = grd.hgrid.mask_rho

# instead of using griddata to find interpolated values, use distance to find the nearest rho point and
# represent the value at (lon_tr, lat_tr).
eta_tr = np.zeros(lat_tr.shape)
xi_tr = np.zeros(lon_tr.shape)

for i in range(len(eta_tr)):
    D2 = (lat-lat_tr[i])**2+(lon-lon_tr[i])**2
    eta_tr[i], xi_tr[i] = np.where(D2==D2.min())

eta_tr = eta_tr.astype(int)
xi_tr = xi_tr.astype(int)

h_tr = h[eta_tr.tolist(), xi_tr.tolist()]
# update lon_tr, lat_tr
lon_tr = lon[eta_tr.tolist(), xi_tr.tolist()]
lat_tr = lat[eta_tr.tolist(), xi_tr.tolist()]

# get rho point depth
Cs = -grd.vgrid.Cs_r
z_tr = np.dot(np.matrix(h_tr).T, np.matrix(Cs))
z_tr = z_tr.T

# calculate distance
dis = np.zeros(lat_tr.size)
for i in range(1, lat_tr.size):
    dis[i] = vincenty(
                      (lat_tr[i-1], lon_tr[i-1]),
                      (lat_tr[i], lon_tr[i])
                     ).meters
dis = np.cumsum(dis)
dis = dis/1000  # [km]
dis = np.tile(dis, (zlev, 1))

# -------------------------------------------------------------------------------
# make plots
try:
    plt.style.use('classic')
except:
    pass

# set the axises
f, (ax1, ax2) = plt.subplots(2, sharex=True)
f.subplots_adjust(hspace=0.05)
ax1.set_xlim(dis[0, 0], dis[0, -1])
ax1.set_ylim(0, depth0)
ax1.set_yticks(range(0, depth0, 10))
ax1.invert_yaxis()
ax2.set_ylim(depth0, depth1)
ax2.set_yticks(range(depth0, depth1, 100))
ax2.invert_yaxis()
# plot bathymetry
ax1.fill_between(dis[0, :], -h_tr, depth1, facecolor='lightgrey')
ax2.fill_between(dis[0, :], -h_tr, depth1, facecolor='lightgrey')

for var in varlist:
    print('For ' + var)

    if var in clim.keys():
        clim_var = clim[var]
        cmap_var = cmap[var]
    else:
        clim_var = [0, 1]
        cmap_var = cmocean.cm.matter

    if var in var_log:
        clevs = np.logspace(clim_var[0], clim_var[1], 11)
    else:
        clevs = np.linspace(clim_var[0], clim_var[1], 21)

    fh = nc.Dataset(file_in)
    t = fh.variables['ocean_time'][:]
    tunit = (fh.variables['ocean_time']).units
    data = fh.variables[var][:]
    fh.close()

    if var in var_omega:
        data = 0.5*(data[:, :-1, :, :] + data[:, 1:, :, :])

    for tt in range(len(t)):
        ttag = nc.num2date(t[tt], tunit).strftime("%Y-%m-%d_%H:%M:%S")
        var_tr = data[tt, :, eta_tr.tolist(), xi_tr.tolist()].T

        # make plot
        if var in var_log:
            pcm1 = ax1.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
            pcm2 = ax2.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
        else:
            pcm1 = ax1.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
            pcm2 = ax2.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
            pcm1.set_clim(clim_var[0], clim_var[1])
            pcm2.set_clim(clim_var[0], clim_var[1])
        # add colorbar axis handle
        cbar_ax = f.add_axes([0.92, 0.10, 0.02, 0.8])
        cb = f.colorbar(pcm1, cax=cbar_ax)

        if plt_contourf==1:
            varc1 = ax1.contour(dis, z_tr, var_tr, clevs, linestyle='--', linewidths=.4, colors='w')
            varc2 = ax2.contour(dis, z_tr, var_tr, clevs, linestyle='--', linewidths=.4, colors='w')
            varc3 = ax1.contour(dis, z_tr, var_tr, clevs[::5], linestyle='--', linewidths=1.0, colors='k')
            varc4 = ax2.contour(dis, z_tr, var_tr, clevs[::5], linestyle='--', linewidths=1.0, colors='k')
            cl3 = plt.clabel(varc3, fontsize=5)
            cl4 = plt.clabel(varc4, fontsize=5)

        ax1.set_title(var + '_' + grd.name + '_' + ttag)
        f.savefig(out_dir + 'figs/init/' + var + '_' + grd.name + '_' + ttag + '.png', dpi=300)

        pcm1.remove()
        pcm2.remove()
        f.delaxes(cbar_ax)
        # cb.remove()
        if plt_contourf==1:
            for cc in varc1.collections:
                cc.remove()
            for cc in varc2.collections:
                cc.remove()
            for cc in varc3.collections:
                cc.remove()
            for cc in varc4.collections:
                cc.remove()

# plt.close()
