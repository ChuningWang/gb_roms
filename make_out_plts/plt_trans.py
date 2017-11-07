'''
Contour transect along Glacier Bay.
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

# my inputs
my_year = 2008
plt_uv = 1
plt_contourf = 1
grd1 = 'GB_lr'
ftype = 'avg'
varlist = ['salt', 'temp', 'dye_01', 'dye_03']
varlist = ['salt', 'temp']
# varlist = ['tke', 'gls']
dd = 3
depth1 = 450
depth0 = 50

# dicts for variable clims, colormaps, and other properties
clim = {'temp': [2, 10],
        'salt': [15, 30],
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

grd = pyroms.grid.get_ROMS_grid(grd1)

if len(sys.argv)>1:
    tag = sys.argv[-1]
else:
    tag = 'GB-ref'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/trans/' + tag +'/'

flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
# flist = flist[-14:]

zlev = grd.vgrid.N
uvar = 'u'
vvar = 'v'
wvar = 'w'

# c0 = np.array([[-137.05000708,   59.05576767],
#                [-137.02431432,   59.03650282],
#                [-136.99075294,   59.02081148],
#                [-136.97549777,   58.99798771],
#                [-136.95566605,   58.99085529],
#                [-136.9251557 ,   58.98086989],
#                [-136.91295157,   58.96375207],
#                [-136.89922191,   58.94806073],
#                [-136.89617088,   58.92666345],
#                [-136.88396674,   58.90954562],
#                [-136.85955846,   58.90526617],
#                [-136.83667571,   58.90098671],
#                [-136.81379295,   58.89528077],
#                [-136.79091019,   58.89528077],
#                [-136.7741295 ,   58.89813374],
#                [-136.74972123,   58.89813374],
#                [-136.72226192,   58.89228077],
#                [-136.69937916,   58.88386889],
#                [-136.67802192,   58.87816295],
#                [-136.66429226,   58.87530998],
#                [-136.63835847,   58.86817755],
#                [-136.60937365,   58.86389809],
#                [-136.58954192,   58.86389809],
#                [-136.5605571 ,   58.8539127 ],
#                [-136.53614882,   58.8439273 ],
#                [-136.50716399,   58.82823596],
#                [-136.49190882,   58.82110354],
#                [-136.47970468,   58.81610354],
#                [-136.45987296,   58.81282408],
#                [-136.44156675,   58.80997111],
#                [-136.42326055,   58.7974122 ],
#                [-136.40953089,   58.7954268 ],
#                [-136.38969917,   58.77366195],
#                [-136.37291848,   58.77250898],
#                [-136.34851021,   58.77260303],
#                [-136.330204  ,   58.76560303],
#                [-136.29816814,   58.76289709],
#                [-136.285964  ,   58.76019115],
#                [-136.2768109 ,   58.74835278],
#                [-136.26003021,   58.73408793],
#                [-136.24630055,   58.71839659],
#                [-136.23714745,   58.70127877],
#                [-136.20968814,   58.68986689],
#                [-136.17002469,   58.67702852],
#                [-136.15171849,   58.67417555],
#                [-136.11358056,   58.66276367],
#                [-136.0906978 ,   58.6399399 ],
#                [-136.07696815,   58.61711614],
#                [-136.06781504,   58.59001292],
#                [-136.05866194,   58.55149781],
#                [-136.05561091,   58.51726216],
#                [-136.0373047 ,   58.47018815],
#                [-136.01899849,   58.44736438],
#                [-135.99916677,   58.42168765],
#                [-135.97323298,   58.38174606],
#                [-135.94272263,   58.35036338],
#                [-135.91373781,   58.33467204],
#                [-135.86949781,   58.3289661 ],
#                [-135.82068126,   58.31898071],
#                [-135.77949229,   58.31898071],
#                [-135.72609919,   58.31042179],
#                [-135.67880816,   58.290451  ],
#                [-135.62999161,   58.27618615],
#                [-135.56286885,   58.26620075],
#                [-135.53083299,   58.26192129],
#                [-135.48049093,   58.25478887],
#                [-135.44082748,   58.25478887],
#                [-135.39658748,   58.24480347],
#                [-135.34777093,   58.21912673],
#                [-135.30810748,   58.21342079],
#                [-135.26844404,   58.21484728],
#                [-135.21505094,   58.19059703],
#                [-135.17386197,   58.16206732],
#                [-135.13419852,   58.13781707],
#                [-135.08995853,   58.12069925],
#                [-135.05792266,   58.10358142]])

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
               [-136.01859570,   58.36218704],
               [-136.06103693,   58.34781386],
               [-136.10532343,   58.32769141],
               [-136.15330047,   58.31906751],
               [-136.19574170,   58.31044360],
               [-136.23818293,   58.30325702],
               [-136.26586199,   58.29894506],
               [-136.30276740,   58.29463311],
               [-136.33782755,   58.29463311],
               [-136.37473296,   58.29032116],
               [-136.39318567,   58.28744652],
               [-136.41901946,   58.29319579],
               [-136.44116271,   58.28025993],
               [-136.48175867,   58.25438822],
               [-136.52604516,   58.21989260],
               [-136.54818841,   58.19977015],
               [-136.56479585,   58.18395966],
               [-136.58693910,   58.16958648]])
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
# lon_tr = lon[eta_tr.tolist(), xi_tr.tolist()]
# lat_tr = lat[eta_tr.tolist(), xi_tr.tolist()]
# get depth as well
z_tr = -grd.vgrid.z_r[:][:, eta_tr.tolist(), xi_tr.tolist()]

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

# initiate vairables
# var_tr = np.zeros((zlev, ct_tr))

# -------------------------------------------------------------------------------
# if plot velocity vector, also calculate and define these variables
if plt_uv==1:
    latu = grd.hgrid.lat_u
    lonu = grd.hgrid.lon_u
    latv = grd.hgrid.lat_v
    lonv = grd.hgrid.lon_v
    msku = grd.hgrid.mask_u
    mskv = grd.hgrid.mask_v

    lonu = lonu[msku==1]
    latu = latu[msku==1]
    lonv = lonv[mskv==1]
    latv = latv[mskv==1]

    eta_tr2 = eta_tr[::dd]
    xi_tr2 = xi_tr[::dd]
    lon_tr2 = lon_tr[::dd]
    lat_tr2 = lat_tr[::dd]
    h_tr2 = h_tr[::dd]
    z_tr2 = z_tr[:, ::dd]
    dis2 = dis[:, ::dd]

    # calculate angle
    ang_tr = ang[eta_tr.tolist(), xi_tr.tolist()]
    ang_tr2 = ang_tr[::dd]
    ang_add = np.zeros(len(ang_tr2))
    dx = 59
    dy = 111
    dvec = np.diff(lon_tr2*dx + 1j*lat_tr2*dy)
    ang_add[:-1] = np.angle(dvec)
    ang_add[-1] = ang_add[-2]
    ang_tr2 = ang_tr2-ang_add

    U_tr2 = np.zeros((zlev, len(lon_tr2)))
    w_tr2_sw = np.zeros((zlev+1, len(lon_tr2)))

# -------------------------------------------------------------------------------
# make plots
plt.switch_backend('Agg')
try:
    plt.style.use('classic')
except:
    pass

# set the axises
f, (ax0, ax1, ax2) = plt.subplots(3, sharex=True, gridspec_kw={'height_ratios':[1, 2, 2]})
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

    for fn in flist:
        # read data
        tag = fn.split('/')[-1].split('.')[0]
        print('processing ' + tag + ' ...')
        fh = nc.Dataset(fn)
        t = fh.variables['ocean_time'][:]
        tunit = (fh.variables['ocean_time']).units
        data = fh.variables[var][:]
        zeta = fh.variables['zeta'][:]
        if plt_uv==1:
            u = fh.variables[uvar][:]
            v = fh.variables[vvar][:]
            w = fh.variables[wvar][:]
            u = 0.5*(u[:, :, 1:, :]+u[:, :, :-1, :])
            v = 0.5*(v[:, :, :, 1:]+v[:, :, :, :-1])
            w = 0.5*(w[:, 1:, :, :]+w[:, :-1, :, :])
        fh.close()

        if var in var_omega:
            data = 0.5*(data[:, :-1, :, :] + data[:, 1:, :, :])

        for tt in range(len(t)):
            ttag = nc.num2date(t[tt], tunit).strftime("%Y-%m-%d_%H:%M:%S")
            var_tr = data[tt, :, eta_tr.tolist(), xi_tr.tolist()].T
            zeta_tr = zeta[tt, eta_tr.tolist(), xi_tr.tolist()].T

            if plt_uv==1:
                u_tr2 = u[tt, :, eta_tr2.tolist(), xi_tr2.tolist()].T
                v_tr2 = u[tt, :, eta_tr2.tolist(), xi_tr2.tolist()].T
                U_tr2 = u_tr2*np.cos(np.tile(ang_tr2, (40, 1)))-v_tr2*np.sin(np.tile(ang_tr2, (40, 1)))
                w_tr2 = w[tt, :, eta_tr2.tolist(), xi_tr2.tolist()].T

            # make plot
            pltz = ax0.plot(dis[0, :], zeta_tr, 'k')
            ax0.set_ylim(-0.25, 0.25)
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
                varc11 = ax1.contour(dis, z_tr, var_tr, clevs[::5], linestyle='--', linewidths=1.0, colors='k')
                varc21 = ax2.contour(dis, z_tr, var_tr, clevs[::5], linestyle='--', linewidths=1.0, colors='k')
                varc12 = ax1.contour(dis, z_tr, var_tr, clevs, linestyle='--', linewidths=.4, colors='w')
                varc22 = ax2.contour(dis, z_tr, var_tr, clevs, linestyle='--', linewidths=.4, colors='w')
                clb11 = plt.clabel(varc11, fontsize=10)
                clb21 = plt.clabel(varc21, fontsize=10)
            if plt_uv==1:
                qv1 = ax1.quiver(dis2, z_tr2, U_tr2, w_tr2, scale=10)
                qv2 = ax2.quiver(dis2, z_tr2, U_tr2, w_tr2, scale=10)

            f.suptitle(var + '_' + grd.name + '_' + ftype + '_' + ttag)
            f.savefig(fig_dir + var + '_' + grd.name + '_' + ftype + '_' + ttag + '.png')

            pltz.pop(0).remove()
            pcm1.remove()
            pcm2.remove()
            f.delaxes(cbar_ax)
            # cb.remove()
            if plt_contourf==1:
                for cc in varc11.collections:
                    cc.remove()
                for cc in varc21.collections:
                    cc.remove()
                for cc in varc12.collections:
                    cc.remove()
                for cc in varc22.collections:
                    cc.remove()

                for cl in clb11:
                    cl.remove()
                for cl in clb21:
                    cl.remove()

            if plt_uv==1:
                qv1.remove()
                qv2.remove()

plt.close()
