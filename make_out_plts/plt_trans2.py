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
plt_contour = 1
grd1 = 'GB_lr'
ftype = 'avg'
varlist = ['salt', 'temp', 'dye_01', 'dye_03']
varlist = ['salt', 'temp', 'v']
# varlist = ['tke', 'gls']
var_contour = 'salt'
dd = 20
depth1 = 450
depth0 = 50

if 'v' in varlist:
    plt_uv = 1
    plt_contour = 1

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
        'v': [-0.5, 0.5],
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
        'v': cmocean.cm.balance,
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
fig_dir = out_dir + 'figs/trans2/' + tag +'/'

flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
# flist = flist[-14:]

zlev = grd.vgrid.N
uvar = 'u'
vvar = 'v'
wvar = 'w'

trs_points = {}

trs_points['c0'] = np.array([[378, 100],
                             [378, 112]])
trs_points['c1'] = np.array([[280, 128],
                             [280, 182]])
trs_points['c2'] = np.array([[210, 127],
                             [210, 145]])
trs_points['c3'] = np.array([[129, 100],
                             [180, 100]])

eta = {}
xi = {}
rdir = {}

for tp in ['c0', 'c1', 'c2', 'c3']:
# for tp in ['c3']:
    cc = trs_points[tp]
    if cc[0, 0] == cc[1, 0]:
        eta = cc[0, 0]*np.ones(abs(cc[1, 1] - cc[0, 1]))
        xi = np.arange(cc[:, 1].min(), cc[:, 1].max())
        rdir = 'ew'
    else:
        eta = np.arange(cc[:, 0].min(), cc[:, 0].max())
        xi = cc[0, 1]*np.ones(abs(cc[1, 0] - cc[0, 0]))
        rdir = 'sn'

    eta = eta.astype(int)
    xi = xi.astype(int)

    # get geo info
    h_tr = grd.vgrid.h[eta, xi]
    lon_tr = grd.hgrid.lon_rho[eta, xi]
    lat_tr = grd.hgrid.lat_rho[eta, xi]
    z_tr = -grd.vgrid.z_r[:][:, eta, xi]

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
    plt.switch_backend('Agg')
    try:
        plt.style.use('classic')
    except:
        pass

    # set the axis
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
            clevs_c = np.logspace(clim[var_contour][0], clim[var_contour][1], 11)
        else:
            clevs = np.linspace(clim_var[0], clim_var[1], 21)
            clevs_c = np.linspace(clim[var_contour][0], clim[var_contour][1], 21)

        for fn in flist:
            # read data
            tag = fn.split('/')[-1].split('.')[0]
            print('processing ' + tag + ' ...')
            fh = nc.Dataset(fn)
            t = fh.variables['ocean_time'][:]
            tunit = (fh.variables['ocean_time']).units
            data = fh.variables[var][:]
            data_c = fh.variables[var_contour][:]
            if plt_uv == 1:
                u = fh.variables[uvar][:]
                v = fh.variables[vvar][:]
                w = fh.variables[wvar][:]
                u = 0.5*(u[:, :, 1:, :]+u[:, :, :-1, :])
                v = 0.5*(v[:, :, :, 1:]+v[:, :, :, :-1])
                w = 0.5*(w[:, 1:, :, :]+w[:, :-1, :, :])
            fh.close()

            if (var == 'v') & (rdir == 'ew'):
                data = v.copy()
            elif (var == 'v') & (rdir == 'sn'):
                data = -u.copy()
                u = v.copy()
                v = data.copy()

            if var in var_omega:
                data = 0.5*(data[:, :-1, :, :] + data[:, 1:, :, :])
            if var_contour in var_omega:
                data_c = 0.5*(data_c[:, :-1, :, :] + data_c[:, 1:, :, :])

            for tt in range(len(t)):
                ttag = nc.num2date(t[tt], tunit).strftime("%Y-%m-%d_%H:%M:%S")
                if plt_uv==1:
                    u_tr = u[tt, :, eta, xi].T
                    v_tr = u[tt, :, eta, xi].T
                    w_tr = w[tt, :, eta, xi].T

                var_tr = data[tt, :, eta, xi].T
                var_c_tr = data_c[tt, :, eta, xi].T

                # make plot
                if var in var_log:
                    # pcm1 = ax1.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
                    # pcm2 = ax2.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
                    pcm1 = ax1.contourf(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
                    pcm2 = ax2.contourf(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
                else:
                    # pcm1 = ax1.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
                    # pcm2 = ax2.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
                    pcm1 = ax1.contourf(dis, z_tr, var_tr, cmap=cmap_var)
                    pcm2 = ax2.contourf(dis, z_tr, var_tr, cmap=cmap_var)
                    pcm1.set_clim(clim_var[0], clim_var[1])
                    pcm2.set_clim(clim_var[0], clim_var[1])

                # add colorbar axis handle
                cbar_ax = f.add_axes([0.90, 0.10, 0.02, 0.8])
                cb = f.colorbar(pcm1, cax=cbar_ax)
                cb.ax.tick_params(labelsize=8)

                if plt_contour==1:
                    varc11 = ax1.contour(dis, z_tr, var_c_tr, clevs_c[::5], linestyle='--', linewidths=1.0, colors='k')
                    varc21 = ax2.contour(dis, z_tr, var_c_tr, clevs_c[::5], linestyle='--', linewidths=1.0, colors='k')
                    varc12 = ax1.contour(dis, z_tr, var_c_tr, clevs_c, linestyle='--', linewidths=.4, colors='w')
                    varc22 = ax2.contour(dis, z_tr, var_c_tr, clevs_c, linestyle='--', linewidths=.4, colors='w')
                    clb11 = plt.clabel(varc11, fontsize=10)
                    clb21 = plt.clabel(varc21, fontsize=10)

                if plt_uv==1:
                    qv1 = ax1.quiver(dis, z_tr, u_tr, w_tr, scale=10)
                    qv2 = ax2.quiver(dis, z_tr, u_tr, w_tr, scale=10)

                f.suptitle(var + '_' + grd.name + '_' + ftype + '_' + tp + '_' + ttag)
                f.savefig(fig_dir + var + '_' + grd.name + '_' + ftype + '_' + tp + '_' + ttag + '.png')

                # pcm1.remove()
                # pcm2.remove()
                for cc in pcm1.collections:
                    cc.remove()
                for cc in pcm2.collections:
                    cc.remove()
                f.delaxes(cbar_ax)
                # cb.remove()
                if plt_contour==1:
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
