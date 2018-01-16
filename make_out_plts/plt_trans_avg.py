"""
Contour transect along Glacier Bay.
2017/10/12
Use nearest neighbor instead griddata to get the transect data.
"""

# --------------------- load modules --------------------------------------
import sys
import csv
import glob

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import path
from matplotlib.gridspec import GridSpec

import cmocean
import pyroms
from geopy.distance import vincenty

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# --------------------- functionals ---------------------------------------
def get_z(h, hc, N, s_rho, Cs_r, zeta, Vtrans):

    z_r = np.empty((N, len(h)), 'd')
    if Vtrans == 1:
        for k in range(N):
            z0 = hc * s_rho[k] + (h - hc) * Cs_r[k]
            z_r[k, :] = z0 + zeta * (1.0 + z0 / h)
    elif Vtrans == 2 or Vtrans == 4 or Vtrans == 5:
        for  k in range(N):
            z0 = (hc * s_rho[k] + h * Cs_r[k]) / (hc + h)
            z_r[k, :] = zeta + (zeta + h) * z0
    return z_r

# --------------------- input arguments -----------------------------------
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

# input for the script
if len(sys.argv)>1:
    tag = sys.argv[-1]
else:
    tag = 'GB-ref'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/trans/' + tag +'/' + ftype + '/'

flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
# flist = flist[90:105]

# dicts for variable clims, colormaps, and other properties
clim = {'temp': [2, 10],
        'salt': [12, 32],
        'zeta': [-0.15, 0.15],
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

long_name = {'temp': r'Temperature [$^{\circ}$C]',
             'salt': r'Salinity [PSU]',
             'zeta': r'Sea Surface Height [m]',
             'wetdry_mask_rho': r'Wet Dry Mask',
             'dye_01': r'Dye 1',
             'dye_02': r'Dye 2',
             'dye_03': r'Dye 3',
             'tke': r'TKE [m$^2$s$^{-2}$]',
             'gls': r'GLS [m$^3$s$^{-2}$]',
            }

var_omega = ['tke', 'gls']
var_log = ['tke', 'gls']

# load grid info
grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho
ang = grd.hgrid.angle_rho
msk = grd.hgrid.mask_rho
h = grd.vgrid.h
zlev = grd.vgrid.N

uvar = 'u'
vvar = 'v'
wvar = 'w'

# river discharge file
river_file = out_dir + 'frc/' + grd.name + '_rivers_' + str(my_year) + '_Hill.nc'

# --------------------- transect coordinates ------------------------------
a0 = []
fh = open('../data/a0.txt')
csvr = csv.reader(fh, delimiter=',')
for line in csvr:
    a0.append(line)
fh.close()
a0 = np.array(a0)
a0 = a0.astype(float)

# --------------------- data preparation ----------------------------------
lon_ct = a0[:, 0]
lat_ct = a0[:, 1]

ct_tr = (len(lon_ct)-1)*dd
lon_tr = np.zeros(ct_tr)
lat_tr = np.zeros(ct_tr)

for i in range(len(lon_ct)-1):
    lon_tr[i*dd:(i+1)*dd] = np.linspace(lon_ct[i], lon_ct[i+1], dd+1)[:-1]
    lat_tr[i*dd:(i+1)*dd] = np.linspace(lat_ct[i], lat_ct[i+1], dd+1)[:-1]

# instead of using griddata to find interpolated values, use distance to find the nearest rho point and
# represent the value at (lon_tr, lat_tr).
eta_tr = np.zeros(lat_tr.shape)
xi_tr = np.zeros(lon_tr.shape)

for i in range(len(eta_tr)):
    D2 = (lat-lat_tr[i])**2+(lon-lon_tr[i])**2
    eta_tr[i], xi_tr[i] = np.where(D2==D2.min())

eta_tr = eta_tr.astype(int)
xi_tr = xi_tr.astype(int)
h_tr = h[eta_tr, xi_tr]

# calculate distance
dis = np.zeros(h_tr.size)
for i in range(1, lat_tr.size):
    dis[i] = vincenty((lat_tr[i-1], lon_tr[i-1]),
                      (lat_tr[i], lon_tr[i])
                     ).meters
dis = np.cumsum(dis)
dis = dis/1000  # [km]
dis = np.tile(dis, (zlev, 1))

# load river file
friver = nc.Dataset(river_file, 'r')
river_time = friver.variables['river_time'][:]
t0 = river_time[0]
river_time = river_time - t0 + 1
river = np.abs(friver.variables['river_transport'][:]).sum(axis=1)

# --------------------- velocity plot preparation -------------------------
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

# --------------------- make plots ----------------------------------------
plt.switch_backend('Agg')
try:
    plt.style.use('classic')
except:
    pass

# set the figure background
# f, (axR, ax0, ax1, ax2) = plt.subplots(4, gridspec_kw={'height_ratios':[1, 1, 4, 4]})
f = plt.figure()
gs = GridSpec(12, 10)
axR = f.add_subplot(gs[0:2, :9])
ax0 = f.add_subplot(gs[2:4, :9])
ax1 = f.add_subplot(gs[4:8, :9], sharex=ax0)
ax2 = f.add_subplot(gs[8:12, :9], sharex=ax0)
f.subplots_adjust(hspace=0.05)

# plot river
axR.xaxis.tick_top()
axR.plot(river_time, river)
axR.text(5, 3500, r'Runoff [m$^3$s$^{-1}$]')
axR.set_xlim(0, 366)
axR.set_ylim(0, 5000)
axR.set_yticks([1000, 3000, 5000])

# set axixes
ax0.plot(dis[0, :], np.zeros(dis.shape[1]), '--k', linewidth=0.3)
ax0.set_ylim(-0.15, 0.15)
ax0.set_yticks([-0.1, 0, 0.1])
ax0.set_xlim(dis[0, 0], dis[0, -1])
ax2.set_xticks(range(10, 151, 20))

ax0.text(5, -0.12, r'SSH [m]')

ax1.set_ylabel(r'Depth [m]')
ax2.set_xlabel(r'Distance [km]')

ax1.set_ylim(-5, depth0)
ax1.set_yticks(range(0, depth0, 10))
ax1.invert_yaxis()

ax2.set_ylim(depth0, depth1)
ax2.set_yticks(range(depth0, depth1, 100))
ax2.invert_yaxis()

# plot bathymetry
ax1.fill_between(dis[0, :], -h_tr, depth1, facecolor='lightgrey')
ax2.fill_between(dis[0, :], -h_tr, depth1, facecolor='lightgrey')

for var in varlist:

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
        tag = fn.split('/')[-1].split('.nc')[0]
        print('processing ' + var + ' for ' + fn + ' ...')
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

            z_tr = -get_z(h_tr, grd.vgrid.hc, grd.vgrid.N, grd.vgrid.s_rho, grd.vgrid.Cs_r, zeta_tr, grd.vgrid.Vtrans)
            z_tr2 = z_tr[:, ::dd]

            if plt_uv==1:
                u_tr2 = u[tt, :, eta_tr2.tolist(), xi_tr2.tolist()].T
                v_tr2 = u[tt, :, eta_tr2.tolist(), xi_tr2.tolist()].T
                U_tr2 = u_tr2*np.cos(np.tile(ang_tr2, (40, 1)))-v_tr2*np.sin(np.tile(ang_tr2, (40, 1)))
                w_tr2 = w[tt, :, eta_tr2.tolist(), xi_tr2.tolist()].T

            # make the plot
            pltr = axR.plot([t[tt]/24/60/60-t0+1, t[tt]/24/60/60-t0+1], [0, 5000], 'k')
            pltz = ax0.plot(dis[0, :], zeta_tr, 'k')
            if var in var_log:
                pcm1 = ax1.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
                pcm2 = ax2.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
            else:
                pcm1 = ax1.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
                pcm2 = ax2.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
                pcm1.set_clim(clim_var[0], clim_var[1])
                pcm2.set_clim(clim_var[0], clim_var[1])
            # add colorbar axis handle
            cbar_ax = f.add_axes([0.85, 0.10, 0.02, 0.8])
            cb = f.colorbar(pcm1, cax=cbar_ax)
            cbar_ax.set_ylabel(long_name[var])

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
                qvkey = ax2.quiverkey(qv2, 0.6, 0.2, 1, r'1 ms$^{-1}$')

            f.suptitle(var + '_' + grd.name + '_' + ftype + '_' + ttag)
            f.savefig(fig_dir + var + '_' + grd.name + '_' + ftype + '_' + ttag + '.png')

            pltr.pop(0).remove()
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
                qvkey.remove()

plt.close()
