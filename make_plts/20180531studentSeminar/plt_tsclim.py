from datetime import datetime, timedelta
import sys
import csv
import glob

import numpy as np
from scipy.signal import filtfilt
import matplotlib.pyplot as plt
from matplotlib import gridspec, style
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc

import pyroms
from geopy.distance import vincenty
from cmocean import cm

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

if grd1=='GB_hr':
    tag = 'GlacierBay_hr'
elif grd1=='GB_lr':
    tag = 'GlacierBay_lr'

mm = 7

# ---------- functionals ---------------------------------------------
def get_zr(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zr = np.empty((ti, vgrid.N) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for n in range(ti):
            for k in range(vgrid.N):
                z0 = vgrid.hc * vgrid.s_rho[k] + (h - vgrid.hc) * vgrid.Cs_r[k]
                zr[n, k, :] = z0 + zeta[n, :] * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for n in range(ti):
            for k in range(vgrid.N):
                z0 = (vgrid.hc * vgrid.s_rho[k] + h * vgrid.Cs_r[k]) / (vgrid.hc + h)
                zr[n, k, :] = zeta[n, :] + (zeta[n, :] + h) * z0

    return zr


def get_zw(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zw = np.empty((ti, vgrid.Np) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for n in range(ti):
            for k in range(vgrid.Np):
                z0 = vgrid.hc * vgrid.s_w[k] + (h - vgrid.hc) * vgrid.Cs_w[k]
                zw[n, k, :] = z0 + zeta[n, :] * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for n in range(ti):
            for k in range(vgrid.Np):
                z0 = (vgrid.hc * vgrid.s_w[k] + h * vgrid.Cs_w[k]) / (vgrid.hc + h)
                zw[n, k, :] = zeta[n, :] + (zeta[n, :] + h) * z0

    return zw


def lfilter(data, dt, dt_filter):
    """ low pass filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


def decomp(phi, dy, dz, t, h, zeta):
    """ decomposate variable """

    nt, nz, ny = phi.shape
    A = np.sum(dy*dz, axis=(1, 2))
    phiint = np.sum(dy*dz*phi, axis=(1, 2))

    dt_hrs = 24*(t[1]-t[0])
    filter_hrs = 24

    Af = lfilter(A, dt_hrs, filter_hrs)
    phiintf = lfilter(phiint, dt_hrs, filter_hrs)
    phi0 = phiintf/Af

    zeta2 = np.tile(np.expand_dims(zeta, 1), (1, nz, 1))
    h2 = np.tile(h, (nt, nz, 1))
    phi02 = phi0[:, np.newaxis, np.newaxis]
    phie = lfilter((zeta2+h2)/h2*phi, dt_hrs, filter_hrs) - \
           np.tile(phi02, (1, nz, ny))
    phit = phi - np.tile(phi02, (1, nz, ny)) - phie

    return phi0, phie, phit

# ---------- load grid info, etc -------------------------------------
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
h = grd.vgrid.h
msk = grd.hgrid.mask_rho == 0
zlev = grd.vgrid.N

pytbase = datetime(2008, 1, 1)
tbase = nc.date2num(pytbase, 'days since 1900-01-01')
pyt0 = datetime(2008, mm, 15)
t0 = nc.date2num(pyt0, 'days since 1900-01-01')

# ---------- get transect data ---------------------------------------
a0 = []
fh = open('../../data/a0.txt')
csvr = csv.reader(fh, delimiter=',')
for line in csvr:
    a0.append(line)
fh.close()
a0 = np.array(a0)
a0 = a0.astype(float)

lon_ct = a0[:, 0]
lat_ct = a0[:, 1]

dd = 3

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

# ---------- get ctd data --------------------------------------------
from ocean_toolbox import ctd

info = {'data_dir': '/glade/p/work/chuning/data/',
        'file_dir': '/glade/p/work/chuning/data/',
        'file_name': 'ctd_clim_all.nc',
        'sl': 'l',
        'var': ['salt', 'temp'],
        'clim_station': range(25),
        'clim_deep_interp': 'no',
        'filter': 'no'}
gb_ctd = ctd.ctd(info)
gb_ctd()
stn_list = [21, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 24]
gb_ctd.get_trans_clim(['salt', 'temp'], stn_list)
lat_ctd = gb_ctd.trans['lat']
lon_ctd = gb_ctd.trans['lon']
s_ctd = gb_ctd.trans['salt']
t_ctd = gb_ctd.trans['temp']
z_ctd = gb_ctd.trans['z']

dis_ctd = np.zeros(len(stn_list))
idx_i = np.zeros(len(lat_ctd), dtype=np.int)
for i in range(len(lat_ctd)):
    dd = (lon_tr-lon_ctd[i])**2 + (lat_tr-lat_ctd[i])**2
    idx_i[i] = np.argmin(dd)
    dis_ctd[i] = dis[0, idx_i[i]]

# ---------- get river data ------------------------------------------
river_file = out_dir + 'frc/GlacierBay_lr_rivers_clim_Hill.nc'
tbase = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01')
fin = nc.Dataset(river_file, 'r')
river_time = fin.variables['river_time'][:] - tbase + 1
river_transport = np.sum(np.abs(fin.variables['river_transport'][:]), axis=-1)
river_temp = fin.variables['river_temp'][:]
fin.close()

# ---------- make plots ----------------------------------------------
style.use('classic')
depth0 = 50
depth1 = 450
dratio = 0.25

h_tr2 = np.zeros(h_tr.shape)
mskh = h_tr <= depth0
h_tr2[mskh] = h_tr[mskh]/depth0*dratio
h_tr2[~mskh] = (h_tr[~mskh]-depth0)/(depth1-depth0)*(1-dratio)+dratio

z_ctd2 = np.zeros(z_ctd.shape)
mskz = z_ctd <= depth0
z_ctd2[mskz] = z_ctd[mskz]/depth0*dratio
z_ctd2[~mskz] = (z_ctd[~mskz]-depth0)/(depth1-depth0)*(1-dratio)+dratio

smin = 12
smax = 32
tmin = 4
tmax = 10

fig, axs = plt.subplots(2, 2)
fig.subplots_adjust(wspace=0.15, hspace=0.10)
for i in range(2):
    for j in range(2):
        axs[i, j].set_xlim(dis[0, 0], dis[0, -1])
        axs[i, j].set_ylim(1, 0)
        axs[i, j].fill_between(dis[0, :], h_tr2, 1, facecolor='lightgrey')
        axs[i, j].plot(dis[0, :], dratio*np.ones(dis[0, :].shape), '--k')
        axs[i, j].set_yticks([0, 0.125, 0.25, 0.4375, 0.625, 0.8125, 1])
        axs[i, j].set_xticks([0, 50, 100, 150])
        axs[i, j].set_yticklabels(['0', '25', '50', '150', '250', '350', ''], fontsize=12)
        axs[i, j].set_xticklabels(['0', '50', '100', '150'], fontsize=12)

axs[0, 0].set_ylabel(r'Depth [m]', fontsize=15)
# axs[0, 0].set_xlabel(r'Along Track Distance [km]', fontsize=15)

axs[1, 0].set_ylabel(r'Depth [m]', fontsize=15)
axs[1, 0].set_xlabel(r'Along Track Distance [km]', fontsize=15)

# axs[0, 1].set_ylabel(r'Depth [m]', fontsize=15)
# axs[0, 1].set_xlabel(r'Along Track Distance [km]', fontsize=15)

# axs[1, 1].set_ylabel(r'Depth [m]', fontsize=15)
axs[1, 1].set_xlabel(r'Along Track Distance [km]', fontsize=15)

cts1 = axs[0, 0].contour(dis_ctd, z_ctd2, s_ctd[:, mm-7, :], np.linspace(30, 32, 11),
                         linestyles='--', linewidths=0.5, colors='k')
cts2 = axs[0, 1].contour(dis_ctd, z_ctd2, s_ctd[:, mm-1, :], np.linspace(30, 32, 11),
                         linestyles='--', linewidths=0.5, colors='k')
cts1.clabel([30, 31, 32], fmt='%d')
cts2.clabel([30, 31, 32], fmt='%d')

ctt1 = axs[1, 0].contour(dis_ctd, z_ctd2, t_ctd[:, mm-7, :], np.linspace(4, 6, 11),
                         linestyles='--', linewidths=0.5, colors='w')
ctt2 = axs[1, 1].contour(dis_ctd, z_ctd2, t_ctd[:, mm-1, :], np.linspace(4, 6, 11),
                         linestyles='--', linewidths=0.5, colors='w')
ctt1.clabel([4, 5, 6], fmt='%d')
ctt2.clabel([4, 5, 6], fmt='%d')

pcms1 = axs[0, 0].contourf(dis_ctd, z_ctd2, s_ctd[:, mm-7, :],
                           np.linspace(smin, smax, 41),
                           vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
pcms2 = axs[0, 1].contourf(dis_ctd, z_ctd2, s_ctd[:, mm-1, :],
                           np.linspace(smin, smax, 41),
                           vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
pcmt1 = axs[1, 0].contourf(dis_ctd, z_ctd2, t_ctd[:, mm-7, :],
                           np.linspace(tmin, tmax, 41),
                           vmin=tmin, vmax=tmax, extend='both', cmap=cm.thermal)
pcmt2 = axs[1, 1].contourf(dis_ctd, z_ctd2, t_ctd[:, mm-1, :],
                           np.linspace(tmin, tmax, 41),
                           vmin=tmin, vmax=tmax, extend='both', cmap=cm.thermal)

axs[0, 0].text(130, 0.125, 'Jan')
axs[1, 0].text(130, 0.125, 'Jan')
axs[0, 1].text(130, 0.125, 'Jul')
axs[1, 1].text(130, 0.125, 'Jul')

# plot colorbar
cbar_ax = fig.add_axes([0.920, 0.53, 0.012, 0.35])
cbar_ax.tick_params(labelsize=12)
cb = fig.colorbar(pcms1, cax=cbar_ax, ticks=np.linspace(smin, smax, 6))
cbar_ax.set_ylabel(r'Salinity [PSU]', fontsize=12)

cbar_ax = fig.add_axes([0.920, 0.11, 0.012, 0.35])
cbar_ax.tick_params(labelsize=12)
cb = fig.colorbar(pcmt1, cax=cbar_ax, ticks=np.linspace(tmin, tmax, 7))
cbar_ax.set_ylabel(r'Temperature [$^{\circ}$C]', fontsize=12)

plt.savefig('TS_clim.png', dpi=600)
plt.close()
