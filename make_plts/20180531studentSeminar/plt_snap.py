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

mm = 7

pytbase = datetime(2008, 1, 1)
tbase = nc.date2num(pytbase, 'days since 1900-01-01')
pyt0 = datetime(2008, mm, 15)
t0 = nc.date2num(pyt0, 'days since 1900-01-01')

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

# # ---------- get river data ------------------------------------------
# river_file = out_dir + 'frc/GlacierBay_lr_rivers_clim_Hill.nc'
# tbase = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01')
# fin = nc.Dataset(river_file, 'r')
# river_time = fin.variables['river_time'][:] - tbase + 1
# river_transport = np.sum(np.abs(fin.variables['river_transport'][:]), axis=-1)
# river_temp = fin.variables['river_temp'][:]
# fin.close()

# # ---------- process data --------------------------------------------
# dsep1 = 10
# dsep2 = 50
# 
# tm = nc.date2num(np.array([datetime(2008, i, 15) for i in range(1, 13)]), 'days since 1900-01-01') - tbase
# s7 = s_ctd[:, :, 6]
# s7s = np.nanmean(s7[:dsep1, :], axis=0)
# s7i = np.nanmean(s7[dsep1:dsep2, :], axis=0)
# s7b = np.nanmean(s7[dsep2:, :], axis=0)
# 
# s7_2 = salt_tr[:, :, idx_i[6]]
# z7 = z_tr[:, :, idx_i[6]].mean(axis=0)
# zw7 = zw_tr[:, :, idx_i[6]].mean(axis=0)
# dz7 = -np.diff(zw7)
# msks = z7 <= 15
# mski = (z7 <= 50) & (z7 > 15)
# mskb = z7 > 50
# s7s_2 = np.zeros(12)
# s7i_2 = np.zeros(12)
# s7b_2 = np.zeros(12)
# for i in range(12):
#     s7s_2[i] = np.sum(s7_2[i, msks]*dz7[msks])/dz7[msks].sum()
#     s7i_2[i] = np.sum(s7_2[i, mski]*dz7[mski])/dz7[mski].sum()
#     s7b_2[i] = np.sum(s7_2[i, mskb]*dz7[mskb])/dz7[mskb].sum()
# 
# s4 = s_ctd[:, :, 9]
# s4s = np.nanmean(s4[:dsep1, :], axis=0)
# s4i = np.nanmean(s4[dsep1:dsep2, :], axis=0)
# s4b = np.nanmean(s4[dsep2:, :], axis=0)
# 
# s4_2 = salt_tr[:, :, idx_i[9]]
# z4 = z_tr[:, :, idx_i[9]].mean(axis=0)
# zw4 = zw_tr[:, :, idx_i[9]].mean(axis=0)
# dz4 = -np.diff(zw4)
# msks = z4 <= 15
# mski = (z4 <= 50) & (z4 > 15)
# mskb = z4 > 50
# s4s_2 = np.zeros(12)
# s4i_2 = np.zeros(12)
# s4b_2 = np.zeros(12)
# for i in range(12):
#     s4s_2[i] = np.sum(s4_2[i, msks]*dz4[msks])/dz4[msks].sum()
#     s4i_2[i] = np.sum(s4_2[i, mski]*dz4[mski])/dz4[mski].sum()
#     s4b_2[i] = np.sum(s4_2[i, mskb]*dz4[mskb])/dz4[mskb].sum()
# 
# # ---------- make plots ----------------------------------------------
# style.use('classic')
# depth0 = 50
# depth1 = 450
# dratio = 0.25
# 
# z_tr2 = np.zeros(z_tr.shape)
# mskz = z_tr <= depth0
# z_tr2[mskz] = z_tr[mskz]/depth0*dratio
# z_tr2[~mskz] = (z_tr[~mskz]-depth0)/(depth1-depth0)*(1-dratio)+dratio
# 
# h_tr2 = np.zeros(h_tr.shape)
# mskh = h_tr <= depth0
# h_tr2[mskh] = h_tr[mskh]/depth0*dratio
# h_tr2[~mskh] = (h_tr[~mskh]-depth0)/(depth1-depth0)*(1-dratio)+dratio
# 
# z_ctd2 = np.zeros(z_ctd.shape)
# mskz = z_ctd <= depth0
# z_ctd2[mskz] = z_ctd[mskz]/depth0*dratio
# z_ctd2[~mskz] = (z_ctd[~mskz]-depth0)/(depth1-depth0)*(1-dratio)+dratio
# 
# lat_min = 57.85
# lat_max = 59.25
# lat_0 = 0.5 * (lat_min + lat_max)
# 
# lon_min = -137.5
# lon_max = -135.0
# lon_0 = 0.5 * (lon_min + lon_max)
# 
# smin = 12
# smax = 32
# 
# # fig, (axz, axt) = plt.subplots(1, 2)
# fig = plt.figure()
# gs = gridspec.GridSpec(10, 10)
# axz = plt.subplot(gs[0:6, 0:4])
# axt = plt.subplot(gs[0:5, 5:10])
# axt2 = plt.subplot(gs[5:10, 5:10])
# axr = plt.subplot(gs[6:8, 0:5])
# axr2 = plt.subplot(gs[8:10, 0:5])
# 
# axt.set_xlim(dis[0, 0], dis[0, -1])
# axt.set_ylim(1, 0)
# axt.fill_between(dis[0, :], h_tr2, 1, facecolor='lightgrey')
# axt.plot(dis[0, :], dratio*np.ones(dis[0, :].shape), '--k')
# axt.set_yticks([0, 0.125, 0.25, 0.4375, 0.625, 0.8125, 1])
# axt.set_yticklabels(['0', '25', '50', '150', '250', '350', ''], fontsize=7)
# axt.yaxis.tick_right()
# axt.set_ylabel(r'Depth [m]', fontsize=7)
# axt.yaxis.set_label_position("right")
# 
# axt2.set_xlim(dis[0, 0], dis[0, -1])
# axt2.set_ylim(1, 0)
# axt2.fill_between(dis[0, :], h_tr2, 1, facecolor='lightgrey')
# axt2.plot(dis[0, :], dratio*np.ones(dis[0, :].shape), '--k')
# axt2.set_yticks([0, 0.125, 0.25, 0.4375, 0.625, 0.8125, 1])
# axt2.set_yticklabels(['0', '25', '50', '150', '250', '350', ''], fontsize=7)
# axt2.yaxis.tick_right()
# axt2.set_ylabel(r'Depth [m]', fontsize=7)
# axt2.yaxis.set_label_position("right")
# 
# axt.set_xticks([0, 50, 100, 150])
# axt.set_xticklabels([''])
# axt2.set_xticks([0, 50, 100, 150])
# axt2.set_xticklabels(['0', '50', '100', '150'], fontsize=7)
# axt2.set_xlabel(r'Along Track Distance [km]', fontsize=7)
# 
# axt.contour(dis, z_tr2[mm-1, :, :], salt_tr[mm-1, :, :], [15, 20 ,25, 30], linestyles='--', linewidths=0.5, colors='k')
# axt2.contour(dis_ctd, z_ctd2, s_ctd[:, mm-1, :], [15, 20 ,25, 30], linestyles='--', linewidths=0.5, colors='k')
# 
# pcm = axt.contourf(dis, z_tr2[mm-1, :, :], salt_tr[mm-1, :, :], np.linspace(smin, smax, 41),
#                    vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# 
# pcm2 = axt2.contourf(dis_ctd, z_ctd2, s_ctd[:, mm-1, :], np.linspace(smin, smax, 41),
#                      vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# 
# m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
#             urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
#             resolution='i', ax=axz)
# 
# m.fillcontinents(color='lightgrey', alpha=0.5)
# mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),
#                      labels=[0,0,0,1], fontsize=6, linewidth=.2)
# pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),
#                      labels=[1,0,0,0], fontsize=6, linewidth=.2)
# 
# x, y = m(lon, lat)
# m.contourf(x, y, salt_s[mm-1, :, :], np.linspace(smin, smax, 41),
#            vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# 
# # plot colorbar
# cbar_ax = fig.add_axes([0.440, 0.45, 0.012, 0.4])
# cbar_ax.tick_params(labelsize=7)
# cb = fig.colorbar(pcm, cax=cbar_ax, ticks=np.linspace(smin, smax, 11))
# cbar_ax.set_ylabel(r'Salinity [PSU]', fontsize=7)
# 
# # axr.plot()
# # axr.set_xlim(river_time[0], river_time[-1])
# # axr.plot(river_time, river_transport, 'b')
# 
# axr.set_xlim(0, 366)
# axr.set_xticks(tm)
# axr.set_xticklabels([''])
# axr.set_ylim(20, 32)
# axr.set_yticks([22, 26, 30])
# axr.set_yticklabels(['22', '26', '30'], fontsize=7)
# axr2.set_xlim(0, 366)
# axr2.set_xticks(tm)
# axr2.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], fontsize=7)
# axr2.set_ylim(20, 32)
# axr2.set_yticks([22, 26, 30])
# axr2.set_yticklabels(['22', '26', '30'], fontsize=7)
# axr2.set_ylabel('                   Salinity [PSU]', fontsize=7)
# 
# axr.plot(tm, s4s, lw=0.75, color='b')
# axr.plot(tm, s4i, lw=0.75, color='g')
# axr.plot(tm, s4b, lw=0.75, color='lawngreen')
# axr.plot(tm, s4s_2, '--', lw=0.75, color='b')
# axr.plot(tm, s4i_2, '--', lw=0.75, color='g')
# axr.plot(tm, s4b_2, '--', lw=0.75, color='lawngreen')
# 
# axr2.plot(tm, s7s, lw=0.75, color='b', label='Surface')
# axr2.plot(tm, s7i, lw=0.75, color='g', label= 'Intermediate')
# axr2.plot(tm, s7b, lw=0.75, color='lawngreen', label='Deep')
# axr2.plot(tm, s7s_2, '--', lw=0.75, color='b')
# axr2.plot(tm, s7i_2, '--', lw=0.75, color='g')
# axr2.plot(tm, s7b_2, '--', lw=0.75, color='lawngreen')
# 
# axr2.legend(loc=3, fontsize=5)
# axr.grid('on')
# axr2.grid('on')
# 
# plt.savefig(out_dir + 'figs/osm2018/snap.png', dpi=600)
# plt.close()

