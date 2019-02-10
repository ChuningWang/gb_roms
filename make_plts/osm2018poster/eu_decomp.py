"""
extract cross transect data from model output history files.
perform eulerian decomposition for U and S.
"""

import sys
from datetime import datetime

import numpy as np
from scipy.signal import filtfilt
import matplotlib as mpl
from matplotlib import gridspec, style
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import netCDF4 as nc
from scipy.signal import filtfilt

from cmocean import cm
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# -------------- functionals --------------------------------
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

    dt_hrs = np.round(24*(t[1]-t[0])*10)/10
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


# -------------- extract data -------------------------------
my_year = 2008
ts = 12
grd1 = 'GB_lr'
# xpos0 = 327
# xpos1 = 337
# ypos0 = 112
# ypos1 = 138
# xpos0 = 211
# xpos1 = 211
# ypos0 = 127
# ypos1 = 144
xpos0 = 184
xpos1 = 175
ypos0 = 126
ypos1 = 145
tbase = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01')
t0 = nc.date2num(datetime(2008, 6, 1), 'days since 1900-01-01')
t1 = nc.date2num(datetime(2008, 8, 1), 'days since 1900-01-01')

grd = pyroms.grid.get_ROMS_grid(grd1)
N = grd.vgrid.N

if len(sys.argv)>1:
    tag = sys.argv[-1]
else:
    tag = 'GB-clim'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
out_file = out_dir + 'tef/trans_' + \
           str(xpos0) + '_' + str(xpos1) + '_' + \
           str(ypos0) + '_' + str(ypos1) + '.nc'

# transect data
fin = nc.Dataset(out_file, 'r')
time = fin.variables['time'][:]/24./3600.
mskt = (time >= t0) & (time < t1)
xx = fin.variables['xx'][:]
yy = fin.variables['yy'][:]
h = fin.variables['h'][:]
ang = fin.variables['ang'][:]
lat = fin.variables['lat'][:]
lon = fin.variables['lon'][:]
dis = fin.variables['dis'][:]
zeta = fin.variables['zeta'][mskt, :]
zr = -fin.variables['zr'][mskt, :]
dz = fin.variables['dz'][mskt, :]
salt = fin.variables['salt'][mskt, :]
temp = fin.variables['temp'][mskt, :]
ubar = fin.variables['ubar'][mskt, :]
vbar = fin.variables['vbar'][mskt, :]
u = fin.variables['u'][mskt, :]
v = fin.variables['v'][mskt, :]
fin.close()

# time convertion
time = time[mskt]
yearday = time - tbase + 1

dis = dis/1000

# -------------- S & V decomposation ------------------------
dy = np.diff(dis).mean()*1000
A = np.sum(dy*dz, axis=(1, 2))
v0, ve, vt = decomp(v, dy, dz, time, h, zeta)
s0, se, st = decomp(salt, dy, dz, time, h, zeta)
q0, qe, qt = decomp(v*salt, dy, dz, time, h, zeta)

vtm = (vt*dz).sum(axis=(1, 2))/(dz.sum(axis=(1,2)))
yearday2 = np.arange(yearday[0], yearday[-1])
vtmax = np.zeros(yearday2.shape)
for i, yd in enumerate(yearday2):
    msk = (yearday > yd) & (yearday < yd+1)
    vtmax[i] = vtm[msk].max()

# -------------- prep for fig 1 -----------------------------
tspr = 157.76 + np.arange(0, 60, 14.76)

smin = 26
smax = 30
slevs = np.linspace(smin, smax, 21)
vmin = -0.4
vmax = 0.4
vlevs = np.linspace(vmin, vmax, 21)
dmin = 0
dmax = 80
# vttmin = -2.5
# vttmax = 2.5
vttmin = -2.
vttmax = 2.

dis2 = np.tile(dis, (40, 1))
sf = lfilter(salt, 1, 24)
vf = lfilter(v, 1, 24)
zf = lfilter(zr, 1, 24)
ss = sf[300, :, :]
vs = vf[300, :, :]
zs = zf[300, :, :]
sn = sf[300+24*7, :, :]
vn = vf[300+24*7, :, :]
zn = zf[300+24*7, :, :]

stm = (sf*dz).sum(axis=(1, 2))/(dz.sum(axis=(1, 2)))

# -------------- make plots ---------------------------------
# style.use('classic')
# mpl.rcParams['font.size'] = 7
# 
# fig = plt.figure()
# gs = gridspec.GridSpec(10, 10)
# axv = plt.subplot(gs[0:2, :])
# axs = axv.twinx()
# axs1 = plt.subplot(gs[2:6, 0:5])
# axs2 = plt.subplot(gs[2:6, 5:10])
# axv1 = plt.subplot(gs[6:10, 0:5])
# axv2 = plt.subplot(gs[6:10, 5:10])
# 
# axv.set_xlim(yearday[0], yearday[-1])
# axv.xaxis.tick_top()
# axv.xaxis.set_label_position('top')
# axv.set_ylim(vttmin, vttmax)
# axv.set_xlabel(r'Yearday')
# axv.set_ylabel(r'$V_T$ [m$\cdot$s$^{-1}$]')
# 
# axs.yaxis.tick_right()
# axs.yaxis.set_label_position('right')
# axs.tick_params(axis='y', colors='r')
# axs.set_ylabel('Salinity [PSU]', color='r')
# axs.set_ylim(29, 31)
# 
# axs1.set_xticklabels([''])
# axs2.set_xticklabels([''])
# axs2.set_yticklabels([''])
# axv2.set_yticklabels([''])
# axs1.set_ylabel(r'Depth [m]')
# axv1.set_xlabel(r'Distance [km]')
# axv1.set_ylabel(r'Depth [m]')
# axv2.set_xlabel(r'Distance [km]')
# 
# axs1.fill_between(dis, h, 100, facecolor='lightgrey')
# axs2.fill_between(dis, h, 100, facecolor='lightgrey')
# axv1.fill_between(dis, h, 100, facecolor='lightgrey')
# axv2.fill_between(dis, h, 100, facecolor='lightgrey')
# 
# for i in tspr:
#     axv.text(i, 1.5, 'N')
#     axv.text(i+7.38, 1.5, 'S')
# 
# axs1.set_xlim(dis[0], dis[-1])
# axs1.set_ylim(dmax, dmin)
# axs2.set_xlim(dis[0], dis[-1])
# axs2.set_ylim(dmax, dmin)
# axv1.set_xlim(dis[0], dis[-1])
# axv1.set_ylim(dmax, dmin)
# axv2.set_xlim(dis[0], dis[-1])
# axv2.set_ylim(dmax, dmin)
# 
# axv.plot(yearday, vtm, lw=1, color='lightgrey')
# axv.plot(yearday2+0.5, vtmax, 'k')
# axs.plot(yearday, stm, 'r')
# 
# pcms = axs1.contourf(dis2, zs, ss, slevs,
#                      vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# pcms = axs2.contourf(dis2, zn, sn, slevs,
#                      vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# pcmv = axv1.contourf(dis2, zs, vs, vlevs,
#                      vmin=vmin, vmax=vmax, extend='both', cmap=cm.balance)
# pcmv = axv2.contourf(dis2, zn, vn, vlevs,
#                      vmin=vmin, vmax=vmax, extend='both', cmap=cm.balance)
# 
# # plot colorbar
# cbar_ax = fig.add_axes([0.915, 0.425, 0.015, 0.31])
# cbar_ax.tick_params(labelsize=7)
# cb = fig.colorbar(pcms, cax=cbar_ax, ticks=np.linspace(smin, smax, 9))
# cbar_ax.set_ylabel(r'Salinity [PSU]', fontsize=7)
# 
# cbar_ax = fig.add_axes([0.915, 0.1, 0.015, 0.31])
# cbar_ax.tick_params(labelsize=7)
# cb = fig.colorbar(pcmv, cax=cbar_ax, ticks=np.linspace(vmin, vmax, 9))
# cbar_ax.set_ylabel(r'V [m$\cdot$s$^{-1}$]', fontsize=7)
# 
# axv.text(yearday[25], -1.8, 'a)')
# axs1.text(0.2, 75, 'b) Spring Tide')
# axs2.text(0.2, 75, 'c) Neap Tide')
# axv1.text(0.2, 75, 'd) Spring Tide')
# axv2.text(0.2, 75, 'e) Neap Tide')
# 
# plt.savefig(out_dir + 'figs/osm2018/decomp.png', dpi=600)
# plt.close()

# -------------- make plots ---------------------------------
# tspr = 157.76 + np.arange(0, 60, 14.76)
# 
# smin = -2
# smax = 2
# slevs = np.linspace(smin, smax, 21)
# vmin = -0.4
# vmax = 0.4
# vlevs = np.linspace(vmin, vmax, 21)
# dmin = 0
# dmax = 80
# vttmin = -2.
# vttmax = 2.
# 
# ss = se[300, :, :]
# vs = ve[300, :, :]
# zs = zf[300, :, :]
# sn = se[300+24*7, :, :]
# vn = ve[300+24*7, :, :]
# zn = zf[300+24*7, :, :]
# 
# stm = (sf*dz).sum(axis=(1, 2))/(dz.sum(axis=(1, 2)))
# 
# style.use('classic')
# mpl.rcParams['font.size'] = 7
# 
# fig = plt.figure()
# gs = gridspec.GridSpec(10, 10)
# axv = plt.subplot(gs[0:2, :])
# axs = axv.twinx()
# axs1 = plt.subplot(gs[2:6, 0:5])
# axs2 = plt.subplot(gs[2:6, 5:10])
# axv1 = plt.subplot(gs[6:10, 0:5])
# axv2 = plt.subplot(gs[6:10, 5:10])
# 
# axv.set_xlim(yearday[0], yearday[-1])
# axv.xaxis.tick_top()
# axv.xaxis.set_label_position('top')
# axv.set_ylim(vttmin, vttmax)
# axv.set_xlabel(r'Yearday')
# axv.set_ylabel(r'$V_T$ [m$\cdot$s$^{-1}$]')
# 
# axs.yaxis.tick_right()
# axs.yaxis.set_label_position('right')
# axs.tick_params(axis='y', colors='r')
# axs.set_ylabel('Salinity [PSU]', color='r')
# axs.set_ylim(29, 31)
# 
# axs1.set_xticklabels([''])
# axs2.set_xticklabels([''])
# axs2.set_yticklabels([''])
# axv2.set_yticklabels([''])
# axs1.set_ylabel(r'Depth [m]')
# axv1.set_xlabel(r'Distance [km]')
# axv1.set_ylabel(r'Depth [m]')
# axv2.set_xlabel(r'Distance [km]')
# 
# axs1.fill_between(dis, h, 100, facecolor='lightgrey')
# axs2.fill_between(dis, h, 100, facecolor='lightgrey')
# axv1.fill_between(dis, h, 100, facecolor='lightgrey')
# axv2.fill_between(dis, h, 100, facecolor='lightgrey')
# 
# for i in tspr:
#     axv.text(i, 1.5, 'N')
#     axv.text(i+7.38, 1.5, 'S')
# 
# axs1.set_xlim(dis[0], dis[-1])
# axs1.set_ylim(dmax, dmin)
# axs2.set_xlim(dis[0], dis[-1])
# axs2.set_ylim(dmax, dmin)
# axv1.set_xlim(dis[0], dis[-1])
# axv1.set_ylim(dmax, dmin)
# axv2.set_xlim(dis[0], dis[-1])
# axv2.set_ylim(dmax, dmin)
# 
# axv.plot(yearday, vtm, lw=1, color='lightgrey')
# axv.plot(yearday2+0.5, vtmax, 'k')
# axs.plot(yearday, stm, 'r')
# 
# axs1.contour(dis2, zs, ss, [0], colors='k', linewidths=0.5)
# axs2.contour(dis2, zn, sn, [0], colors='k', linewidths=0.5)
# axv1.contour(dis2, zs, vs, [0], colors='k', linewidths=0.5)
# axv2.contour(dis2, zn, vn, [0], colors='k', linewidths=0.5)
# 
# pcms = axs1.contourf(dis2, zs, ss, slevs,
#                      vmin=smin, vmax=smax, extend='both', cmap=cm.balance)
# pcms = axs2.contourf(dis2, zn, sn, slevs,
#                      vmin=smin, vmax=smax, extend='both', cmap=cm.balance)
# pcmv = axv1.contourf(dis2, zs, vs, vlevs,
#                      vmin=vmin, vmax=vmax, extend='both', cmap=cm.balance)
# pcmv = axv2.contourf(dis2, zn, vn, vlevs,
#                      vmin=vmin, vmax=vmax, extend='both', cmap=cm.balance)
# 
# # plot colorbar
# cbar_ax = fig.add_axes([0.915, 0.425, 0.015, 0.31])
# cbar_ax.tick_params(labelsize=7)
# cb = fig.colorbar(pcms, cax=cbar_ax, ticks=np.linspace(smin, smax, 9))
# cbar_ax.set_ylabel(r'Salinity [PSU]', fontsize=7)
# 
# cbar_ax = fig.add_axes([0.915, 0.1, 0.015, 0.31])
# cbar_ax.tick_params(labelsize=7)
# cb = fig.colorbar(pcmv, cax=cbar_ax, ticks=np.linspace(vmin, vmax, 9))
# cbar_ax.set_ylabel(r'V [m$\cdot$s$^{-1}$]', fontsize=7)
# 
# axv.text(yearday[25], -1.8, 'a)')
# axs1.text(0.2, 75, 'b) Spring Tide')
# axs2.text(0.2, 75, 'c) Neap Tide')
# axv1.text(0.2, 75, 'd) Spring Tide')
# axv2.text(0.2, 75, 'e) Neap Tide')
# 
# plt.savefig(out_dir + 'figs/osm2018/decomp5.png', dpi=600)
# plt.close()

# -------------- make plots 2 -------------------------------
fig, axarr = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(right=0.75, hspace=0.05)
axarr[0].set_xlim(yearday[0], yearday[-1])
axarr[0].set_ylim(0, 35)
axarr[1].set_ylim(-100, 100)
axarr[2].set_ylim(-1100, 0)
axarr[0].set_ylabel(r'$F_e$ [$10^3$ kg$\cdot$s$^{-1}$]')
axarr[1].set_ylabel(r'$F_t$ [$10^3$ kg$\cdot$s$^{-1}$]')
axarr[2].set_ylabel(r'$F_0$ [$10^3$ kg$\cdot$s$^{-1}$]')
axarr[2].set_xlabel('Yearday')

axarr[0].plot(yearday, (dy*dz*se*ve).sum(axis=(1, 2))/1e3, 'k', lw=1)
axarr[1].plot(yearday, np.zeros(yearday.shape), '--k', lw=0.5)
axarr[1].plot(yearday, (dy*dz*st*vt).sum(axis=(1, 2))/1e3, 'lightgrey', lw=1)
axarr[1].plot(yearday, lfilter((dy*dz*st*vt).sum(axis=(1, 2)), 1, 24)/1e3, 'k', lw=1)
axarr[2].plot(yearday, A*v0*s0/1e3, 'k', lw=1, label=r'$Q_0S_0$')
axarr[2].plot(yearday, A*q0/1e3, 'r', lw=1, label=r'$F_S$')

axarr[2].legend(loc=3)
axarr[0].grid('on')
axarr[1].grid('on')
axarr[2].grid('on')

# plot spring-neap time
for i in tspr:
    axarr[0].text(i, 30, 'N')
    axarr[0].text(i+7.38, 30, 'S')
    axarr[1].text(i, 75, 'N')
    axarr[1].text(i+7.38, 75, 'S')
    axarr[2].text(i, -100, 'N')
    axarr[2].text(i+7.38, -100, 'S')

axarr[0].text(yearday[24], 30, 'a)')
axarr[1].text(yearday[24], 75, 'b)')
axarr[2].text(yearday[24], -100, 'c)')

fig.savefig(out_dir + 'figs/osm2018//tef_fig3.png', dpi=500)
plt.close()

# # -------------- S & V decomposation ------------------------
# dy = np.diff(dis).mean()
# A = np.sum(dy*dz, axis=(1, 2))
# v0, ve, vt = decomp(v, dy, dz, time, h, zeta)
# s0, se, st = decomp(salt, dy, dz, time, h, zeta)
# q0, qe, qt = decomp(v*salt, dy, dz, time, h, zeta)
# 
# vtm = vt.mean(axis=(1, 2))
# yearday2 = np.arange(yearday[0], yearday[-1])
# vtmax = np.zeros(yearday2.shape)
# for i, yd in enumerate(yearday2):
#     msk = (yearday > yd) & (yearday < yd+1)
#     vtmax[i] = vtm[msk].max()
# 
# fig, axarr = plt.subplots(3, 1, sharex=True)
# axarr[0].set_xlim(yearday[0], yearday[-1])
# axarr[0].set_ylim(-1.5, 1.5)
# axarr[1].set_ylim(-4, 1)
# axarr[2].set_ylim(28, 32)
# axarr[0].set_ylabel(r'Tidal Amplitude [m$\cdot$s$^{-1}$]')
# axarr[1].set_ylabel(r'$Q_f$ [$\times10^4$ m$^3$s$^{-1}$]')
# axarr[2].set_ylabel(r'$S_0$ [PSU]')
# 
# axarr[0].plot(yearday, vt.mean(axis=(1, 2)), 'k', lw=0.3)
# axarr[0].plot(yearday2+0.5, vtmax, 'b')
# axarr[1].plot(yearday, v0*A/1e4)
# axarr[2].plot(yearday, s0)
# 
# # plot spring-neap time
# spr = np.arange(0, 365, 14.76) + 107
# spr = np.vstack((spr, spr))
# ver = np.ones(spr.shape[-1])
# axarr[0].plot(spr, np.vstack((-100*ver, 100*ver)), '--k', lw=0.3)
# axarr[1].plot(spr, np.vstack((-100*ver, 100*ver)), '--k', lw=0.3)
# axarr[2].plot(spr, np.vstack((-100*ver, 100*ver)), '--k', lw=0.3)
# 
# fig.savefig(out_dir + 'figs/tef/tef_fig1.png', dpi=500)
# plt.close()
# 
# # -------------- spring/neap transect -----------------------
# idxs = np.argmin(np.abs(yearday-166.04))
# idxn = np.argmin(np.abs(yearday-173.42))
# 
# zrs = zr[idxs, :, :]
# zrn = zr[idxn, :, :]
# ves = ve[idxs, :, :]
# ven = ve[idxn, :, :]
# ses = se[idxs, :, :]
# sen = se[idxn, :, :]
# dis2 = np.tile(dis, (40, 1))
# cmin = -0.5
# cmax = 0.5
# smin = -1.
# smax = 1.
# clevs = np.linspace(cmin, cmax, 101)
# clevs2 = clevs[::10]
# slevs = np.linspace(smin, smax, 101)
# 
# fig, axarr = plt.subplots(2, 1, sharex=True, sharey=True)
# fig.subplots_adjust(right=0.85)
# axarr[0].set_xlim(dis[0], dis[-1])
# axarr[0].set_ylim(80, -5)
# axarr[0].plot(dis, h, 'k')
# axarr[1].plot(dis, h, 'k')
# cf0 = axarr[0].contourf(dis2, zrs, ses, slevs,
#                         vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# cf1 = axarr[1].contourf(dis2, zrn, sen, slevs,
#                         vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
# cs0 = axarr[0].contour(dis2, zrs, ves, clevs2,
#                        colors='k', extend='both', linewidths=1)
# cs1 = axarr[1].contour(dis2, zrn, ven, clevs2,
#                        colors='k', extend='both', linewidths=1)
# cs0.clabel()
# cs1.clabel()
# 
# # add colorbar axis handle
# cbar_ax = fig.add_axes([0.87, 0.1, 0.015, 0.8])
# cb = fig.colorbar(cf0, cax=cbar_ax, ticks=np.linspace(smin, smax, 11))
# cbar_ax.set_ylabel(r'$S_e$ [PSU]')
# 
# fig.savefig(out_dir + 'figs/tef/tef_fig2.png', dpi=500)
# plt.close()
# 
# # -------------- decomposed salt flux -----------------------
# fig, axarr = plt.subplots(3, 1, sharex=True)
# axarr[0].set_xlim(yearday[0], yearday[-1])
# axarr[0].set_ylim(-10, 0)
# # axarr[1].set_ylim(-200, 200)
# axarr[1].set_ylim(-2, 2)
# axarr[2].set_ylim(-15, 5)
# axarr[0].set_ylabel(r'$F_e$ [$10^5$ kg$\cdot$s$^{-1}$]')
# axarr[1].set_ylabel(r'$F_t$ [$10^5$ kg$\cdot$s$^{-1}$]')
# axarr[2].set_ylabel(r'$F_0$ [$10^5$ kg$\cdot$s$^{-1}$]')
# 
# axarr[0].plot(yearday, A*qe.mean(axis=(1, 2))/1e5, 'k', lw=0.5)
# # axarr[1].plot(yearday, A*qt.mean(axis=(1, 2))/1e5, 'grey', lw=0.5)
# axarr[1].plot(yearday, lfilter(A*qt.mean(axis=(1, 2))/1e5, dt_hrs, filter_hrs), 'k', lw=0.5)
# axarr[2].plot(yearday, A*q0/1e5, 'k', lw=0.5)
# axarr[2].plot(yearday, A*v0*s0/1e5, 'r', lw=0.5)
# 
# # plot spring-neap time
# axarr[0].plot(spr, np.vstack((-1000*ver, 1000*ver)), '--k', lw=0.3)
# axarr[1].plot(spr, np.vstack((-1000*ver, 1000*ver)), '--k', lw=0.3)
# axarr[2].plot(spr, np.vstack((-1000*ver, 1000*ver)), '--k', lw=0.3)
# 
# fig.savefig(out_dir + 'figs/tef/tef_fig3.png', dpi=500)
# plt.close()

# # -------------- area weighted average ----------------------
# Ut = np.zeros(tt)
# Vt = np.zeros(tt)
# Ue = np.zeros(tt)
# Ve = np.zeros(tt)
# 
# for i in range(tt):
#     Ut[i] = np.sum((h-zeta[i, :])*pts*ut[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
#     Vt[i] = np.sum((h-zeta[i, :])*pts*vt[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
#     Ue[i] = np.sum((h-zeta[i, :])*pts*ue[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
#     Ve[i] = np.sum((h-zeta[i, :])*pts*ve[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)



# # -------------- harmonic analysis --------------------------
# # perform harmonic analysis
# ut = np.zeros(zeta.shape)
# vt = np.zeros(zeta.shape)
# for i in range(pts):
#     tfit = ttide.t_tide(ubar[:, i] + 1j*vbar[:, i], dt=2,
#                         stime=pytime[0], out_style=None, errcalc='wboot')
#     Ut = tfit(pytime)
#     ut[:, i] = Ut.real
#     vt[:, i] = Ut.imag
# 
# ue = ubar - ut
# ve = vbar - vt
# 
# urf = filter(ur, 2, 30)
# vrf = filter(vr, 2, 30)
# urr = ur - urf
# vrr = vr - vrf
# 
# # only use a piece of the data
# time = time[t0:t1]
# pytime = pytime[t0:t1]
# yearday = yearday[t0:t1]
# zeta = zeta[t0:t1, :]
# zr = zr[t0:t1, :]
# dz = dz[t0:t1, :]
# 
# u = u[t0:t1, :, :]
# v = v[t0:t1, :, :]
# ubar = ubar[t0:t1, :]
# vbar = vbar[t0:t1, :]
# 
# ut = ut[t0:t1, :]
# vt = vt[t0:t1, :]
# ue = ue[t0:t1, :]
# ve = ve[t0:t1, :]
# ur = ur[t0:t1, :, :]
# vr = vr[t0:t1, :, :]
# urf = urf[t0:t1, :, :]
# vrf = vrf[t0:t1, :, :]
# urr = urr[t0:t1, :, :]
# vrr = vrr[t0:t1, :, :]
# 
# # -------------- area weighted average ----------------------
# Ut = np.zeros(tt)
# Vt = np.zeros(tt)
# Ue = np.zeros(tt)
# Ve = np.zeros(tt)
# 
# for i in range(tt):
#     Ut[i] = np.sum((h-zeta[i, :])*pts*ut[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
#     Vt[i] = np.sum((h-zeta[i, :])*pts*vt[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
#     Ue[i] = np.sum((h-zeta[i, :])*pts*ue[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
#     Ve[i] = np.sum((h-zeta[i, :])*pts*ve[i, :]) / \
#         np.sum((h-zeta[i, :])*pts)
# 
# # -------------- make plots ---------------------------------
# fig, (axt, axe, ax, ax2) = plt.subplots(4, gridspec_kw={'height_ratios': [2, 2, 4, 4]})
# fig.subplots_adjust(right=0.83, hspace=.05)
# # ax.set_position([0.15, 0.1, 0.7, 0.4])
# # axe.set_position([0.15, 0.6, 0.7, 0.15])
# # axt.set_position([0.15, 0.76, 0.7, 0.15])
# 
# # ax.set_ylabel(r'Depth [m]')
# ax.set_xlim(lon.min(), lon.max())
# ax.set_ylim(-5, h.max())
# ax.invert_yaxis()
# ax.fill_between(lon, -5, h.max(), facecolor='lightgrey')
# 
# ax2.set_xlim(lon.min(), lon.max())
# ax2.set_ylim(-5, h.max())
# ax2.invert_yaxis()
# ax2.fill_between(lon, -5, h.max(), facecolor='lightgrey')
# 
# ax2.set_ylabel(r'Depth [m]')
# ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax2.set_xticks([-136.08, -136.06, -136.04, -136.02])
# ax2.set_xticklabels(['136.08', '136.06', '136.04', '136.02'])
# ax2.set_xlabel(r'Longitude [$^{\circ}$W]')
# 
# axt.xaxis.tick_top()
# axt.set_xlabel(r'Yearday')
# # axe.set_ylabel(r'$U_e$ [ms$^{-1}$]')
# # axt.set_ylabel(r'$U_f$ [ms$^{-1}$]')
# 
# axt.set_xlim(yearday[0], yearday[-1])
# axe.set_xlim(yearday[0], yearday[-1])
# axt.set_ylim(-2.5, 2.5)
# axe.set_ylim(-0.05, 0.05)
# axt.set_yticks([-2., 0, 2.])
# axe.set_yticks([-0.05, 0, 0.05])
# 
# axt.plot(yearday, Vt, 'k')
# axe.plot(yearday, Ve, '--', color='grey', linewidth=.5)
# axe.plot(yearday, filter(Ve, 2, 30), 'k')
# 
# for i in range(len(time)):
# 
#     ttag = pytime[i].strftime("%Y-%m-%d_%H:%M:%S")
#     print('processing ' + ttag + ' ...')
#     pltue = axe.plot([yearday[i], yearday[i]], [-0.05, 0.05], 'r')
#     pltut = axt.plot([yearday[i], yearday[i]], [-2.5, 2.5], 'r')
# 
#     pct = ax.contour(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], vrf[i, :, :], [0],
#                      colors='k')
#     pctf = ax.contourf(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], vrf[i, :, :],
#                        np.linspace(-0.5, 0.5, 51), extend='both',
#                        cmap=cmocean.cm.balance)
# 
#     pct2 = ax2.contour(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], vrr[i, :, :], [0],
#                        colors='k')
#     pctf2 = ax2.contourf(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], vrr[i, :, :],
#                          np.linspace(-0.5, 0.5, 51), extend='both',
#                          cmap=cmocean.cm.balance)
# 
#     cbar_ax = fig.add_axes([0.87, 0.10, 0.02, 0.8])
#     cb = fig.colorbar(pctf, cax=cbar_ax, ticks=np.linspace(-0.5, 0.5, 11))
#     cbar_ax.set_ylabel(r'Velocity [ms$^{-1}$]')
# 
#     fig.suptitle(ttag)
#     # save fig
#     fig.savefig(out_dir + 'figs/tef/trans_' +
#                 str(xpos0) + '_' + str(xpos1) + '_' +
#                 str(ypos0) + '_' + str(ypos1) + '/' + 
#                 'trans_' + ttag + '.png')
# 
#     # remove plot objects
#     pltue.pop(0).remove()
#     pltut.pop(0).remove()
#     for cc in pct.collections:
#         cc.remove()
#     for cc in pctf.collections:
#         cc.remove()
#     for cc in pct2.collections:
#         cc.remove()
#     for cc in pctf2.collections:
#         cc.remove()
# 
#     fig.delaxes(cbar_ax)
# 
# plt.close()
