"""
extract cross transect data from model output history files.
perform eulerian decomposition for U and S.
"""

import glob
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.signal import filtfilt

import cmocean
import pyroms
import ttide

from ocean_toolbox import noaa_adcp
import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# -------------- functionals --------------------------------
def get_z(h, hc, N, s, Cs, zeta, Vtrans):
    """ get z points from grid and zeta info. """

    ti = zeta.shape[0]
    z = np.empty((ti, N, len(h)), 'd')
    if Vtrans == 1:
        for n in range(ti):
            for k in range(N):
                z0 = hc * s[k] + (h - hc) * Cs[k]
                z[n, k, :] = z0 + zeta[n, :] * (1.0 + z0 / h)
    elif Vtrans == 2 or Vtrans == 4 or Vtrans == 5:
        for n in range(ti):
            for  k in range(N):
                z0 = (hc * s[k] + h * Cs[k]) / (hc + h)
                z[n, k, :] = zeta[n, :] + (zeta[n, :] + h) * z0

    return z


def filter(data, dt, dt_filter):
    """ filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


def extr_itrans(out_file, flist, grd, ts, xpos, ypos0, ypos1):
    """ load 3d S, U, V from history files and save it to out_file """

    tt = len(flist)*ts
    dx = ypos1-ypos0
    N = grd.vgrid.N

    # variables
    h = grd.vgrid.h[xpos, ypos0:ypos1]
    time = np.zeros(tt)
    zeta = np.zeros((tt, dx))
    zr = np.zeros((tt, N, dx))
    dz = np.zeros((tt, N, dx))

    salt = np.zeros((tt, N, dx))
    temp = np.zeros((tt, N, dx))

    ubar = np.zeros((tt, dx))
    vbar = np.zeros((tt, dx))
    u = np.zeros((tt, N, dx))
    v = np.zeros((tt, N, dx))

    for i, fname in enumerate(flist):
        print('Loading file ' + fname)
        fin = nc.Dataset(fname, 'r')
        time[i*ts:(i+1)*ts] = fin.variables['ocean_time'][:]

        # zeta, z, dz
        zeta[i*ts:(i+1)*ts, :] = fin.variables['zeta'][:, xpos, ypos0:ypos1]
        zr[i*ts:(i+1)*ts, :] = get_z(h, grd.vgrid.hc, N, grd.vgrid.s_rho,
                                     grd.vgrid.Cs_r, zeta[i*ts:(i+1)*ts, :], grd.vgrid.Vtrans)
        zw = get_z(h, grd.vgrid.hc, N+1, grd.vgrid.s_w,
                   grd.vgrid.Cs_w, zeta[i*ts:(i+1)*ts, :], grd.vgrid.Vtrans)
        dz[i*ts:(i+1)*ts, :] = np.diff(zw, axis=1)

        # salt, temp
        salt[i*ts:(i+1)*ts, :, :] = fin.variables['salt'][:, :, xpos, ypos0:ypos1]
        temp[i*ts:(i+1)*ts, :, :] = fin.variables['salt'][:, :, xpos, ypos0:ypos1]

        # u, v, ubar, vbar
        ubari = fin.variables['ubar'][:]
        vbari = fin.variables['vbar'][:]
        ui = fin.variables['u'][:]
        vi = fin.variables['v'][:]
        ubari[ubari.mask] = 0
        vbari[vbari.mask] = 0
        ui[ui.mask] = 0
        vi[vi.mask] = 0
        ubar[i*ts:(i+1)*ts, :] = 0.5*(ubari[:, xpos, ypos0:ypos1] + ubari[:, xpos, ypos0-1:ypos1-1])
        vbar[i*ts:(i+1)*ts, :] = 0.5*(vbari[:, xpos, ypos0:ypos1] + vbari[:, xpos-1, ypos0:ypos1])
        u[i*ts:(i+1)*ts, :, :] = 0.5*(ui[:, :, xpos, ypos0:ypos1] + ui[:, :, xpos, ypos0-1:ypos1-1])
        v[i*ts:(i+1)*ts, :, :] = 0.5*(vi[:, :, xpos, ypos0:ypos1] + vi[:, :, xpos-1, ypos0:ypos1])
        fin.close()

    fout = nc.Dataset(out_file, 'w')
    fout.createDimension('time')
    fout.createDimension('N', N)
    fout.createDimension('x', dx)
    fout.createVariable('time', 'd', ('time'))
    fout.createVariable('zeta', 'd', ('time', 'x'))
    fout.createVariable('zr', 'd', ('time', 'N', 'x'))
    fout.createVariable('dz', 'd', ('time', 'N', 'x'))
    fout.createVariable('salt', 'd', ('time', 'N', 'x'))
    fout.createVariable('temp', 'd', ('time', 'N', 'x'))
    fout.createVariable('ubar', 'd', ('time', 'x'))
    fout.createVariable('vbar', 'd', ('time', 'x'))
    fout.createVariable('u', 'd', ('time', 'N', 'x'))
    fout.createVariable('v', 'd', ('time', 'N', 'x'))

    fout.variables['time'][:] = time
    fout.variables['zeta'][:] = zeta
    fout.variables['zr'][:] = zr
    fout.variables['dz'][:] = dz
    fout.variables['salt'][:] = salt
    fout.variables['temp'][:] = temp
    fout.variables['ubar'][:] = ubar
    fout.variables['vbar'][:] = vbar
    fout.variables['u'][:] = u
    fout.variables['v'][:] = v
    fout.close()

    return None

# -------------- extract data -------------------------------
my_year = 2008
ts = 12
grd1 = 'GB_lr'
ftype = 'his'
itrans_sl = 's'
xpos = 211
ypos0 = 127
ypos1 = 145
t0 = 720
t1 = 1080

xx = ypos1 - ypos0
tt = t1 - t0
grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho[xpos, ypos0:ypos1]
lon = grd.hgrid.lon_rho[xpos, ypos0:ypos1]
h = grd.vgrid.h[xpos, ypos0:ypos1]
dx = grd.hgrid.dx[xpos, ypos0:ypos1]

if len(sys.argv)>1:
    tag = sys.argv[-1]
else:
    tag = 'GB-ref'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
out_file = out_dir + 'tef/i_trans_' + \
           str(xpos) + '_' + str(ypos0) + '_' + str(ypos1) + '.nc'
flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
flist = flist[1:]
# flist = flist[1:5]

# transect data
if itrans_sl == 's':
    extr_itrans(out_file, flist, grd, ts, xpos, ypos0, ypos1)

fin = nc.Dataset(out_file, 'r')
time = fin.variables['time'][:]
zeta = fin.variables['zeta'][:]
zr = -fin.variables['zr'][:]
dz = fin.variables['dz'][:]
salt = fin.variables['salt'][:]
temp = fin.variables['temp'][:]
ubar = fin.variables['ubar'][:]
vbar = fin.variables['vbar'][:]
u = fin.variables['u'][:]
v = fin.variables['v'][:]
fin.close()

ur = u - np.tile(np.expand_dims(ubar, axis=1), (1, 40, 1))
vr = v - np.tile(np.expand_dims(vbar, axis=1), (1, 40, 1))

# time convertion
time = time/24./3600.
pytime = nc.num2date(time, 'days since 1900-01-01')
yearday = time - \
    nc.date2num(datetime(pytime[0].year, 1, 1), 'days since 1900-01-01') + \
    1

# -------------- harmonic analysis --------------------------
# perform harmonic analysis
ut = np.zeros(zeta.shape)
vt = np.zeros(zeta.shape)
for i in range(xx):
    tfit = ttide.t_tide(ubar[:, i] + 1j*vbar[:, i], dt=2,
                        stime=pytime[0], out_style=None, errcalc='wboot')
    Ut = tfit(pytime)
    ut[:, i] = Ut.real
    vt[:, i] = Ut.imag

ue = ubar - ut
ve = vbar - vt

# only use a piece of the data
time = time[t0:t1]
pytime = pytime[t0:t1]
yearday = yearday[t0:t1]
zeta = zeta[t0:t1, :]
zr = zr[t0:t1, :]
dz = dz[t0:t1, :]

u = u[t0:t1, :, :]
v = v[t0:t1, :, :]
ubar = ubar[t0:t1, :]
vbar = vbar[t0:t1, :]

ut = ut[t0:t1, :]
vt = vt[t0:t1, :]
ue = ue[t0:t1, :]
ve = ve[t0:t1, :]
ur = ur[t0:t1, :, :]
vr = vr[t0:t1, :, :]

# -------------- area weighted average ----------------------
Ut = np.zeros(tt)
Vt = np.zeros(tt)
Ue = np.zeros(tt)
Ve = np.zeros(tt)

for i in range(tt):
    Ut[i] = np.sum((h-zeta[i, :])*dx*ut[i, :]) / \
        np.sum((h-zeta[i, :])*dx)
    Vt[i] = np.sum((h-zeta[i, :])*dx*vt[i, :]) / \
        np.sum((h-zeta[i, :])*dx)
    Ue[i] = np.sum((h-zeta[i, :])*dx*ue[i, :]) / \
        np.sum((h-zeta[i, :])*dx)
    Ve[i] = np.sum((h-zeta[i, :])*dx*ve[i, :]) / \
        np.sum((h-zeta[i, :])*dx)

# -------------- read in station data -----------------------
info = {'stn' : 'SEA0847',
        'file_dir': out_dir + 'tef/',
        'sl': 'l',
        'Wp_hrs': 2}

crt = noaa_adcp.get_noaa_current(info)
crt()

# -------------- make plots ---------------------------------
# fig, (axt, axe, ax) = plt.subplots(3, gridspec_kw={'height_ratios': [2, 2, 6]})
# ax.set_position([0.15, 0.1, 0.7, 0.4])
# axe.set_position([0.15, 0.6, 0.7, 0.15])
# axt.set_position([0.15, 0.76, 0.7, 0.15])
# 
# ax.set_xlabel(r'Longitude [$^{\circ}$W]')
# ax.set_ylabel(r'Depth [m]')
# ax.set_xlim(lon.min(), lon.max())
# ax.set_ylim(-5, 60)
# ax.set_xticks([-136.08, -136.04, -136.00])
# ax.set_xticklabels(['136.08', '136.04', '136.00'])
# ax.invert_yaxis()
# ax.fill_between(lon, -5, 60, facecolor='lightgrey')
# 
# axe.set_xlabel(r'Yearday')
# axe.set_ylabel(r'$U_e$ [ms$^{-1}$]')
# axt.set_ylabel(r'$U_f$ [ms$^{-1}$]')
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
# # for i in range(2):
#     pltue = axe.plot([yearday[i], yearday[i]], [-0.05, 0.05], 'r')
#     pltut = axt.plot([yearday[i], yearday[i]], [-2.5, 2.5], 'r')
# 
#     pct = ax.contour(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], ur[i, :, :], [0],
#                      colors='k')
#     pctf = ax.contourf(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], ur[i, :, :],
#                        np.linspace(-0.5, 0.5, 51), extend='both',
#                        cmap=cmocean.cm.balance)
#     cbar_ax = fig.add_axes([0.87, 0.10, 0.02, 0.8])
#     cb = fig.colorbar(pctf, cax=cbar_ax, ticks=np.linspace(-0.5, 0.5, 11))
#     cbar_ax.set_ylabel(r'$U$ [ms$^{-1}$]')
# 
#     ttag = pytime[i].strftime("%Y-%m-%d_%H:%M:%S")
#     fig.suptitle(ttag)
#     # save fig
#     fig.savefig(out_dir + 'figs/tef/itrans_' + ttag + '.png')
# 
#     # remove plot objects
#     pltue.pop(0).remove()
#     pltut.pop(0).remove()
#     for cc in pct.collections:
#         cc.remove()
#     for cc in pctf.collections:
#         cc.remove()
# 
#     fig.delaxes(cbar_ax)
# 
# plt.close()
