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
        temp[i*ts:(i+1)*ts, :, :] = fin.variables['temp'][:, :, xpos, ypos0:ypos1]

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
itrans_sl = 'l'
xpos = 211
ypos0 = 127
ypos1 = 145
t0 = 720 + 52*12
t1 = 1080 + 60*12
yidx = 10

xx = ypos1 - ypos0
tt = t1 - t0
grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho[xpos, ypos0:ypos1]
lon = grd.hgrid.lon_rho[xpos, ypos0:ypos1]
ang = grd.hgrid.angle_rho[xpos, ypos0:ypos1]
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
zeta = fin.variables['zeta'][:, yidx]
zr = -fin.variables['zr'][:, :, yidx]
dz = fin.variables['dz'][:, :, yidx]
salt = fin.variables['salt'][:, :, yidx]
temp = fin.variables['temp'][:, :, yidx]
ubar = fin.variables['ubar'][:, yidx]
vbar = fin.variables['vbar'][:, yidx]
u = fin.variables['u'][:, :, yidx]
v = fin.variables['v'][:, :, yidx]
fin.close()

ur = u - np.tile(np.expand_dims(ubar, axis=1), (1, 40))
vr = v - np.tile(np.expand_dims(vbar, axis=1), (1, 40))

# time convertion
time = time/24./3600.
pytime = nc.num2date(time, 'days since 1900-01-01')
yearday = time - \
    nc.date2num(datetime(pytime[0].year, 1, 1), 'days since 1900-01-01') + \
    1

# -------------- harmonic analysis --------------------------
# perform harmonic analysis
tfit = ttide.t_tide(ubar + 1j*vbar, dt=2,
                    stime=pytime[0], out_style=None, errcalc='wboot')
Ut = tfit(pytime)
ut = Ut.real
vt = Ut.imag

ue = ubar - ut
ve = vbar - vt

# only use a piece of the data
time = time[t0:t1]
pytime = pytime[t0:t1]
yearday = yearday[t0:t1]
zeta = zeta[t0:t1]
zr = zr[t0:t1]
dz = dz[t0:t1]

u = u[t0:t1, :]
v = v[t0:t1, :]
ubar = ubar[t0:t1]
vbar = vbar[t0:t1]

ut = ut[t0:t1]
vt = vt[t0:t1]
ue = ue[t0:t1]
ve = ve[t0:t1]
ur = ur[t0:t1, :]
vr = vr[t0:t1, :]

# -------------- read in station data -----------------------
info = {'stn' : 'SEA0847',
        'file_dir': out_dir + 'tef/',
        'sl': 'l',
        'Wp_hrs': 2}

crt = noaa_adcp.get_noaa_current(info)
crt()

time2 = crt.ctime
yearday2 = time2 - \
    nc.date2num(datetime(pytime[0].year, 1, 1), 'days since 1900-01-01') + \
    1
z2 = crt.z
U2 = crt.u.T + 1j*crt.v.T
U2 = U2*np.exp(-ang[10]*1j)

u2 = U2.real
v2 = U2.imag
ubar2 = u2.mean(axis=1)
vbar2 = v2.mean(axis=1)

ur2 = u2 - np.tile(np.expand_dims(ubar2, axis=1), (1, len(z2)))
vr2 = v2 - np.tile(np.expand_dims(vbar2, axis=1), (1, len(z2)))

# -------------- harmonic analysis --------------------------
tfit = ttide.t_tide(ubar2 + 1j*vbar2, dt=0.1,
                    stime=nc.num2date(time2[0], 'days since 1900-01-01'),
                    out_style=None)
Ut = tfit(nc.num2date(time2, 'days since 1900-01-01'))
ut2 = Ut.real
vt2 = Ut.imag

ue2 = ubar2 - ut2
ve2 = vbar2 - vt2

# -------------- filter bcl ---------------------------------
vef = filter(ve, 2, 30)
ve2f = filter(ve2, 0.1, 30)
vrf = filter(vr, 2, 30)
vr2f = filter(vr2, 0.1, 30)

saltf = filter(salt, 2, 30)
tempf = filter(temp, 2, 30)

# -------------- make plots ---------------------------------
fig, (axt, axe, axbcl1, axbcl2, axbcl3, axbcl4) = \
        plt.subplots(6, 1, sharex=True,
                     gridspec_kw={'height_ratios':[1, 1, 1, 1, 1, 1]})
fig.subplots_adjust(right=0.83, hspace=.05)
axt.set_xlim(yearday[0], yearday[-1])

axt.set_ylim(-2.7, 2.7)
axe.set_ylim(-0.15, 0.15)

axt.set_yticks([-2., 0., 2.])
axe.set_yticks([-0.1, 0., 0.1])

axbcl1.set_ylim(-5, 65)
axbcl2.set_ylim(-5, 65)
axbcl3.set_ylim(-5, 65)
axbcl4.set_ylim(-5, 65)

axbcl1.invert_yaxis()
axbcl2.invert_yaxis()
axbcl3.invert_yaxis()
axbcl4.invert_yaxis()

axbcl4.set_xlabel(r'Yearday')
axbcl1.set_ylabel(r'Depth [m]')

# axt.set_ylabel(r'Tides [ms$^{-1}$]')
# axe.set_ylabel(r'2D Res [ms$^{-1}$]')

axt.plot(yearday, vt, 'r', linewidth=.7)
axt.plot(yearday2, vt2, 'b', linewidth=.7)

axe.plot(yearday, vef, 'r', linewidth=1)
axe.plot(yearday2, ve2f, 'b', linewidth=1)
axe.plot(yearday, ve, '--r', linewidth=.5)
axe.plot(yearday2, ve2, '--b', linewidth=.5)
axe.legend(['Model', 'Obs'])

axt.plot(yearday, np.zeros(yearday.shape), '--k', linewidth=.3)
axe.plot(yearday, np.zeros(yearday.shape), '--k', linewidth=.3)

pctf1 = axbcl1.contourf(np.tile(np.expand_dims(yearday, axis=1),
                        (1, grd.vgrid.N)),
                        zr, vrf,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

pctf2 = axbcl2.contourf(yearday2, z2, vr2f.T,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

pctf3 = axbcl3.contourf(np.tile(np.expand_dims(yearday, axis=1),
                        (1, grd.vgrid.N)),
                        zr, vr - vrf,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

pctf4 = axbcl4.contourf(yearday2, z2, vr2.T - vr2f.T,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

cbar_ax = fig.add_axes([0.85, 0.10, 0.02, 0.8])
cb = fig.colorbar(pctf1, cax=cbar_ax, ticks=np.linspace(-0.3, 0.3, 13))
cbar_ax.set_ylabel(r'Baroclinic [ms$^{-1}$]')
plt.savefig(out_dir + 'figs/tef/stn_cmp.png', dpi=300)
plt.close()


fig, (axbcl1, axbcl2, axbcl3, axbcl4) = \
        plt.subplots(4, 1, sharex=True, sharey=True)
fig.subplots_adjust(right=0.83, hspace=.05)
axbcl1.set_xlim(yearday[0], yearday[-1])
axbcl1.set_ylim(-5, 65)
axbcl1.invert_yaxis()
axbcl4.set_xlabel(r'Yearday')
axbcl4.set_ylabel(r'Depth [m]')

pctf1 = axbcl1.contourf(np.tile(np.expand_dims(yearday, axis=1),
                        (1, grd.vgrid.N)),
                        zr, saltf,
                        np.linspace(-29., 31., 21.),
                        extend='both', cmap=cmocean.cm.balance)

pctf2 = axbcl2.contourf(np.tile(np.expand_dims(yearday, axis=1),
                        (1, grd.vgrid.N)),
                        zr, salt - saltf,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

pctf3 = axbcl3.contourf(yearday2, z2, tempf.T,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

pctf4 = axbcl4.contourf(yearday2, z2, temp.T - tempf.T,
                        np.linspace(-.3, .3, 13),
                        extend='both', cmap=cmocean.cm.balance)

cbar_ax1 = fig.add_axes([0.85, 0.10, 0.02, 0.8])
cb1 = fig.colorbar(pctf1, cax=cbar_ax1, ticks=np.linspace(-0.3, 0.3, 13))
cbar_ax1.set_ylabel(r'Salinity [PSU]')
plt.savefig(out_dir + 'figs/tef/stn_st_cmp.png', dpi=300)
plt.close()
