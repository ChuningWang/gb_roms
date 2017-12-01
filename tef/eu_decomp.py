"""
extract cross transect data from model output history files.
perform eulerian decomposition for U and S.
"""

import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import netCDF4 as nc
from scipy.signal import filtfilt

import cmocean
import pyroms
import ttide

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# -------------- functionals --------------------------------
def filter(data, dt, dt_filter):
    """ filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


# -------------- extract data -------------------------------
my_year = 2008
ts = 12
grd1 = 'GB_lr'
xpos0 = 327
xpos1 = 337
ypos0 = 112
ypos1 = 138
# xpos0 = 211
# xpos1 = 211
# ypos0 = 127
# ypos1 = 144
# xpos0 = 184
# xpos1 = 175
# ypos0 = 126
# ypos1 = 145
t0 = 720 + 52*12
t1 = 1080 + 60*12

tt = t1 - t0
grd = pyroms.grid.get_ROMS_grid(grd1)
N = grd.vgrid.N

if len(sys.argv)>1:
    tag = sys.argv[-1]
else:
    tag = 'GB-ref'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
out_file = out_dir + 'tef/trans_' + \
           str(xpos0) + '_' + str(xpos1) + '_' + \
           str(ypos0) + '_' + str(ypos1) + '.nc'

# transect data
fin = nc.Dataset(out_file, 'r')
time = fin.variables['time'][:]
xx = fin.variables['xx'][:]
yy = fin.variables['yy'][:]
h = fin.variables['h'][:]
ang = fin.variables['ang'][:]
lat = fin.variables['lat'][:]
lon = fin.variables['lon'][:]
dis = fin.variables['dis'][:]
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

pts = len(h)

ur = u - np.tile(np.expand_dims(ubar, axis=1), (1, N, 1))
vr = v - np.tile(np.expand_dims(vbar, axis=1), (1, N, 1))

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
for i in range(pts):
    tfit = ttide.t_tide(ubar[:, i] + 1j*vbar[:, i], dt=2,
                        stime=pytime[0], out_style=None, errcalc='wboot')
    Ut = tfit(pytime)
    ut[:, i] = Ut.real
    vt[:, i] = Ut.imag

ue = ubar - ut
ve = vbar - vt

urf = filter(ur, 2, 30)
vrf = filter(ur, 2, 30)
urr = ur - urf
vrr = vr - vrf

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
urf = urf[t0:t1, :, :]
vrf = vrf[t0:t1, :, :]
urr = urr[t0:t1, :, :]
vrr = vrr[t0:t1, :, :]

# -------------- area weighted average ----------------------
Ut = np.zeros(tt)
Vt = np.zeros(tt)
Ue = np.zeros(tt)
Ve = np.zeros(tt)

for i in range(tt):
    Ut[i] = np.sum((h-zeta[i, :])*pts*ut[i, :]) / \
        np.sum((h-zeta[i, :])*pts)
    Vt[i] = np.sum((h-zeta[i, :])*pts*vt[i, :]) / \
        np.sum((h-zeta[i, :])*pts)
    Ue[i] = np.sum((h-zeta[i, :])*pts*ue[i, :]) / \
        np.sum((h-zeta[i, :])*pts)
    Ve[i] = np.sum((h-zeta[i, :])*pts*ve[i, :]) / \
        np.sum((h-zeta[i, :])*pts)

# -------------- make plots ---------------------------------
fig, (axt, axe, ax, ax2) = plt.subplots(4, gridspec_kw={'height_ratios': [2, 2, 4, 4]})
fig.subplots_adjust(right=0.83, hspace=.05)
# ax.set_position([0.15, 0.1, 0.7, 0.4])
# axe.set_position([0.15, 0.6, 0.7, 0.15])
# axt.set_position([0.15, 0.76, 0.7, 0.15])

# ax.set_ylabel(r'Depth [m]')
ax.set_xlim(lon.min(), lon.max())
ax.set_ylim(-5, h.max())
ax.invert_yaxis()
ax.fill_between(lon, -5, h.max(), facecolor='lightgrey')

ax2.set_xlim(lon.min(), lon.max())
ax2.set_ylim(-5, h.max())
ax2.invert_yaxis()
ax2.fill_between(lon, -5, h.max(), facecolor='lightgrey')

ax2.set_ylabel(r'Depth [m]')
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax2.set_xlabel(r'Longitude [$^{\circ}$W]')

axt.xaxis.tick_top()
axt.set_xlabel(r'Yearday')
axe.set_ylabel(r'$U_e$ [ms$^{-1}$]')
axt.set_ylabel(r'$U_f$ [ms$^{-1}$]')

axt.set_xlim(yearday[0], yearday[-1])
axe.set_xlim(yearday[0], yearday[-1])
axt.set_ylim(-2.5, 2.5)
axe.set_ylim(-0.05, 0.05)
axt.set_yticks([-2., 0, 2.])
axe.set_yticks([-0.05, 0, 0.05])

axt.plot(yearday, Vt, 'k')
axe.plot(yearday, Ve, '--', color='grey', linewidth=.5)
axe.plot(yearday, filter(Ve, 2, 30), 'k')

for i in range(len(time)):

    ttag = pytime[i].strftime("%Y-%m-%d_%H:%M:%S")
    print('processing ' + ttag + ' ...')
    pltue = axe.plot([yearday[i], yearday[i]], [-0.05, 0.05], 'r')
    pltut = axt.plot([yearday[i], yearday[i]], [-2.5, 2.5], 'r')

    pct = ax.contour(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], urf[i, :, :], [0],
                     colors='k')
    pctf = ax.contourf(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], urf[i, :, :],
                       np.linspace(-0.5, 0.5, 51), extend='both',
                       cmap=cmocean.cm.balance)

    pct2 = ax2.contour(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], urr[i, :, :], [0],
                       colors='k')
    pctf2 = ax2.contourf(np.tile(lon, (grd.vgrid.N, 1)), zr[i, :, :], urr[i, :, :],
                         np.linspace(-0.5, 0.5, 51), extend='both',
                         cmap=cmocean.cm.balance)

    cbar_ax = fig.add_axes([0.87, 0.10, 0.02, 0.8])
    cb = fig.colorbar(pctf, cax=cbar_ax, ticks=np.linspace(-0.5, 0.5, 11))
    cbar_ax.set_ylabel(r'$U$ [ms$^{-1}$]')

    fig.suptitle(ttag)
    # save fig
    fig.savefig(out_dir + 'figs/tef/trans_' +
                str(xpos0) + '_' + str(xpos1) + '_' +
                str(ypos0) + '_' + str(ypos1) + '/' + 
                'trans_' + ttag + '.png')

    # remove plot objects
    pltue.pop(0).remove()
    pltut.pop(0).remove()
    for cc in pct.collections:
        cc.remove()
    for cc in pctf.collections:
        cc.remove()
    for cc in pct2.collections:
        cc.remove()
    for cc in pctf2.collections:
        cc.remove()

    fig.delaxes(cbar_ax)

plt.close()
