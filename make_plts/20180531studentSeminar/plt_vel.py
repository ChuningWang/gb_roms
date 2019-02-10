from datetime import datetime

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import style
import netCDF4 as nc
from scipy.signal import filtfilt

import ttide
import pyroms
from ocean_toolbox import noaa_adcp

import read_host_info
sv = read_host_info.read_host_info()
frc_dir = sv['out_dir']
data_dir = sv['in_dir']

# ---------------- functionals ---------------------------------------
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


def lfilter(data, dt, dt_filter):
    """ low pass filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


# ---------------- some inputs ---------------------------------------
tbase = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01')
time2 = np.array([datetime(2008, mm+1, 15) for mm in range(12)])
time_tick = nc.date2num(time2, 'days since 1900-01-01') - tbase

grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)

# ---------------- freshwater runoff ---------------------------------
river_file = frc_dir + 'frc/GlacierBay_lr_rivers_2008_Hill.nc'
fin = nc.Dataset(river_file, 'r')
river_time = fin.variables['river_time'][:]-tbase+1
river_transport = np.sum(np.abs(fin.variables['river_transport'][:]), axis=-1)
river_temp = fin.variables['river_temp'][:]
fin.close()

# ---------------- ADCP data -----------------------------------------
info = {'stn' : 'SEA0850',
        'file_dir': data_dir + 'NOAA_ADCP/',
        'sl': 'l',
        'Wp_hrs': -1}
crt = noaa_adcp.get_noaa_current(info)
crt()
crt_time = crt.ctime-tbase+1
pytime = nc.num2date(crt.ctime, 'days since 1900-01-01')
ubar = crt.u.mean(axis=0)
vbar = crt.v.mean(axis=0)

tfit = ttide.t_tide(ubar + 1j*vbar, dt=0.1,
                    stime=pytime[0])

Ut = tfit(pytime)
ut = Ut.real
vt = Ut.imag

ue = ubar - ut
ve = vbar - vt
uef = lfilter(ue, 0.1, 24)
vef = lfilter(ve, 0.1, 24)

fig, axs = plt.subplots(3, 1)
fig.subplots_adjust(right=0.9)
axs1c = axs[1].twinx()
axs[0].set_xlim(0, 366)
axs[0].set_ylim(0, 5000)
axs1c.set_ylim(0, 5000)
for i in [1, 2]:
    axs[i].set_xlim(crt_time[0], crt_time[-1])

axs[0].plot([crt_time[0], crt_time[0]], [0, 5000], 'k')
axs[0].plot([crt_time[-1], crt_time[-1]], [0, 5000], 'k')
axs[0].plot(river_time, river_transport, label='Freshwater Runoff')
axs[0].legend()
axs[0].set_ylabel(r'[m$\cdot$s$^{-3}$]')

axs1c.plot(river_time, river_transport, label='Runoff')
axs1c.set_ylabel(r'[m$\cdot$s$^{-3}$]')
axs[1].plot(crt_time, ue, 'r', alpha=0.3, linewidth=0.5)
axs[1].plot(crt_time, ve, 'y', alpha=0.3, linewidth=0.5)
axs[1].plot(crt_time, uef, 'r', label='$U_{mean}$')
axs[1].plot(crt_time, vef, 'y', label='$V_{mean}$')
axs[1].plot(crt_time, np.zeros(crt_time.shape), '--k', linewidth=0.5)
axs[1].legend(loc=3)
axs[1].set_ylabel(r'[m$\cdot$s$^{-1}$]')

axs[2].plot(crt_time, ut, 'r', linewidth=0.75, label='$U_{tide}$')
axs[2].plot(crt_time, vt, 'y', linewidth=0.75, label='$V_{tide}$')
axs[2].legend(loc=3)
axs[2].set_ylabel(r'[m$\cdot$s$^{-1}$]')
axs[2].plot(crt_time, np.zeros(crt_time.shape), '--k', linewidth=0.5)
axs[2].set_xlabel('Yearday in 2008')

plt.savefig('vel0850.png', dpi=600)
plt.close()
