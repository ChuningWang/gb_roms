"""
extract station data.
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
def filter(data, dt, dt_filter):
    """ filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


# -------------- read in grid info --------------------------
grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)
ang = grd.hgrid.angle_rho
ang = ang[332, 125]

# -------------- read in station data -----------------------
info = {'stn' : 'SEA0850',
        'file_dir': out_dir + 'tef/',
        'sl': 'l',
        'Wp_hrs': 2}

crt = noaa_adcp.get_noaa_current(info)
crt()

time2 = crt.ctime
yearday2 = time2 - \
    nc.date2num(time2[0], 'days since 1900-01-01') + \
    1
z2 = crt.z
U2 = crt.u.T + 1j*crt.v.T
U2 = U2*np.exp(-ang*1j)

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


