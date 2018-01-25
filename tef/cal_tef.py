""" calculate total exchange flow. """

import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.signal import filtfilt

import cmocean
import ttide

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# -------------- functionals --------------------------------
def lfilter(data, dt, dt_filter):
    """ low pass filter. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


def tef(v, salt, slev, dy, dz, t, h, zeta):
    """ calculate total exchange flow. """

    # salinity levels
    ds = np.mean(np.diff(slev))
    Q = np.zeros((len(t), len(slev)))
    dQds = np.zeros((len(t), len(slev)))

    for i, sl in enumerate(slev):
        Qi = np.zeros(salt.shape)
        msk = salt >= sl
        Qi[msk] = dz[msk]*v[msk]*dy
        Q[:, i] = np.sum(Qi, axis=(-2, -1))

    # calculate dQds, use second order
    dQds[:, 2:-2] = (1/ds) * (- (1./12.)*Q[:, 4:] + (8./12.)*Q[:, 3:-1]
                              - (8./12.)*Q[:, 1:-3] + (1./12.)*Q[:, :-4])

    dQds = lfilter(dQds, 2, 336)

    Qini = np.zeros(dQds.shape)
    Qouti = np.zeros(dQds.shape)

    msk = -dQds > 0

    Qini[msk] = -dQds[msk]*ds
    Qouti[~msk] = -dQds[~msk]*ds
    Qin = np.sum(Qini, axis=-1)
    Qout = np.sum(Qouti, axis=-1)

    return Q, dQds, Qin, Qout

# -------------- transect coords ----------------------------
xpos0 = np.array([184, 199, 211, 232, 246, 280, 317, 327, 362])
xpos1 = np.array([175, 190, 211, 232, 246, 280, 317, 337, 372])
ypos0 = np.array([126, 128, 127, 134, 131, 128, 115, 112, 106])
ypos1 = np.array([145, 155, 144, 160, 167, 180, 179, 138, 113])

# -------------- data prep ----------------------------------
t0 = 0
t1 = 1000
ds = .1
slev = np.arange(0., 35., ds)
# spring-neap time
spt = np.arange(0, 365, 14.76) + 107
npt = np.arange(0, 365, 14.76) + 107 + 7.38

Qa = []
dQdsa = []
Qina = []
Qouta = []

out_file = out_dir + 'tef/trans_' + \
           str(xpos0[0]) + '_' + str(xpos1[0]) + '_' + \
           str(ypos0[0]) + '_' + str(ypos1[0]) + '.nc'

# load basic info
fin = nc.Dataset(out_file, 'r')
time = fin.variables['time'][t0:t1]
fin.close()

# seconds to days
time = time/24./3600.
pytime = nc.num2date(time, 'days since 1900-01-01')
yearday = time - \
    nc.date2num(datetime(pytime[0].year, 1, 1), 'days since 1900-01-01') + \
    1

# -------------- loop through to calculate tef --------------
for i in range(len(xpos0)):
    out_file = out_dir + 'tef/trans_' + \
               str(xpos0[i]) + '_' + str(xpos1[i]) + '_' + \
               str(ypos0[i]) + '_' + str(ypos1[i]) + '.nc'

    # transect data
    fin = nc.Dataset(out_file, 'r')
    xx = fin.variables['xx'][:]
    yy = fin.variables['yy'][:]
    h = fin.variables['h'][:]
    ang = fin.variables['ang'][:]
    lat = fin.variables['lat'][:]
    lon = fin.variables['lon'][:]
    dis = fin.variables['dis'][:]
    zeta = fin.variables['zeta'][t0:t1, :]
    zr = -fin.variables['zr'][t0:t1, :, :]
    dz = fin.variables['dz'][t0:t1, :, :]
    salt = fin.variables['salt'][t0:t1, :, :]
    temp = fin.variables['temp'][t0:t1, :, :]
    ubar = fin.variables['ubar'][t0:t1, :]
    vbar = fin.variables['vbar'][t0:t1, :]
    u = fin.variables['u'][t0:t1, :, :]
    v = fin.variables['v'][t0:t1, :, :]
    fin.close()

    # -------------- tef analysis -------------------------------
    dy = np.mean(np.diff(dis))
    Q, dQds, Qin, Qout = tef(v, salt, slev, dy, dz, time, h, zeta)
    Q = lfilter(Q, 2, 336)

    Qa.append(Q)
    dQdsa.append(dQds)
    Qina.append(Qin)
    Qouta.append(Qout)

# -------------- analyze transect ---------------------------
sptime = 348
nptime = 348 + 12*7

Qs = np.zeros((len(xpos0), len(slev)))
Qn = np.zeros((len(xpos0), len(slev)))
for i in range(len(xpos0)):
    Qs[i, :] = Qa[i][sptime, :]
    Qn[i, :] = Qa[i][nptime, :]
