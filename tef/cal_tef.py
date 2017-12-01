""" calculate total exchange flow. """

import sys

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
def filter(data, dt, dt_filter):
    """ filter. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


# -------------- load data ----------------------------------
xpos0 = 327
xpos1 = 337
ypos0 = 112
ypos1 = 138
# xpos0 = 211
# xpos1 = 211
# ypos0 = 337
# ypos1 = 144
# xpos0 = 184
# xpos1 = 175
# ypos0 = 126
# ypos1 = 145
# t0 = 720 + 52*12
# t1 = 1080 + 60*12
t0 = 0
t1 = -1
ds = .5

out_file = out_dir + 'tef/trans_' + \
           str(xpos0) + '_' + str(xpos1) + '_' + \
           str(ypos0) + '_' + str(ypos1) + '.nc'

# transect data
fin = nc.Dataset(out_file, 'r')
time = fin.variables['time'][t0:t1]
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

[tt, N, pts] = zr.shape

# salinity levels
slev = np.arange(1., 35., ds)
slev2 = 0.5*(slev[1:] + slev[:-1])
Q = np.zeros((tt, len(slev)))
Qin = np.zeros(tt)
Qout = np.zeros(tt)
fin = np.zeros(tt)
fout = np.zeros(tt)

dx = np.zeros(pts)
for i in range(1, pts-1):
    dx[i] = 0.5*(dis[i+1] - dis[i-1])

dx[0] = dis[1] - dis[0]
dx[-1] = dis[-1] - dis[-2]
dx = np.tile(dx, (tt, N, 1))

for t in range(tt):
    for i, sl in enumerate(slev):
        msk = salt[t, :, :] >= sl
        Q[t, i] = np.sum(dz[t, msk]*dx[t, msk]*u[t, msk])

Q = filter(Q, 2, 30)

dQds = np.diff(Q, axis=-1)/ds
for t in range(tt):
    msk = -dQds[t, :] > 0
    Qin[t] = np.sum(-dQds[t, msk]*ds)
    Qout[t] = np.sum(-dQds[t, ~msk]*ds)
    fin[t] = np.sum(slev2[msk]*-dQds[t, msk]*ds)
    fout[t] = np.sum(slev2[~msk]*-dQds[t, ~msk]*ds)
