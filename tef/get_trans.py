""" extract transect data from history files and save to netcdf file. """

import glob
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.signal import filtfilt

import cmocean
import pyroms
from geopy.distance import vincenty
import ttide

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


# -------------- data preparation ---------------------------
my_year = 2008
ts = 12
grd1 = 'GB_lr'
ftype = 'his'
# xpos0 = 184
# xpos1 = 175
# ypos0 = 126
# ypos1 = 145
xpos0 = 211
xpos1 = 211
ypos0 = 127
ypos1 = 145

if len(sys.argv)>1:
    tag = sys.argv[-1]
else:
    tag = 'GB-ref'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
out_file = out_dir + 'tef/trans_' + \
           str(xpos0) + '_' + str(xpos1) + '_' + \
           str(ypos0) + '_' + str(ypos1) + '.nc'
flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
flist = flist[1:]
# flist = flist[0:5]

xp = abs(xpos1 - xpos0)
yp = abs(ypos1 - ypos0)
pts = max(xp, yp) + 1

xx = np.linspace(xpos0, xpos1, pts)
yy = np.linspace(ypos0, ypos1, pts)
xx = xx.round().astype(int)
yy = yy.round().astype(int)

grd = pyroms.grid.get_ROMS_grid(grd1)
lat_raw = grd.hgrid.lat_rho[xx, yy]
lon_raw = grd.hgrid.lon_rho[xx, yy]
dx_raw = grd.hgrid.dx[xx, yy]
dy_raw = grd.hgrid.dy[xx, yy]

msk = grd.hgrid.mask_rho[xx, yy]
ang = grd.hgrid.angle_rho[xx, yy]
h = grd.vgrid.h[xx, yy]
N = grd.vgrid.N

lat0 = lat_raw[0]
lat1 = lat_raw[-1]
lon0 = lon_raw[0]
lon1 = lon_raw[-1]

lat = np.linspace(lat0, lat1, pts)
lon = np.linspace(lon0, lon1, pts)
dis = np.zeros(pts)

# calculate distance
for i in range(pts):
    dis[i] = vincenty((lat[i], lon[i]),
                      (lat0, lon0)).meters

if xpos0 == xpos1:
    ang_corr = 0
elif ypos0 == ypos1:
    ang_corr = 0.5*np.pi
else:
    # correct the angle
    disx = np.sign(lat1-lat0)*vincenty((lat0, lon0),
                                       (lat1, lon0)).meters
    disy = np.sign(lon1-lon0)*vincenty((lat0, lon0),
                                       (lat0, lon1)).meters
    ang_corr = np.arctan(disx/disy)

ang = ang - ang_corr

tt = len(flist)*ts
N = grd.vgrid.N

# define variables
time = np.zeros(tt)
zeta = np.zeros((tt, pts))
zr = np.zeros((tt, N, pts))
dz = np.zeros((tt, N, pts))

salt = np.zeros((tt, N, pts))
temp = np.zeros((tt, N, pts))

ubar = np.zeros((tt, pts))
vbar = np.zeros((tt, pts))
u = np.zeros((tt, N, pts))
v = np.zeros((tt, N, pts))

# -------------- load data ----------------------------------
for i, fname in enumerate(flist):

    print('Loading file ' + fname)

    fin = nc.Dataset(fname, 'r')
    time[i*ts:(i+1)*ts] = fin.variables['ocean_time'][:]
    ubari = fin.variables['ubar'][:]
    vbari = fin.variables['vbar'][:]
    ui = fin.variables['u'][:]
    vi = fin.variables['v'][:]

    ubari[ubari.mask] = 0
    vbari[vbari.mask] = 0
    ui[ui.mask] = 0
    vi[vi.mask] = 0

    if xpos0 == xpos1:
        zeta[i*ts:(i+1)*ts, :] = fin.variables['zeta'][:, xx[0], yy]
        salt[i*ts:(i+1)*ts, :] = fin.variables['salt'][:, :, xx[0], yy]
        temp[i*ts:(i+1)*ts, :] = fin.variables['temp'][:, :, xx[0], yy]

        ubar[i*ts:(i+1)*ts, :] = 0.5*(ubari[:, xx[0], yy] + ubari[:, xx[0], yy-1])
        vbar[i*ts:(i+1)*ts, :] = 0.5*(vbari[:, xx[0], yy] + vbari[:, xx[0]-1, yy])
        u[i*ts:(i+1)*ts, :, :] = 0.5*(ui[:, :, xx[0], yy] + ui[:, :, xx[0], yy-1])
        v[i*ts:(i+1)*ts, :, :] = 0.5*(vi[:, :, xx[0], yy] + vi[:, :, xx[0]-1, yy])

    elif ypos0 == ypos1:
        zeta[i*ts:(i+1)*ts, :] = fin.variables['zeta'][:, xx, yy[0]]
        salt[i*ts:(i+1)*ts, :] = fin.variables['salt'][:, :, xx, yy[0]]
        temp[i*ts:(i+1)*ts, :] = fin.variables['temp'][:, :, xx, yy[0]]

        ubar[i*ts:(i+1)*ts, :] = 0.5*(ubari[:, xx, yy[0]] + ubari[:, xx, yy[0]-1])
        vbar[i*ts:(i+1)*ts, :] = 0.5*(vbari[:, xx, yy[0]] + vbari[:, xx-1, yy[0]])
        u[i*ts:(i+1)*ts, :, :] = 0.5*(ui[:, :, xx, yy[0]] + ui[:, :, xx, yy[0]-1])
        v[i*ts:(i+1)*ts, :, :] = 0.5*(vi[:, :, xx, yy[0]] + vi[:, :, xx-1, yy[0]])

    else:
        # load vars for this file
        zetai = fin.variables['zeta'][:]
        salti = fin.variables['salt'][:]
        tempi = fin.variables['temp'][:]

        for p in range(pts):
            # zeta
            zeta[i*ts:(i+1)*ts, p] = zetai[:, xx[p], yy[p]]

            # salt, temp
            salt[i*ts:(i+1)*ts, :, p] = salti[:, :, xx[p], yy[p]]
            temp[i*ts:(i+1)*ts, :, p] = tempi[:, :, xx[p], yy[p]]

            # u, v, ubar, vbar
            ubar[i*ts:(i+1)*ts, p] = 0.5*(ubari[:, xx[p], yy[p]] + ubari[:, xx[p], yy[p]-1])
            vbar[i*ts:(i+1)*ts, p] = 0.5*(vbari[:, xx[p], yy[p]] + vbari[:, xx[p]-1, yy[p]])
            u[i*ts:(i+1)*ts, :, p] = 0.5*(ui[:, :, xx[p], yy[p]] + ui[:, :, xx[p], yy[p]-1])
            v[i*ts:(i+1)*ts, :, p] = 0.5*(vi[:, :, xx[p], yy[p]] + vi[:, :, xx[p]-1, yy[p]])

    fin.close()

zr = get_z(h, grd.vgrid.hc, N, grd.vgrid.s_rho,
           grd.vgrid.Cs_r, zeta, grd.vgrid.Vtrans)
zw = get_z(h, grd.vgrid.hc, N+1, grd.vgrid.s_w,
           grd.vgrid.Cs_w, zeta, grd.vgrid.Vtrans)
dz = np.diff(zw, axis=1)

fout = nc.Dataset(out_file, 'w')
fout.createDimension('time')
fout.createDimension('N', N)
fout.createDimension('x', pts)
fout.createVariable('time', 'd', ('time'))
fout.createVariable('xx', 'd', ('x'))
fout.createVariable('yy', 'd', ('x'))
fout.createVariable('h', 'd', ('x'))
fout.createVariable('ang', 'd', ('x'))
fout.createVariable('lat', 'd', ('x'))
fout.createVariable('lon', 'd', ('x'))
fout.createVariable('dis', 'd', ('x'))
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
fout.variables['xx'][:] = xx
fout.variables['yy'][:] = yy
fout.variables['h'][:] = h
fout.variables['ang'][:] = ang
fout.variables['lat'][:] = lat
fout.variables['lon'][:] = lon
fout.variables['dis'][:] = dis
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
