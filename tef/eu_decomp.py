"""
extract cross transect data from model output history files.
perform eulerian decomposition for U and S.
"""

import glob
import sys

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import pyroms
import ttide as tt

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
    fout.createVariable('time', 'f', ('time'))
    fout.createVariable('zeta', 'f', ('time', 'x'))
    fout.createVariable('zr', 'f', ('time', 'N', 'x'))
    fout.createVariable('dz', 'f', ('time', 'N', 'x'))
    fout.createVariable('salt', 'f', ('time', 'N', 'x'))
    fout.createVariable('temp', 'f', ('time', 'N', 'x'))
    fout.createVariable('ubar', 'f', ('time', 'x'))
    fout.createVariable('vbar', 'f', ('time', 'x'))
    fout.createVariable('u', 'f', ('time', 'N', 'x'))
    fout.createVariable('v', 'f', ('time', 'N', 'x'))

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
xpos = 210
ypos0 = 127
ypos1 = 145
t0 = 720
t1 = 2160

dx = ypos1 - ypos0
grd = pyroms.grid.get_ROMS_grid(grd1)

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

# transect data
if itrans_sl == 's':
    extr_itrans(out_file, flist, grd, ts, xpos, ypos0, ypos1)
elif itrans_sl == 'l':
    fin = nc.Dataset(out_file, 'r')
    time = fin.variables['time'][:]
    zeta = fin.variables['zeta'][:]
    zr = fin.variables['zr'][:]
    dz = fin.variables['dz'][:]
    salt = fin.variables['salt'][:]
    temp = fin.variables['temp'][:]
    ubar = fin.variables['ubar'][:]
    vbar = fin.variables['vbar'][:]
    u = fin.variables['u'][:]
    v = fin.variables['v'][:]
    fin.close()

# harmonic analysis
utide = np.zeros(zeta.shape)
vtide = np.zeros(zeta.shape)
for i in range(dx):
    tfit = tt.t_tide(ubar[:, i] + 1j*vbar[:, i], dt=2,
                     out_style=None, errcalc='wboot')
    utide[:, i] = tfit(time).real
    vtide[:, i] = tfit(time).imag

ures = ubar - utide
vres = vbar - vtide

ubar = ubar[t0:t1, :]
vbar = vbar[t0:t1, :]
utide = utide[t0:t1, :]
vtide = vtide[t0:t1, :]
ures = ures[t0:t1, :]
vres = vres[t0:t1, :]
u = u[t0:t1, :, :]
v = v[t0:t1, :, :]

