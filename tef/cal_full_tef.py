""" calculate total exchange flow. """

import sys
import glob
from datetime import datetime
import pdb

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.signal import filtfilt

import pyroms
# import cmocean
# import ttide

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


def get_zw(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zw = np.empty((ti, vgrid.Np) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for n in range(ti):
            for k in range(vgrid.Np):
                z0 = vgrid.hc * vgrid.s_w[k] + (h - vgrid.hc) * vgrid.Cs_w[k]
                zw[n, k, :] = z0 + zeta[n, :] * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for n in range(ti):
            for k in range(vgrid.Np):
                z0 = (vgrid.hc * vgrid.s_w[k] + h * vgrid.Cs_w[k]) / (vgrid.hc + h)
                zw[n, k, :] = zeta[n, :] + (zeta[n, :] + h) * z0

    return zw


def tef(v, salt, dy, dz, slev, t):
    """ calculate total exchange flow. """

    # salinity levels
    ds = np.mean(np.diff(slev))
    Q = np.zeros((len(t), len(slev)))
    dQds = np.zeros((len(t), len(slev)))

    for i, sl in enumerate(slev):
        Qi = np.zeros(salt.shape)
        msks = (salt >= sl).data & ~salt.mask
        Qi[msks] = dz[msks]*v[msks]*dy[msks]
        Q[:, i] = np.sum(Qi, axis=(-2, -1))

    # calculate dQds, use second order
    dQds[:, 2:-2] = (1/ds) * (- (1./12.)*Q[:, 4:] + (8./12.)*Q[:, 3:-1]
                              - (8./12.)*Q[:, 1:-3] + (1./12.)*Q[:, :-4])

    Q = lfilter(Q, 1, 24)
    dQds = lfilter(dQds, 1, 24)

    Qini = np.zeros(dQds.shape)
    Qouti = np.zeros(dQds.shape)

    msk = -dQds > 0

    Qini[msk] = -dQds[msk]*ds
    Qouti[~msk] = -dQds[~msk]*ds
    Qin = np.sum(Qini, axis=-1)
    Qout = np.sum(Qouti, axis=-1)

    return Q, dQds, Qin, Qout


def cal_trans(u, v, grd, istart, iend, jstart, jend):

    i0=istart; j0=jstart; i1=iend;  j1=jend
    istart = float(istart); iend = float(iend)
    jstart = float(jstart); jend = float(jend)

    # Compute equation:  j = aj i + bj
    if istart != iend:
        aj = (jend - jstart ) / (iend - istart)
        bj = jstart - aj * istart
    else:
        aj=10000.
        bj=0.

    # Compute equation:  i = ai j + bi
    if jstart != jend:
        ai = (iend - istart ) / ( jend - jstart )
        bi = istart - ai * jstart
    else:
        ai=10000.
        bi=0.

    # Compute the integer pathway:
    # Chose the strait line with the smallest slope
    if (abs(aj) <=  1 ):
        # Here, the best line is y(x)
        print 'Here, the best line is y(x)'
        # If i1 < i0 swap points and remember it has been swapped
        if i1 <  i0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if j1 >= j0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 1; jst = 0
            norm_u = -1; norm_v = -1

        near = []
        # compute the nearest j point on the line crossing at i
        for i in range(i0,i1+1):
            j = aj*i + bj
            near.append(i + round(j)*1j)

    else:
        # Here, the best line is x(y)
        print 'Here, the best line is x(y)'
        # If j1 < j0 swap points and remember it has been swapped
        if j1 <  j0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if i1 >= i0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 0; jst = 1
            norm_u = 1; norm_v = 1

        near = []
        # compute the nearest i point on the line crossing at j
        for j in range(j0,j1+1):
            i = ai*j + bi
            near.append(round(i) + j*1j)

    inear = np.copy(near)

    n = len(near)
    nn=1

    for k in range(1,n):
        # distance between 2 neighbour points
        d = abs(inear[k] - inear[k-1])

        if ( d > 1 ):
            # intermediate points required if d>1
            neari = interm_pt(inear, k, ai, bi, aj, bj)
            near.insert(nn,neari)
            nn=nn+1

        nn=nn+1

    #get metrics
    dx = grd.hgrid.dx
    dy = grd.hgrid.dy
    z_w = grd.vgrid.z_w[:]
    # average z_w at Arakawa-C u points
    zu = 0.5 * (z_w[:, :, :, :-1] + z_w[:, :, :, 1:])
    dzu = zu[:, 1:, :, :] - zu[:, :-1, :, :]
    # average z_w at Arakawa-C v points
    zv = 0.5 * (z_w[:, :, :-1, :] + z_w[:, :, 1:, :])
    dzv = zv[:, 1:, :, :] - zv[:, :-1, :, :]

    #set u and v to zero where u and v are masked for the sum
    for k in range(u.shape[0]):
        u[k,:] = np.where(grd.hgrid.mask_u == 1, u[k,:], 0)
        v[k,:] = np.where(grd.hgrid.mask_v == 1, v[k,:], 0)

    n = len(near)
    transpu = np.zeros(tt)
    transpv = np.zeros(tt)

    for l in range(0,n-1):
        ii = int(np.real(near[l])); jj = int(np.imag(near[l]))
        for k in range(0, dzu.shape[1]):
            if np.real(near[l]) == np.real(near[l+1]):
                trans = u[:, k, jj+jst, ii] * dy[jj+jst, ii] * \
                        dzu[:, k, jj+jst, ii] * norm_u * norm
                transpu = transpu + trans

            elif np.imag(near[l]) == np.imag(near[l+1]):
                trans = v[:, k, jj, ii+ist] * dx[jj, ii+ist] * \
                        dzv[:, k, jj, ii+ist] * norm_v * norm
                transpv = transpv + trans

    return transpu, transpv

# -------------- my inputs ----------------------------------
grd1 = 'GB_lr'
outputs_dir = model_dir + 'tmpdir_GB-clim/outputs/2008/'
ftype = 'his'
eta = 502
xi = 252
N = 40
ts = 24
slev = np.arange(0, 32, 0.1)

read_data = True

# -------------- data prep ----------------------------------
flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
# flist = flist[150:200]
# flist = flist[150:180]
flist = flist[150:155]
tt = len(flist)*ts
Ns = len(slev)

# -------------- load data ----------------------------------
if read_data:
    # define variables
    time = np.zeros(tt)

    zeta = np.zeros((tt, eta, xi))
    salt = np.zeros((tt, N, eta, xi))
    ubar = np.zeros((tt, eta, xi-1))
    vbar = np.zeros((tt, eta-1, xi))
    u = np.zeros((tt, N, eta, xi-1))
    v = np.zeros((tt, N, eta-1, xi))

    for i, fname in enumerate(flist):

        print('Loading file ' + fname)

        fin = nc.Dataset(fname, 'r')
        time[i*ts:(i+1)*ts] = fin.variables['ocean_time'][:]/24/60/60
        zeta[i*ts:(i+1)*ts, :, :] = fin.variables['zeta'][:]
        salt[i*ts:(i+1)*ts, :, :, :] = fin.variables['salt'][:]
        u[i*ts:(i+1)*ts, :, :, :] = fin.variables['u'][:]
        v[i*ts:(i+1)*ts, :, :, :] = fin.variables['v'][:]
        fin.close()

# -------------- load grid info -----------------------------
grd = pyroms.grid.get_ROMS_grid(grd1, zeta=zeta)
transpu = np.zeros((tt, 30))
transpv = np.zeros((tt, 30))

istart = 40
iend = 65
# jstart = 450
# jend = 450
jstart = np.arange(450, 480)
jend = np.arange(450, 480)

for i in range(len(jstart)):
    transpu[:, i], transpv[:, i] = \
        cal_trans(u, v, grd, istart, iend, jstart[i], jend[i])
