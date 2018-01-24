from datetime import datetime
import sys
import csv

import numpy as np
from scipy.signal import filtfilt
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc

import pyroms
from geopy.distance import vincenty
from cmocean import cm

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

if grd1=='GB_hr':
    tag = 'GlacierBay_hr'
elif grd1=='GB_lr':
    tag = 'GlacierBay_lr'

# ---------- functionals ---------------------------------------------
def get_z(h, hc, N, s_rho, Cs_r, zeta, Vtrans):

    z_r = np.empty((N, len(h)), 'd')
    if Vtrans == 1:
        for k in range(N):
            z0 = hc * s_rho[k] + (h - hc) * Cs_r[k]
            z_r[k, :] = z0 + zeta * (1.0 + z0 / h)
    elif Vtrans == 2 or Vtrans == 4 or Vtrans == 5:
        for  k in range(N):
            z0 = (hc * s_rho[k] + h * Cs_r[k]) / (hc + h)
            z_r[k, :] = zeta + (zeta + h) * z0
    return z_r


def lfilter(data, dt, dt_filter):
    """ low pass filter residual flow. """

    wp = int(float(dt_filter)/dt)
    b = np.ones(wp)/float(wp)
    a = 1
    data_filtered = filtfilt(b, a, data, axis=0)

    return data_filtered


def decomp(phi, dy, dz, t, h, zeta):
    """ decomposate variable """

    nt, nz, ny = phi.shape
    A = np.sum(dy*dz, axis=(1, 2))
    phiint = np.sum(dy*dz*phi, axis=(1, 2))

    dt_hrs = 24*(t[1]-t[0])
    filter_hrs = 24

    Af = lfilter(A, dt_hrs, filter_hrs)
    phiintf = lfilter(phiint, dt_hrs, filter_hrs)
    phi0 = phiintf/Af

    zeta2 = np.tile(np.expand_dims(zeta, 1), (1, nz, 1))
    h2 = np.tile(h, (nt, nz, 1))
    phi02 = phi0[:, np.newaxis, np.newaxis]
    phie = lfilter((zeta2+h2)/h2*phi, dt_hrs, filter_hrs) - \
           np.tile(phi02, (1, nz, ny))
    phit = phi - np.tile(phi02, (1, nz, ny)) - phie

    return phi0, phie, phit


# ---------- load grid info, etc -------------------------------------
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
h = grd.vgrid.h
msk = grd.hgrid.mask_rho
zlev = grd.vgrid.N

tbase = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01')
pyt0 = datetime(2008, 7, 1)
t0 = nc.date2num(pyt0, 'days since 1900-01-01')

# ---------- load river discharge ------------------------------------
fh = nc.Dataset(out_dir + 'frc/' + tag + '_rivers_clim_Hill.nc', 'r')
time = fh.variables['river_time'][:]
trs = fh.variables['river_transport'][:]
fh.close()

time = time-tbase+1
runoff = np.abs(trs).sum(axis=1)

# ---------- load salinity -------------------------------------------
ttag = pyt0.strftime('%Y-%m-%d') + 'T12:00:00'
fname = model_dir + 'tmpdir_GB-ref/outputs/2008/GB-ref_avg_' + ttag + '.nc'
fin = nc.Dataset(fname, 'r')
salt = fin.variables['salt'][:]
temp = fin.variables['temp'][:]
zeta = fin.variables['zeta'][:]
fin.close()

salt_s = salt[0, -1, :, :]

# ---------- get transect data ---------------------------------------
a0 = []
fh = open('../../data/a0.txt')
csvr = csv.reader(fh, delimiter=',')
for line in csvr:
    a0.append(line)
fh.close()
a0 = np.array(a0)
a0 = a0.astype(float)

lon_ct = a0[:, 0]
lat_ct = a0[:, 1]

dd = 3

ct_tr = (len(lon_ct)-1)*dd
lon_tr = np.zeros(ct_tr)
lat_tr = np.zeros(ct_tr)

for i in range(len(lon_ct)-1):
    lon_tr[i*dd:(i+1)*dd] = np.linspace(lon_ct[i], lon_ct[i+1], dd+1)[:-1]
    lat_tr[i*dd:(i+1)*dd] = np.linspace(lat_ct[i], lat_ct[i+1], dd+1)[:-1]

# instead of using griddata to find interpolated values, use distance to find the nearest rho point and
# represent the value at (lon_tr, lat_tr).
eta_tr = np.zeros(lat_tr.shape)
xi_tr = np.zeros(lon_tr.shape)

for i in range(len(eta_tr)):
    D2 = (lat-lat_tr[i])**2+(lon-lon_tr[i])**2
    eta_tr[i], xi_tr[i] = np.where(D2==D2.min())

eta_tr = eta_tr.astype(int)
xi_tr = xi_tr.astype(int)
h_tr = h[eta_tr, xi_tr]

# calculate distance
dis = np.zeros(h_tr.size)
for i in range(1, lat_tr.size):
    dis[i] = vincenty((lat_tr[i-1], lon_tr[i-1]),
                      (lat_tr[i], lon_tr[i])
                     ).meters
dis = np.cumsum(dis)
dis = dis/1000  # [km]
dis = np.tile(dis, (zlev, 1))

zeta_tr = zeta[0, eta_tr.tolist(), xi_tr.tolist()].T
z_tr = -get_z(h_tr, grd.vgrid.hc, grd.vgrid.N, grd.vgrid.s_rho, grd.vgrid.Cs_r, zeta_tr, grd.vgrid.Vtrans)
salt_tr = salt[0, :, eta_tr.tolist(), xi_tr.tolist()].T

# ---------- get tidal velocity --------------------------------------
xpos0 = 184
xpos1 = 175
ypos0 = 126
ypos1 = 145
file_in = out_dir + 'tef/trans_' + \
          str(xpos0) + '_' + str(xpos1) + '_' + \
          str(ypos0) + '_' + str(ypos1) + '.nc'

# transect data
fin = nc.Dataset(file_in, 'r')
time_tr2 = fin.variables['time'][:]
xx_tr2 = fin.variables['xx'][:]
yy_tr2 = fin.variables['yy'][:]
h_tr2 = fin.variables['h'][:]
ang_tr2 = fin.variables['ang'][:]
lat_tr2 = fin.variables['lat'][:]
lon_tr2 = fin.variables['lon'][:]
dis_tr2 = fin.variables['dis'][:]
zeta_tr2 = fin.variables['zeta'][:]
zr_tr2 = -fin.variables['zr'][:]
dz_tr2 = fin.variables['dz'][:]
salt_tr2 = fin.variables['salt'][:]
temp_tr2 = fin.variables['temp'][:]
ubar_tr2 = fin.variables['ubar'][:]
vbar_tr2 = fin.variables['vbar'][:]
u_tr2 = fin.variables['u'][:]
v_tr2 = fin.variables['v'][:]
fin.close()

time_tr2 = time_tr2/24./3600.
yearday_tr2 = time_tr2-tbase+1
dy_tr2 = np.diff(dis_tr2).mean()
v0, ve, vt = decomp(v_tr2, dy_tr2, dz_tr2, time_tr2, h_tr2, zeta_tr2)
vt2 = (vt*dz_tr2*dy_tr2).mean(axis=(1, 2))
yearday = np.arange(1, 367)
vtmax = np.zeros(yearday.shape)*np.nan
for i, yd in enumerate(yearday):
    msk = (yearday_tr2 > yd) & (yearday_tr2 < yd+1)
    if msk.sum()>0:
        vtmax[i] = vt2[msk].max()

# ---------- make plots ----------------------------------------------

lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

smin = 12
smax = 32

fig = plt.figure()
gs = gridspec.GridSpec(3, 2)
axz = plt.subplot(gs[:2, 0])
axr = plt.subplot(gs[-1, 0])
axt = plt.subplot(gs[-1, -1])
axtr1 = plt.subplot(gs[0, 1])
axtr2 = plt.subplot(gs[1, 1])

axtr1.set_xlim(dis[0, 0], dis[0, -1])
axtr1.set_ylim(50, 0)
axtr2.set_xlim(dis[0, 0], dis[0, -1])
axtr2.set_ylim(450, 50)
axtr1.fill_between(dis[0, :], h_tr, 450, facecolor='lightgrey')
axtr2.fill_between(dis[0, :], h_tr, 450, facecolor='lightgrey')

axr.plot(time, runoff)
axt.plot(yearday, vtmax)

pcm1 = axtr1.contourf(dis, z_tr, salt_tr, np.linspace(smin, smax, 41),
                      vmin=smin, vmax=smax, extend='both', cmap=cm.haline)
pcm2 = axtr2.contourf(dis, z_tr, salt_tr, np.linspace(smin, smax, 41),
                      vmin=smin, vmax=smax, extend='both', cmap=cm.haline)

# add colorbar
cbar_ax = fig.add_axes([0.92, 0.10, 0.02, 0.8])
cb = fig.colorbar(pcm1, cax=cbar_ax, ticks=np.linspace(smin, smax, 11))
cbar_ax.set_ylabel(r'Salinity [PSU]', fontsize=7)
cb.ax.tick_params(labelsize=7)

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='i', ax=axz)

m.fillcontinents(color='lightgrey', alpha=0.5)
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),
                     labels=[0,0,0,1], fontsize=6, linewidth=.2)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),
                     labels=[1,0,0,0], fontsize=6, linewidth=.2)

x, y = m(lon, lat)
m.contourf(x, y, salt_s, np.linspace(smin, smax, 41),
           vmin=smin, vmax=smax, extend='both', cmap=cm.haline)

plt.savefig(out_dir + 'figs/osm2018/snap.png')
plt.close()

