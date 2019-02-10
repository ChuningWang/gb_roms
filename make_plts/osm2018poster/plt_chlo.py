from datetime import datetime, timedelta
import sys
import csv
import glob

import numpy as np
from scipy.signal import filtfilt
import matplotlib.pyplot as plt
from matplotlib import gridspec, style
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


# ---------- load grid info, etc -------------------------------------
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
h = grd.vgrid.h
msk = grd.hgrid.mask_rho == 0
zlev = grd.vgrid.N

mm = 7

pytbase = datetime(2008, 1, 1)
tbase = nc.date2num(pytbase, 'days since 1900-01-01')
pyt0 = datetime(2008, mm, 15)
t0 = nc.date2num(pyt0, 'days since 1900-01-01')

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

# ---------- get ctd data --------------------------------------------
from ocean_toolbox import ctd

info = {'data_dir': '/glade/p/work/chuning/data/',
        'file_dir': '/glade/p/work/chuning/data/',
        'file_name': 'ctd_clim_all.nc',
        'sl': 'l',
        'var': ['salt', 'temp'],
        'clim_station': range(25),
        'clim_deep_interp': 'no',
        'filter': 'no'}
gb_ctd = ctd.ctd(info)
gb_ctd()
stn_list = [21, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 24]
gb_ctd.get_trans_clim(['salt', 'temp'], stn_list)
lat_ctd = gb_ctd.trans['lat']
lon_ctd = gb_ctd.trans['lon']
s_ctd = gb_ctd.trans['salt']
t_ctd = gb_ctd.trans['temp']
z_ctd = gb_ctd.trans['z']

dis_ctd = np.zeros(len(stn_list))
idx_i = np.zeros(len(lat_ctd), dtype=np.int)
for i in range(len(lat_ctd)):
    dd = (lon_tr-lon_ctd[i])**2 + (lat_tr-lat_ctd[i])**2
    idx_i[i] = np.argmin(dd)
    dis_ctd[i] = dis[0, idx_i[i]]
