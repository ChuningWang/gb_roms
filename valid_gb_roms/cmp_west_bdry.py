import numpy as np
import netCDF4 as nc
from datetime import datetime
from gb_toolbox import gb_current
from scipy.signal import buttord, butter, filtfilt

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

ctime2008 = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01 00:00:00')
ctime2010 = nc.date2num(datetime(2010, 1, 1), 'days since 1900-01-01 00:00:00')

# -----------------------------------------------------------------------
# read ADCP data
stn_list = ['SEA1008', 'SEA1009', 'SEA1010']
bdate_list = ['20100525', '20100625', '20100725']
edate_list = ['20100625', '20100725', '20100815']

stn_list = ['SEA0845', 'SEA0846', 'SEA0847', 'SEA0848', 'SEA0849', 'SEA0850']
bdate_list = ['20080810']
edate_list = ['20080910']

t = {}
u = {}
v = {}
z = {}
lat = {}
lon = {}

for stn in stn_list:
    for nct in range(len(bdate_list)):
        bdate = bdate_list[nct]
        edate = edate_list[nct]
        info = {'stn' : stn,
                'bdate' : bdate,
                'edate' : edate,
                'filename': '/glade/p/work/chuning/data/NOAA_ADCP/'+stn+'_'+bdate+'_'+edate+'.nc',
                'sl': 'l',
                'Wp_hrs': -1}

        crt = gb_current.get_noaa_current(info)
        crt()

        if nct==0:
            t[stn] = crt.ctime
            u[stn] = crt.u
            v[stn] = crt.v
            z[stn] = crt.z
            lat[stn] = crt.info['lat']
            lon[stn] = crt.info['lon']
        else:
            t[stn] = np.concatenate((t[stn], crt.ctime))
            u[stn] = np.concatenate((u[stn], crt.u), axis=1)
            v[stn] = np.concatenate((v[stn], crt.v), axis=1)

# -----------------------------------------------------------------------
# read bc data
fh = nc.Dataset(out_dir+'bc_ic/GlacierBay_usgs_bdry_2008_SODA3.3.1_0.25.nc', 'r')
u_west = fh.variables['u_west'][:]
v_west = fh.variables['v_west'][:]
zeta_west = fh.variables['zeta_west'][:]
ocean_time = fh.variables['ocean_time'][:]
lat_west = fh.variables['lat_psi'][:]
lon_west = fh.variables['lon_psi'][:]
h = fh.variables['h'][:]
ang = fh.variables['angle'][:]
s_rho = fh.variables['s_rho'][:]
fh.close()

lat_west = lat_west[:, 0]
lon_west = lon_west[:, 0]
h = h[:, 0]
ang = ang[:, 0]

u_west = 0.5*(u_west[:, :, 1:]+u_west[:, :, :-1])
zeta_west = 0.5*(zeta_west[:, 1:]+zeta_west[:, :-1])
h = 0.5*(h[1:]+h[:-1])
ang = 0.5*(ang[1:]+ang[:-1])

# -----------------------------------------------------------------------
# station SEA1009
z2n = -z['SEA1009']
u2n = u['SEA1009']
v2n = v['SEA1009']
t2n = t['SEA1009']-ctime2010

def filter(ctime, uraw, vraw, Wp_hrs):
    dt = (ctime[1]-ctime[0])*24.  # hours
    samplefreq = 24./dt  # rad per day
    stopfreq = 24./Wp_hrs
    Wp = stopfreq*2./samplefreq
    Ws = 2.*Wp

    n, Wn = buttord(Wp, Ws, 3, 60)
    b, a = butter(n, Wn)

    us = filtfilt(b, a, uraw)
    vs = filtfilt(b, a, vraw)
    return us, vs

u2n, v2n = filter(t2n, u2n, v2n, 72)
u2n, v2n = u2n.T, v2n.T

u2 = u_west[:, :, 279]
v2 = v_west[:, :, 279]
zeta2 = zeta_west[:, 279]
h2 = h[279]
z2 = s_rho*h2
ang2 = ang[279]

# rotate velocity vector
u2r = u2*np.cos(ang2)+v2*np.sin(ang2)
v2r = -u2*np.sin(ang2)+v2*np.cos(ang2)

# interpolate
u2s = np.zeros((len(ocean_time), len(z2n)))
v2s = np.zeros((len(ocean_time), len(z2n)))

for i in range(len(ocean_time)):
    u2s[i, :] = np.interp(z2n, z2, u2r[i, :])
    v2s[i, :] = np.interp(z2n, z2, v2r[i, :])
t2s = ocean_time-ctime2008
