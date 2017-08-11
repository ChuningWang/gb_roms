import numpy as np
import netCDF4 as nc
from gb_toolbox import gb_current

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

stn_list = ['SEA1008', 'SEA1009', 'SEA1010']
bdate = '20100601'
edate = '20100701'

t = {}
u = {}
v = {}
lat = {}
lon = {}

for stn in stn_list:
    info = {'stn' : stn,
            'bdate' : bdate,
            'edate' : edate,
            'filename': '/glade/p/work/chuning/data/NOAA_ADCP/'+stn+'.nc',
            'sl': 'l',
            'Wp_hrs': 2}

    crt = gb_current.get_noaa_current(info)
    crt()
    t[stn] = crt.ctime
    u[stn] = crt.u
    v[stn] = crt.v
    lat[stn] = crt.info['lat']
    lon[stn] = crt.info['lon']

fh = nc.Dataset(out_dir+'bc_ic/GlacierBay_usgs_bdry_2008_SODA3.3.1_0.25.nc', 'r')
u_west = fh.variables['u_west'][:]
v_west = fh.variables['v_west'][:]
zeta_west = fh.variables['zeta_west'][:]
ocean_time = fh.variables['ocean_time'][:]
lat_west = fh.variables['lat_psi'][:]
lon_west = fh.variables['lon_psi'][:]
h = fh.variables['h'][:]
s_rho = fh.variables['s_rho'][:]
fh.close()

lat_west = lat_west[:, 0]
lon_west = lon_west[:, 0]
h = h[:, 0]

u_west = 0.5*(u_west[:, :, 1:]+u_west[:, :, :-1])
zeta_west = 0.5*(zeta_west[:, 1:]+zeta_west[:, :-1])
h = 0.5*(h[1:]+h[:-1])

u2 = u_west[:, :, 279]
v2 = v_west[:, :, 279]
