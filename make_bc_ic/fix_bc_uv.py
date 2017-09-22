import pyroms
import netCDF4 as nc
import numpy as np
import sys
from gb_toolbox import gb_current
from scipy.interpolate import interp1d
from scipy.signal import filtfilt

import read_host_info
sv = read_host_info.read_host_info()
dst_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

my_year = '2008'
tag = 'SODA3.3.1_0.25'

# load bc data
grd = pyroms.grid.get_ROMS_grid(grd1)
fh = nc.Dataset(dst_dir+'bc_ic/'+grd.name+'_bdry_'+my_year+'_'+tag+'.nc')
h = fh.variables['h'][:]
lat_u = fh.variables['lat_u'][:]
lon_u = fh.variables['lon_u'][:]
lat_v = fh.variables['lat_v'][:]
lon_v = fh.variables['lon_v'][:]
msk_u = fh.variables['mask_u'][:]
msk_v = fh.variables['mask_v'][:]
u_w = fh.variables['u_west'][:]
v_w = fh.variables['v_west'][:]
Cs_r = fh.variables['Cs_r'][:]
Cs_w = fh.variables['Cs_w'][:]
ang = fh.variables['angle'][:]
fh.close()

h_u = 0.5*(h[:,1:]+h[:,:-1])
h_v = 0.5*(h[:,1:]+h[:,:-1])

# only extract the west boundary
lat_u = lat_u[:, 0]
lat_v = lat_v[:, 0]
lon_u = lon_u[:, 0]
lon_v = lon_v[:, 0]
msk_u = msk_u[:, 0]
msk_v = msk_v[:, 0]
h_u = h_u[:, 0]
h_v = h_v[:, 0]

ang = ang.mean()

# load ADCP data near the grid boundary
stn_list = ['SEA1008', 'SEA1009', 'SEA1010']
bdate_list = ['20100525', '20100625', '20100725']
edate_list = ['20100625', '20100725', '20100815']

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

    # sort the ADCP data. get rid of overlap periods.
    t[stn], t_idx = np.unique(t[stn], return_index=True)
    u[stn] = u[stn][:, t_idx]
    v[stn] = v[stn][:, t_idx]

# process, interpolate and smooth data
lat08 = lat['SEA1008']
lon08 = lon['SEA1008']
z08 = z['SEA1008']
u08 = u['SEA1008'].mean(axis=1)
v08 = v['SEA1008'].mean(axis=1)

lat09 = lat['SEA1009']
lon09 = lon['SEA1009']
z09 = z['SEA1009']
u09 = u['SEA1009'].mean(axis=1)
v09 = v['SEA1009'].mean(axis=1)

lat10 = lat['SEA1010']
lon10 = lon['SEA1010']
z10 = z['SEA1010']
u10 = u['SEA1010'].mean(axis=1)
v10 = v['SEA1010'].mean(axis=1)

zz = np.arange(325)
uu08 = np.nan*np.zeros(zz.shape)
vv08 = np.nan*np.zeros(zz.shape)
uu09 = np.nan*np.zeros(zz.shape)
vv09 = np.nan*np.zeros(zz.shape)
uu10 = np.nan*np.zeros(zz.shape)
vv10 = np.nan*np.zeros(zz.shape)

msk08 = (zz<z08[0]) & (zz>z08[-1])
uu08[msk08] = interp1d(z08, u08)(zz[msk08])
vv08[msk08] = interp1d(z08, v08)(zz[msk08])
msk09 = (zz<z09[0]) & (zz>z09[-1])
uu09[msk09] = interp1d(z09, u09)(zz[msk09])
vv09[msk09] = interp1d(z09, v09)(zz[msk09])
msk10 = (zz<z10[0]) & (zz>z10[-1])
uu10[msk10] = interp1d(z10, u10)(zz[msk10])
vv10[msk10] = interp1d(z10, v10)(zz[msk10])

uu = np.nanmean(np.array([uu08, uu09, uu10]), axis=0)
vv = np.nanmean(np.array([vv08, vv09, vv10]), axis=0)

uu[0] = -0.15
uu[-1] = 0.
vv[0] = -0.25
vv[-1] = 0.

msk = ~np.isnan(uu)
uu = interp1d(zz[msk], uu[msk])(zz)
vv = interp1d(zz[msk], vv[msk])(zz)

uu = filtfilt(np.ones(10)/10, 1, uu)
vv = filtfilt(np.ones(10)/10, 1, vv)

# rotate the velocity vector
UU = uu+1j*vv
UU = UU*np.exp(ang*1j)
ux = UU.real
uy = UU.imag

# calculate distance from station to boundary. find the nearest point index
dis = np.sqrt((lon_u-lon09)**2+(lat_u-lat09)**2)
idx = np.argmin(dis)
hi = h_u[idx]
zi = -hi*Cs_r
ux = interp1d(zz, ux)(zi)

dis = np.sqrt((lon_v-lon09)**2+(lat_v-lat09)**2)
idx = np.argmin(dis)
hi = h_v[idx]
zi = -hi*Cs_r
uy = interp1d(zz, uy)(zi)

u_w_new = np.zeros(u_w.shape)
v_w_new = np.zeros(v_w.shape)
for i in range(u_w_new.shape[-1]):
    u_w_new[:, :, i] = np.tile(ux, (u_w.shape[0], 1))
for i in range(v_w_new.shape[-1]):
    v_w_new[:, :, i] = np.tile(uy, (v_w.shape[0], 1))

u_w_new = np.ma.masked_where(u_w.mask, u_w_new)
v_w_new = np.ma.masked_where(v_w.mask, v_w_new)

# update u_w, v_w to bdry file
fh = nc.Dataset(dst_dir+'bc_ic/'+grd.name+'_bdry_'+my_year+'_'+tag+'.nc', 'a')
fh.variables['u_west'][:] = u_w_new
fh.variables['v_west'][:] = v_w_new
fh.close()

