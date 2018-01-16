'''
This script is used to generate ROMS boundary file from WOD CTD stations and NOAA ADCP measurements.
This is a messy script, since it involves a lot of subjective choices of profiles, values.
'''

import sys
import glob
import numpy as np
from datetime import datetime
from scipy.signal import filtfilt
from scipy.interpolate import interp1d
import netCDF4 as nc
import pyroms
import pyroms_toolbox
from matplotlib import path
from wodpy import wod
from ocean_toolbox import noaa_adcp

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
dst_dir = sv['out_dir']
soda_dir = sv['soda_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

my_year = 2009
zlev_ts = 500
grd = pyroms.grid.get_ROMS_grid(grd1)

class nctime(object):
    pass

nctime.long_name = 'time'
nctime.units = 'days since 1900-01-01 00:00:00'
bc_file = dst_dir + 'bc_ic/' + grd.name + '_bc_' + str(my_year) + '_CTD_ADCP' + '.nc'
pyroms_toolbox.nc_create_roms_bdry_file(bc_file, grd, nctime)

cal_ts = 1
if cal_ts == 1:
    # first, process salinity and temperature
    # load ctd casts
    c_east = [(-134.50, 57.95),
              (-135.50, 57.95),
              (-135.50, 58.50),
              (-134.50, 58.50)]
    p_east = path.Path(c_east)

    c_west = [(-136.30, 57.50),
              (-137.50, 57.50),
              (-137.50, 58.50),
              (-136.30, 58.50)]
    p_west = path.Path(c_west)

    varlist = ['time', 'yearday', 'lat', 'lon', 'salt', 'temp']
    for var in varlist:
        exec(var + '_east = []')
        exec(var + '_west = []')

    z_ts = np.arange(zlev_ts)

    fid = open(in_dir + 'CTDS7513')
    n = 0  # counter
    t = True
    while t:
        pr = wod.WodProfile(fid)
        lon0, lat0 = pr.longitude(), pr.latitude()
        if p_east.contains_point((lon0, lat0)):
            zpr = pr.z()
            if len(zpr)>1:
                zmin, zmax = zpr[0], zpr[-1]
                zmsk = (z_ts>=zmin) & (z_ts<=zmax)
                spr = np.zeros(zlev_ts)*np.NaN
                spr[zmsk] = np.interp(z_ts[zmsk], pr.z(), pr.s())
                tpr = np.zeros(zlev_ts)*np.NaN
                tpr[zmsk] = np.interp(z_ts[zmsk], pr.z(), pr.t())
                salt_east.append(spr)
                temp_east.append(tpr)
                time_east.append(nc.date2num(pr.datetime(), 'days since 1900-01-01'))
                yearday_east.append(pr.datetime().timetuple().tm_yday)
                lat_east.append(lat0)
                lon_east.append(lon0)
        elif p_west.contains_point((lon0, lat0)):
            zpr = pr.z()
            if len(zpr)>1:
                zmin, zmax = zpr[0], zpr[-1]
                zmsk = (z_ts>=zmin) & (z_ts<=zmax)
                spr = np.zeros(zlev_ts)*np.NaN
                spr[zmsk] = np.interp(z_ts[zmsk], pr.z(), pr.s())
                tpr = np.zeros(zlev_ts)*np.NaN
                tpr[zmsk] = np.interp(z_ts[zmsk], pr.z(), pr.t())
                salt_west.append(spr)
                temp_west.append(tpr)
                time_west.append(nc.date2num(pr.datetime(), 'days since 1900-01-01'))
                yearday_west.append(pr.datetime().timetuple().tm_yday)
                lat_west.append(lat0)
                lon_west.append(lon0)
        t = not pr.is_last_profile_in_file(fid)
        n += 1

    for var in varlist:
        exec(var + '_east = ' + 'np.array(' + var + '_east)')
        exec(var + '_west = ' + 'np.array(' + var + '_west)')

    # get average profile
    for var in ['salt', 'temp']:
        exec(var + '_avg_east = np.nanmean(' + var + '_east, axis=0)')
        exec(var + '_avg_west = np.nanmean(' + var + '_west, axis=0)')

    temp_avg_east[z_ts>200] = np.NaN

    # fill up NaNs
    temp_avg_east[-1] = 5.0
    temp_avg_west[-1] = 4.3
    salt_avg_east[-1] = 34.0
    salt_avg_west[-1] = 34.0

    msk = ~np.isnan(temp_avg_east)
    temp_avg_east = np.interp(z_ts, z_ts[msk], temp_avg_east[msk])
    salt_avg_east = np.interp(z_ts, z_ts[msk], salt_avg_east[msk])
    msk = ~np.isnan(temp_avg_west)
    temp_avg_west = np.interp(z_ts, z_ts[msk], temp_avg_west[msk])
    salt_avg_west = np.interp(z_ts, z_ts[msk], salt_avg_west[msk])

    # filter
    b = np.ones(20)/20.
    a = 1
    temp_avg_east = filtfilt(b, a, temp_avg_east)
    salt_avg_east = filtfilt(b, a, salt_avg_east)
    temp_avg_west = filtfilt(b, a, temp_avg_west)
    salt_avg_west = filtfilt(b, a, salt_avg_west)

    # # deep water filter
    # b = np.ones(20)/20.
    # a = 1
    # temp_avg_east[100:] = filtfilt(b, a, temp_avg_east[100:])
    # salt_avg_east[100:] = filtfilt(b, a, salt_avg_east[100:])
    # temp_avg_west[100:] = filtfilt(b, a, temp_avg_west[100:])
    # salt_avg_west[100:] = filtfilt(b, a, salt_avg_west[100:])

    # give salinity some seasonal variability
    t_clim = np.array([0., 31., 59., 90., 120., 151., 181., 212., 243., 273., 304., 334.])+14
    salt_clim_s_east = np.array([31, 31, 30.5, 30,   29,   29,   28, 28,   27, 26.5, 25, 26])
    salt_clim_s_west = np.array([31, 31, 30.5, 30.5, 30.5, 30.5, 30, 29.5, 29, 30,   31, 31])
    temp_clim_s_east = np.array([8, 8, 8, 9, 10, 10, 10, 11, 11, 11, 10, 9])
    temp_clim_s_west = np.array([8, 8, 8, 9, 10, 10, 10, 11, 11, 11, 10, 9])

    salt_clim_east = np.zeros((zlev_ts, 12))
    salt_clim_west = np.zeros((zlev_ts, 12))
    temp_clim_east = np.zeros((zlev_ts, 12))
    temp_clim_west = np.zeros((zlev_ts, 12))

    for i in range(zlev_ts):
        salt_clim_east[i, :] = (salt_avg_east[i]-34)/(salt_avg_east[0]-34)*(salt_clim_s_east-34)+34
        salt_clim_west[i, :] = (salt_avg_west[i]-34)/(salt_avg_west[0]-34)*(salt_clim_s_west-34)+34
        temp_clim_east[i, :] = (temp_avg_east[i]-5)/(temp_avg_east[0]-5)*(temp_clim_s_east-5)+5
        temp_clim_west[i, :] = (temp_avg_west[i]-4.3)/(temp_avg_west[0]-4.3)*(temp_clim_s_west-4.3)+4.3

# --------------------------------------------------------------------------------
# second process velocity
# load ADCP data near the grid boundary
stn_list = ['SEA1008', 'SEA1009', 'SEA1010', 'SEA0839']
zlev_uv = 450

t = {}
u = {}
v = {}
z = {}
lat = {}
lon = {}

for stn in stn_list:

    # info = {'stn' : stn,
    #         'file_dir': in_dir + 'NOAA_ADCP/',
    #         'sl': 'l',
    #         'Wp_hrs': -1}

    # crt = noaa_adcp.get_noaa_current(info)
    # crt()

    # t[stn] = crt.ctime
    # u[stn] = crt.u
    # v[stn] = crt.v
    # z[stn] = crt.z
    # lat[stn] = crt.info['lat']
    # lon[stn] = crt.info['lon']

    fname = glob.glob(in_dir + 'NOAA_ADCP/' + stn + '*.nc')[0]

    fin = nc.Dataset(fname, 'r')
    t[stn] = fin.variables['time'][:]
    u[stn] = fin.variables['u'][:]
    v[stn] = fin.variables['v'][:]
    z[stn] = fin.variables['z'][:]
    lat[stn] = fin.lat
    lon[stn] = fin.lon

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

z_uv = np.arange(zlev_uv)
uu08 = np.nan*np.zeros(z_uv.shape)
vv08 = np.nan*np.zeros(z_uv.shape)
uu09 = np.nan*np.zeros(z_uv.shape)
vv09 = np.nan*np.zeros(z_uv.shape)
uu10 = np.nan*np.zeros(z_uv.shape)
vv10 = np.nan*np.zeros(z_uv.shape)

msk08 = (z_uv<z08[0]) & (z_uv>z08[-1])
uu08[msk08] = interp1d(z08, u08)(z_uv[msk08])
vv08[msk08] = interp1d(z08, v08)(z_uv[msk08])
msk09 = (z_uv<z09[0]) & (z_uv>z09[-1])
uu09[msk09] = interp1d(z09, u09)(z_uv[msk09])
vv09[msk09] = interp1d(z09, v09)(z_uv[msk09])
msk10 = (z_uv<z10[0]) & (z_uv>z10[-1])
uu10[msk10] = interp1d(z10, u10)(z_uv[msk10])
vv10[msk10] = interp1d(z10, v10)(z_uv[msk10])

# uu = np.nanmean(np.array([uu08, uu09, uu10]), axis=0)
# vv = np.nanmean(np.array([vv08, vv09, vv10]), axis=0)
uu = np.nanmean(np.array([uu08, uu10]), axis=0)
vv = np.nanmean(np.array([vv08, vv10]), axis=0)

uu[0] = -0.15
uu[350] = 0.
uu[-1] = 0.
vv[0] = -0.25
vv[350] = 0.
vv[-1] = 0.

msk = ~np.isnan(uu)
uu = interp1d(z_uv[msk], uu[msk])(z_uv)
vv = interp1d(z_uv[msk], vv[msk])(z_uv)

uu = filtfilt(np.ones(10)/10, 1, uu)
vv = filtfilt(np.ones(10)/10, 1, vv)

UU_west = uu+1j*vv

lat39 = lat['SEA0839']
lon39 = lon['SEA0839']
z39 = z['SEA0839']
u39 = u['SEA0839'].mean(axis=1)
v39 = v['SEA0839'].mean(axis=1)

uu39 = np.nan*np.zeros(z_uv.shape)
vv39 = np.nan*np.zeros(z_uv.shape)

msk39 = (z_uv<z39[0]) & (z_uv>z39[-1])
uu39[msk39] = interp1d(z39, u39)(z_uv[msk39])
vv39[msk39] = interp1d(z39, v39)(z_uv[msk39])

uu39[0] = 0.02
uu39[-1] =0
vv39[0] = -0.025
vv39[-1] =0

msk = ~np.isnan(uu39)
uu39 = interp1d(z_uv[msk], uu39[msk])(z_uv)
vv39 = interp1d(z_uv[msk], vv39[msk])(z_uv)

uu39 = filtfilt(np.ones(10)/10, 1, uu39)
vv39 = filtfilt(np.ones(10)/10, 1, vv39)

UU_east = uu39+1j*vv39

# --------------------------------------------------------------------------------
# remap TS and write into nc file
Cs_r = grd.vgrid.Cs_r
zlev = len(Cs_r)
h_ts = grd.vgrid.h
h_u = 0.5*(h_ts[:, 1:] + h_ts[:, :-1])
h_v = 0.5*(h_ts[1:, :] + h_ts[:-1, :])

dt = (datetime(my_year, 01, 01) - datetime(1900, 01, 01)).days

fh = nc.Dataset(bc_file, 'r+')
fh.variables['ocean_time'][:] = t_clim + dt
spval = -1.0e20

# --------------------------------------------------------------------------------
varlist = ['temp', 'salt', 'u', 'v']
for var in varlist:
    var_name_east = var + '_east'
    var_name_west = var + '_west'
    # write var info
    if var in['temp']:
        h = h_ts.copy()
        msk = grd.hgrid.mask_rho.copy()
        d3 = 'eta_rho'
        long_name_east = 'potential temperature east boundary condition'
        field_east = 'temp_east, scalar, series'
        long_name_west = 'potential temperature west boundary condition'
        field_west = 'temp_west, scalar, series'
        units = 'Celsius'
    elif var in['salt']:
        h = h_ts.copy()
        msk = grd.hgrid.mask_rho.copy()
        d3 = 'eta_rho'
        long_name_east = 'salinity east boundary condition'
        field_east = 'salt_east, scalar, series'
        long_name_west = 'salinity west boundary condition'
        field_west = 'salt_west, scalar, series'
        units = 'PSU'
    elif var in ['u']:
        h = h_u.copy()
        msk = grd.hgrid.mask_u.copy()
        d3 = 'eta_u'
        long_name_east = '3D u-momentum east boundary condition'
        field_east = 'u_east, scalar, series'
        long_name_west = '3D u-momentum west boundary condition'
        field_west = 'u_west, scalar, series'
        units = 'meter second-1'

        long_name_east2 = '2D u-momentum east boundary condition'
        field_east2 = 'ubar_east, scalar, series'
        long_name_west2 = '2D u-momentum west boundary condition'
        field_west2 = 'ubar_west, scalar, series'
        units2 = 'meter second-1'
    elif var in ['v']:
        h = h_v.copy()
        msk = grd.hgrid.mask_v.copy()
        d3 = 'eta_v'
        long_name_east = '3D v-momentum east boundary condition'
        field_east = 'v_east, scalar, series'
        long_name_west = '3D v-momentum west boundary condition'
        field_west = 'v_west, scalar, series'
        units = 'meter second-1'

        long_name_east2 = '2D v-momentum east boundary condition'
        field_east2 = 'vbar_east, scalar, series'
        long_name_west2 = '2D v-momentum west boundary condition'
        field_west2 = 'vbar_west, scalar, series'
        units2 = 'meter second-1'

    # some grid info
    eta, xi = h.shape
    h_east = h[:, -1]
    h_west = h[:, 0]
    msk_east = msk[:, -1]
    msk_west = msk[:, 0]
    h_east[msk_east==0] = np.NaN
    h_west[msk_west==0] = np.NaN
    z_east = -grd.vgrid.z_r[:][:, :, -1]
    z_west = -grd.vgrid.z_r[:][:, :, 0]
    ang_east = grd.hgrid.angle_rho[: ,-1]
    ang_west = grd.hgrid.angle_rho[: ,0]

    if var in ['v']:
        z_east = 0.5*(z_east[:, 1:] + z_east[:, :-1])
        z_west = 0.5*(z_west[:, 1:] + z_west[:, :-1])
        ang_east = 0.5*(ang_east[1:] + ang_east[:-1])
        ang_west = 0.5*(ang_west[1:] + ang_west[:-1])

    data_east = np.NaN*np.zeros((len(t_clim), zlev, eta))
    data_west = np.NaN*np.zeros((len(t_clim), zlev, eta))

    for i in range(z_east.shape[1]):
        for tt in range(len(t_clim)):
            if var == 'salt':
                data_east[tt, :, i] = interp1d(z_ts, salt_clim_east[:, tt])(z_east[:, i])
                data_west[tt, :, i] = interp1d(z_ts, salt_clim_west[:, tt])(z_west[:, i])
            elif var == 'temp':
                data_east[tt, :, i] = interp1d(z_ts, temp_clim_east[:, tt])(z_east[:, i])
                data_west[tt, :, i] = interp1d(z_ts, temp_clim_west[:, tt])(z_west[:, i])
            elif var == 'u':
                u_east = (UU_east * np.exp(-ang_east[i]*1j)).real
                data_east[tt, :, i] = interp1d(z_uv, u_east)(z_east[:, i])
                u_west = (UU_west * np.exp(-ang_west[i]*1j)).real
                data_west[tt, :, i] = interp1d(z_uv, u_west)(z_west[:, i])
            elif var == 'v':
                v_east = (UU_east * np.exp(-ang_east[i]*1j)).imag
                data_east[tt, :, i] = interp1d(z_uv, v_east)(z_east[:, i])
                v_west = (UU_west * np.exp(-ang_west[i]*1j)).imag
                data_west[tt, :, i] = interp1d(z_uv, v_west)(z_west[:, i])

    data_east[:, :, msk_east == 0] = np.NaN
    data_west[:, :, msk_west == 0] = np.NaN
    data_east = np.ma.masked_invalid(data_east)
    data_west = np.ma.masked_invalid(data_west)

    fh.createVariable(var_name_east, 'f8', ('ocean_time', 's_rho', d3), fill_value=spval)
    fh.variables[var_name_east].long_name = long_name_east
    fh.variables[var_name_east].units = units
    fh.variables[var_name_east].field = field_east
    fh.variables[var_name_east][:] = data_east

    fh.createVariable(var_name_west, 'f8', ('ocean_time', 's_rho', d3), fill_value=spval)
    fh.variables[var_name_west].long_name = long_name_west
    fh.variables[var_name_west].units = units
    fh.variables[var_name_west].field = field_west
    fh.variables[var_name_west][:] = data_west

    if var == 'temp':
        fh.createVariable('zeta_east', 'f8', ('ocean_time', d3), fill_value=spval)
        fh.variables['zeta_east'].long_name = 'free-surface east boundary condition'
        fh.variables['zeta_east'].units = 'meter'
        fh.variables['zeta_east'].field = 'zeta_east, scaler, series'
        zeta_east = h_east*0.
        zeta_east = np.ma.masked_invalid(zeta_east)
        for i in range(len(t_clim)):
            fh.variables['zeta_east'][i, :] = zeta_east

        fh.createVariable('zeta_west', 'f8', ('ocean_time', d3), fill_value=spval)
        fh.variables['zeta_west'].long_name = 'free-surface west boundary condition'
        fh.variables['zeta_west'].units = 'meter'
        fh.variables['zeta_west'].field = 'zeta_west, scaler, series'
        zeta_west = h_west*0.
        zeta_west = np.ma.masked_invalid(zeta_west)
        for i in range(len(t_clim)):
            fh.variables['zeta_west'][i, :] = zeta_west

    if var == 'u':
        ubar_east = np.NaN*np.zeros((len(t_clim), eta))
        ubar_west = np.NaN*np.zeros((len(t_clim), eta))
        for i in range(z_east.shape[1]):
            mskz_east = z_uv < h_east[i]
            mskz_west = z_uv < h_west[i]
            for tt in range(len(t_clim)):
                ubar_east[tt, i] = UU_east[mskz_east].mean().real
                ubar_west[tt, i] = UU_west[mskz_west].mean().real

        ubar_east = np.ma.masked_invalid(ubar_east)
        ubar_west = np.ma.masked_invalid(ubar_west)

        fh.createVariable('ubar_east', 'f8', ('ocean_time', d3), fill_value=spval)
        fh.variables['ubar_east'].long_name = long_name_east2
        fh.variables['ubar_east'].units = units2
        fh.variables['ubar_east'].field = field_east2
        fh.variables['ubar_east'][:] = ubar_east

        fh.createVariable('ubar_west', 'f8', ('ocean_time', d3), fill_value=spval)
        fh.variables['ubar_west'].long_name = long_name_west2
        fh.variables['ubar_west'].units = units2
        fh.variables['ubar_west'].field = field_west2
        fh.variables['ubar_west'][:] = ubar_west

    if var == 'v':
        vbar_east = np.NaN*np.zeros((len(t_clim), eta))
        vbar_west = np.NaN*np.zeros((len(t_clim), eta))
        for i in range(z_east.shape[1]):
            mskz_east = z_uv < h_east[i]
            mskz_west = z_uv < h_west[i]
            for tt in range(len(t_clim)):
                vbar_east[tt, i] = UU_east[mskz_east].mean().imag
                vbar_west[tt, i] = UU_west[mskz_west].mean().imag

        vbar_east = np.ma.masked_invalid(vbar_east)
        vbar_west = np.ma.masked_invalid(vbar_west)

        fh.createVariable('vbar_east', 'f8', ('ocean_time', d3), fill_value=spval)
        fh.variables['vbar_east'].long_name = long_name_east2
        fh.variables['vbar_east'].units = units2
        fh.variables['vbar_east'].field = field_east2
        fh.variables['vbar_east'][:] = vbar_east

        fh.createVariable('vbar_west', 'f8', ('ocean_time', d3), fill_value=spval)
        fh.variables['vbar_west'].long_name = long_name_west2
        fh.variables['vbar_west'].units = units2
        fh.variables['vbar_west'].field = field_west2
        fh.variables['vbar_west'][:] = vbar_west

fh.close()
