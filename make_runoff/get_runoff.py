# Read freshwater runoff data (hydrology)
# data are stored in /Volumes/R1/ROMS/hydrology/GOA/

import netCDF4 as nc
import glob
from datetime import datetime, timedelta
import pdb
import matplotlib.dates as mdates
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import path
from geopy.distance import vincenty

def get_discharge(pth, savepth, latlim, lonlim):
    # Reading lat lon
    print 'Loading Lat & lon...'
    fh = nc.Dataset(pth+'lat_lon.nc')
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    lon = lon-360

    # load coast cell mask
    coast = np.squeeze(fh.variables['coast_cells'][:])
    fh.close()

    # Cut out useful portion
    msk_gb = (lat>latlim[0]) & (lat<latlim[1]) & (lon>lonlim[0]) & (lon<lonlim[1])
    msk_lat = np.any(msk_gb,axis=1)
    msk_lon = np.any(msk_gb,axis=0)

    lat = lat[msk_lat, :][:, msk_lon]
    lon = lon[msk_lat, :][:, msk_lon]
    coast = coast[msk_lat, :][:, msk_lon]

    # get grid corner coordinates
    Mp, Lp = lon.shape
    lon_corner = np.zeros([Mp+1,Lp+1])
    lat_corner = np.zeros([Mp+1,Lp+1])
    lon_corner[1:Mp,1:Lp] = 0.25*(lon[:Mp-1,:Lp-1] + lon[1:Mp,:Lp-1] + \
                                  lon[:Mp-1,1:Lp] + lon[1:Mp,1:Lp])
    lat_corner[1:Mp,1:Lp] = 0.25*(lat[:Mp-1,:Lp-1] + lat[1:Mp,:Lp-1] + \
                                  lat[:Mp-1,1:Lp] + lat[1:Mp,1:Lp])
    lon_corner[0,1:Lp] = 2*lon_corner[1,1:Lp] - lon_corner[2,1:Lp]
    lon_corner[Mp,1:Lp] = 2*lon_corner[Mp-1,1:Lp] - lon_corner[Mp-2,1:Lp]
    lon_corner[:,0] = 2*lon_corner[:,1] - lon_corner[:,2]
    lon_corner[:,Lp] = 2*lon_corner[:,Lp-1] - lon_corner[:,Lp-2]
    lat_corner[0,1:Lp] = 2*lat_corner[1,1:Lp] - lat_corner[2,1:Lp]
    lat_corner[Mp,1:Lp] = 2*lat_corner[Mp-1,1:Lp] - lat_corner[Mp-2,1:Lp]
    lat_corner[:,0] = 2*lat_corner[:,1] - lat_corner[:,2]
    lat_corner[:,Lp] = 2*lat_corner[:,Lp-1] - lat_corner[:,Lp-2]
     
    dx = np.diff(lon_corner, axis=1)
    dy = np.diff(lat_corner, axis=0)
    dx = 0.5*(dx[:-1, :]+dx[1:, :])
    dy = 0.5*(dy[:, :-1]+dy[:, 1:])

    # calculate grid area in [m^2]
    area = np.zeros(lon.shape)
    for i in range(Mp):
        for j in range(Lp):
            dxi = vincenty((lat[i, j], lon[i, j]-0.5*dx[i, j]), (lat[i, j], lon[i, j]+0.5*dx[i, j])).meters
            dyi = vincenty((lat[i, j]-0.5*dy[i, j], lon[i, j]), (lat[i, j]+0.5*dy[i, j], lon[i, j])).meters
            area[i, j] = dxi*dyi
           
    # base time
    t0 = datetime(1900, 01, 01)

    # load data
    print 'Getting file names...'
    flist = glob.glob(pth+'discharge'+'*.nc')

    rddata = 1
    if rddata == 1:
        # initiate
        print 'Loading '+flist[0]+'...'
        fh = nc.Dataset(flist[0], mode='r')
        t = fh.variables['time'][:]
        t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
        t = t/24. + (t_ini-t0).days
        d = fh.variables['discharge'][:]
        fh.close()
        d = d[:, :, msk_lat, :][:, :, :, msk_lon]

        for ff in flist[1:]:
            print 'Loading '+ff+'...'
            fh = nc.Dataset(ff,mode='r')
            t_in = fh.variables['time'][:]
            t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
            t_in = t_in/24. + (t_ini-t0).days
            d_in = fh.variables['discharge'][:]
            fh.close()
            d_in = d_in[:, :, msk_lat, :][:, :, :, msk_lon]
            
            d = np.concatenate([d, d_in], axis=0)
            t = np.concatenate([t, t_in], axis=0)


    # mask out invalid data
    print 'Setting invalid data to NaN...'
    d = np.squeeze(d)
    d[d<-1000] = np.nan

    lon_dim = d.shape[1]
    lat_dim = d.shape[2]

    # unit coversion
    # first, convert discharge unit from m3day-1 to kgs-1
    d = d*1000./24./60./60.

    # second, divide by grid area to get kg/s/m^2
    for i in range(d.shape[0]):
        d[i, :, :] = d[i, :, :]/area

    svdata = 1
    if svdata == 1:
        # write data into netCDF file
        print 'Saving data to ' + savepth + '...'
        f = nc.Dataset(savepth, 'w', format='NETCDF4')
        f.description = 'Glacier Bay river discharge and deglaciation'

        f.createDimension('time', None)
        f.createDimension('lat', lat_dim)
        f.createDimension('lon', lon_dim)

        t_nc = f.createVariable('t', 'f8', ('time'))
        t_nc.long_name = 'days since 1900-01-01'
        t_nc.units = 'days'
        lat_nc = f.createVariable('lat', 'f8', ('lon', 'lat'))
        lat_nc.long_name = 'latitude'
        lon_nc = f.createVariable('lon', 'f8', ('lon', 'lat'))
        lon_nc.long_name = 'longitude'
        area_nc = f.createVariable('area', 'f8', ('lon', 'lat'))
        area_nc.long_name = 'grid cell area'
        area_nc.units = 'm^2'
        coast_nc = f.createVariable('coast', 'f8', ('lon', 'lat'))
        coast_nc.long_name = 'coast cell mask'
        d_nc = f.createVariable('discharge', 'f8', ('time', 'lon', 'lat'))
        d_nc.long_name = 'fast dicharge'
        d_nc.units = 'm3s-1'

        t_nc[:] = t
        lat_nc[:, :] = lat
        lon_nc[:, :] = lon
        area_nc[:, :] = area
        coast_nc[:, :] = coast
        d_nc[:, :, :] = d

        f.close()

pth = '/Volumes/R1/ROMS/hydrology/GOA/'
savepth = '/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/gb_discharge.nc'
latlim = [57.5, 60.5]
lonlim = [-138.5, -134.5]

# ------------------------------------------------------------------------------------------------------
get_discharge(pth, savepth, latlim, lonlim)

# ------------------------------------------------------------------------------------------------------

