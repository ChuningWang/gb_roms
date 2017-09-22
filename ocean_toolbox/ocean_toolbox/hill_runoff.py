import numpy as np
import pdb

def get_discharge(pth, savepth, latlim, lonlim):
    import netCDF4 as nc
    import glob
    from datetime import datetime, timedelta
    import matplotlib.dates as mdates

    # Reading lat lon
    print 'Loading Lat & lon...'
    fh = nc.Dataset(pth+'lat_lon.nc')
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]
    lon = lon-360

    # Cut out useful portion
    msk_gb = (lat>latlim[0]) & (lat<latlim[1]) & (lon>lonlim[0]) & (lon<lonlim[1])
    msk_lat = np.any(msk_gb,axis=1)
    msk_lon = np.any(msk_gb,axis=0)

    lat = lat[msk_lat, :][:, msk_lon]
    lon = lon[msk_lat, :][:, msk_lon]

    # Load data
    print 'Getting file names...'
    flist = glob.glob(pth+'discharge'+'*.nc')

    rddata = 1
    if rddata == 1:
        # Initiate
        print 'Loading '+flist[0]+'...'
        fh = nc.Dataset(flist[0],mode='r')
        d = fh.variables['discharge'][:, :, msk_lat, :][:, :, :, msk_lon]
        t = fh.variables['time'][:]
        t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
        pyt = np.array([t_ini+timedelta(hours=t[i]) for i in range(t.size)])

        for ff in flist[1:]:
            print 'Loading '+ff+'...'
            fh = nc.Dataset(ff,mode='r')
            d_in = fh.variables['discharge'][:]
            d_in = d_in[:, :, msk_lat, :][:, :, :, msk_lon]
            t = fh.variables['time'][:]
            t_ini = datetime.strptime(fh.variables['time'].units[12:], '%Y-%m-%d %H:%M:%S')
            pyt_in = np.array([t_ini+timedelta(hours=t[i]) for i in range(t.size)])
            
            d = np.concatenate([d, d_in], axis=0)
            pyt = np.concatenate([pyt, pyt_in])

    # mask out invalid data
    print 'Setting invalid data to NaN...'
    d = np.squeeze(d)
    d[d<-1000] = np.nan

    lon_dim = d.shape[1]
    lat_dim = d.shape[2]

    mt = mdates.date2num(pyt)
    
    # Load coast cell file
    print 'Getting coast cells...'
    fh = nc.Dataset(pth+'lat_lon.nc')
    coast_cells = np.squeeze(fh.variables['coast_cells'][:, msk_lat, :][:, :, msk_lon])

    # d1 = np.squeeze(d[0,:,:,:])
    # plt.close()
    # plt.contourf(lon,lat,d1)
    # plt.show(block=False)

    svdata = 1
    if svdata == 1:
        # write data into netCDF file
        print 'Saving data as netCDF4 file...'
        f = nc.Dataset(savepth+'discharge_gb.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay river discharge and deglaciation'

        f.createDimension('time', None)
        f.createDimension('lat', lat_dim)
        f.createDimension('lon', lon_dim)

        t_nc = f.createVariable('t', 'f8', ('time'))
        lat_nc = f.createVariable('lat', 'f8', ('lon', 'lat'))
        lon_nc = f.createVariable('lon', 'f8', ('lon', 'lat'))
        coast_nc = f.createVariable('coast_cells', 'f8', ('lon', 'lat'))
        d_nc = f.createVariable('discharge', 'f8', ('time', 'lon', 'lat'))

        t_nc[:] = mt
        lat_nc[:, :] = lat
        lon_nc[:, :] = lon
        coast_nc[:, :] = coast_cells
        d_nc[:, :, :] = d

        t_nc.units = 'days since 0001-01-01'
        d_nc.units = 'm^3day^-1'

        f.close()

def get_avgbox(boxMethod=1):

    if boxMethod==1:
        print 'Generating Box Type 1...'
        box0 = np.array([[  -137.7,     59.1  ],
                         [  -137.7,     59.25 ],
                         [  -136.6,     59.25 ],
                         [  -136.15,    58.725],
                         [  -136.65,    58.725],
                         [  -137.15,    58.65 ]])

        box1 = np.array([[  -136.6,     59.25 ],
                         [  -136.3,     59.25 ],
                         [  -135.7,     59.0  ],
                         [  -135.7,     58.725],
                         [  -136.,      58.825],
                         [  -136.15,    58.725],
                         [  -136.3,     58.9  ]])

        box2 = np.array([[  -136.15,    58.725],
                         [  -136.,      58.825],
                         [  -135.7,     58.725],
                         [  -135.65,    58.55 ],
                         [  -136.65,    58.55 ],
                         [  -136.65,    58.725]])

        box3 = np.array([[  -136.65,    58.55 ],
                         [  -135.65,    58.55 ],
                         [  -135.65,    58.45 ],
                         [  -136.0,     58.375]])

        box = {'box0': box0,
               'box1': box1,
               'box2': box2,
               'box3': box3
              }


    elif boxMethod==2:
        print 'Generating Box Type 2...'
        box0 = np.array([[  -137.7,     59.1  ],
                         [  -137.7,     59.25 ],
                         [  -136.6,     59.25 ],
                         [  -136.15,    58.725],
                         [  -136.65,    58.725],
                         [  -137.15,    58.65 ]])

        box1 = np.array([[  -136.6,     59.25 ],
                         [  -136.3,     59.25 ],
                         [  -135.7,     59.0  ],
                         [  -135.7,     58.725],
                         [  -136.,      58.825],
                         [  -136.15,    58.725],
                         [  -136.3,     58.9  ]])

        box2 = np.array([[  -136.15,    58.725],
                         [  -136.,      58.825],
                         [  -135.7,     58.725],
                         [  -135.6,     58.425],
                         [  -135.975,   58.375],
                         [  -136.0,     58.575]])

        box3 = np.array([[  -136.65,    58.725],
                         [  -136.15,    58.725],
                         [  -136.0,     58.575],
                         [  -135.975,   58.375],
                         [  -136.65,    58.575]])

        box = {'box0': box0,
               'box1': box1,
               'box2': box2,
               'box3': box3
              }

    return box

def get_discharge_avgbox(pth, savepth=-1, boxMethod=1):
    import netCDF4 as nc
    from matplotlib import path

    print 'Loading Boxes...'
    boxes = get_avgbox(boxMethod=boxMethod)
    box0 = boxes['box0']
    box1 = boxes['box1']
    box2 = boxes['box2']
    box3 = boxes['box3']
    # ------------------------------------------------------------------------------------------------------
    # Load discharge data
    print 'Loading freshwater discharge...'
    fh = nc.Dataset(pth,'r')

    t = fh.variables['t'][:]
    lat = fh.variables['lat'][:]
    lon = fh.variables['lon'][:]

    # Get points in boxes
    hydro_box = np.ones(lon.shape)*(-1)
    p0 = path.Path(box0)
    p1 = path.Path(box1)
    p2 = path.Path(box2)
    p3 = path.Path(box3)

    for i in range(lon.shape[0]):
        for j in range(lon.shape[1]):
            if p0.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 0
            elif p1.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 1
            elif p2.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 2
            elif p3.contains_points([(lon[i, j], lat[i, j])]):
                hydro_box[i, j] = 3

    # Find coastal cells
    coast = np.squeeze(fh.variables['coast_cells'][:])
    
    # gridcoast = np.zeros(d0.shape)
    # 
    # for i in range(1, d0.shape[0]-1):
    #     for j in range(1, d0.shape[1]-1):
    #         # if np.any(np.isnan([d0[i-1, j], d0[i+1, j], d0[i, j-1], d0[i, j+1]])):
    #         if np.any(np.isnan(d0[i-1:i+2, j-1:j+2])):
    #             gridcoast[i, j] = 1
    # 
    # hydro_box[gridcoast==0] = -1

    hydro_box[coast.mask] = -1

    import matplotlib.pyplot as plt
    plt.contourf(lon, lat, hydro_box)
    plt.savefig('/Users/chuning/projects/glacierbay/figs/hydro_box.eps', format='eps')

    # pdb.set_trace()
    # ------------------------------------------------------------------------------------------------------
    # Divide GB into several hydro regions
    print 'Averaging data within each box...'
    avg_runoff = 1
    if avg_runoff == 1:
        d = np.empty((t.size, 4))
        d[:] = np.NaN
        for i in range(t.size):
            d0 = np.squeeze(fh.variables['discharge'][i, :, :])
            d0[d0<=0] = np.NaN
            # d0[gridcoast==0] = np.NaN
            d[i, 0] = np.nansum(d0[hydro_box==0])
            d[i, 1] = np.nansum(d0[hydro_box==1])
            d[i, 2] = np.nansum(d0[hydro_box==2])
            d[i, 3] = np.nansum(d0[hydro_box==3])

        # plt.close()
        # plt.figure()
        # plt.plot(pyt, d)
        # plt.legend(('1','2','3','4'))
        # plt.savefig('../figs/runoff.eps', format='eps')

        discharge = {'mt':          t,
                     'box0':        box0,
                     'box1':        box1,
                     'box2':        box2,
                     'box3':        box3,
                     'discharge':   d
                    }

    if savepth!=-1:
        print 'Saving data to '+savepth+'discharge_gb_box.nc...'
        f = nc.Dataset(savepth+'discharge_gb_box.nc', 'w', format='NETCDF4')
        f.description = 'Glacier Bay freshwater discharge and deglaciation, sum of each box'

        f.createDimension('time', None)
        f.createDimension('box', None)

        t_nc = f.createVariable('t', 'f8', ('time'))
        d_nc = f.createVariable('discharge', 'f8', ('time', 'box'))

        t_nc[:] = t
        d_nc[:, :] = d
        t_nc.units = 'days since 0001-01-01'
        d_nc.units = 'm^3day^-1'

        f.close()

    return discharge


#     # ------------------------------------------------------------------------------------------------------
#     # Load discharge data
#     fh = nc.Dataset(pth,'r')
# 
#     t = fh.variables['t'][:]
#     lat = fh.variables['lat'][:]
#     lon = fh.variables['lon'][:]
# 
#     # Get points in boxes
#     hydro_box = np.zeros(lon.shape)
#     p1 = path.Path(box1)
#     p2 = path.Path(box2)
#     p3 = path.Path(box3)
#     p4 = path.Path(box4)
# 
#     for i in range(lon.shape[0]):
#         for j in range(lon.shape[1]):
#             if p1.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 1
#             elif p2.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 2
#             elif p3.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 3
#             elif p4.contains_points([(lon[i, j], lat[i, j])]):
#                 hydro_box[i, j] = 4


