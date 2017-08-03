import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from datetime import datetime
import pyroms

def get_avgbox(boxMethod=1):
    if boxMethod==1:
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

    elif boxMethod==2:
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
                         [  -135.6,     58.425],
                         [  -136.0,     58.375]])

        box = {'box0': box0,
               'box1': box1,
               'box2': box2,
               'box3': box3
              }

    return box

def get_discharge_avgbox(t, lat, lon, discharge, coast, savepth=-1, boxMethod=1):
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
    hydro_box[coast.mask] = -1

    # pdb.set_trace()
    # ------------------------------------------------------------------------------------------------------
    # Divide GB into several hydro regions
    print 'Averaging data within each box...'
    avg_runoff = 1
    if avg_runoff == 1:
        d = np.empty((t.size, 4))
        d[:] = np.NaN
        for i in range(t.size):
            d0 = np.squeeze(discharge[i, :, :])
            d0[d0<=0] = np.NaN
            d[i, 0] = np.nansum(d0[hydro_box==0])
            d[i, 1] = np.nansum(d0[hydro_box==1])
            d[i, 2] = np.nansum(d0[hydro_box==2])
            d[i, 3] = np.nansum(d0[hydro_box==3])

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

t_base = datetime(1900, 01, 01)
t_ini = datetime(1999, 12, 31)
t_end = datetime(2001, 01, 01)

fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/gb_discharge.nc', 'r')
time = fh.variables['t'][:]
t1 = int((t_ini-t_base).days-time[0])
t2 = int((t_end-t_base).days-time[0])
time = time[t1:t2]
lon = fh.variables['lon'][:]
lat = fh.variables['lat'][:]
area = fh.variables['area'][:]
r = fh.variables['discharge'][t1:t2, :, :]
coast = fh.variables['coast'][:]
fh.close()

for i in range(5):
    r[i, :, :] = r[i, :, :]*area

rbox = get_discharge_avgbox(time, lat, lon, r, coast)

runoff_file = '/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/runoff_GB_hill.nc'
fh = nc.Dataset(runoff_file)
lon2 = fh.variables['lon_rho'][:]
lat2 = fh.variables['lat_rho'][:]
time2 = fh.variables['time'][:]
r2 = fh.variables['Runoff'][:]
fh.close()

grd = pyroms.grid.get_ROMS_grid('GB')
area2 = grd.hgrid.dx*grd.hgrid.dy
coast2 = np.zeros(lon2.shape) 
# coast2 = np.ma.masked_array(coast2, mask=1-grd.hgrid.mask)
coast2 = np.ma.masked_array(coast2, mask=0)

for i in range(5):
    r2[i, :, :] = r2[i, :, :]*area2

rbox2 = get_discharge_avgbox(time2, lat2, lon2, r2, coast2)

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(rbox['mt'][:5], rbox['discharge'][:5])
axarr[1].plot(rbox2['mt'], rbox2['discharge'])
axarr[2].plot(rbox2['mt'], rbox2['discharge']/rbox['discharge'])
plt.show(block=False)

