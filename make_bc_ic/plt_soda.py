import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# load topography
bathydir = '/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/ARDEMv2.0.nc'
fh = nc.Dataset(bathydir, mode='r')
topo = fh.variables['z'][:]
lons = fh.variables['lon'][:] 
lons[lons>180] = lons[lons>180]-360  # lons from -180 to 180
lats = fh.variables['lat'][:] 
fh.close()

xr = (lons>=-140) & (lons<=-132)
yr = (lats>=52) & (lats<=60)
topo = topo[yr, :][:, xr]
lons = lons[xr]
lats = lats[yr]

soda_dir = '/Volumes/P1/Data/SODA/SODA_3.3.1/monthly/'
tag = '2000_01'
soda_filein = soda_dir+'soda3.3.1_monthly_ocean_reg_'+tag+'.nc'

x0 = 446
y0 = 266 
xt = np.arange(440, 455)
yt = np.arange(255, 270)

fh = nc.Dataset(soda_filein, 'r')
lont = fh.variables['xt_ocean'][xt]
lont = lont-360
latt = fh.variables['yt_ocean'][yt]
zt = fh.variables['st_ocean'][:]
t = fh.variables['temp'][:, :, yt, xt].squeeze()
s = fh.variables['salt'][:, :, yt, xt].squeeze()
h = fh.variables['ssh'][:, yt, xt].squeeze()

x0 = 446
y0 = 266 
xu = np.arange(440, 455)
yu = np.arange(255, 270)

lonu = fh.variables['xu_ocean'][xu]
lonu = lonu-360
latu = fh.variables['yu_ocean'][yu]
zu = fh.variables['sw_ocean'][:]
u = fh.variables['u'][:, :, yu, xu].squeeze()
v = fh.variables['v'][:, :, yu, xu].squeeze()

# make plots
for i in range(25):
    plt.figure()
    plt.pcolormesh(lont, latt, t[i, :, :])
    plt.clim([4, 7])
    plt.colorbar()
    plt.quiver(lonu, latu, u[i, :, :], v[i, :, :,])
    plt.contour(lons, lats, topo, [0, 0], colors='k')
    plt.xlim([-139, -133])
    plt.ylim([54, 60])
    plt.title('Depth= %02d' % zt[i])
    plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/sync/soda/t_2000Jan%02d.tiff' % i, format='tiff')
    plt.close()

    plt.figure()
    plt.pcolormesh(lont, latt, s[i, :, :])
    plt.clim([28, 34])
    plt.colorbar()
    plt.quiver(lonu, latu, u[i, :, :], v[i, :, :,])
    plt.contour(lons, lats, topo, [0, 0], colors='k')
    plt.xlim([-139, -133])
    plt.ylim([54, 60])
    plt.title('Depth= %02d' % zt[i])
    plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/sync/soda/s_2000Jan%02d.tiff' % i, format='tiff')
    plt.close()
