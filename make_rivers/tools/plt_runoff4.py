import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc

fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/runoff_GB_hill.nc', 'r')
lon = fh.variables['lon_rho'][:]
lat = fh.variables['lat_rho'][:]
msk = fh.variables['Runoff'][:].mask
fh.close()

fh = nc.Dataset('../Glacier_Bay_rivers.nc', 'r')
xi = fh.variables['river_Xposition'][:]
yi = fh.variables['river_Eposition'][:]
t = fh.variables['river_time'][:]
t = t-t[0]  # yearday
runoff = fh.variables['river_transport'][:]
fh.close()

r = np.ma.zeros((len(t), lon.shape[0], lon.shape[1]))
r.mask = True
for i in range(len(t)):
    for j in range(len(xi)):
        r.mask[i, yi[j], xi[j]] = False
        r[i, yi[j], xi[j]] = runoff[i, j]

fh = nc.Dataset('gb_runoff_grid.nc', 'w')
fh.description = 'Glacier Bay river discharge and deglaciation'

fh.createDimension('time', None)
fh.createDimension('lat', lon.shape[1])
fh.createDimension('lon', lon.shape[0])

t_nc = fh.createVariable('t', 'f8', ('time'))
t_nc.long_name = 'yearday'
t_nc.units = 'days'
lat_nc = fh.createVariable('lat', 'f8', ('lon', 'lat'))
lat_nc.long_name = 'latitude'
lon_nc = fh.createVariable('lon', 'f8', ('lon', 'lat'))
lon_nc.long_name = 'longitude'
d_nc = fh.createVariable('runoff', 'f8', ('time', 'lon', 'lat'))
d_nc.long_name = 'runoff'

t_nc[:] = t
lat_nc[:, :] = lat
lon_nc[:, :] = lon
d_nc[:, :, :] = r

fh.close()

# r = np.ma.masked_where(msk, r)

for i in range(len(t)):
    plt.figure()
    # plt.scatter(lon, lat, c=np.log10(runoff[i, :]), edgecolors='none', alpha=0.5)
    plt.pcolor(lon, lat, np.log10(r[i, :, :]))
    # plt.xlim(0, 600)
    # plt.ylim(0, 1000)
    plt.clim(-4, 4)
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'Discharge (log scale) [m$^3$s$^{-1}$]')
    ttl = "%03d" % t[i]
    plt.title('Day '+ttl)
    plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/runoff4/'+ttl+'.tiff', format='tiff')
    plt.close()
