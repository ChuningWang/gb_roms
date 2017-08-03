import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']

data_dir = out_dir + 'frc/'
fig_dir = out_dir + 'figs/frc/'
var = 'Tair'

# switch backend
plt.switch_backend('Agg') 

fh = nc.Dataset(data_dir + 'Uwind_2000_JRA55v1.1.nc')
lats = fh.variables['lat'][:]
lons = fh.variables['lon'][:]
time = fh.variables['wind_time'][:]
u10 = fh.variables['Uwind'][:]
fh.close()

fh = nc.Dataset(data_dir + 'Vwind_2000_JRA55v1.1.nc')
v10 = fh.variables['Vwind'][:]
fh.close()

fh = nc.Dataset(data_dir + var + '_2000_JRA55v1.1.nc')
v = fh.variables[var][:]
fh.close()

lon, lat = np.meshgrid(lons, lats)

m = Basemap(projection='merc', llcrnrlon=lon.min(), llcrnrlat=lat.min(),
            urcrnrlon=lon.max(), urcrnrlat=lat.max(), lat_0=0.5*(lat.min()+lat.max()), lon_0=0.5*(lon.min()+lon.max()),
            resolution='i')

x, y = m(lon, lat)

# draw background
m.drawcoastlines(linewidth=0.5)
m.drawmeridians(np.arange(lon.min(), lon.max(), 1.),labels=[0,0,0,1],fontsize=6, linewidth=.2)
m.drawparallels(np.arange(lat.min(), lat.max(), 1.),labels=[1,0,0,0],fontsize=6, linewidth=.2)

pc = m.pcolormesh(x, y, v[0, :, :], cmap='Reds')
cb = m.colorbar()
uproj, vproj = m.transform_vector(u10[0, :, :], v10[0, :, :], lons, lats, len(lons), len(lats))
Q = m.quiver(x, y, uproj, vproj, scale=100)
qk = plt.quiverkey(Q, 0.35, 0.8, 10, r'10 m$\cdot$s$^{-1}$', labelpos='E',
                   coordinates='figure')

for i in range(len(time[:240])):

    uproj, vproj = m.transform_vector(u10[i, :, :], v10[i, :, :], lons, lats, len(lons), len(lats))
    pc = m.pcolormesh(x, y, v[i, :, :], cmap='Reds')
    cb.set_clim([v[i, :, :].min(), v[i, :, :].max()])
    Q = m.quiver(x, y, uproj, vproj, scale=100)

    plt.savefig(fig_dir + 'Wind_' + var + '_' + "%.03f" % time[i] + '.png')
    pc.remove()
    Q.remove()
    # plt.close()

