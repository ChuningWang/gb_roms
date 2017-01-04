# Plot Initial Condition
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# ------------------------------------------------------------------------------------------------------
# Read grid
fh = nc.Dataset('../data/GlacierBay_ic_1980_01_08_SODA3.3.1.nc', mode='r')
lon = fh.variables['lon_rho'][:]
# lon = lon-360
lat = fh.variables['lat_rho'][:]
h = fh.variables['h'][:]
s = fh.variables['salt'][:]
msk = fh.variables['mask_rho'][:]
fh.close()

plt.close()
fig = plt.figure()

lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -134.5
lon_0 = 0.5 * (lon_min + lon_max)

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f')

# xstn, ystn = m(ctd['lon_stn'],ctd['lat_stn'])
# m.plot(xstn,ystn,'ok',ms=2)

# for i in stn:
#     plt.text(xstn[i]+1000, ystn[i], "%02d"%i, fontsize=5, va='center')

m.drawcoastlines(linewidth=.2)
m.fillcontinents(color='lightgrey')
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),labels=[0,0,0,1],fontsize=6, linewidth=.2)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),labels=[1,0,0,0],fontsize=6, linewidth=.2)
# setlabelrot(mr,-90)

xh, yh = m(lon, lat)
m.contourf(xh, yh, s[:, 0, :, :].squeeze())
m.colorbar()

plt.savefig('../figs/map_ic_s.eps',format='eps')
plt.close()

