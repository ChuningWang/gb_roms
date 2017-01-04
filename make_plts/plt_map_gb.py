# Plot Glacier Bay map
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
# clevs_land = 4*10**np.linspace(0, 3, 101)
# clevs_land = np.linspace(0, 2000, 11)
# clevs_land = []
# clevs_sea = np.linspace(-200, 0, 11)
# clevs_sea = []
# clevs = np.concatenate((clevs_sea, clevs_land))
clevs = np.linspace(-200, 200, 21)
# ------------------------------------------------------------------------------------------------------
# ctd = rd_ctd('../data/ctd.nc')

# Read topography
fh = nc.Dataset('../data/ARDEMv2.0.nc', mode='r')
lont = fh.variables['lon'][:]
lont = lont-360
latt = fh.variables['lat'][:]
z = fh.variables['z'][:]
fh.close()

# stn = np.arange(24)

# Read grid
fh = nc.Dataset('../data/GB_grd_nocurv.nc', mode='r')
lon = fh.variables['lon_rho'][:]
# lon = lon-360
lat = fh.variables['lat_rho'][:]
h = fh.variables['h'][:]
hraw = np.squeeze(fh.variables['hraw'][:])
s = fh.variables['s_rho'][:]
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

x, y = m(lon, lat)
x[msk==0] = np.nan
y[msk==0] = np.nan
m.plot(x, y, 'k', lw=0.1)
m.plot(x.T, y.T, 'k', lw=0.1)

# # plot topography

msk1 = (lont>lon_min) & (lont<lon_max)
msk2 = (latt>lat_min) & (latt<lat_max)

lont = lont[msk1]
latt = latt[msk2]
z = z[msk2,:][:,msk1]

z = -z
z[z>500] = 500
z[z<0] = 0
clevs = np.arange(0, 501, 10)

lont, latt = np.meshgrid(lont, latt)
xt, yt = m(lont, latt)

# m.contourf(xt, yt, z, clev)
# m.colorbar()

xh, yh = m(lon, lat)
m.contourf(xh, yh, h, clevs)
m.colorbar()

plt.savefig('../figs/map_h.eps',format='eps')
plt.close()

