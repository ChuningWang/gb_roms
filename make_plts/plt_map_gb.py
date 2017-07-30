# Plot Glacier Bay map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc
import pyroms

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_USGS'
plt_map = 0

ctd = rd_ctd(in_dir + 'ctd.nc')
lat_ctd = ctd['lat_stn']
lon_ctd = ctd['lon_stn']

# Read grid
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
z = grd.vgrid.h
msk = grd.hgrid.mask_rho

plt.close()
fig = plt.figure()

lat_min = 58.00
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

if plt_map == 1:

    m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
                urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
                resolution='f')

    # m.drawcoastlines(linewidth=0.01)
    mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),labels=[0,0,0,1],fontsize=6, linewidth=.2)
    pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),labels=[1,0,0,0],fontsize=6, linewidth=.2)
    # setlabelrot(mr,-90)

    x, y = m(lon, lat)
    x2, y2 = m(lon_ctd, lat_ctd)
    m.pcolor(x, y, z, cmap='Greens')
    plt.clim(0, 400)
    plt.colorbar()
    m.contour(x, y, msk, np.array([0.5, 0.5]), linewidths=0.05, colors='k')
    m.plot(x2, y2, '.k', ms=2)

    plt.savefig(out_dir + 'figs/map_grd.tiff', format='tiff', dpi=600)

elif plt_map == 0:

    plt.pcolormesh(z, cmap='Greens')
    plt.clim(400, 450)
    plt.colorbar()
    plt.contour(msk, np.array([0.5, 0.5]), linewidths=0.05, colors='k')

    plt.savefig(out_dir + 'figs/map_grd_noproj.tiff', format='tiff', dpi=600)

plt.close()

