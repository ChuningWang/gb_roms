# Plot Glacier Bay map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc
import pyroms
import sys

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

plt_proj = 0
plt_hill_coast = 1

# ctd = rd_ctd(in_dir + 'ctd.nc')
# lat_ctd = ctd['lat_stn']
# lon_ctd = ctd['lon_stn']

# Read grid
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
z = grd.vgrid.h
msk = grd.hgrid.mask_rho
Mp, Np = msk.shape

plt.close()
fig = plt.figure()

lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

if plt_proj == 1:

    m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
                urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
                resolution='f')

    m.drawcoastlines(linewidth=0.01)
    mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),labels=[0,0,0,1],fontsize=6, linewidth=.2)
    pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),labels=[1,0,0,0],fontsize=6, linewidth=.2)
    # setlabelrot(mr,-90)

    x, y = m(lon, lat)
    m.pcolor(x, y, z, cmap='Greens')
    plt.clim(-1, 20)
    plt.colorbar()
    m.contour(x, y, msk, [0.5], linewidths=0.05, colors='k')
    # x2, y2 = m(lon_ctd, lat_ctd)
    # m.plot(x2, y2, '.k', ms=2)

    if plt_hill_coast==1:
        # Overlay Hill discharge point on the map
        fh = nc.Dataset(in_dir+'gb_discharge.nc', 'r')
        lonh = fh.variables['lon'][:]
        lath = fh.variables['lat'][:]
        coast = fh.variables['coast'][:]
        fh.close
        lonh = lonh[~coast.mask]
        lath = lath[~coast.mask]
        lonh = lonh-0.010
        lath = lath-0.005
        xh, yh = m(lonh, lath)
        m.plot(xh, yh, '.k', ms=2)

    plt.savefig(out_dir + 'figs/'+grd1+'_grd.tiff', format='tiff',dpi=600)

elif plt_proj == 0:

    plt.pcolormesh(z, cmap='Greens')
    plt.clim(-1, 20)
    plt.colorbar()
    plt.contour(msk, [0.5], linewidths=0.05, colors='k')

    plt.xticks(np.arange(0, Np, 10), rotation='vertical')
    plt.yticks(np.arange(0, Mp, 20))
    plt.tick_params(axis='both', which='major', labelsize=5)
    plt.grid()

    plt.savefig(out_dir + 'figs/'+grd1+'_grd_noproj.tiff', format='tiff',dpi=600)

plt.close()
