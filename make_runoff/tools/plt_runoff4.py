import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc

fh = nc.Dataset('../Glacier_Bay_rivers.nc', 'r')
lon = fh.variables['river_Xposition'][:]
lat = fh.variables['river_Eposition'][:]
t = fh.variables['river_time'][:]
t = t-t[0]  # yearday
runoff = fh.variables['river_transport'][:]
fh.close()

for i in range(len(t)):
    plt.figure()
    plt.scatter(lon, lat, c=np.log10(runoff[i, :]), edgecolors='none', alpha=0.5)
    plt.xlim(0, 600)
    plt.ylim(0, 1000)
    plt.clim(-4, 4)
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'Discharge (log scale) [m$^3$s$^{-1}$]')
    ttl = "%03d" % t[i]
    plt.title('Day '+ttl)
    plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/runoff4/'+ttl+'.tiff', format='tiff')
    plt.close()
