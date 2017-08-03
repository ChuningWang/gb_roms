import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc

fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/runoff_GB_hill.nc', 'r')
lon = fh.variables['lon_rho'][:]
lat = fh.variables['lat_rho'][:]
t = fh.variables['time'][:]
t = t-t[0]  # yearday
runoff = fh.variables['Runoff'][:]
fh.close()

for i in range(len(t)):
    plt.figure()
    plt.pcolor(lon, lat, np.log10(runoff[i, :, :]))
    plt.clim(-4, 4)
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'Discharge (log scale) [m$^3$s$^{-1}$]')
    ttl = "%03d" % t[i]
    plt.title('Day '+ttl)
    plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/runoff2/'+ttl+'.tiff', format='tiff')
    plt.close()
