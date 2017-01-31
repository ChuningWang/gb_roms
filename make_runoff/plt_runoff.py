import numpy as np
import matplotlib.pyplot as plt
import pyroms
import netCDF4 as nc

fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/GOA_RUNOFF.2000.nc', 'r')
t = fh.variables['time'][:]
t = t-t[0]+1  # year day
lon = fh.variables['lon'][:]
lat = fh.variables['lat'][:]
strm = fh.variables['strm'][:]
r = fh.variables['q'][:]
topo = fh.variables['topo'][:]
fh.close()

plt.figure()
plt.pcolor(lon, lat, topo)
# plt.clim(-4, 4)
cb = plt.colorbar()
cb.ax.set_ylabel(r'Height [m]')
ttl = 'Topography'
plt.title(ttl)
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/runoff/'+ttl+'.tiff', format='tiff')
plt.close()

plt.figure()
plt.pcolor(lon, lat, strm)
# plt.clim(-4, 4)
cb = plt.colorbar()
cb.ax.set_ylabel(r'Stream Line (log scale)')
ttl = 'StreamLine'
plt.title(ttl)
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/runoff/'+ttl+'.tiff', format='tiff')
plt.close()

# for i in range(len(t)):
#     plt.figure()
#     plt.pcolor(lon, lat, np.log10(r[i, :, :]))
#     plt.clim(-4, 4)
#     cb = plt.colorbar()
#     cb.ax.set_ylabel(r'Discharge (log scale) [m$^3$s$^{-1}$]')
#     ttl = "%03d" % t[i]
#     plt.title('Day '+ttl)
#     plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/runoff/'+ttl+'.tiff', format='tiff')
#     plt.close()
