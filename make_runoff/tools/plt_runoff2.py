import numpy as np
from datetime import datetime
import netCDF4 as netCDF
import matplotlib.pyplot as plt

nc_data = netCDF.Dataset('/Volumes/R1/ROMS/hydrology/GOA/lat_lon.nc', 'r')
lon = nc_data.variables['lon'][:]
lon = lon-360
lat = nc_data.variables['lat'][:]
coast = nc_data.variables['coast_cells'][:]
coast = np.squeeze(coast)
Mp, Lp = lon.shape
nc_data.close()

nc_data = netCDF.Dataset('/Volumes/R1/ROMS/hydrology/GOA/discharge_1980_1981.nc', 'r')
r = nc_data.variables['discharge'][180, :, :]
r = np.squeeze(r)
r[r<0] = np.nan
nc_data.close()

msk1 = (lon[0, :]>=-139) & (lon[0, :]<=-135)
msk2 = (lat[:, 0]>=58)   & (lat[:, 0]<=61)

lon = lon[msk2, :][:, msk1]
lat = lat[msk2, :][:, msk1]
coast = coast[msk2, :][:, msk1]
r = r[msk2, :][:, msk1]

plt.figure()
plt.pcolor(lon, lat, coast)
plt.colorbar()
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/coast_old.eps', format='eps')
plt.close()

plt.figure()
plt.pcolor(lon, lat, r)
plt.colorbar()
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/discharge_old.eps', format='eps')
plt.close()

# plt.figure()
# plt.pcolor(lon, lat, topo)
# plt.colorbar()
# plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/topo_old.eps', format='eps')
# plt.close()
