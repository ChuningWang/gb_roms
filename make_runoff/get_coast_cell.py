import numpy as np
import netCDF4 as netCDF
import matplotlib.pyplot as plt

print 'Load lat_lon'
fh = netCDF.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/hydrology/GOA_RUNOFF.2000.nc', 'r')
lon = fh.variables['lon'][:]
lat = fh.variables['lat'][:]
topo = fh.variables['topo'][:]
strm = fh.variables['strm'][:]
Mp, Lp = lon.shape
fh.close()

msk_w = topo.mask*1 
msk = msk_w.copy()
# stack the mask over surrounding grid cells
msk[:-1, :] = msk[:-1, :]+msk_w[1:, :]
msk[1:, :] = msk[1:, :]+msk_w[:-1, :]
msk[:, :-1] = msk[:, :-1]+msk_w[:, 1:]
msk[:, 1:] = msk[:, 1:]+msk_w[:, :-1]
# diagonal
msk[:-1, :-1] = msk[:-1, :-1]+msk_w[1:, 1:]
msk[:-1, 1:] = msk[:-1, 1:]+msk_w[1:, :-1]
msk[1:, 1:] = msk[1:, 1:]+msk_w[:-1, :-1]
msk[1:, :-1] = msk[1:, :-1]+msk_w[:-1, 1:]

msk_coast = ((msk_w==0) & (msk>0))*1

plt.figure()
plt.pcolor(lon, lat, strm)
plt.show()

