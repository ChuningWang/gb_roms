import netCDF4 as nc
import numpy as np
from geopy.distance import vincenty

pth = '/Volumes/R1/ROMS/hydrology/GOA/'
fh = nc.Dataset(pth+'lat_lon.nc')
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]
lon = lon-360

# load coast cell mask
coast = np.squeeze(fh.variables['coast_cells'][:])
fh.close()

r = 6371000

# get grid corner coordinates
Mp, Lp = lon.shape
lon_corner = np.zeros([Mp+1,Lp+1])
lat_corner = np.zeros([Mp+1,Lp+1])
lon_corner[1:Mp,1:Lp] = 0.25*(lon[:Mp-1,:Lp-1] + lon[1:Mp,:Lp-1] + \
                              lon[:Mp-1,1:Lp] + lon[1:Mp,1:Lp])
lat_corner[1:Mp,1:Lp] = 0.25*(lat[:Mp-1,:Lp-1] + lat[1:Mp,:Lp-1] + \
                              lat[:Mp-1,1:Lp] + lat[1:Mp,1:Lp])
lon_corner[0,1:Lp] = 2*lon_corner[1,1:Lp] - lon_corner[2,1:Lp]
lon_corner[Mp,1:Lp] = 2*lon_corner[Mp-1,1:Lp] - lon_corner[Mp-2,1:Lp]
lon_corner[:,0] = 2*lon_corner[:,1] - lon_corner[:,2]
lon_corner[:,Lp] = 2*lon_corner[:,Lp-1] - lon_corner[:,Lp-2]
lat_corner[0,1:Lp] = 2*lat_corner[1,1:Lp] - lat_corner[2,1:Lp]
lat_corner[Mp,1:Lp] = 2*lat_corner[Mp-1,1:Lp] - lat_corner[Mp-2,1:Lp]
lat_corner[:,0] = 2*lat_corner[:,1] - lat_corner[:,2]
lat_corner[:,Lp] = 2*lat_corner[:,Lp-1] - lat_corner[:,Lp-2]
 
dphi_h = np.diff(lat_corner, axis=0)*np.pi/180.
dphi_v = np.diff(lat_corner, axis=1)*np.pi/180.

sint = np.sin(lon_corner*np.pi/180.)

ws =  r*dphi_h[:, :-1]*0.5*(sint[:-1, :-1] + sint[1:, :-1])
wn = -r*dphi_h[:, :1 ]*0.5*(sint[:-1, 1: ] + sint[1:, 1: ])
ww =  r*dphi_v[:-1, :]*0.5*(sint[:-1, :-1] + sint[:-1, 1:])
we = -r*dphi_v[1: , :]*0.5*(sint[1: , :-1] + sint[1: , 1:])

weight = ws+wn+ww+we


dx = np.diff(lon_corner, axis=1)
dy = np.diff(lat_corner, axis=0)
dx = 0.5*(dx[:-1, :]+dx[1:, :])
dy = 0.5*(dy[:, :-1]+dy[:, 1:])

# calculate grid area in [m^2]
area = np.zeros(lon.shape)
for i in range(Mp):
    for j in range(Lp):
        dxi = vincenty((lat[i, j], lon[i, j]-0.5*dx[i, j]), (lat[i, j], lon[i, j]+0.5*dx[i, j])).meters
        dyi = vincenty((lat[i, j]-0.5*dy[i, j], lon[i, j]), (lat[i, j]+0.5*dy[i, j], lon[i, j])).meters
        area[i, j] = dxi*dyi
 
