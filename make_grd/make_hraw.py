import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import pyroms

from scipy.interpolate import griddata
from matplotlib import path

import read_host_info
sv = read_host_info.read_host_info()
bathy_dir = sv['in_dir']
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_300m_orig'
hgrd = pyroms.grid.get_ROMS_hgrid(grd1)
lat = hgrd.lat_rho
lon = hgrd.lon_rho

# ------------------------------------------------------------------------
# defind the boundary of mapping domain
lat_min = 57.
lat_max = 60.
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -138.
lon_max = -134.
lon_0 = 0.5 * (lon_min + lon_max)

# ------------------------------------------------------------------------
# generate the base bathymetry
fh = nc.Dataset(bathy_dir + 'ARDEMv2.0.nc', mode='r')
h0 = fh.variables['z'][:]
lon0 = fh.variables['lon'][:] 
lon0[lon0>180] = lon0[lon0>180]-360  # lons from -180 to 180
lat0 = fh.variables['lat'][:] 
fh.close()

msk1 = (lon0>lon_min) & (lon0<lon_max)
msk2 = (lat0>lat_min) & (lat0<lat_max)

lon0 = lon0[msk1]
lat0 = lat0[msk2]
h0 = h0[msk2,:][:,msk1]

# depth positive
h0 = -h0

# interpolate new bathymetry
lon0, lat0 = np.meshgrid(lon0, lat0)
h = griddata((lon0.flatten(), lat0.flatten()), h0.flatten(), (lon, lat), method='linear')
print 'griddata from ARDEM done...'

# ------------------------------------------------------------------------
# fix bathymetry with USGS and NOAA data
fh = nc.Dataset(bathy_dir + 'bathy_noaa.nc', 'r')
lon1 = fh.variables['lon'][:]
lat1 = fh.variables['lat'][:]
h1 = fh.variables['z'][:]
fh.close()

bdry = np.loadtxt('bdry_usgs.txt')
p0 = [(bdry[i, 0], bdry[i, 1]) for i in range(len(bdry[:, 0]))]
p = path.Path(p0)
pc = ~p.contains_points(np.array([lon1, lat1]).T) 

lon1 = lon1[pc]
lat1 = lat1[pc]
h1 = h1[pc]

fh = nc.Dataset(bathy_dir + 'bathy_usgs.nc', 'r')
lon2 = fh.variables['lon'][:][1::3, 1::3]
lat2 = fh.variables['lat'][:][1::3, 1::3]
h2 = fh.variables['z'][:][1::3, 1::3]
fh.close()

msk = ~h2.mask
lon2 = lon2[msk]
lat2 = lat2[msk]
h2 = h2[msk]

lon0 = np.concatenate((lon1, lon2))
lat0 = np.concatenate((lat1, lat2))
h0 = np.concatenate((h1, h2))

# load grid boundary
bdry = np.loadtxt('bdry.txt')
x = bdry[:, 0]
y = bdry[:, 1]

p0 = [(x[i], y[i]) for i in range(len(x))]

p = path.Path(p0)
pc = p.contains_points(np.array([lon.flatten(), lat.flatten()]).T).reshape(h.shape)
 
# interpolate new bathymetry
h[pc] = griddata((lon0.flatten(), lat0.flatten()), h0.flatten(),
            (lon[pc], lat[pc]), method='linear')

# copy raw bathymetry
hraw = h.copy()

# save hraw to netcdf
fout = nc.Dataset(in_dir + 'hraw.nc', 'w')
fout.createDimension('x', h.shape[0])
fout.createDimension('y', h.shape[1])

fout.createVariable('lat', 'f8', ('x', 'y'))
fout.createVariable('lon', 'f8', ('x', 'y'))
fout.createVariable('hraw', 'f8', ('x', 'y'))
fout.variables['lat'][:] = lat
fout.variables['lon'][:] = lon
fout.variables['hraw'][:] = hraw
fout.close()
