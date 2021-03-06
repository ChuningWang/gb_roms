import numpy as np
from datetime import datetime
import netCDF4 as netCDF
import sys
from matplotlib import path

import pyroms
import pyroms_toolbox

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
home_dir = sv['home_dir']

grd1 = in_dir + 'gb_discharge.nc'
if len(sys.argv)>1:
    grd2 = sys.argv[-1]
else:
    grd2 = 'GB_lr'

# load 2-dimentional interannual discharge data 
print 'Load lat_lon'
nc_data = netCDF.Dataset(grd1, 'r')
lon = nc_data.variables['lon'][:]
lat = nc_data.variables['lat'][:]
# the coordinate of Hill's grid is slightly off - correct it here
lon = lon-0.010
lat = lat-0.010
mask = nc_data.variables['coast'][:]
# here use polygon to mask data out of Glacier Bay
box = np.array([[-137.40, 59.10],
                [-137.00, 58.50],
                [-136.55, 58.30],
                [-136.40, 58.15],
                [-136.00, 57.95],
                [-135.00, 58.05],
                [-136.10, 59.35]])
p0 = path.Path(box)
coords = np.array([lon.flatten(), lat.flatten()]).T
msk2 = p0.contains_points(coords)
msk2 = np.reshape(msk2, mask.shape)
mask = np.ma.masked_where(~msk2, mask)
mask = np.where(mask < 0, 0, mask)

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

##  create data remap file for scrip
print 'Create remap grid file for Hill and Beamer runoff'
remap_filename = 'remap_grid_runoff.nc'
nc = netCDF.Dataset(remap_filename, 'w', format='NETCDF3_CLASSIC')
nc.Description = 'remap grid file for Hill and Beamer runoff data'
nc.Author = 'build_runoff'
nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
nc.title = 'Hill and Beamer runoff'

grid_center_lon = lon.flatten()
grid_center_lat = lat.flatten()
Mp, Lp = lon.shape
grid_imask = mask.flatten()
grid_size = Lp * Mp
grid_corner_lon = np.zeros((grid_size, 4))
grid_corner_lat = np.zeros((grid_size, 4))
k = 0
for j in range(Mp):
    for i in range(Lp):
        grid_corner_lon[k,0] = lon_corner[j,i]
        grid_corner_lat[k,0] = lat_corner[j,i]
        grid_corner_lon[k,1] = lon_corner[j,i+1]
        grid_corner_lat[k,1] = lat_corner[j,i+1]
        grid_corner_lon[k,2] = lon_corner[j+1,i+1]
        grid_corner_lat[k,2] = lat_corner[j+1,i+1]
        grid_corner_lon[k,3] = lon_corner[j+1,i]
        grid_corner_lat[k,3] = lat_corner[j+1,i]
        k = k + 1

nc.createDimension('grid_size', grid_size)
nc.createDimension('grid_corners', 4)
nc.createDimension('grid_rank', 2)

nc.createVariable('grid_dims', 'i4', ('grid_rank'))
nc.variables['grid_dims'].long_name = 'grid size along x and y axis'
nc.variables['grid_dims'].units = 'None'
nc.variables['grid_dims'][:] = [(Lp, Mp)]

nc.createVariable('grid_center_lon', 'f8', ('grid_size'))
nc.variables['grid_center_lon'].long_name = 'longitude of cell center'
nc.variables['grid_center_lon'].units = 'degrees'
nc.variables['grid_center_lon'][:] = grid_center_lon

nc.createVariable('grid_center_lat', 'f8', ('grid_size'))
nc.variables['grid_center_lat'].long_name = 'latitude of cell center'
nc.variables['grid_center_lat'].units = 'degrees'
nc.variables['grid_center_lat'][:] = grid_center_lat

nc.createVariable('grid_imask', 'i4', ('grid_size'))
nc.variables['grid_imask'].long_name = 'mask'
nc.variables['grid_imask'].units = 'None'
nc.variables['grid_imask'][:] = grid_imask

nc.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
nc.variables['grid_corner_lon'].long_name = 'longitude of cell corner'
nc.variables['grid_corner_lon'].units = 'degrees'
nc.variables['grid_corner_lon'][:] = grid_corner_lon

nc.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
nc.variables['grid_corner_lat'].long_name = 'latitude of cell corner'
nc.variables['grid_corner_lat'].units = 'degrees'
nc.variables['grid_corner_lat'][:] = grid_corner_lat

nc.close()


#  create remap file for scrip
print 'Create remap grid file for target grid'
dstgrd = pyroms.grid.get_ROMS_grid(grd2)
dstgrd.hgrid.mask_rho = np.ones(dstgrd.hgrid.mask_rho.shape)
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')

## compute remap weights
print 'compute remap weights using scrip'
# input namelist variables for conservative remapping at rho points
grid1_file = remap_filename
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_runoff_to_' + dstgrd.name + '_conservative_nomask.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_runoff_conservative_nomask.nc'
map1_name = 'runoff to ' + dstgrd.name + ' conservative Mapping'
map2_name = dstgrd.name + ' to runoff conservative Mapping'
num_maps = 1
map_method = 'conservative'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)
