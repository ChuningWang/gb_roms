import pyroms
import pyroms_toolbox
import sys

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
soda_dir = sv['soda_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

# load the grid
# srcgrd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_dir + 'grid/SODA3_0.5deg_grid.nc', name='SODA3.3.1', xrange=(400, 500), yrange=(180, 280))
srcgrd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL(soda_dir + 'grid/SODA3_0.25deg_grid.nc', 
                                                     name='SODA3.3.1_0.25', xrange=(550, 600), yrange=(750, 800))
dstgrd = pyroms.grid.get_ROMS_grid(grd1)

# make remap grid file for scrip
pyroms_toolbox.BGrid_GFDL.make_remap_grid_file(srcgrd, Bpos='t')
pyroms_toolbox.BGrid_GFDL.make_remap_grid_file(srcgrd, Bpos='uv')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_' + srcgrd.name + '_t.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_t.nc'
map1_name = srcgrd.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name = dstgrd.name + ' to ' + srcgrd.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_' + srcgrd.name + '_uv.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_uv_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_uv.nc'
map1_name = srcgrd.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name = dstgrd.name + ' to ' + srcgrd.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_' + srcgrd.name + '_t.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_u.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_u_to_t.nc'
map1_name = srcgrd.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name = dstgrd.name + ' to ' + srcgrd.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_' + srcgrd.name + '_t.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_v.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_v_to_t.nc'
map1_name = srcgrd.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name = dstgrd.name + ' to ' + srcgrd.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)

