import numpy as np
import netCDF4 as nc
from matplotlib import path

import pyroms
import sys

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
data_dir = sv['soda_dir']
model_dir = sv['model_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

grd = pyroms.grid.get_ROMS_grid(grd1)
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho
Cs_r = grd.vgrid.Cs_r
h = grd.vgrid.h
eta, xi = h.shape
N = grd.vgrid.N
zz = np.tile(Cs_r, (eta, xi, 1)).transpose((2, 0, 1))*np.tile(h, (N, 1, 1))
mskz = zz<=-50

ic_file = model_dir + 'tmpdir_GB-CIRC/GB-CIRC_rst.nc'
ic_file = out_dir + 'bc_ic/' + grd.name + '_ic_2008_04_15_CTD_floodFill.nc'
fh = nc.Dataset(ic_file, 'a')
salt = fh.variables['salt'][:]
dims = len(salt.shape)
if dims == 5:
    mskd = ~salt.mask[0, 0, :, :, :]
else:
    mskd = ~salt.mask[0, :, :, :]

dye_1 = np.zeros(salt.shape)
dye_2 = np.zeros(salt.shape)
dye_3 = np.zeros(salt.shape)

# box for tracers
box1 = np.array([[-136.10, 58.40],
                 [-135.975, 58.3125],
                 [-135.925, 58.3125],
                 [-135.875, 58.40],
                 [-135.50, 58.90],
                 [-136.00, 59.20],
                 [-137.25, 59.20],
                 [-137.25, 58.75]])

box2 = np.array([[-136.10, 58.40],
                 [-135.975, 58.3125],
                 [-135.925, 58.3125],
                 [-135.875, 58.40],
                 [-135.50, 58.55],
                 [-135.00, 58.25],
                 [-135.00, 57.90],
                 [-136.75, 57.90],
                 [-136.75, 58.60]])

p1 = path.Path(box1)
pc1 = p1.contains_points(np.array([lon.flatten(), lat.flatten()]).T).reshape((eta, xi))
msk1 = np.tile(pc1, (N, 1, 1)) & mskz & mskd
if dims == 5:
    dye_1[:, :, msk1] = 1
else:
    dye_1[:, msk1] = 1

p2 = path.Path(box2)
pc2 = p2.contains_points(np.array([lon.flatten(), lat.flatten()]).T).reshape((eta, xi))
msk2 = np.tile(pc2, (N, 1, 1)) & mskz & mskd
if dims == 5:
    dye_2[:, :, msk1] = 1
else:
    dye_2[:, msk1] = 1

dye_1 = np.ma.masked_where(salt.mask, dye_1)
dye_2 = np.ma.masked_where(salt.mask, dye_2)
dye_3 = np.ma.masked_where(salt.mask, dye_3)

# add dye tracers to ic file
if dims == 5:
    fh.createVariable('dye_01', 'f8', ('ocean_time', 'two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1.e37)
else:
    fh.createVariable('dye_01', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1.e37)
fh.variables['dye_01'].long_name = "Glacier Bay deep water dye"
fh.variables['dye_01'].units = "kilogram meter-3"
fh.variables['dye_01'].time = "ocean_time"
fh.variables['dye_01'].grid = "grid"
fh.variables['dye_01'].location = "face"
fh.variables['dye_01'].coordinates = "lon_rho lat_rho s_rho ocean_time"
fh.variables['dye_01'].field = "dye_01, scalar, series"
fh.variables['dye_01'][:] = dye_1

if dims == 5:
    fh.createVariable('dye_02', 'f8', ('ocean_time', 'two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1.e37)
else:
    fh.createVariable('dye_02', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1.e37)
fh.variables['dye_02'].long_name = "Shelf deep water dye"
fh.variables['dye_02'].units = "kilogram meter-3"
fh.variables['dye_02'].time = "ocean_time"
fh.variables['dye_02'].grid = "grid"
fh.variables['dye_02'].location = "face"
fh.variables['dye_02'].coordinates = "lon_rho lat_rho s_rho ocean_time"
fh.variables['dye_02'].field = "dye_02, scalar, series"
fh.variables['dye_02'][:] = dye_2

if dims == 5:
    fh.createVariable('dye_03', 'f8', ('ocean_time', 'two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1.e37)
else:
    fh.createVariable('dye_03', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=1.e37)
fh.variables['dye_03'].long_name = "River plume dye"
fh.variables['dye_03'].units = "kilogram meter-3"
fh.variables['dye_03'].time = "ocean_time"
fh.variables['dye_03'].grid = "grid"
fh.variables['dye_03'].location = "face"
fh.variables['dye_03'].coordinates = "lon_rho lat_rho s_rho ocean_time"
fh.variables['dye_03'].field = "dye_03, scalar, series"
fh.variables['dye_03'][:] = dye_3

fh.close()
