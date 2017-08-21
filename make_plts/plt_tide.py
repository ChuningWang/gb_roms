import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import pyroms
import pyroms_toolbox

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

grd1 = 'GB_USGS'
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
z = grd.vgrid.h
msk = grd.hgrid.mask_rho

fh = nc.Dataset(out_dir+'frc/'+grd.name+'_tides_otps.nc')
name = fh.variables['tide_name'][:]
prd = fh.variables['tide_period'][:]
Eamp = fh.variables['tide_Eamp'][:]
Epha = fh.variables['tide_Ephase'][:]
Cmax = fh.variables['tide_Cmax'][:]
Cmin = fh.variables['tide_Cmin'][:]
Cang = fh.variables['tide_Cangle'][:]
Cpha = fh.variables['tide_Cphase'][:]
fh.close()

consts_num = prd.shape[0]

for i in range(consts_num):

    plt.contourf(Cmax[i, :, :])
    plt.colorbar()
    # plt.contour(Epha[i, :, :], colors='k')
    # plt.contour(msk, np.array([0.5, 0.5]), linewidths=0.5, colors='k')
    plt.title(name[i])

    plt.savefig(out_dir+'figs/tides/tide' + str(i) + 'cmax.tiff',format='tiff')
    plt.close()

