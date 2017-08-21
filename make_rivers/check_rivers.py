import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['out_dir']

grd1 = 'GB_USGS'
my_year = 2008

grd = pyroms.grid.get_ROMS_grid(grd1)
msk = grd.hgrid.mask_rho

fh = nc.Dataset(data_dir + 'frc/GlacierBay_usgs_rivers_'+str(my_year)+'_Hill.nc', 'r')
rsign = fh.variables['river_sign'][:]
rdir = fh.variables['river_direction'][:]
xpos = fh.variables['river_Xposition'][:]
epos = fh.variables['river_Eposition'][:]
rtrans = fh.variables['river_transport'][:]
fh.close()

# find rho index
xrho = xpos.copy()
erho = epos.copy()
for i in range(len(rsign)):
    if rsign[i]==-1:
        if rdir[i]==0:
            xrho[i] = xrho[i]-1
        if rdir[i]==1:
            erho[i] = erho[i]-1
            
rho_msk = np.array([msk[erho[i], xrho[i]] for i in range(len(rsign))])
