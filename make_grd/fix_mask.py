import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import netCDF4

import pyroms
import pyroms_toolbox
import bathy_smoother

# ----------------------------------------------------------------------------------------------------------
# defind the boundary of mapping domain
lat_min = 57.
lat_max = 60.
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -138.
lon_max = -134.
lon_0 = 0.5 * (lon_min + lon_max)

# load grid
hgrd = pyroms.grid.get_ROMS_hgrid('GB')

# use GUI masking tool
m = Basemap(projection='lcc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f')
# plt.ion()
# plt.figure()
# m.drawcoastlines()
# plt.show()
# pyroms.grid.edit_mask_mesh(hgrd, proj=m)
# or (this is faster 'I guess')
coast = pyroms.utility.get_coast_from_map(m)
# pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

