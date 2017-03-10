# Plot Glacier Bay map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc
import pyroms

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# Read grid
grd = pyroms.grid.get_ROMS_grid('GB')
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
msk = grd.hgrid.mask_rho

plt.close()
fig = plt.figure()
plt.pcolormesh(msk, lw=0.1)
plt.show()

# plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/map_grd.eps',format='eps')
# plt.close()

