import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import pyroms

# blowup point (446, 0269, 13)

fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/outputs/UPWELLING_ICE_his_00004.nc')

u = fh.variables['u'][:]
v = fh.variables['v'][:]
w = fh.variables['w'][:]

fh.close()

# find maximum w
ind = np.where(w == w.max())
zlevel = ind[1].data[0]
xx = ind[2].data[0]
yy = ind[3].data[0]

grd = pyroms.grid.get_ROMS_grid('GB')

plt.figure()
plt.pcolor(w[0, zlevel, :, :].squeeze())
plt.colorbar()
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/blowup_w.png', format='png', dpi=900)
plt.close()

plt.figure()
plt.pcolor(grd.vgrid.h)
plt.colorbar()
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/depth.png', format='png', dpi=900)
plt.close()

