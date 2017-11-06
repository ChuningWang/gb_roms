""" hypsometric plot of glacier bay """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
import netCDF4 as nc
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']

grd1 = 'GB_lr'
grd1_dir = '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_hr_grd.nc'
grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
h = grd.vgrid.hraw[0, :, :]

MP, NP = h.shape
coords = []
for i in range(MP):
    for j in range(NP):
        coords.append((lon[i, j], lat[i, j]))

fin = nc.Dataset(grd1_dir, 'r')
pm = fin.variables['pm'][:]
pn = fin.variables['pn'][:]
fin.close()

area = ((1./pm)/1000.)*((1./pn)/1000.)

# boundary of Glacier Bay
box = np.array([[-137.30, 58.75],
                [-137.30, 59.15],
                [-135.70, 59.15],
                [-135.70, 58.45],
                [-135.98, 58.37],
                [-136.58, 58.51]])

path_points = path.Path(box)
inbox = path_points.contains_points(coords)
inbox = np.reshape(inbox, (MP, NP))

h[~inbox] = -9999.
zbins = np.arange(450)
hypso = np.zeros(zbins.shape)
for z in zbins:
    msk = h >= z
    hypso[z] = area[msk].sum()

plt.figure()
plt.fill_between(hypso, zbins)
plt.ylim(zbins.min(), zbins.max())
plt.xlim(0, 400)
plt.gca().invert_yaxis()
plt.xlabel('Area [km$^2$]')
plt.ylabel('Depth [m]')
plt.title('Hypsometric Curve in Glacier Bay')
plt.savefig(out_dir + 'figs/hypso.png', dpi=300)
