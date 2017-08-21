import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import tidal_ellipse

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']


fh = nc.Dataset(in_dir+'tpxo8nc/grid_tpxo8_atlas6.nc')
latz = fh.variables['lat_z'][:]
lonz = fh.variables['lon_z'][:]
latu = fh.variables['lat_u'][:]
lonu = fh.variables['lon_u'][:]
latv = fh.variables['lat_v'][:]
lonv = fh.variables['lon_v'][:]
hz = fh.variables['hz'][:]
hu = fh.variables['hu'][:]
hv = fh.variables['hv'][:]
fh.close()

fh = nc.Dataset(in_dir+'tpxo8nc/hf.mf_tpxo8_atlas_6.nc')
hre = fh.variables['hRe'][:]
him = fh.variables['hIm'][:]
fh.close()

fh = nc.Dataset(in_dir+'tpxo8nc/uv.mf_tpxo8_atlas_6.nc')
ure = fh.variables['uRe'][:]
uim = fh.variables['uIm'][:]
vre = fh.variables['vRe'][:]
vim = fh.variables['vIm'][:]
fh.close()

ure = 0.5*(ure[:-1, :]+ure[1:, :])
uim = 0.5*(uim[:-1, :]+uim[1:, :])
vre = 0.5*(vre[:, :-1]+vre[:, 1:])
vim = 0.5*(vim[:, :-1]+vim[:, 1:])

latv = 0.5*(latv[:-1]+latv[1:])
lonu = 0.5*(lonu[:-1]+lonu[1:])

# ------------------------------------------------------------------------------------------
# extract data near Glacier Bay
lonmin = 220
lonmax = 230
latmin = 55
latmax = 60

mskz1 = (lonz>=lonmin) & (lonz<=lonmax)
mskz2 = (latz>=latmin) & (latz<=latmax)

lonz = lonz[mskz1]
latz = latz[mskz2]

hz = hz[mskz1, :][:, mskz2]
hz[hz<=0] = np.nan
hre = hre[mskz1, :][:, mskz2]
him = him[mskz1, :][:, mskz2]

msku1 = (lonu>=lonmin) & (lonu<=lonmax)
msku2 = (latu>=latmin) & (latu<=latmax)

ure = ure[msku1, :][:, msku2]
uim = uim[msku1, :][:, msku2]

mskv1 = (lonv>=lonmin) & (lonv<=lonmax)
mskv2 = (latv>=latmin) & (latv<=latmax)

vre = vre[mskv1, :][:, mskv2]
vim = vim[mskv1, :][:, mskv2]

# ------------------------------------------------------------------------------------------
# convert from h,u,v to tidal consts
hamp = (hre**2+him**2)**0.5  # mm
hpha = np.arctan(-him/hre)/np.pi*180.  # deg
uamp = (ure**2+uim**2)**0.5  # mm
upha = np.arctan(-uim/ure)/np.pi*180.  # deg
vamp = (vre**2+vim**2)**0.5  # mm
vpha = np.arctan(-vim/vre)/np.pi*180.  # deg

# convert ap to ep
cmax, ecc, cang, cpha, w = tidal_ellipse.ap2ep(uamp, upha, vamp, vpha)
cmin = cmax*ecc

# ------------------------------------------------------------------------------------------
# plotting
plt.figure()
m = Basemap(projection='merc', llcrnrlon=-137.5+360, llcrnrlat=57.75,
            urcrnrlon=-134.5+360, urcrnrlat=59.25,
            resolution='h')

m.drawcoastlines()
mr = m.drawmeridians(np.arange(-137.5, 57.75, 0.5),labels=[0,0,0,1],fontsize=6, linewidth=.2)
pr = m.drawparallels(np.arange(57.75, 59.25, 0.25),labels=[1,0,0,0],fontsize=6, linewidth=.2)

lat, lon = np.meshgrid(latz, lonz)
xh, yh = m(lon, lat)

CS = m.contour(xh, yh, hamp, colors='k')
# m.contour(xh, yh, hpha, colors='k')
m.pcolor(xh, yh, hz)
plt.clim(0,400)
plt.clabel(CS, fontsize=9, inline=1)
m.colorbar()

plt.savefig(out_dir+'figs/tides/tide_mf_orig.tiff',format='tiff')
plt.close()

