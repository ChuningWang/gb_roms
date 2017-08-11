import subprocess
import os
import numpy as np
import pyroms

import netCDF4 as nc

import scipy as sp
import scipy.ndimage

from gb_toolbox import gb_ctd

'''
The goal of this script is to reconstruct January T and S at all depth as initial condition
for Glacier Bay ROMS model. The purpose of using CTD instead of SODA data is that SODA data is
Too coarse for detailed mapping of T and S inside the bay.

The whole script takes about 1 hour to run.
'''

def floodFill(c, r, mask): 
    """
    produce a mask array that identifies the flood of ctd stations.
    """

    # number of seeds
    seednum = c.size
    # cells already filled
    filled = set()
    # cells to fill
    fill = [set() for _ in range(seednum)] 
    for i in range(seednum):
        fill[i].add((c[i], r[i]))

    eta = mask.shape[0]
    xi = mask.shape[1]

    # Our output inundation array
    flood = np.ones_like(mask, dtype=np.int8)*-1
    # Loop through and modify the cells which
    # need to be checked.

    while any(fill):
        fill_p1 = [set() for _ in range(seednum)]

        # loop through the stations
        for i in range(seednum):

            while fill[i]:
                # Grab a cell
                x, y = fill[i].pop()
                if x == eta or y == xi or x < 0 or y < 0:
                    # Boundary, don't fill
                    continue
                if mask[x][y] == 1:
                    # Do fill
                    flood[x][y] = i
                    filled.add((x, y))
                    # Check neighbors for 1 values
                    west = (x-1, y)
                    east = (x+1, y)
                    north = (x, y-1)
                    south = (x, y+1)
                    if west not in filled:
                        fill_p1[i].add(west)
                    if east not in filled:
                        fill_p1[i].add(east)
                    if north not in filled:
                        fill_p1[i].add(north)
                    if south not in filled:
                        fill_p1[i].add(south)

        fill = fill_p1[:]

    return flood

def nan_gaussian_filter(U, sigma):
    '''
    http://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    '''
    V = U.copy()
    msk = U!=U
    V[msk] = 0
    VV = sp.ndimage.gaussian_filter(V, sigma = sigma)

    W = 0*U.copy()+1
    W[msk] = 0
    WW = sp.ndimage.gaussian_filter(W, sigma = sigma)

    Z=VV/WW
    Z[msk] = np.NaN
    return Z

# ----------------------------------------------------------------
# preparation
# output dir
import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
dst_dir = sv['out_dir']
soda_dir = sv['soda_dir']

grd1 = 'GB_USGS'
month = 7

# Load target grid and land mask
dst_grd = pyroms.grid.get_ROMS_grid(grd1)
msk = dst_grd.hgrid.mask
h = dst_grd.vgrid.h.squeeze()
lat_grd = dst_grd.hgrid.lat_rho
lon_grd = dst_grd.hgrid.lon_rho
lat_grd[msk==0] = np.NaN
lon_grd[msk==0] = np.NaN

eta = msk.shape[0]
xi = msk.shape[1]

# ----------------------------------------------------------------
# Load CTD climatology
ctd_dir = in_dir + 'ctd.nc'
ctd = gb_ctd.rd_ctd(ctd_dir)
ctd_clim = gb_ctd.cal_ctd_climatology(ctd)

# Extract data
lat = ctd_clim['lat']
lon = ctd_clim['lon']
t = ctd_clim['t'][month-1, :, :].squeeze()
s = ctd_clim['s'][month-1, :, :].squeeze()

# Clear NaN profiles in data
si = np.where(np.isnan(t[:, 0]))[0]
lat = np.delete(lat, si, axis=0)
lon = np.delete(lon, si, axis=0)
t = np.delete(t, si, axis=0)
s = np.delete(s, si, axis=0)

hctd = ctd_clim['d']
zlev = hctd.shape[0]

# ----------------------------------------------------------------
# Load SODA data
tag = '2000_01'
soda_filein = soda_dir + 'monthly/soda3.3.1_monthly_ocean_reg_' + tag + '.nc'

xt = 446
yt = 265

fh = nc.Dataset(soda_filein, 'r')
lont_soda = fh.variables['xt_ocean'][xt]
lont_soda = lont_soda-360
latt_soda = fh.variables['yt_ocean'][yt]
zt_soda = fh.variables['st_ocean'][:].squeeze()
t_soda = fh.variables['temp'][:, :, yt, xt].squeeze()
s_soda = fh.variables['salt'][:, :, yt, xt].squeeze()
spval = fh.variables['temp'].missing_value
fh.close()

t_soda[t_soda.mask] = np.NaN
s_soda[s_soda.mask] = np.NaN

zt_soda = np.insert(zt_soda, 0, 0)
t_soda = np.insert(t_soda, 0, t_soda[0])
s_soda = np.insert(s_soda, 0, s_soda[0])

t_soda = np.interp(hctd, zt_soda, t_soda)
s_soda = np.interp(hctd, zt_soda, s_soda)

# add the profile as a station in Glacier Bay
lat = np.append(lat, 58.1)
lon = np.append(lon, -136.7)
t = np.vstack((t, t_soda))
s = np.vstack((s, s_soda))

xt = 449
yt = 266

fh = nc.Dataset(soda_filein, 'r')
lont_soda = fh.variables['xt_ocean'][xt]
lont_soda = lont_soda-360
latt_soda = fh.variables['yt_ocean'][yt]
zt_soda = fh.variables['st_ocean'][:].squeeze()
t_soda = fh.variables['temp'][:, :, yt, xt].squeeze()
s_soda = fh.variables['salt'][:, :, yt, xt].squeeze()
fh.close()

t_soda[t_soda.mask] = np.NaN
s_soda[s_soda.mask] = np.NaN

zt_soda = np.insert(zt_soda, 0, 0)
t_soda = np.insert(t_soda, 0, t_soda[0])
s_soda = np.insert(s_soda, 0, s_soda[0])

t_soda = np.interp(hctd, zt_soda, t_soda)
s_soda = np.interp(hctd, zt_soda, s_soda)

# add the profile as a station in Glacier Bay
lat = np.append(lat, 58.12)
lon = np.append(lon, -135.07)
t = np.vstack((t, t_soda))
s = np.vstack((s, s_soda))

# t[t<-1000] = np.NaN
# s[s<-1000] = np.NaN

# ----------------------------------------------------------------
# flood CTD stations to the whole grid
seednum = lat.size

c = np.zeros(seednum)
r = np.zeros(seednum)

for i in range(seednum):
    dis = (lat_grd-lat[i])**2+(lon_grd-lon[i])**2
    c[i], r[i] = np.where(dis==np.nanmin(dis))

fl = np.zeros((msk.shape[0], msk.shape[1], zlev))
t_ini = np.zeros((1, zlev, eta, xi))*np.NaN
s_ini = np.zeros((1, zlev, eta, xi))*np.NaN

for i in range(zlev):
    mask = (msk==1) & (h>=i)
    fl[:, :, i] = floodFill(c, r, mask)

    t_sl = np.zeros((eta, xi))*np.NaN
    s_sl = np.zeros((eta, xi))*np.NaN

    for j in range(seednum):
        msk_fl = fl[:, :, i]==j
        t_sl[msk_fl] = t[j, i]
        s_sl[msk_fl] = s[j, i]

    t_ini[:, i, :, :] = t_sl
    s_ini[:, i, :, :] = s_sl

# flood deep water value
for i in range(eta):
    for j in range(xi):
        if msk[i, j]==1:
            bi = np.where(np.isnan(t_ini[:, :, i, j].squeeze()))[0][0]-1
            tb = t_ini[:, bi, i, j]
            sb = s_ini[:, bi, i, j]
            if hctd[bi]<=h[i, j]:
                t_ini[:, bi:int(h[i, j])+1, i, j] = tb
                s_ini[:, bi:int(h[i, j])+1, i, j] = sb

# ----------------------------------------------------------------
# post-process
# filter the results
for i in range(zlev):
    t_ini[:, i, :, :] = nan_gaussian_filter(t_ini[:, i, :, :], 10.0)
    s_ini[:, i, :, :] = nan_gaussian_filter(s_ini[:, i, :, :], 10.0)

# mask data below water depth
# t_ini = np.ma.masked_invalid(t_ini)
# s_ini = np.ma.masked_invalid(s_ini)

# interpolate to S-coord 
# let 1 m extend to surface
hctd[0] = 1
t_ini2 = np.ones((1, dst_grd.vgrid.N, eta, xi ))*spval
s_ini2 = np.ones((1, dst_grd.vgrid.N, eta, xi ))*spval
zs = dst_grd.vgrid.z_r[:]

for i in range(eta):
    for j in range(xi):
        if msk[i, j]==1:
            t_ini2[:, :, i, j] = np.interp(-zs[:, i, j], hctd, t_ini[:, :, i, j].squeeze())
            s_ini2[:, :, i, j] = np.interp(-zs[:, i, j], hctd, s_ini[:, :, i, j].squeeze())

t_ini2 = np.ma.masked_where(t_ini2==spval, t_ini2)
s_ini2 = np.ma.masked_where(s_ini2==spval, s_ini2)

# ----------------------------------------------------------------
# write into nc file
src_grd_name = '2000_01_03_SODA3.3.1'
ic_file = dst_dir + 'bc_ic/' + dst_grd.name + '_ic_' + src_grd_name + '.nc'
fh = nc.Dataset(ic_file, 'r+') 
fh.variables['temp'][:, :, :, :] = t_ini2
fh.variables['salt'][:, :, :, :] = s_ini2
fh.close()

