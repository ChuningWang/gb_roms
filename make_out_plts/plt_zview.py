import numpy as np
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cmocean
from mpl_toolkits.basemap import Basemap
import pyroms
import pyroms_toolbox as prt
import netCDF4 as nc
import sys

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# my inputs
my_year = 2008
varlist = ['zeta', 'temp', 'salt', 'dye_03']
varlist = ['salt', 'dye_03']
depth = 1
dd = 10
uscale = 20

depth = -abs(depth)
lat_min = 57.75
lat_max = 59.25
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -137.5
lon_max = -135.0
lon_0 = 0.5 * (lon_min + lon_max)

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

grd = pyroms.grid.get_ROMS_grid(grd1)
lon = grd.hgrid.lon_rho
lat = grd.hgrid.lat_rho
mask = grd.hgrid.mask_rho

if grd1=='GB_lr':
    tag = 'GB-CIRC'
if grd1=='GB_hr':
    tag = 'GB-TIDE'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/zview/' + tag + '/' + str(my_year) + '/'

flist = sorted(glob.glob(outputs_dir+'*his*.nc'))
flist = flist[:7]

plt.switch_backend('Agg')

m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
            urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
            resolution='f')

# draw background
m.fillcontinents(color='lightgrey', alpha=0.5)
mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),labels=[0,0,0,1],fontsize=6, linewidth=.2)
pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),labels=[1,0,0,0],fontsize=6, linewidth=.2)

for var in varlist:
    print('For ' + var)
    if var in ['temp']:
        clim = [2, 10]
        cmap_var = cmocean.cm.thermal
    elif var in ['salt']:
        clim = [25, 35]
        cmap_var = cmocean.cm.haline
    elif var in ['zeta']:
        clim = [-3, 3]
        cmap_var = cmocean.cm.balance
    elif var in ['wetdry_mask_rho']:
        clim = [-1, 1]
        cmap_var = cmocean.cm.balance
    elif var in ['dye_03']:
        clim = [0, 0.1]
        cmap_var = cmocean.cm.matter

    if var in ['zeta', 'wetdry_mask_rho']:
        uvar = 'ubar'
        vvar = 'vbar'
    else:
        uvar = 'u'
        vvar = 'v'

    x, y = m(lon, lat)

    for fn in flist:

        tag = fn.split('/')[-1].split('.')[0]
        print('  processing ' + tag + ' ...')

        fh = nc.Dataset(fn)
        t = fh.variables['ocean_time'][:]
        tunit = (fh.variables['ocean_time']).units
        data = fh.variables[var][:]
        u = fh.variables[uvar][:]
        v = fh.variables[vvar][:]
        fh.close()

        for tindex in range(len(t)):
            ttag = nc.num2date(t[tindex], tunit).strftime("%Y-%m-%d_%H:%M:%S")

            if var not in ['zeta', 'wetdry_mask_rho']:
                # get z slice
                zslice, _, _ = pyroms.tools.zslice(data[tindex, :, :, :], depth, grd, \
                                          Cpos='rho', vert=True)
                # get u and v slice at requested depth
                zsliceu, _, _ = pyroms.tools.zslice(u[tindex, :, :, :], depth, grd, Cpos='u')
                zslicev, _, _ = pyroms.tools.zslice(v[tindex, :, :, :], depth, grd, Cpos='v')

            else:
                zslice = data[tindex, :, :]
                zsliceu = u[tindex, :, :]
                zslicev = v[tindex, :, :]

            zslice = np.ma.masked_where(mask == 0, zslice)
            # average field at rho point position
            zsliceu = 0.5 * (zsliceu[:,:-1] + zsliceu[:,1:])
            zsliceu = zsliceu[:,np.r_[0,:np.size(zsliceu,1),-1]]
            zsliceu = np.ma.masked_where(mask == 0, zsliceu)
            zsliceu = np.ma.masked_where(zsliceu >= 1000, zsliceu)
            zslicev = 0.5 * (zslicev[:-1,:] + zslicev[1:,:])
            zslicev = zslicev[np.r_[0,:np.size(zslicev,0),-1],:]
            zslicev = np.ma.masked_where(mask == 0, zslicev)
            zslicev = np.ma.masked_where(zslicev >= 1000, zslicev)
            U = zsliceu + 1j * zslicev
            # rotate velocity vector according to grid angle
            U = U * np.exp(1j * grd.hgrid.angle_rho)

            pc = m.pcolormesh(x, y, zslice, cmap=cmap_var)
            plt.clim(clim[0], clim[1])
            cb = m.colorbar()
            qv = m.quiver(x[::dd,::dd], y[::dd,::dd], \
                          np.real(U[::dd,::dd]), np.imag(U[::dd,::dd]), \
                          scale = uscale, linewidths=0.01)

            if var not in ['zeta', 'wetdry_mask_rho']:
                ttl = plt.title(var + '_' + str(int(-depth)) + 'm_' + ttag + '.png')
                plt.savefig(fig_dir + var + '_' + str(int(-depth)) + 'm_' + ttag + '.png')
            else:
                ttl = plt.title(var + '_' + ttag)
                plt.savefig(fig_dir + var + '_' + ttag + '.png')

            pc.remove()
            cb.remove()
            qv.remove()
            # ttl.remove()

plt.close()
