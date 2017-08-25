# plot ROMS history file
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms
from pyroms_toolbox import rx0

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

grd1 = 'GB_lr'

ex = 'abs'
rmap = 'full'
outpath = out_dir + 'figs/rst/'

grd = pyroms.grid.get_ROMS_grid(grd1)
msk = grd.hgrid.mask_rho
msku = grd.hgrid.mask_u
mskv = grd.hgrid.mask_v

h = grd.vgrid.h
h = np.ma.masked_where(grd.hgrid.mask==0, h)
hraw = grd.vgrid.hraw.squeeze()
hraw = np.ma.masked_where(grd.hgrid.mask==0, hraw)

rx = rx0(h, grd.hgrid.mask)

tag = 'GB-CIRC'

varlist = ['u', 'v', 'temp', 'salt', 'zeta', 'wetdry_mask_rho']
# varlist = ['wetdry_mask_rho']
for var in varlist:
    # var = 'u'

    fh = nc.Dataset(model_dir + 'tmpdir_GB-CIRC/' + tag + '_rst.nc')
    ocean_time = fh.variables['ocean_time'][:]
    tindex = np.argmax(ocean_time)

    if var=='zeta':
        data = fh.variables[var][tindex, 0, :, :].squeeze()
    elif var=='wetdry_mask_rho':
        data = fh.variables[var][tindex, :, :].squeeze()
    else:
        data = fh.variables[var][tindex, 0, :, :, :].squeeze()

    if var=='u':
        data = np.ma.masked_where(np.tile(msku==0, (grd.vgrid.N, 1, 1)), data)
    elif var=='v':
        data = np.ma.masked_where(np.tile(mskv==0, (grd.vgrid.N, 1, 1)), data)
    elif var=='wetdry_mask_rho':
        data = np.ma.masked_where(msk==0, data)

    fh.close()

    x = np.arange(h.shape[0])
    y = np.arange(h.shape[1])

    if ex=='abs':
        idx = np.where(np.abs(data) == np.nanmax(np.abs(data)))
    elif ex=='max':
        idx = np.where(data == np.nanmax(data))
    elif ex=='min':
        idx = np.where(data == np.nanmin(data))

    if var=='zeta':
        zlevel = 0
        xx = idx[0][0]
        yy = idx[1][0]
    elif var=='wetdry_mask_rho':
        zlevel = 0
        xx = 0
        yy = 0
    else:
        zlevel = idx[0].data[0]
        data = data[zlevel, :, :]
        # zlevel = 40  # surface
        yy = idx[1].data[0]
        xx = idx[2].data[0]

    if rmap=='full':
        xmin = 0
        xmax = 500
        ymin = 0
        ymax = 250
    elif rmap=='zoom':
        xmin = xx-20
        xmax = xx+20
        ymin = yy-20
        ymax = yy+20

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], data[xmin:xmax, ymin:ymax].squeeze())
    plt.plot([ymin, ymax], [yy, yy], '--k', lw=0.01)
    plt.plot([xx, xx], [xmin, xmax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    if var=='temp':
        plt.clim(0, 20)
    elif var=='salt':
        plt.clim(0, 35)
    plt.colorbar()
    plt.title(str(xx)+ ',' + str(yy) + ',' + str(zlevel))
    plt.savefig(outpath + tag + '_rst_' + var + '_' + ex + '.png', format='png', dpi=900)
    plt.close()

    if var=='zeta':
        plt.figure()
        plt.pcolor(y[ymin:ymax], x[xmin:xmax], h[xmin:xmax, ymin:ymax]-data[xmin:xmax, ymin:ymax].squeeze())
        plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
        plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
        plt.xlim(ymin, ymax)
        plt.ylim(xmin, xmax)
        plt.clim(-5, 5)
        plt.colorbar()
        plt.title(str(xx)+ ',' + str(yy) + ',' + str(zlevel))
        plt.savefig(outpath + tag + '_rst_h-' + var + '_' + ex + '.png', format='png', dpi=900)
        plt.close()
