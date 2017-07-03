# plot ROMS history file

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import subprocess
import pyroms
from pyroms_toolbox import rx0

import glob

grd = pyroms.grid.get_ROMS_grid('GB2')
msku = grd.hgrid.mask_u
mskv = grd.hgrid.mask_v

ex = 'abs'
rmap = 'full'
plth = 1
outpath = '/Volumes/R1/scratch/chuning/gb_spinup_roms/figs/diagnostic/'

h = grd.vgrid.h
h = np.ma.masked_where(grd.hgrid.mask==0, h)
hraw = grd.vgrid.hraw.squeeze()
hraw = np.ma.masked_where(grd.hgrid.mask==0, hraw)

inpath = '/Volumes/R1/scratch/chuning/'
flist = glob.glob(inpath + '*.nc')

rx = rx0(h, grd.hgrid.mask)

tag = 'GB-TIDE'

varlist = ['u', 'v', 'temp', 'salt']
for var in varlist:
    # var = 'u'

    fh = nc.Dataset(inpath + tag + '_rst.nc')

    data = fh.variables[var][:, 0, :, :, :].squeeze()
    if var=='u':
        data = np.ma.masked_where(np.tile(msku==0, (40, 1, 1)), data)
    elif var=='v':
        data = np.ma.masked_where(np.tile(mskv==0, (40, 1, 1)), data)

    fh.close()

    x = np.arange(h.shape[0])
    y = np.arange(h.shape[1])

    if ex=='abs':
        idx = np.where(np.abs(data) == np.nanmax(np.abs(data)))
    elif ex=='max':
        idx = np.where(data == np.nanmax(data))
    elif ex=='min':
        idx = np.where(data == np.nanmin(data))
    zlevel = idx[0].data[0]
    # zlevel = 40  # surface
    yy = idx[1].data[0]
    xx = idx[2].data[0]

    if rmap=='full':
        xmin = 0
        xmax = 1000
        ymin = 0
        ymax = 500
    elif rmap=='zoom':
        xmin = xx-20
        xmax = xx+20
        ymin = yy-20
        ymax = yy+20

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], data[zlevel, xmin:xmax, ymin:ymax].squeeze())
    plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.colorbar()
    plt.title(str(xx)+ ',' + str(yy) + ',' + str(zlevel))
    plt.savefig(outpath + tag + '_rst_' + var + '_' + ex + '.png', format='png', dpi=900)
    plt.close()

if plth==1:
    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], h[xmin:xmax, ymin:ymax], cmap='OrRd')
    plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.clim(0, 500)
    # plt.clim(-5, 5)
    plt.colorbar()
    plt.title('h [m]')
    plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + tag + '_h_' + '_' + ex + '.png', format='png', dpi=900)
    plt.close()

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], hraw[xmin:xmax, ymin:ymax], cmap='OrRd')
    plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.clim(0, 500)
    # plt.clim(-5, 5)
    plt.colorbar()
    plt.title('hraw [m]')
    plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + tag + '_hraw_' + '_' + ex + '.png', format='png', dpi=900)
    plt.close()

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], h[xmin:xmax, ymin:ymax]-hraw[xmin:xmax, ymin:ymax], cmap='OrRd')
    plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.colorbar()
    plt.title('hdiff [m]')
    plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + tag + '_hdiff_' + '_' + ex + '.png', format='png', dpi=900)
    plt.close()

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], cmap='OrRd')
    plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    # plt.clim(0, 0.3)
    plt.clim(0, 1)
    plt.colorbar()
    plt.title('rx0')
    plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + tag + '_rx0_' + ex + '.png', format='png', dpi=900)
    plt.close()

