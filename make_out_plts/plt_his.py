# plot ROMS history file

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import subprocess
import pyroms
from pyroms_toolbox import rx0

import glob

grd = pyroms.grid.get_ROMS_grid('GB')
msku = grd.hgrid.mask_u
mskv = grd.hgrid.mask_v

plth = 0
outpath = '/Volumes/R1/scratch/chuning/gb_roms/sync/'

h = grd.vgrid.h
h = np.ma.masked_where(grd.hgrid.mask==0, h)
hraw = grd.vgrid.hraw.squeeze()
hraw = np.ma.masked_where(grd.hgrid.mask==0, hraw)

inpath = '/Volumes/R1/scratch/chuning/gb_roms/outputs/01/his/'
flist = glob.glob(inpath + '*.nc')

# fh = nc.Dataset(inpath + flist[0])
# wdmsk = fh.variables['wetdry_mask_psi'][:].squeeze()
# fh.close()

rx = rx0(h, grd.hgrid.mask)

zlevel = 39  # surface
zlevel2 = 0  # bottom

xmin = 0
xmax = 1000
ymin = 0
ymax = 500

x = np.arange(h.shape[0])
y = np.arange(h.shape[1])
x2 = 0.5*(x[1:]+x[:-1])
y2 = 0.5*(y[1:]+y[:-1])

varlist = ['temp', 'salt']

for fl in flist:
    # fl = flist[0]
    fname = fl.split('/')[-1]

    for var in varlist:

        fh = nc.Dataset(fl)
        data = fh.variables[var][:, :, :, :].squeeze()
        fh.close()

        plt.figure()
        plt.pcolor(y, x, data[zlevel, :, :].squeeze())
        plt.xlim(ymin, ymax)
        plt.ylim(xmin, xmax)
        plt.colorbar()
        plt.title(var)
        plt.savefig(outpath + 'his/' + var + '_surf_' + fname + '_his' + '.png', format='png', dpi=900)
        plt.close()

        plt.figure()
        plt.pcolor(y, x, data[zlevel2, :, :].squeeze())
        plt.xlim(ymin, ymax)
        plt.ylim(xmin, xmax)
        plt.colorbar()
        plt.title(var)
        plt.savefig(outpath + 'his/' + var + '_bott_' + fname + '_his' + '.png', format='png', dpi=900)
        plt.close()

    # plot zeta & velocity vector
    fh = nc.Dataset(fl)
    zeta = fh.variables['zeta'][:, :, :].squeeze()
    ubar = fh.variables['ubar'][:, :, :].squeeze()
    vbar = fh.variables['vbar'][:, :, :].squeeze()
    fh.close()

    # mask land
    # ubar = np.ma.masked_where(np.tile(msku==0, (40, 1, 1)), ubar)
    # vbar = np.ma.masked_where(np.tile(mskv==0, (40, 1, 1)), vbar)
    ubar = np.ma.masked_where(msku==0, ubar)
    vbar = np.ma.masked_where(mskv==0, vbar)

    ubar = 0.5*(ubar[:-1, :] + ubar[1:, :])
    vbar = 0.5*(vbar[:, :-1] + vbar[:, 1:])

    plt.figure()
    plt.pcolor(y, x, zeta.squeeze())
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.colorbar()
    plt.quiver(y2[::10], x2[::10], ubar[::10, ::10], vbar[::10, ::10], 
               headwidth=1, headaxislength=2, color='k')
    plt.title('zeta')
    plt.savefig(outpath + 'his/' + 'zeta_' + fname + '_his_surf' + '.png', format='png', dpi=900)
    plt.close()


if plth==1:
    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], h[xmin:xmax, ymin:ymax], cmap='OrRd')
    # plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    # plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.clim(0, 500)
    # plt.clim(-5, 5)
    plt.colorbar()
    plt.title('h [m]')
    # plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + fname + '_h' + '.png', format='png', dpi=900)
    plt.close()

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], hraw[xmin:xmax, ymin:ymax], cmap='OrRd')
    # plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    # plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.clim(0, 500)
    # plt.clim(-5, 5)
    plt.colorbar()
    plt.title('hraw [m]')
    # plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + fname + '_hraw' + '.png', format='png', dpi=900)
    plt.close()

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], h[xmin:xmax, ymin:ymax]-hraw[xmin:xmax, ymin:ymax], cmap='OrRd')
    # plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    # plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.colorbar()
    plt.title('hdiff [m]')
    # plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + fname + '_hdiff' + '.png', format='png', dpi=900)
    plt.close()

    plt.figure()
    plt.pcolor(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], cmap='OrRd')
    # plt.plot([xmin, xmax], [yy, yy], '--k', lw=0.01)
    # plt.plot([xx, xx], [ymin, ymax], '--k', lw=0.01)
    plt.xlim(ymin, ymax)
    plt.ylim(xmin, xmax)
    plt.clim(0, 0.3)
    # plt.clim(0, 1)
    plt.colorbar()
    plt.title('rx0')
    # plt.contour(y[ymin:ymax], x[xmin:xmax], rx[xmin:xmax, ymin:ymax], levels=[0.35], colors='k', linewidths=0.1)
    plt.savefig(outpath + fname + '_rx0' + '.png', format='png', dpi=900)
    plt.close()

