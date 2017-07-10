import numpy as np

import pyroms
import pyroms_toolbox
import bathy_smoother

grd1 = 'GB3'
grd_name = 'GlacierBay_Kate'
tag = 'smooth'
out_file = '/glade/p/work/chuning/gb_roms/grd/' + grd_name + '_grd_' + tag + '.nc'

grd = pyroms.grid.get_ROMS_grid(grd1)
h0 = grd.vgrid.h
msk0 = grd.hgrid.mask_rho

# ------------------------------------------------
# mask out small channels
msk0[:230, :75] = 0
msk0[:150, :150] = 0
grd.hgrid.mask_rho = msk0

# ------------------------------------------------
# shapiro filter
h0 = pyroms_toolbox.shapiro_filter.shapiro2(h0, 32)

# ------------------------------------------------
# first blow-up point
x0 = 339
y0 = 112
dy = 28
dx = 28

h = h0[x0-dx:x0+dx, y0-dy:y0+dy]
msk = msk0[x0-dx:x0+dx, y0-dy:y0+dy]

# smooth bathymetry
rx0_max = 0.3
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, msk)
print 'Max Roughness value is: ', RoughMat.max()
hs = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(msk, h, rx0_max)
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hs, msk)
print 'Max Roughness value is: ', RoughMat.max()

grd.vgrid.h[x0-dx:x0+dx, y0-dy:y0+dy] = hs

# ------------------------------------------------
# # second blow-up point
# x0 = 248
# y0 = 61
# dy = 10
# dx = 15
# 
# h = h0[x0-dx:x0+dx, y0-dy:y0+dy]
# msk = msk0[x0-dx:x0+dx, y0-dy:y0+dy]
# 
# # smooth bathymetry
# rx0_max = 0.3
# RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, msk)
# print 'Max Roughness value is: ', RoughMat.max()
# hs = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(msk, h, rx0_max)
# RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(hs, msk)
# print 'Max Roughness value is: ', RoughMat.max()
# 
# grd.vgrid.h[x0-dx:x0+dx, y0-dy:y0+dy] = hs

# ------------------------------------------------
# # second blow-up point
# # fix point by point
# x0 = 248
# y0 = 61
# dy = 10
# dx = 15
# 
# h = h0[x0-dx:x0+dx, y0-dy:y0+dy]
# msk = msk0[x0-dx:x0+dx, y0-dy:y0+dy]
# 
# h[10:16, 5:10][h[10:16, 5:10]<50] = 50
# h[16:19, 6:10][h[16:19, 6:10]<50] = 50
# 
# grd.vgrid.h[x0-dx:x0+dx, y0-dy:y0+dy] = h

# apply filter again
rx0_max = 0.35
grd.vgrid.h = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(grd.hgrid.mask_rho, grd.vgrid.h, rx0_max)

pyroms.grid.write_ROMS_grid(grd, filename=out_file)
