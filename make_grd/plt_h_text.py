import numpy as np

import pyroms
import pyroms_toolbox
import bathy_smoother

import matplotlib.pyplot as plt

x0 = 248
y0 = 61
dy = 10
dx = 10

grd = pyroms.grid.get_ROMS_grid('GB')
h0 = grd.vgrid.h
msk0 = grd.hgrid.mask_rho

# ------------------------------------------------
# shapiro filter
h0 = pyroms_toolbox.shapiro_filter.shapiro2(h0, 32)

h = h0[x0-dx:x0+dx, y0-dy:y0+dy]
msk = msk0[x0-dx:x0+dx, y0-dy:y0+dy]

h[10:16, 5:10][h[10:16, 5:10]<50] = 50
h[16:19, 6:10][h[16:19, 6:10]<50] = 50

rx0_max = 0.35
hs = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(msk, h, rx0_max)
# ------------------------------------------------
plt.pcolor(hs)
plt.colorbar
for i in range(2*dx):
    for j in range(2*dy):
        plt.text(j, i, "%03d" % hs[i, j])

plt.show(block=False)
