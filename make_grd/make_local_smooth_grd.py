import numpy as np

import pyroms
import pyroms_toolbox
import bathy_smoother

grd = pyroms.grid.get_ROMS_grid('GB')
h = grd.vgrid.h
msk = grd.hgrd.mask_rho

smth = 10

# blowup coord 1 (443, 271, 16)

h2 = h[443-smth:443+smth, 271-smth:271+smth]
msk2 = msk[443-smth:443+smth, 271-smth:271+smth]

# smooth bathymetry
rx0_max = 0.3
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()
h2 = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h2, rx0_max)
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print 'Max Roughness value is: ', RoughMat.max()
