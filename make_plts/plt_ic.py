import numpy as np
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pyroms
import pyroms_toolbox as prt

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
grd1 = 'GB_USGS'
grd = pyroms.grid.get_ROMS_grid(grd1)

in_file = out_dir + 'bc_ic/' + grd.name + '_ic_2000_01_03_SODA3.3.1.nc'

depth = 5
tindex = 0
var = 'salt'
# var = 'zeta'
uvar = 'u'
vvar = 'v'
clim = [28, 32]
# clim = [0, 0.15]

print 'processing ' + in_file + ' ...'
if var == 'zeta':
    prj = prt.twoDview(var, tindex, grd,
                       filename=in_file, title='IC_' + var, cmin=clim[0], cmax=clim[1],
                       outfile = out_dir + 'figs/' + grd.name + '_' + var + '.png'
                      )

else:
    prj = prt.zview(var, tindex, depth, grd,
                    filename=in_file, title='IC_' + var, cmin=clim[0], cmax=clim[1]
                   )
    prt.quiver(uvar, vvar, tindex, depth, grd,
               filename = in_file, proj=prj, d=10, uscale=10,
               outfile = out_dir + 'figs/' + grd.name + '_' + var + '.png'
              )

plt.close()
