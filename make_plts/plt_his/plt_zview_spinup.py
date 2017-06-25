import numpy as np
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pyroms
import pyroms_toolbox as prt

outputs_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/outputs/spinup/'
fig_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/figs/zview/spinup/'
depth = 200
tindex = 0
# var = 'salt'
var = 'zeta'
uvar = 'u'
vvar = 'v'
# clim = [28, 32]
clim = [0, 0.15]

grd = pyroms.grid.get_ROMS_grid('GB')

flist = glob.glob(outputs_dir+'*.nc')
# flist = flist[-10:]

for fn in flist:
    tag = fn.split('/')[-1].split('.')[0]
    print 'processing ' + tag + ' ...'
    if var == 'zeta':
        prj = prt.twoDview(var, tindex, grd,
                           filename=fn, cmin=clim[0], cmax=clim[1], title=tag,
                           outfile = fig_dir + var + '_' + tag + '.png'
                          )

    else:
        prj = prt.zview(var, tindex, depth, grd,
                        filename=fn, cmin=clim[0], cmax=clim[1], title=tag
                       )
        prt.quiver(uvar, vvar, tindex, depth, grd,
                   filename = fn, proj=prj, d=10,
                   outfile = fig_dir + var + '_' + str(int(depth)) + 'm_' + tag + '.png'
                  )

    plt.close()
