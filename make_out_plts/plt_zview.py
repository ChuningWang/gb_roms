import numpy as np
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pyroms
import pyroms_toolbox as prt

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

grd1 = 'GB_USGS'
model = 'tmpdir_GB-TIDE/outputs/2000/'
# model = 'tmpdir_GB-TIDE/'

outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/zview/2000/'
depth = 150
tindex = 0
var = 'salt'
# var = 'zeta'
uvar = 'u'
vvar = 'v'
clim = [28, 32]
# clim = [0, 0.15]

grd = pyroms.grid.get_ROMS_grid(grd1)

flist = glob.glob(outputs_dir+'*his*.nc')
flist = flist[-48:]

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
                   filename = fn, proj=prj, d=10, uscale=20,
                   outfile = fig_dir + var + '_' + str(int(depth)) + 'm_' + tag + '.png'
                  )

    plt.close()
