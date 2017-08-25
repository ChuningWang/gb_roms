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
my_year = 2008
model = 'tmpdir_GB-TIDE/outputs/2008/'

varlist = ['temp', 'salt', 'zeta']

outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/zview/GB-TIDE/'

grd = pyroms.grid.get_ROMS_grid(grd1)

flist = sorted(glob.glob(outputs_dir+'*his*.nc'))
flist = flist[-30:]

depth = 5
tindex = 0
# var = 'temp'
# var = 'salt'
# var = 'zeta'

plt.switch_backend('Agg')

for var in varlist:
    if var=='temp':
        clim = [0, 8]
    elif var=='salt':
        clim = [28, 32]
    elif var=='zeta':
        clim = [-1, 1]

    if var=='zeta':
        uvar = 'ubar'
        vvar = 'vbar'
    else:
        uvar = 'u'
        vvar = 'v'

    for fn in flist:
        tag = fn.split('/')[-1].split('.')[0]
        print 'processing ' + tag + ' ...'
        if var == 'zeta':
            prj = prt.twoDview(var, tindex, grd,
                               filename=fn, cmin=clim[0], cmax=clim[1], title=tag
                              )
            prt.quiver2D(uvar, vvar, tindex, grd,
                         filename = fn, proj=prj, d=10, uscale=20,
                         outfile = fig_dir + str(my_year) + '/' + var + '_' + tag + '.png'
                        )

        else:
            prj = prt.zview(var, tindex, depth, grd,
                            filename=fn, cmin=clim[0], cmax=clim[1], title=tag
                           )
            prt.quiver(uvar, vvar, tindex, depth, grd,
                       filename = fn, proj=prj, d=10, uscale=20,
                       outfile = fig_dir + str(my_year) + '/' + var + '_' + str(int(depth)) + 'm_' + tag + '.png'
                      )

        plt.close()
