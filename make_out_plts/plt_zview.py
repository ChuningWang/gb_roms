import numpy as np
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pyroms
import pyroms_toolbox as prt
import netCDF4 as nc
import sys

import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# my inputs
my_year = 2008
varlist = ['zeta', 'temp', 'salt', 'dye_03']
varlist = ['temp']
depth = 1
# tindex = 0

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

grd = pyroms.grid.get_ROMS_grid(grd1)

if grd1=='GB_lr':
    tag = 'GB-CIRC'
if grd1=='GB_hr':
    tag = 'GB-TIDE'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/zview/' + tag + '/' + str(my_year) + '/'

flist = sorted(glob.glob(outputs_dir+'*his*.nc'))
flist = flist[-30:]

# var = 'temp'
# var = 'salt'
# var = 'zeta'

plt.switch_backend('Agg')

for var in varlist:
    if var=='temp':
        clim = [5, 10]
    elif var=='salt':
        clim = [30, 35]
    elif var=='zeta':
        clim = [-3, 3]
    elif var=='dye_03':
        clim = [0, 1]

    if var=='zeta':
        uvar = 'ubar'
        vvar = 'vbar'
    else:
        uvar = 'u'
        vvar = 'v'

    for fn in flist:
        tag = fn.split('/')[-1].split('.')[0]
        print 'processing ' + tag + ' ...'
        fh = nc.Dataset(fn)
        t = fh.variables['ocean_time'][:]
        tunit = (fh.variables['ocean_time']).units
        fh.close()
        for tindex in range(len(t)):
            ttag = nc.num2date(t[tindex], tunit).strftime("%Y-%m-%d_%H:%M:%S")
            if var == 'zeta':
                prj = prt.twoDview(var, tindex, grd,
                                   filename=fn, cmin=clim[0], cmax=clim[1], title=tag
                                  )
                prt.quiver2D(uvar, vvar, tindex, grd,
                             filename = fn, proj=prj, d=10, uscale=20,
                             outfile = fig_dir + var + '_' + ttag + '.png'
                            )

            else:
                prj = prt.zview(var, tindex, depth, grd,
                                filename=fn, cmin=clim[0], cmax=clim[1], title=tag
                               )
                prt.quiver(uvar, vvar, tindex, depth, grd,
                           filename = fn, proj=prj, d=10, uscale=20,
                           outfile = fig_dir + var + '_' + str(int(depth)) + 'm_' + ttag + '.png'
                          )

            plt.close()
