import subprocess
import os
import numpy as np

import pyroms
import pyroms_toolbox

my_year = 2001
tag = 'JRA55v0.8'
vlist = ['q_10', 'rain', 'rlds', 'rsds', 'slp', 'snow', 't_10', 'u_10', 'v_10']

dst_grd = pyroms.grid.get_ROMS_grid('GB')

data_dir = '/Volumes/R1/scratch/chuning/data/jra/jra55do/'
data_dir_year = data_dir + str(my_year) + '/'
dst_dir='/Users/chuning/projects/gb_roms/data/'

# filelst = subprocess.check_output(['ls', data_dir_year]).replace('/n',' ').split()

# for filein in filelst:
for var in vlist:
    filein = var + '*.nc'
    fileout = var + '.' + str(my_year) + '.' + tag + '.nc'
    command1 = 'ncks -d latitude,255,270 -d longitude,390,410 ' \
                + data_dir_year + filein + ' -O ' + dst_dir + 'temp/' + fileout
    subprocess.call(command1, shell=True)

# merge file
atmos_file = dst_dir + dst_grd.name + '_atmos_' + str(my_year) + '_' + tag + '.nc'

kk = 0
for var in vlist:
    filein = dst_dir + 'temp/' + var + '*.nc'
    if kk==0:
        subprocess.call('mv ' + filein + ' ' + atmos_file, shell=True)
    else:
        subprocess.call('ncks -A ' + filein + ' -o ' + atmos_file ,shell=True)
    kk = kk+1

# clean up
subprocess.call('rm ' + dst_dir + 'temp/*', shell=True)

