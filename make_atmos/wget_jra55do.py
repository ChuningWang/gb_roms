import subprocess
import os
import numpy as np
# import wget

url = 'http://amaterasu.ees.hokudai.ac.jp/~tsujino/JRA55-do-v1.1/'

my_year = 2000
tag = 'JRA55do'
vlist = ['q_10', 'rain', 'rlds', 'rsds', 'slp', 'snow', 't_10', 'u_10', 'v_10']
vlist2 = ['Qair', 'rain', 'lwrad_down', 'swrad', 'Pair', 'snow', 'Tair', 'Uwind', 'Vwind']
# vlist = ['rain', 'snow']
ttag1 = '30Jun2016'
ttag2 = '15Dec2016'
save_dir = '/Volumes/R1/scratch/chuning/gb_roms/data/jra55do/' + str(my_year) + '/'

if not os.path.exists(save_dir):
    subprocess.call('mkdir ' + save_dir, shell=True)

cwd = os.getcwd()
os.chdir(save_dir)

kk = 0
for v in vlist:
    if v in ['q_10', 'rlds', 'rsds', 'slp', 't_10', 'u_10', 'v_10']:
        urlvar = url+v+'.'+str(my_year)+'.'+ttag1+'.nc'
    else:
        urlvar = url+v+'.'+str(my_year)+'.'+ttag2+'.nc'
    subprocess.call('wget -O '+ save_dir + vlist2[kk] + '.' + str(my_year) + '.nc' + ' ' + urlvar, shell=True)
    kk = kk+1

os.chdir(cwd)

# http://amaterasu.ees.hokudai.ac.jp/~tsujino/JRA55_v0.8/rsds.2001.30Jun2016.nc
