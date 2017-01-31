import subprocess
import os
import numpy as np
# import wget

url = 'http://amaterasu.ees.hokudai.ac.jp/~tsujino/JRA55_v1.1/'

my_year = 2000
tag = 'JRA55do'
# vlist = ['q_10', 'rain', 'rlds', 'rsds', 'slp', 'snow', 't_10', 'u_10', 'v_10']
vlist = ['rain', 'snow']
ttag1 = '30Jun2016'
ttag2 = '15Dec2016'
save_dir = '/Volumes/R1/scratch/chuning/data/jra/jra55do/' + str(my_year) + '/'

if not os.path.exists(save_dir):
    subprocess.call('mkdir ' + save_dir, shell=True)

cwd = os.getcwd()
os.chdir(save_dir)

for v in vlist:
    urlvar = url+v+'.'+str(my_year)+'.'+ttag1+'.nc'
    subprocess.call('wget -O '+ save_dir + v + '.' + str(my_year) + '.nc' + ' ' + urlvar, shell=True)
    urlvar = url+v+'.'+str(my_year)+'.'+ttag2+'.nc'
    subprocess.call('wget -O '+ save_dir + v + '.' + str(my_year) + '.nc' + ' ' + urlvar, shell=True)

os.chdir(cwd)

# http://amaterasu.ees.hokudai.ac.jp/~tsujino/JRA55_v0.8/rsds.2001.30Jun2016.nc
