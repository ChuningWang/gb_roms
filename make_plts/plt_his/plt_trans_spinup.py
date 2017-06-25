import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms
import glob
from gb_toolbox import gb_ctd

# load data
outputs_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/outputs/spinup/'
fig_dir = '/Volumes/R1/scratch/chuning/gb_spinup_roms/figs/trans/spinup/'
depth = 200
tindex = 0
var = 'salt'
uvar = 'u'
vvar = 'v'
clim = [28, 32]

# load geological coordinates
ctd = gb_ctd.rd_ctd('/Volumes/R1/scratch/chuning/gb_roms/data/ctd/ctd.nc')
lat_ctd = ctd['lat']
lon_ctd = ctd['lon']

station = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 21]
lat_ctd = lat_ctd[station]
lon_ctd = lon_ctd[station]

lat_t = np.zeros((lat_ctd.shape[0]-1)*10)
lon_t = np.zeros((lon_ctd.shape[0]-1)*10)

for i in range(lat_ctd.shape[0]-1):
    # lat_t = np.linspace(lat_ctd)

grd = pyroms.grid.get_ROMS_grid('GB')
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho


flist = glob.glob(outputs_dir+'*.nc')
flist = flist[-1:]

for fn in flist:
    tag = fn.split('/')[-1].split('.')[0]
    print 'processing ' + tag + ' ...'
    fh = nc.Dataset(fn)
    s = fh.variables['salt'][:].squeeze()
