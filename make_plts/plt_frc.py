import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from gb_toolbox import buoy
import pandas

fh = nc.Dataset('../../data/GlacierBay_atmos_2000_JRA55v0.8.nc')

lat = fh.variables['latitude'][:]
lon = fh.variables['longitude'][:]
t = fh.variables['time'][:]

rl = fh.variables['rlds'][:]
rs = fh.variables['rsds'][:]
hu10 = fh.variables['huss_10m'][:]
prcp = fh.variables['prrn'][:]
pr = fh.variables['psl'][:]
t10 = fh.variables['tas_10m'][:]
u10 = fh.variables['uas_10m'][:]
v10 = fh.variables['vas_10m'][:]

fh.close()

noc_dir = '/Users/chuning/projects/gb_roms/data/noc/surface/nocs_v2_0_lw_2000.nc'
fh = nc.Dataset(noc_dir, 'r')
t_noc = fh.variables['time'][:]
lat_noc = fh.variables['lat'][:]
lon_noc = fh.variables['lon'][:]
var_noc = fh.variables['lw'][:]
fh.close()


jundir = '/Users/chuning/projects/gb_roms/data/cdo_stn/juneau.csv'
jun = pandas.read_csv(jundir).values

# make plots
pltfig = 1
if pltfig==1:
    for i in range(12):
        plt.figure()
        plt.pcolor(lon_noc, lat_noc, var_noc[i, :, :])
        plt.colorbar()
        plt.savefig('../../figs/noc_sw' + str(i+1) + '.tiff', format='tiff')
        plt.close()

    # plt.figure()
    # plt.plot(t, rs[:, 8, 7])
    # plt.show(block=False)
