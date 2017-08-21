import numpy as np
import netCDF4 as nc
import datetime
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

# -------------------------------------------------------------------------------
stn_file = in_dir+'station/juneau2000.csv'
from station import read_station
data = read_station(stn_file)
lat_stn = data['lat']
lon_stn = data['lon']
t_stn = data['ctime']
tmax_stn = data['tmax']
tmin_stn = data['tmin']

# -------------------------------------------------------------------------------
jra_dir = out_dir+'frc/Tair_2000_JRA55v1.1.nc'
fh = nc.Dataset(jra_dir, 'r')
t_jra = fh.variables['tair_time'][:]
lat_jra = fh.variables['lat'][:]
lon_jra = fh.variables['lon'][:]
tair_jra = fh.variables['Tair'][:]-273.16
fh.close()

noc_dir = in_dir+'noc/surface/nocs_v2_0_at_2000.nc'
fh = nc.Dataset(noc_dir, 'r')
t_noc = fh.variables['time'][:]
lat_noc = fh.variables['lat'][:]
lon_noc = fh.variables['lon'][:]
tair_noc = fh.variables['at'][:]
# if v=='hair':
#     var_noc = fh.variables[vnoc][:]/1000
# elif v=='pre':
#     var_noc = fh.variables[vnoc][:]*100
# else:
#     var_noc = fh.variables[vnoc][:]
fh.close()

t0 = datetime.datetime(1900, 1, 1)
t1 = datetime.datetime(2000, 1, 1)
dt = (t1-t0).days

t_noc = t_noc+dt

# -------------------------------------------------------------------------------
# take a small region
msk1 = (lat_jra>=56.) & (lat_jra<=60.)
msk2 = (lon_jra>=-140.) & (lon_jra<=-130.)
lat_jra = lat_jra[msk1]
lon_jra = lon_jra[msk2]
tair_jra = tair_jra[:, msk1, :][:, :, msk2]

msk1 = (lat_noc>=56.) & (lat_noc<=60.)
msk2 = (lon_noc>=-140.) & (lon_noc<=-130.)
lat_noc = lat_noc[msk1]
lon_noc = lon_noc[msk2]
tair_noc = tair_noc[:, msk1, :][:, :, msk2]

# -------------------------------------------------------------------------------
tair_jra = np.array([interp2d(lon_jra, lat_jra, tair_jra[i, :, :])(lon_stn, lat_stn) for i in range(len(t_jra))])
tair_noc = np.array([interp2d(lon_noc, lat_noc, tair_noc[i, :, :])(lon_stn, lat_stn) for i in range(len(t_noc))])

# monthly bin average
mm_jra = np.array([nc.num2date(i, 'days since 1900-01-01 00:00:00').month for i in t_jra])
# t_jra_clim = np.zeros(12)
# var_jra_clim = np.zeros(12)
# for i in range(12):
#     msk = mm_jra==i+1
#     t_jra_clim[i] = np.mean(t_jra[msk])
#     var_jra_clim[i] = np.mean(var_jra[msk])

# plt.plot(t_jra, tair_jra)
# plt.plot(t_noc, tair_noc, '.k')
# plt.plot(t_stn, tmax_stn, 'r')
# plt.plot(t_stn, tmin_stn, 'r')

# shp_jra = var_jra.shape
# shp_noc = var_noc.shape
# 
# # -------------------------------------------------------------------------------
# var_jra = np.reshape(var_jra, (shp_jra[0], shp_jra[1]*shp_jra[2]))
# var_noc = np.reshape(var_noc, (shp_noc[0], shp_noc[1]*shp_noc[2]))
# 
# var_jra = np.mean(var_jra, axis=1)
# var_noc = np.mean(var_noc, axis=1)
# 
# # monthly bin average
# mm_jra = np.array([(t1+datetime.timedelta(i)).month for i in t_jra])
# t_jra_clim = np.zeros(12)
# var_jra_clim = np.zeros(12)
# for i in range(12):
#     msk = mm_jra==i+1
#     t_jra_clim[i] = np.mean(t_jra[msk])
#     var_jra_clim[i] = np.mean(var_jra[msk])
# 
# corr = var_noc/var_jra_clim
# # corr = np.concatenate((corr[-1], corr, corr[0]), axis=0) 
# corr = np.insert(corr, 0,  corr[-1])
# corr = np.append(corr, corr[1])
# t = np.insert(t_noc, 0, -14.5)
# t = np.append(t, 380.5)
# 
# -------------------------------------------------------------------------------
# jra_dir = '/Users/chuning/projects/gb_roms/data/GlacierBay_atmos_2000_JRA55v0.8.nc'
# fh = nc.Dataset(jra_dir, 'r')
# t_jra = fh.variables['time'][:]
# lat_jra = fh.variables['latitude'][:]
# lon_jra = fh.variables['longitude'][:]
# if v=='t10':
#     var_jra = fh.variables[vjra][:]-273.15
# else:
#     var_jra = fh.variables[vjra][:]
# fh.close()

# noc_dir = '/Users/chuning/projects/gb_roms/data/noc/' + cnoc + '/nocs_v2_0_' + vnoc + '_2000.nc'
# fh = nc.Dataset(noc_dir, 'r')
# t_noc = fh.variables['time'][:]
# lat_noc = fh.variables['lat'][:]
# lon_noc = fh.variables['lon'][:]
# if v=='hair':
#     var_noc = fh.variables[vnoc][:]/1000
# elif v=='pre':
#     var_noc = fh.variables[vnoc][:]*100
# else:
#     var_noc = fh.variables[vnoc][:]
# fh.close()

# mer_dir = '/Volumes/P1/Data/MERRA/drowned/drowned_MERRA_' + vmer + '_3hours_2000.nc'
# fh = nc.Dataset(mer_dir, 'r')
# t_mer = fh.variables['time'][:]
# lat_mer = fh.variables['lat'][:]
# lon_mer = fh.variables['lon'][:]
# var_mer = fh.variables[vmer][:]
# fh.close()

# cor_dir = '/Volumes/P1/Data/CORE/CORE2/drowned/drowned_' + vcor + '.1948-2007.nc'
# fh = nc.Dataset(cor_dir, 'r')
# t_cor = fh.variables[vcor3][:]
# lat_cor = fh.variables['lat'][:]
# lon_cor = fh.variables['lon'][:]
# if v=='pre':
#     var_cor = fh.variables[vcor2][:]
# else:
#     var_cor = fh.variables[vcor2][:]
# fh.close()

# var_mer_cp = var_mer.copy()

# t0 = datetime.datetime(1900, 1, 1)
# t1 = datetime.datetime(2000, 1, 1)
# dt = (t1-t0).days

# t_jra = t_jra-dt
# t_mer = t_mer-dt
# t_cor = t_cor-dt

# # take CORE data in 2000
# msk = (t_cor>=0) & (t_cor<=366)
# t_cor = t_cor[msk]
# var_cor = var_cor[msk, :, :]

# -------------------------------------------------------------------------------
# # take a small region
# msk1 = (lat_jra>=53.5) & (lat_jra<=55.5)
# msk2 = (lon_jra>=219.5) & (lon_jra<=221.5)
# var_jra = var_jra[:, msk1, :][:, :, msk2]

# msk1 = (lat_noc>=53.5) & (lat_noc<=55.5)
# msk2 = (lon_noc>=219.5-360) & (lon_noc<=221.5-360)
# var_noc = var_noc[:, msk1, :][:, :, msk2]

# msk1 = (lat_mer>=53.5) & (lat_mer<=55.5)
# msk2 = (lon_mer>=219.5) & (lon_mer<=221.5)
# var_mer = var_mer[:, msk1, :][:, :, msk2]

# msk1 = (lat_cor>=53.5) & (lat_cor<=55.5)
# msk2 = (lon_cor>=219.5) & (lon_cor<=221.5)
# var_cor = var_cor[:, msk1, :][:, :, msk2]

# shp_jra = var_jra.shape
# shp_noc = var_noc.shape
# shp_mer = var_mer.shape
# shp_cor = var_cor.shape

# var_jra = np.reshape(var_jra, (shp_jra[0], shp_jra[1]*shp_jra[2]))
# var_noc = np.reshape(var_noc, (shp_noc[0], shp_noc[1]*shp_noc[2]))
# var_mer = np.reshape(var_mer, (shp_mer[0], shp_mer[1]*shp_mer[2]))
# var_cor = np.reshape(var_cor, (shp_cor[0], shp_cor[1]*shp_cor[2]))

# var_jra = np.mean(var_jra, axis=1)
# var_noc = np.mean(var_noc, axis=1)
# var_mer = np.mean(var_mer, axis=1)
# var_cor = np.mean(var_cor, axis=1)

# # monthly bin average
# mm_jra = np.array([(t1+datetime.timedelta(i)).month for i in t_jra])
# t_jra_clim = np.zeros(12)
# var_jra_clim = np.zeros(12)
# for i in range(12):
#     msk = mm_jra==i+1
#     t_jra_clim[i] = np.mean(t_jra[msk])
#     var_jra_clim[i] = np.mean(var_jra[msk])

# mm_mer = np.array([(t1+datetime.timedelta(i)).month for i in t_mer])
# t_mer_clim = np.zeros(12)
# var_mer_clim = np.zeros(12)
# for i in range(12):
#     msk = mm_mer==i+1
#     t_mer_clim[i] = np.mean(t_mer[msk])
#     var_mer_clim[i] = np.mean(var_mer[msk])

# mm_cor = np.array([(t1+datetime.timedelta(i)).month for i in t_cor])
# t_cor_clim = np.zeros(12)
# var_cor_clim = np.zeros(12)
# for i in range(12):
#     msk = mm_cor==i+1
#     t_cor_clim[i] = np.mean(t_cor[msk])
#     var_cor_clim[i] = np.mean(var_cor[msk])


# pltfig = 1
# if pltfig==1:
#     import matplotlib.pyplot as plt
#     plt.figure()
#     plt.plot(t_jra_clim, var_jra_clim, 'r')
#     plt.plot(t_mer_clim, var_mer_clim, 'k')
#     plt.plot(t_cor_clim, var_cor_clim, 'g')
#     plt.plot(t_noc, var_noc, 'b')
#     # plt.show(block=False)
#     plt.savefig('/Users/chuning/projects/gb_roms/figs/' + v + '_compare.tiff', format='tiff')
#     plt.close()
