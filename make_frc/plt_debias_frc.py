import numpy as np
import netCDF4 as nc
import datetime
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

lon0 = -136.0
lat0 = 58.5

vlist1 = ['at', 'lw', 'qair', 'slp', 'sw', 'wspd']
vlist2 = ['Tair', 'lwrad_down', 'Qair', 'Pair', 'swrad', 'wspd']
vlist3 = ['tair', 'lrf', 'qair', 'pair', 'srf', 'wind']
vlist1 = ['at']
vlist2 = ['Tair']
vlist3 = ['tair']

# -------------------------------------------------------------------------------
stn_file = in_dir+'station/juneau2000.csv'
from station import read_station
data = read_station(stn_file)
lat_stn = data['lat']
lon_stn = data['lon']
t_stn = data['ctime']
tmax_stn = data['tmax']
tmin_stn = data['tmin']

for i in range(len(vlist1)):
    v1 = vlist1[i]
    v2 = vlist2[i]
    v3 = vlist3[i]

    # -------------------------------------------------------------------------------
    noc_dir = in_dir+'noc/surface/nocs_v2_0_'+v1+'_2000.nc'
    fh = nc.Dataset(noc_dir, 'r')
    t_noc = fh.variables['time'][:]
    lat_noc = fh.variables['lat'][:]
    lon_noc = fh.variables['lon'][:]
    data_noc = fh.variables[v1][:]
    if v1=='qair':
        data_noc = fh.variables[v1][:]/1000
    elif v1=='slp':
        data_noc = fh.variables[v1][:]*100
    else:
        data_noc = fh.variables[v1][:]
    fh.close()

    if v2=='wspd':
        u_jra_dir = out_dir+'frc/Uwind_2000_JRA55v1.1.nc'
        v_jra_dir = out_dir+'frc/Vwind_2000_JRA55v1.1.nc'
        fh = nc.Dataset(u_jra_dir, 'r')
        t_jra = fh.variables[v3+'_time'][:]
        lat_jra = fh.variables['lat'][:]
        lon_jra = fh.variables['lon'][:]
        u_jra = fh.variables['Uwind'][:]
        fh.close()
        fh = nc.Dataset(v_jra_dir, 'r')
        v_jra = fh.variables['Vwind'][:]
        fh.close()
        data_jra = np.sqrt(u_jra**2+v_jra**2)
    else:
        jra_dir = out_dir+'frc/'+v2+'_2000_JRA55v1.1.nc'
        fh = nc.Dataset(jra_dir, 'r')
        t_jra = fh.variables[v3+'_time'][:]
        lat_jra = fh.variables['lat'][:]
        lon_jra = fh.variables['lon'][:]
        data_jra = fh.variables[v2][:]
        fh.close()

    if v2=='Tair':
        data_jra += -273.16

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
    data_jra = data_jra[:, msk1, :][:, :, msk2]

    msk1 = (lat_noc>=56.) & (lat_noc<=60.)
    msk2 = (lon_noc>=-140.) & (lon_noc<=-130.)
    lat_noc = lat_noc[msk1]
    lon_noc = lon_noc[msk2]
    data_noc = data_noc[:, msk1, :][:, :, msk2]

    # -------------------------------------------------------------------------------
    # interpolate to Glacier Bay
    data_jra = np.array([interp2d(lon_jra, lat_jra, data_jra[i, :, :])(lon0, lat0) for i in range(len(t_jra))])
    data_noc = np.array([interp2d(lon_noc, lat_noc, data_noc[i, :, :])(lon0, lat0) for i in range(len(t_noc))])

    # monthly bin average
    mm_jra = np.array([nc.num2date(i, 'days since 1900-01-01 00:00:00').month for i in t_jra])
    t_jra_clim = np.zeros(12)
    data_jra_clim = np.zeros(12)
    for i in range(12):
        msk = mm_jra==i+1
        t_jra_clim[i] = np.mean(t_jra[msk])
        data_jra_clim[i] = np.mean(data_jra[msk])

    plt.plot(t_jra_clim, data_jra_clim, 'r', lw=3)
    plt.plot(t_noc, data_noc, 'k', lw=3)
    plt.plot(t_jra, data_jra, '--r', lw=0.1)
    if v1=='wspd':
        plt.plot(t_stn, data['awnd'], lw=0.1)
    elif v1=='at':
        plt.plot(t_stn, data['tair'], lw=0.1)
    plt.legend(['jra', 'noc'])
    plt.savefig(out_dir+'figs/debias/'+v1+'.tiff', format='tiff')
    plt.close()

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
