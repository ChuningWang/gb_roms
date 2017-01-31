import numpy as np
import netCDF4 as nc
import datetime

vlist = ['sw', 'lw', 't10', 'pre', 'hair', 'u10', 'v10', 'rain']
vlist = ['sw']

for v in vlist:
    if v=='sw':
        vjra = 'rsds'
        vnoc = 'sw'
        cnoc = 'surface'
        vmer = 'swrad'
        vcor = 'swrad'
        vcor2 = 'swrad'
        vcor3 = 'srf_time'
    elif v=='lw':
        vjra = 'rlds'
        vnoc = 'lw'
        cnoc = 'surface'
        vmer = 'lwrad_down'
        vcor = 'lwrad'
        vcor2 = 'lwrad_down'
        vcor3 = 'lrf_time'
    elif v=='t10':
        vjra = 'tas_10m'
        vnoc = 'at'
        cnoc = 'surface'
        vmer = 'Tair'
        vcor = 'Tair'
        vcor2 = 'Tair'
        vcor3 = 'tair_time'
    elif v=='pre':
        vjra = 'psl'
        vnoc = 'slp'
        cnoc = 'surface'
        vmer = 'Pair'
        vcor = 'Pair'
        vcor2 = 'Pair'
        vcor3 = 'pair_time'
    elif v=='hair':
        vjra = 'huss_10m'
        vnoc = 'qair'
        cnoc = 'surface'
        vmer = 'Qair'
        vcor = 'Qair'
        vcor2 = 'Qair'
        vcor3 = 'qair_time'
    elif v=='u10':
        vjra = 'uas_10m'
        vnoc = 'wspd'
        cnoc = 'surface'
        vmer = 'Uwind'
        vcor = 'Uwind'
        vcor2 = 'Uwind'
        vcor3 = 'wind_time'
    elif v=='v10':
        vjra = 'vas_10m'
        vnoc = 'wspd'
        cnoc = 'surface'
        vmer = 'Vwind'
        vcor = 'Vwind'
        vcor2 = 'Vwind'
        vcor3 = 'wind_time'
    elif v=='rain':
        vjra = 'prrn'
        vnoc = 'cldc'
        cnoc = 'surface'
        vmer = 'rain'
        vcor = 'rain'
        vcor2 = 'rain'
        vcor3 = 'rain_time'

    # -------------------------------------------------------------------------------
    jra_dir = '/Users/chuning/projects/gb_roms/data/GlacierBay_atmos_2000_JRA55v0.8.nc'
    fh = nc.Dataset(jra_dir, 'r')
    t_jra = fh.variables['time'][:]
    lat_jra = fh.variables['latitude'][:]
    lon_jra = fh.variables['longitude'][:]
    if v=='t10':
        var_jra = fh.variables[vjra][:]-273.15
    else:
        var_jra = fh.variables[vjra][:]
    fh.close()

    noc_dir = '/Users/chuning/projects/gb_roms/data/noc/' + cnoc + '/nocs_v2_0_' + vnoc + '_2000.nc'
    fh = nc.Dataset(noc_dir, 'r')
    t_noc = fh.variables['time'][:]
    lat_noc = fh.variables['lat'][:]
    lon_noc = fh.variables['lon'][:]
    if v=='hair':
        var_noc = fh.variables[vnoc][:]/1000
    elif v=='pre':
        var_noc = fh.variables[vnoc][:]*100
    else:
        var_noc = fh.variables[vnoc][:]
    fh.close()

    t0 = datetime.datetime(1900, 1, 1)
    t1 = datetime.datetime(2000, 1, 1)
    dt = (t1-t0).days

    t_jra = t_jra-dt

    # -------------------------------------------------------------------------------
    # take a small region
    msk1 = (lat_jra>=53.5) & (lat_jra<=55.5)
    msk2 = (lon_jra>=219.5) & (lon_jra<=221.5)
    var_jra = var_jra[:, msk1, :][:, :, msk2]

    msk1 = (lat_noc>=53.5) & (lat_noc<=55.5)
    msk2 = (lon_noc>=219.5-360) & (lon_noc<=221.5-360)
    var_noc = var_noc[:, msk1, :][:, :, msk2]

    shp_jra = var_jra.shape
    shp_noc = var_noc.shape

    # -------------------------------------------------------------------------------
    var_jra = np.reshape(var_jra, (shp_jra[0], shp_jra[1]*shp_jra[2]))
    var_noc = np.reshape(var_noc, (shp_noc[0], shp_noc[1]*shp_noc[2]))

    var_jra = np.mean(var_jra, axis=1)
    var_noc = np.mean(var_noc, axis=1)

    # monthly bin average
    mm_jra = np.array([(t1+datetime.timedelta(i)).month for i in t_jra])
    t_jra_clim = np.zeros(12)
    var_jra_clim = np.zeros(12)
    for i in range(12):
        msk = mm_jra==i+1
        t_jra_clim[i] = np.mean(t_jra[msk])
        var_jra_clim[i] = np.mean(var_jra[msk])

    corr = var_noc/var_jra_clim
    # corr = np.concatenate((corr[-1], corr, corr[0]), axis=0) 
    corr = np.insert(corr, 0,  corr[-1])
    corr = np.append(corr, corr[1])
    t = np.insert(t_noc, 0, -14.5)
    t = np.append(t, 380.5)

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
