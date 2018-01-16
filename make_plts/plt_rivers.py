import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

if grd1=='GB_hr':
    tag = 'GlacierBay_hr'
elif grd1=='GB_lr':
    tag = 'GlacierBay_lr'

dtype = 3

grd = pyroms.grid.get_ROMS_grid(grd1)
lon_grd = grd.hgrid.lon_rho
lat_grd = grd.hgrid.lat_rho
hraw = grd.vgrid.h
msk = grd.hgrid.mask_rho

if dtype==1:
    # load remapped data
    fh = nc.Dataset(data_dir+tag+'_runoff_' + str(rspread) + '_2008_Hill.nc', 'r')
    time = fh.variables['time'][:]
    lon = fh.variables['lon_rho'][:]
    lat = fh.variables['lat_rho'][:]
    r = fh.variables['Runoff'][:]
    fh.close
elif dtype==2:
    # load original data
    fh = nc.Dataset(data_dir+'gb_discharge.nc', 'r')
    time = fh.variables['t'][:]
    lon = fh.variables['lon'][:]
    lat = fh.variables['lat'][:]
    r = fh.variables['discharge'][:]
    coast = fh.variables['coast'][:]
    tmsk = (time>=39446) & (time <=39812)
    time = time[tmsk]
    r = r[tmsk, :, :]
    for i in range(len(time)):
        r[i, :, :] = np.ma.masked_where(coast.mask, r[i, :, :])
    fh.close()
elif dtype==3:
    # load original data
    fh = nc.Dataset(out_dir + 'frc/' + tag + '_rivers_2008_Hill_ana.nc', 'r')
    time = fh.variables['river_time'][:]
    eta = fh.variables['river_Eposition'][:]
    xi = fh.variables['river_Xposition'][:].astype(int)
    trs = fh.variables['river_transport'][:].astype(int)
    rdir = fh.variables['river_direction'][:]
    sign = fh.variables['river_sign'][:]
    fh.close()
    lon = lon_grd
    lat = lat_grd
    r = np.NaN*np.zeros((len(time), lon.shape[0], lon.shape[1]))
    h = np.NaN*np.zeros(lon.shape)
    for i in range(len(time)):
        for j in range(len(eta)):
            if sign[j] == 1:
                r[i, eta[j], xi[j]] = trs[i, j]
                h[eta[j], xi[j]] = hraw[eta[j], xi[j]]
            else:
                if rdir[j] == 0:
                    r[i, eta[j]-1, xi[j]] = trs[i, j]
                    h[eta[j]-1, xi[j]] = hraw[eta[j]-1, xi[j]]
                elif rdir[j] == 1:
                    r[i, eta[j], xi[j]-1] = trs[i, j]
                    h[eta[j], xi[j]-1] = hraw[eta[j], xi[j]-1]

    r = np.ma.masked_invalid(r)
    h = np.ma.masked_invalid(h)
    plt.pcolormesh(lon, lat, h, cmap='Greens')
    plt.savefig(out_dir+'figs/rivers/h.png')
    plt.close()

tp, xp, yp = r.shape

fig, ax1 = plt.subplots()

for i, t in enumerate(time):
    data = r[i, :, :]
    # data = np.ma.masked_where(data, msk == 1)
    ttag = nc.num2date(t, 'days since 1900-01-01').strftime("%Y-%m-%d_%H:%M:%S")
    pcm = plt.pcolor(lon, lat, abs(data), cmap='Greens', edgecolors='k', linewidth=0.005)
    # pcm = plt.pcolor(abs(data), cmap='Greens', edgecolors='k', linewidth=0.005)
    plt.clim(0, 10)

    if i == 0:
        plt.xlim(-137.5, -135)
        plt.ylim(58., 59.25)
        plt.colorbar()

    plt.savefig(out_dir+'figs/rivers/runoff_proj_' + ttag + '.png', dpi=300)
    pcm.remove()

plt.close()

# def plt_rivers(t, r, out_dir):
#     plt.pcolormesh(lon, lat, r, cmap='Greens')
#     plt.xlim(-137.5, -135)
#     plt.ylim(58., 59.25)
#     plt.clim(0, 100)
#     plt.colorbar()
#     plt.contour(lon_grd, lat_grd, msk, np.array([0.5, 0.5]), linewidths=0.05, colors='k')
#     plt.savefig(out_dir+'figs/rivers/runoff_'+str(int(t))+'.png')
#     plt.close()
#     return None

# from joblib import Parallel, delayed
# Parallel(n_jobs=16)(delayed(plt_rivers)(time[i], r[i, :, :], out_dir) for i in range(len(time)))
# Parallel(n_jobs=16)(delayed(plt_rivers)(time[i], r[i, :, :], out_dir) for i in range(150, 165))
