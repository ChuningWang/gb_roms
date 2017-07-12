# Plot Glacier Bay map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from gb_toolbox.gb_ctd import rd_ctd
from matplotlib.mlab import griddata
import numpy as np
import netCDF4 as nc

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# ------------------------------------------------------------------------------------------------------
# clevs_land = 4*10**np.linspace(0, 3, 101)
# clevs_land = np.linspace(0, 2000, 11)
# clevs_land = []
# clevs_sea = np.linspace(-200, 0, 11)
# clevs_sea = []
# clevs = np.concatenate((clevs_sea, clevs_land))
clevs = np.linspace(-200, 200, 21)
# ------------------------------------------------------------------------------------------------------
# ctd = rd_ctd('../data/ctd.nc')

# Read topography
fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/ARDEMv2.0.nc', mode='r')
lont = fh.variables['lon'][:]
lont = lont-360
latt = fh.variables['lat'][:]
z = fh.variables['z'][:]
fh.close()

# stn = np.arange(24)

# Read grid
fh = nc.Dataset('/Volumes/R1/scratch/chuning/gb_roms/data/roms_prep/GlacierBay_usgs_bdry_2000_SODA3.3.1.nc', mode='r')
lon_psi = fh.variables['lon_psi'][:]
lat_psi = fh.variables['lat_psi'][:]
lon_rho = fh.variables['lon_rho'][:]
lat_rho = fh.variables['lat_rho'][:]
lon_u = fh.variables['lon_u'][:]
lat_u = fh.variables['lat_u'][:]
lon_v = fh.variables['lon_v'][:]
lat_v = fh.variables['lat_v'][:]

lon_psi = lon_psi-360
lon_rho = lon_rho-360
lon_u = lon_u-360
lon_v = lon_v-360

h = fh.variables['h'][:]
# hraw = np.squeeze(fh.variables['hraw'][:])
s_rho = fh.variables['s_rho'][:]
s_w = fh.variables['s_w'][:]
msk = fh.variables['mask_psi'][:]

s_e = fh.variables['salt_east'][0, :].squeeze()
s_w = fh.variables['salt_west'][0, :].squeeze()

t_e = fh.variables['temp_east'][0, :].squeeze()
t_w = fh.variables['temp_west'][0, :].squeeze()

u_e = fh.variables['u_east'][0, :].squeeze()
u_w = fh.variables['u_west'][0, :].squeeze()

v_e = fh.variables['v_east'][0, :].squeeze()
v_w = fh.variables['v_west'][0, :].squeeze()

fh.close()

lyr = np.size(s_rho)

lat_st_e = np.tile(lat_rho[:, -1].squeeze(), (40, 1))
lat_st_w = np.tile(lat_rho[:, 0].squeeze(), (40, 1))

lon_st_e = np.tile(lon_rho[:, -1].squeeze(), (40, 1))
lon_st_w = np.tile(lon_rho[:, 0].squeeze(), (40, 1))

lat_u_e = np.tile(lat_u[:, -1].squeeze(), (40, 1))
lat_u_w = np.tile(lat_u[:, 0].squeeze(), (40, 1))

lon_u_e = np.tile(lon_u[:, -1].squeeze(), (40, 1))
lon_u_w = np.tile(lon_u[:, 0].squeeze(), (40, 1))

lat_v_e = np.tile(lat_v[:, -1].squeeze(), (40, 1))
lat_v_w = np.tile(lat_v[:, 0].squeeze(), (40, 1))

lon_v_e = np.tile(lon_v[:, -1].squeeze(), (40, 1))
lon_v_w = np.tile(lon_v[:, 0].squeeze(), (40, 1))

s_rho = np.matrix(s_rho).T

h_st_e = s_rho*np.matrix(h[:, -1])
h_st_w = s_rho*np.matrix(h[:, 0])

h_u_e = h_st_e
h_u_w = h_st_w

h_v_e = 0.5*(h_st_e[:-1]+h_st_e[1:])
h_v_w = 0.5*(h_st_w[:-1]+h_st_w[1:])

plt.close()
fig = plt.figure()

clevs = np.arange(30., 35., 0.05)

plt.contourf(lat_st_e, h_st_e, s_e)
plt.colorbar()
plt.title('Salinity [PSU] at East Boundary')
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/bdry_e_s.eps',format='eps')
plt.close()

plt.contourf(lat_st_w, h_st_w, s_w)
plt.colorbar()
plt.title('Salinity [PSU] at West Boundary')
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')
plt.savefig('/Volumes/R1/scratch/chuning/gb_roms/figs/bdry_w_s.eps',format='eps')
plt.close()


# plt.close()
# fig = plt.figure()
# 
# lat_min = 57.75
# lat_max = 59.25
# lat_0 = 0.5 * (lat_min + lat_max)
# 
# lon_min = -137.5
# lon_max = -134.5
# lon_0 = 0.5 * (lon_min + lon_max)
# 
# m = Basemap(projection='merc', llcrnrlon=lon_min, llcrnrlat=lat_min,
#             urcrnrlon=lon_max, urcrnrlat=lat_max, lat_0=lat_0, lon_0=lon_0,
#             resolution='f')
# 
# # xstn, ystn = m(ctd['lon_stn'],ctd['lat_stn'])
# # m.plot(xstn,ystn,'ok',ms=2)
# 
# # for i in stn:
# #     plt.text(xstn[i]+1000, ystn[i], "%02d"%i, fontsize=5, va='center')
# 
# m.drawcoastlines(linewidth=.2)
# m.fillcontinents(color='lightgrey')
# mr = m.drawmeridians(np.arange(lon_min, lon_max, 0.5),labels=[0,0,0,1],fontsize=6, linewidth=.2)
# pr = m.drawparallels(np.arange(lat_min, lat_max, 0.25),labels=[1,0,0,0],fontsize=6, linewidth=.2)
# # setlabelrot(mr,-90)
# 
# x, y = m(lon, lat)
# x[msk==0] = np.nan
# y[msk==0] = np.nan
# m.plot(x, y, 'k', lw=0.1)
# m.plot(x.T, y.T, 'k', lw=0.1)
# 
# # # plot topography
# 
# msk1 = (lont>lon_min) & (lont<lon_max)
# msk2 = (latt>lat_min) & (latt<lat_max)
# 
# lont = lont[msk1]
# latt = latt[msk2]
# z = z[msk2,:][:,msk1]
# 
# z = -z
# z[z>500] = 500
# z[z<0] = 0
# clevs = np.arange(0, 501, 10)
# 
# lont, latt = np.meshgrid(lont, latt)
# xt, yt = m(lont, latt)
# 
# # m.contourf(xt, yt, z, clev)
# # m.colorbar()
# 
# xh, yh = m(lon, lat)
# m.contourf(xh, yh, h, clevs)
# m.colorbar()
# 
# plt.savefig('../figs/map_bdry.eps',format='eps')
# plt.close()

