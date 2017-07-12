# from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# from gb_toolbox.gb_ctd import rd_ctd
import numpy as np
import netCDF4 as nc

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# ------------------------------------------------------------------------------------------------------
import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']

# Read grid
fh = nc.Dataset(out_dir + 'bc_ic/GlacierBay_usgs_bdry_spinup_SODA3.3.1.nc', mode='r')
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

# ------------------------------------------------------------------------------------------------------
plt.pcolormesh(lat_st_e, h_st_e, s_e)
plt.colorbar()
plt.title('Salinity [PSU] at East Boundary')
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')
plt.savefig(out_dir + 'figs/bdry_e_s.tiff',format='tiff')
plt.close()

plt.pcolormesh(lat_st_w, h_st_w, s_w)
plt.colorbar()
plt.title('Salinity [PSU] at West Boundary')
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')
plt.savefig(out_dir + 'figs/bdry_w_s.tiff',format='tiff')

plt.close()
plt.pcolormesh(lat_st_e, h_st_e, t_e)
plt.colorbar()
plt.title('Temperature [degC] at East Boundary')
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')
plt.savefig(out_dir + 'figs/bdry_e_t.tiff',format='tiff')
plt.close()

plt.pcolormesh(lat_st_w, h_st_w, t_w)
plt.colorbar()
plt.title('Temperature [degC] at West Boundary')
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')
plt.savefig(out_dir + 'figs/bdry_w_t.tiff',format='tiff')
plt.close()

