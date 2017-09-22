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
fh = nc.Dataset(out_dir + 'bc_ic/GlacierBay_lr_bdry_2008_SODA3.3.1_0.25.nc', mode='r')
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
Cs_r = fh.variables['Cs_r'][:]
Cs_w = fh.variables['Cs_w'][:]
msk = fh.variables['mask_psi'][:]
t = fh.variables['ocean_time'][:]

s_e = fh.variables['salt_east'][:].squeeze()
s_w = fh.variables['salt_west'][:].squeeze()

t_e = fh.variables['temp_east'][:].squeeze()
t_w = fh.variables['temp_west'][:].squeeze()

u_e = fh.variables['u_east'][:].squeeze()
u_w = fh.variables['u_west'][:].squeeze()

v_e = fh.variables['v_east'][:].squeeze()
v_w = fh.variables['v_west'][:].squeeze()

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

Cs_r = np.matrix(Cs_r).T

h_st_e = Cs_r*np.matrix(h[:, -1])
h_st_w = Cs_r*np.matrix(h[:, 0])

h_u_e = h_st_e
h_u_w = h_st_w

h_v_e = 0.5*(h_st_e[:, :-1]+h_st_e[:, 1:])
h_v_w = 0.5*(h_st_w[:, :-1]+h_st_w[:, 1:])

# ------------------------------------------------------------------------------------------------------
# for i in range(len(t)):
for i in range(5):
    # plt.pcolormesh(lat_st_e, h_st_e, s_e[i, :, :])
    # plt.colorbar()
    # plt.title('Salinity [PSU] at East Boundary')
    # plt.xlabel('Latitude')
    # plt.ylabel('Depth [m]')
    # plt.savefig(out_dir + 'figs/bdry_e_s_' + str(t[i]) + '.tiff',format='tiff')
    # plt.close()

    # plt.pcolormesh(lat_st_w, h_st_w, s_w[i, :, :])
    # plt.colorbar()
    # plt.title('Salinity [PSU] at West Boundary')
    # plt.xlabel('Latitude')
    # plt.ylabel('Depth [m]')
    # plt.savefig(out_dir + 'figs/bdry_w_s_' + str(t[i]) + '.tiff',format='tiff')

    # plt.close()
    # plt.pcolormesh(lat_st_e, h_st_e, t_e[i, :, :])
    # plt.colorbar()
    # plt.title('Temperature [degC] at East Boundary')
    # plt.xlabel('Latitude')
    # plt.ylabel('Depth [m]')
    # plt.savefig(out_dir + 'figs/bdry_e_t_' + str(t[i]) + '.tiff',format='tiff')
    # plt.close()

    # plt.pcolormesh(lat_st_w, h_st_w, t_w[i, :, :])
    # plt.colorbar()
    # plt.title('Temperature [degC] at West Boundary')
    # plt.xlabel('Latitude')
    # plt.ylabel('Depth [m]')
    # plt.savefig(out_dir + 'figs/bdry_w_t_' + str(t[i]) + '.tiff',format='tiff')
    # plt.close()

    plt.pcolormesh(lat_u_w, h_u_w, u_w[i, :, :])
    plt.colorbar()
    plt.title('U at West Boundary')
    plt.xlabel('Latitude')
    plt.ylabel('Depth [m]')
    plt.savefig(out_dir + 'figs/bdry/bdry_w_u_' + str(t[i]) + '.tiff',format='tiff')
    plt.close()

    plt.pcolormesh(lat_v_w, h_v_w, v_w[i, :, :])
    plt.colorbar()
    plt.title('V at West Boundary')
    plt.xlabel('Latitude')
    plt.ylabel('Depth [m]')
    plt.savefig(out_dir + 'figs/bdry/bdry_w_v_' + str(t[i]) + '.tiff',format='tiff')
    plt.close()
