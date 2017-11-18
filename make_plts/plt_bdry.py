import sys
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pyroms

# ------------------------------------------------------------------------------------------------------
def setlabelrot(x, rot):
    for m in x:
        for t in x[m][1]:
            t.set_rotation(rot)

# ------------------------------------------------------------------------------------------------------
import read_host_info
sv = read_host_info.read_host_info()
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

# Read grid
grd = pyroms.grid.get_ROMS_grid(grd1)
lon_rho = grd.hgrid.lon_rho
lat_rho = grd.hgrid.lat_rho
lon_u = grd.hgrid.lon_u
lat_u = grd.hgrid.lat_u
lon_v = grd.hgrid.lon_v
lat_v = grd.hgrid.lat_v

h = grd.vgrid.h
Cs_r = grd.vgrid.Cs_r

fh = nc.Dataset(out_dir + 'bc_ic/GlacierBay_lr_bc_2008_CTD_ADCP.nc', mode='r')
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

fig, ax = plt.subplots()
plt.xlabel('Latitude')
plt.ylabel('Depth [m]')

# ------------------------------------------------------------------------------------------------------
for i in range(len(t)):
    pcm = ax.pcolormesh(lat_st_e, h_st_e, s_e[i, :, :])
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    # plt.title('Salinity [PSU] at East Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_e_s_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)

    pcm = ax.pcolormesh(lat_st_w, h_st_w, s_w[i, :, :])
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    # plt.title('Salinity [PSU] at West Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_w_s_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)

    pcm = ax.pcolormesh(lat_st_e, h_st_e, t_e[i, :, :])
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    # plt.title('Temperature [degC] at East Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_e_t_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)

    pcm = ax.pcolormesh(lat_st_w, h_st_w, t_w[i, :, :])
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    # plt.title('Temperature [degC] at West Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_w_t_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)

    pcm = ax.pcolormesh(lat_u_e, h_u_e, u_e[i, :, :])
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    pct = ax.contour(lat_u_e, h_u_e, u_e[i, :, :], [0.], linewidths=0.5, colors='k')
    # plt.title('U at East Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_e_u_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)
    for cc in pct.collections:
        cc.remove()

    pcm = ax.pcolormesh(lat_u_w, h_u_w, u_w[i, :, :])
    pcm.set_clim(-0.2, 0.2)
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    pct = ax.contour(lat_u_w, h_u_w, u_w[i, :, :], [0.], linewidths=0.5, colors='k')
    # plt.title('U at West Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_w_u_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)
    for cc in pct.collections:
        cc.remove()

    pcm = ax.pcolormesh(lat_v_e, h_v_e, v_e[i, :, :])
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    pct = ax.contour(lat_v_e, h_v_e, v_e[i, :, :], [0.], linewidths=0.5, colors='k')
    # plt.title('V at East Boundary')
    fig.savefig(out_dir + 'figs/bdry/bdry_e_v_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)
    for cc in pct.collections:
        cc.remove()

    pcm = ax.pcolormesh(lat_v_w, h_v_w, v_w[i, :, :])
    pcm.set_clim(-0.2, 0.2)
    cbar_ax = fig.add_axes([0.9, 0.10, 0.02, 0.8])
    cb = plt.colorbar(pcm, cax=cbar_ax)
    pct = ax.contour(lat_v_w, h_v_w, v_w[i, :, :], [0.], linewidths=0.5, colors='k')
    # plt.title('V at West Boundary')
    plt.savefig(out_dir + 'figs/bdry/bdry_w_v_' + str(t[i]) + '.png',format='png')
    pcm.remove()
    fig.delaxes(cbar_ax)
    for cc in pct.collections:
        cc.remove()
