from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import style
import netCDF4 as nc

import pyroms

import read_host_info
sv = read_host_info.read_host_info()
frc_dir = sv['out_dir']

# ---------------- functionals ---------------------------------------
def get_zr(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zr = np.empty((ti, vgrid.N) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for n in range(ti):
            for k in range(vgrid.N):
                z0 = vgrid.hc * vgrid.s_rho[k] + (h - vgrid.hc) * vgrid.Cs_r[k]
                zr[n, k, :] = z0 + zeta[n, :] * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for n in range(ti):
            for k in range(vgrid.N):
                z0 = (vgrid.hc * vgrid.s_rho[k] + h * vgrid.Cs_r[k]) / (vgrid.hc + h)
                zr[n, k, :] = zeta[n, :] + (zeta[n, :] + h) * z0

    return zr

# ---------------- some inputs ---------------------------------------
tbase = nc.date2num(datetime(2008, 1, 1), 'days since 1900-01-01')
time2 = np.array([datetime(2008, mm+1, 15) for mm in range(12)])
time_tick = nc.date2num(time2, 'days since 1900-01-01') - tbase

grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)

# ---------------- freshwater runoff ---------------------------------
river_file = frc_dir + 'frc/GlacierBay_lr_rivers_clim_Hill.nc'
fin = nc.Dataset(river_file, 'r')
river_time = fin.variables['river_time'][:] - tbase + 1
river_transport = np.sum(np.abs(fin.variables['river_transport'][:]), axis=-1)
river_temp = fin.variables['river_temp'][:]
fin.close()

# ---------------- atmospheric forcing -------------------------------
swrad_file = frc_dir + 'frc_clim/swrad_clim_JRA55v1.3.nc'
fin = nc.Dataset(swrad_file, 'r')
srf_time = fin.variables['srf_time'][:]-tbase
swrad = fin.variables['swrad'][:].mean(axis=(-2, -1))
fin.close()

tair_file = frc_dir + 'frc_clim/Tair_clim_JRA55v1.3.nc'
fin = nc.Dataset(tair_file, 'r')
tair_time = fin.variables['tair_time'][:]-tbase
tair = fin.variables['Tair'][:].mean(axis=(-2, -1))
fin.close()

frc_yd = np.arange(366)+0.5
swrad2 = np.zeros(frc_yd.shape)
tair2 = np.zeros(frc_yd.shape)
for i, yd in enumerate(frc_yd):
    msk = (srf_time <= yd+0.5) & (srf_time > yd-0.5)
    swrad2[i] = swrad[msk].mean()
    msk = (tair_time <= yd+0.5) & (tair_time > yd-0.5)
    tair2[i] = tair[msk].mean()

# ---------------- bc ------------------------------------------------
bc_file = frc_dir + 'bc_ic/GlacierBay_lr_bc_2008_CTD_ADCP.nc'
fin = nc.Dataset(bc_file, 'r')
h_east = fin.variables['h'][:, -1]
h_west = fin.variables['h'][:, 0]
bc_time = fin.variables['ocean_time'][:]
se = fin.variables['salt_east'][:]
te = fin.variables['temp_east'][:]
ue = fin.variables['u_east'][:]
ve = fin.variables['v_east'][:]
sw = fin.variables['salt_west'][:]
tw = fin.variables['temp_west'][:]
uw = fin.variables['u_west'][:]
vw = fin.variables['v_west'][:]
fin.close()

i_east = np.argmax(h_east)
i_west = np.argmax(h_west)
se = se[:, :, i_east].mean(axis=0)
te = te[:, :, i_east].mean(axis=0)
ue = ue[:, :, i_east].mean(axis=0)
ve = ve[:, :, i_east].mean(axis=0)
sw = sw[:, :, i_west].mean(axis=0)
tw = tw[:, :, i_west].mean(axis=0)
uw = uw[:, :, i_west].mean(axis=0)
vw = vw[:, :, i_west].mean(axis=0)

ze = -grd.vgrid.z_r[:][:, i_east, -1]
zw = -grd.vgrid.z_r[:][:, i_west, 0]
ange = grd.hgrid.angle_rho[i_east, -1]
angw = grd.hgrid.angle_rho[i_west, -1]

Ue = ue + 1j*ve
Uw = uw + 1j*vw
Ue2 = Ue*np.exp(1j*ange)
Uw2 = Uw*np.exp(1j*angw)
ue = Ue2.real
ve = Ue2.imag
uw = Uw2.real
vw = Uw2.imag

# ---------------- make plots ----------------------------------------
style.use('classic')

gs = gridspec.GridSpec(12, 9)
axr1 = plt.subplot(gs[0:3, :])
axs = plt.subplot(gs[3:6, :])
axtse = plt.subplot(gs[7:, 5:7])
axuve = plt.subplot(gs[7:, 7:9])
axtsw = plt.subplot(gs[7:, 0:2])
axuvw = plt.subplot(gs[7:, 2:4])

# fig, (axr1, axs) = plt.subplots(2, 1)
axr2 = axr1.twinx()
axt = axs.twinx()

axr1.tick_params(axis='y', colors='b')
axr2.tick_params(axis='y', colors='r')
axr1.set_ylabel(r'$R$ [m$^3$s$^{-1}$]', color='b', fontsize=7)
axr2.set_ylabel(r'$T_R$ [$^{\circ}$C]', color='r', fontsize=7)
axr1.set_xticks(time_tick)
axr1.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], fontsize=7)
axr1.xaxis.tick_top()
axr2.set_xticks(time_tick)
# axr2.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])

axr1.grid('on')
axr1.set_ylim(0, 3000)
axr1.set_yticks([0, 1000, 2000, 3000])
axr2.set_ylim(3, 7)
axr2.set_yticks([3, 5, 7])

axs.tick_params(axis='y', colors='b')
axt.tick_params(axis='y', colors='r')
axs.set_ylabel(r'SW Rad [W$\cdot$m$^2$]', color='b')
axt.set_ylabel(r'$T_{air}$ [$^{\circ}$C]', color='r')
axs.set_xticks(time_tick)
axs.set_xticklabels([])
# axs.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
axt.set_xticks(time_tick)
# axt.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])

axs.grid('on')
axs.set_ylim(0, 250)
axs.set_yticks([0, 100, 200])
axt.set_ylim(-10, 15)
axt.set_yticks([-10, 0, 10])

axr1.set_xlim(river_time[0], river_time[-1])
axs.set_xlim(river_time[0], river_time[-1])

axr1.plot(river_time, river_transport, 'b')
axr2.plot(river_time, river_temp, 'r')

axs.plot(frc_yd, swrad2, 'b')
axt.plot(frc_yd, tair2-273.15, 'r')

# ---------------- pt II ---------------------------------------------
axtse2 = axtse.twiny()
axtsw2 = axtsw.twiny()
axuve2 = axuve.twiny()
axuvw2 = axuvw.twiny()

axtse.set_ylim(450, 0)
axuve.set_ylim(450, 0)
axtsw.set_ylim(350, 0)
axuvw.set_ylim(350, 0)

axtse.set_xlim(28, 34)
axtse.set_xticks([29, 31, 33])
axtsw.set_xlim(28, 34)
axtsw.set_xticks([29, 31, 33])

axtse2.set_xlim(4, 10)
axtse2.set_xticks([5, 7, 9])
axtsw2.set_xlim(4, 10)
axtsw2.set_xticks([5, 7, 9])

axuve.set_xlim(-0.3, 0.3)
axuve.set_xticks([-0.2, 0, 0.2])
axuvw.set_xlim(-0.3, 0.3)
axuvw.set_xticks([-0.2, 0, 0.2])

axuve2.set_xlim(-0.3, 0.3)
axuve2.set_xticks([-0.2, 0, 0.2])
axuvw2.set_xlim(-0.3, 0.3)
axuvw2.set_xticks([-0.2, 0, 0.2])

axtse.tick_params(axis='x', colors='b')
axtse2.tick_params(axis='x', colors='r')
axuve.tick_params(axis='x', colors='b')
axuve2.tick_params(axis='x', colors='r')
axtsw.tick_params(axis='x', colors='b')
axtsw2.tick_params(axis='x', colors='r')
axuvw.tick_params(axis='x', colors='b')
axuvw2.tick_params(axis='x', colors='r')

axtse2.xaxis.tick_top()
axtsw2.xaxis.tick_top()
axuve2.xaxis.tick_top()
axuvw2.xaxis.tick_top()

axtse.grid('on')
axuve.grid('on')
axtsw.grid('on')
axuvw.grid('on')

axtse.set_xlabel(r'Salt [PSU]', color='b', fontsize=7)
axtse2.set_xlabel(r'Temp [$^{\circ}$C]', color='r', fontsize=7)
axtsw.set_xlabel(r'Salt [PSU]', color='b', fontsize=7)
axtsw2.set_xlabel(r'Temp [$^{\circ}$C]', color='r', fontsize=7)

axuve.set_xlabel(r'U [m$\cdot$s$^{-1}$]', color='b', fontsize=7)
axuve2.set_xlabel(r'V [m$\cdot$s$^{-1}$]', color='r', fontsize=7)
axuvw.set_xlabel(r'U [m$\cdot$s$^{-1}$]', color='b', fontsize=7)
axuvw2.set_xlabel(r'V [m$\cdot$s$^{-1}$]', color='r', fontsize=7)

axtse.set_ylabel(r'Depth [m]', fontsize=7)
axtsw.set_ylabel(r'Depth [m]', fontsize=7)

axuve.set_yticklabels([])
axuvw.set_yticklabels([])

axtse.plot(se, ze, 'b')
axtse2.plot(te, ze, 'r')
axuve.plot(ue, ze, 'b')
axuve2.plot(ve, ze, 'r')

axtsw.plot(sw, zw, 'b')
axtsw2.plot(tw, zw, 'r')
axuvw.plot(uw, zw, 'b')
axuvw2.plot(vw, zw, 'r')

axr1.text(20, 2500, 'a)')
axs.text(20, 200, 'b)')
axtsw.text(29, 330, 'c)')
axuvw.text(-0.2, 330, 'd)')
axtse.text(29, 430, 'e)')
axuve.text(-0.2, 430, 'f)')

plt.savefig(frc_dir + 'figs/osm2018/forcing.png', dpi=600)
plt.close()
