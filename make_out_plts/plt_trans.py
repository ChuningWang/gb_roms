import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyroms
import glob
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm
import cmocean
from geopy.distance import vincenty
from matplotlib import path
import sys

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

# my inputs
my_year = 2008
pltuv = 1
ftype = 'his'
varlist = ['salt', 'temp', 'dye_01', 'dye_03']
varlist = ['tke', 'gls']
dd = 3

# dicts for variable clims, colormaps, and other properties
clim = {'temp': [2, 10],
        'salt': [15, 35],
        'zeta': [-3, 3],
        'wetdry_mask_rho': [0, 1],
        'dye_01': [0, 1],
        'dye_02': [0, 1],
        'dye_03': [0, 1],
        'tke': [1e-5, 1e0],
        'gls': [1e-6, 1e-4],
       }

cmap = {'temp': cmocean.cm.thermal,
        'salt': cmocean.cm.haline,
        'zeta': cmocean.cm.balance,
        'wetdry_mask_rho': cmocean.cm.balance,
        'dye_01': cmocean.cm.matter,
        'dye_02': cmocean.cm.matter,
        'dye_03': cmocean.cm.matter,
        'tke': cmocean.cm.matter,
        'gls': cmocean.cm.matter,
       }

var_omega = ['tke', 'gls']
var_log = ['tke', 'gls']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

grd = pyroms.grid.get_ROMS_grid(grd1)

if grd1=='GB_lr':
    tag = 'GB-CIRC'
if grd1=='GB_hr':
    tag = 'GB-TIDE'

model = 'tmpdir_' + tag + '/outputs/' + str(my_year) + '/'
outputs_dir = model_dir + model
fig_dir = out_dir + 'figs/trans/' + tag +'/' + str(my_year) + '/'

flist = sorted(glob.glob(outputs_dir + '*' + ftype + '*.nc'))
flist = flist[-1:]

zlev = grd.vgrid.N
uvar = 'u'
vvar = 'v'
wvar = 'w'

c0 = np.array([[-137.05000708,   59.05576767],
               [-137.02431432,   59.03650282],
               [-136.99075294,   59.02081148],
               [-136.97549777,   58.99798771],
               [-136.95566605,   58.99085529],
               [-136.9251557 ,   58.98086989],
               [-136.91295157,   58.96375207],
               [-136.89922191,   58.94806073],
               [-136.89617088,   58.92666345],
               [-136.88396674,   58.90954562],
               [-136.85955846,   58.90526617],
               [-136.83667571,   58.90098671],
               [-136.81379295,   58.89528077],
               [-136.79091019,   58.89528077],
               [-136.7741295 ,   58.89813374],
               [-136.74972123,   58.89813374],
               [-136.72226192,   58.89228077],
               [-136.69937916,   58.88386889],
               [-136.67802192,   58.87816295],
               [-136.66429226,   58.87530998],
               [-136.63835847,   58.86817755],
               [-136.60937365,   58.86389809],
               [-136.58954192,   58.86389809],
               [-136.5605571 ,   58.8539127 ],
               [-136.53614882,   58.8439273 ],
               [-136.50716399,   58.82823596],
               [-136.49190882,   58.82110354],
               [-136.47970468,   58.81610354],
               [-136.45987296,   58.81282408],
               [-136.44156675,   58.80997111],
               [-136.42326055,   58.7974122 ],
               [-136.40953089,   58.7954268 ],
               [-136.38969917,   58.77366195],
               [-136.37291848,   58.77250898],
               [-136.34851021,   58.77260303],
               [-136.330204  ,   58.76560303],
               [-136.29816814,   58.76289709],
               [-136.285964  ,   58.76019115],
               [-136.2768109 ,   58.74835278],
               [-136.26003021,   58.73408793],
               [-136.24630055,   58.71839659],
               [-136.23714745,   58.70127877],
               [-136.20968814,   58.68986689],
               [-136.17002469,   58.67702852],
               [-136.15171849,   58.67417555],
               [-136.11358056,   58.66276367],
               [-136.0906978 ,   58.6399399 ],
               [-136.07696815,   58.61711614],
               [-136.06781504,   58.59001292],
               [-136.05866194,   58.55149781],
               [-136.05561091,   58.51726216],
               [-136.0373047 ,   58.47018815],
               [-136.01899849,   58.44736438],
               [-135.99916677,   58.42168765],
               [-135.97323298,   58.38174606],
               [-135.94272263,   58.35036338],
               [-135.91373781,   58.33467204],
               [-135.86949781,   58.3289661 ],
               [-135.82068126,   58.31898071],
               [-135.77949229,   58.31898071],
               [-135.72609919,   58.31042179],
               [-135.67880816,   58.290451  ],
               [-135.62999161,   58.27618615],
               [-135.56286885,   58.26620075],
               [-135.53083299,   58.26192129],
               [-135.48049093,   58.25478887],
               [-135.44082748,   58.25478887],
               [-135.39658748,   58.24480347],
               [-135.34777093,   58.21912673],
               [-135.30810748,   58.21342079],
               [-135.26844404,   58.21484728],
               [-135.21505094,   58.19059703],
               [-135.17386197,   58.16206732],
               [-135.13419852,   58.13781707],
               [-135.08995853,   58.12069925],
               [-135.05792266,   58.10358142]])

# c0 = c0[::2, :]

lon_ct = c0[:, 0]
lat_ct = c0[:, 1]

ct_tr = (len(lon_ct)-1)*dd
lon_t = np.zeros(ct_tr)
lat_t = np.zeros(ct_tr)

lon_b0 = lon_ct[:1]-0.01
lat_b0 = lat_ct[:1]+0.01
lon_b1 = lon_ct-0.01
lat_b1 = lat_ct-0.01
lon_b2 = lon_ct[-1:]+0.01
lat_b2 = lat_ct[-1:]-0.01
lon_b3 = lon_ct[::-1]+0.01
lat_b3 = lat_ct[::-1]+0.01
lon_bdry = np.concatenate((lon_b0, lon_b1, lon_b2, lon_b3))
lat_bdry = np.concatenate((lat_b0, lat_b1, lat_b2, lat_b3))

for i in range(len(lon_ct)-1):
    lon_t[i*dd:(i+1)*dd] = np.linspace(lon_ct[i], lon_ct[i+1], dd+1)[:-1]
    lat_t[i*dd:(i+1)*dd] = np.linspace(lat_ct[i], lat_ct[i+1], dd+1)[:-1]

# read grid information
lat = grd.hgrid.lat_rho
lon = grd.hgrid.lon_rho

# interpolate topography
h = grd.vgrid.h
ang = grd.hgrid.angle_rho
msk = grd.hgrid.mask_rho

lon = lon[msk==1]
lat = lat[msk==1]
h = h[msk==1]
ang = ang[msk==1]

p0 = [(lon_bdry[i], lat_bdry[i]) for i in range(len(lon_bdry))]
p = path.Path(p0)
pc = p.contains_points(np.array([lon, lat]).T)
lon = lon[pc]
lat = lat[pc]
h = h[pc]
ang = ang[pc]

h_tr = griddata(lon, lat, h, lon_t, lat_t, interp='linear').diagonal()

# get rho point depth
Cs = grd.vgrid.Cs_r
z_tr = np.dot(np.matrix(h_tr).T, np.matrix(Cs))
z_tr = z_tr.T

# calculate distance
dis = np.zeros(lat_t.size)
for i in range(1, lat_t.size):
    dis[i] = vincenty(
                      (lat_t[i-1], lon_t[i-1]),
                      (lat_t[i], lon_t[i])
                     ).meters
dis = np.cumsum(dis)
dis = dis/1000  # [km]
dis = np.tile(dis, (zlev, 1))

# initiate vairables
var_tr = np.zeros((zlev, ct_tr))

# -------------------------------------------------------------------------------
# if plot velocity vector, also calculate and define these variables
if pltuv==1:
    latu = grd.hgrid.lat_u
    lonu = grd.hgrid.lon_u
    latv = grd.hgrid.lat_v
    lonv = grd.hgrid.lon_v
    msku = grd.hgrid.mask_u
    mskv = grd.hgrid.mask_v

    lonu = lonu[msku==1]
    latu = latu[msku==1]
    lonv = lonv[mskv==1]
    latv = latv[mskv==1]

    pcu = p.contains_points(np.array([lonu, latu]).T)
    pcv = p.contains_points(np.array([lonv, latv]).T)
    lonu = lonu[pcu]
    latu = latu[pcu]
    lonv = lonv[pcv]
    latv = latv[pcv]

    lon_t2 = lon_t[::dd]
    lat_t2 = lat_t[::dd]
    h_tr2 = h_tr[::dd]
    z_tr2 = z_tr[:, ::dd]
    dis2 = dis[:, ::dd]

    # calculate angle
    ang_tr = griddata(lon, lat, ang, lon_t, lat_t, interp='linear').diagonal()
    ang_tr2 = ang_tr[::dd]
    ang_add = np.zeros(len(ang_tr2))
    dx = 59
    dy = 111
    dvec = np.diff(lon_t2*dx + 1j*lat_t2*dy)
    ang_add[:-1] = np.angle(dvec)
    ang_add[-1] = ang_add[-2]
    ang_tr2 = ang_tr2-ang_add

    U_tr2 = np.zeros((zlev, len(lon_t2)))
    w_tr2_sw = np.zeros((zlev+1, len(lon_t2)))

# -------------------------------------------------------------------------------
# make plots
plt.switch_backend('Agg')

for var in varlist:
    print('For ' + var)
    # if var in ['temp']:
    #     clim = [2, 10]
    #     cmap_var = cmocean.cm.thermal
    # elif var in ['salt']:
    #     clim = [15, 35]
    #     cmap_var = cmocean.cm.haline
    # elif var in ['dye_01', 'dye_02', 'dye_03']:
    #     clim = [0, 1]
    #     cmap_var = cmocean.cm.matter
    # else:
    #     clim = [0, 1]
    #     cmap_var = cmocean.cm.matter

    if var in clim.keys():
        clim_var = clim[var]
        cmap_var = cmap[var]
    else:
        clim_var = [0, 1]
        cmap_var = cmocean.cm.matter

    for fn in flist:
        # read data
        tag = fn.split('/')[-1].split('.')[0]
        print('processing ' + tag + ' ...')
        fh = nc.Dataset(fn)
        t = fh.variables['ocean_time'][:]
        tunit = (fh.variables['ocean_time']).units
        data = fh.variables[var][:]
        if pltuv==1:
            u = fh.variables[uvar][:]
            v = fh.variables[vvar][:]
            w = fh.variables[wvar][:]
        fh.close()

        if var in var_omega:
            data = 0.5*(data[:, :-1, :, :] + data[:, 1:, :, :])

        for tt in range(len(t)):
            ttag = nc.num2date(t[tt], tunit).strftime("%Y-%m-%d_%H:%M:%S")
            for i in range(zlev):
                dslice = data[tt, i, :, :][msk==1][pc]
                var_tr[i, :] = griddata(lon, lat, dslice, lon_t, lat_t, interp='linear').diagonal()

            if pltuv==1:
                for i in range(zlev):
                    uslice = u[tt, i, :, :][msku==1][pcu]
                    vslice = v[tt, i, :, :][mskv==1][pcv]
                    u_tr2 = griddata(lonu, latu, uslice, lon_t2, lat_t2, interp='linear').diagonal()
                    v_tr2 = griddata(lonv, latv, vslice, lon_t2, lat_t2, interp='linear').diagonal()

                    U_tr2[i, :] = u_tr2*np.cos(ang_tr2)-v_tr2*np.sin(ang_tr2)

                for i in range(zlev+1):
                    wslice = w[tt, i, :, :][msk==1][pc]
                    w_tr2_sw[i, :] = griddata(lon, lat, wslice, lon_t2, lat_t2, interp='linear').diagonal()

                w_tr2 = 0.5*(w_tr2_sw[1:, :]+w_tr2_sw[:-1, :])

            # make plot
            if var in var_log:
                pcm = plt.pcolormesh(dis, z_tr, var_tr, norm=LogNorm(vmin=clim_var[0], vmax=clim_var[1]), cmap=cmap_var)
            else:
                pcm = plt.pcolormesh(dis, z_tr, var_tr, cmap=cmap_var)
            plt.clim(clim_var[0], clim_var[1])
            cb = plt.colorbar()

            if pltuv==1:
                qv = plt.quiver(dis2, z_tr2, U_tr2, w_tr2, scale=100)

            plt.title(var + '_' + grd.name + '_' + ftype + '_' + ttag)
            plt.savefig(fig_dir + var + '_' + grd.name + '_' + ftype + '_' + ttag + '.png')

            pcm.remove()
            cb.remove()
            if pltuv==1:
                qv.remove()

plt.close()
