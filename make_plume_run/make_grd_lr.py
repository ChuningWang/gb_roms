import numpy as np
import netCDF4 as nc

import pyroms
import pyroms_toolbox
import bathy_smoother

import read_host_info
sv = read_host_info.read_host_info()
bathy_dir = sv['in_dir']
out_dir = sv['out_dir']

# ------------------------------------------------------------------------
grd1 = 'GB_lr'
grd2 = 'GB_300m_orig'
grd_name = 'GlacierBay_lr'

hgrd = pyroms.grid.get_ROMS_hgrid(grd2)
water = hgrd.mask_rho

# ------------------------------------------------------------------------
# mask out small channels
water[190:260, :10] = 0
water[:100, :55] = 0
water[90:125, 10:25] = 0
water[:45, 95:120] = 0
water[:15, 115:135] = 0
water[330:360, 210:] = 0

# ------------------------------------------------------------------------
msk_c = np.loadtxt('mask_change_GB_lr.txt')
for i in range(len(msk_c)):
    hgrd.mask_rho[int(msk_c[i, 1]), int(msk_c[i, 0])] = msk_c[i, 2]
print 'mask done...'
 
# ------------------------------------------------------------------------
# defind the boundary of mapping domain
lat_min = 57.
lat_max = 60.
lat_0 = 0.5 * (lat_min + lat_max)

lon_min = -138.
lon_max = -134.
lon_0 = 0.5 * (lon_min + lon_max)

# ------------------------------------------------------------------------
# load bathymetry
fin = nc.Dataset('hraw.nc', 'r')
lat_raw = fin.variables['lat'][:]
lon_raw = fin.variables['lon'][:]
hraw = fin.variables['hraw'][:]
fin.close()

h = hraw.copy()

# shapiro filter
h = pyroms_toolbox.shapiro_filter.shapiro2(h, 32)

# ------------------------------------------------------------------------
# locally constrain hmin at some location
hmin0 = 20  # m
h1 = h[160:180, 139:144]
h1[h1<hmin0] = hmin0
h[160:180, 139:144] = h1

# ------------------------------------------------------------------------
# insure that depth is always deeper than hmin
hmin = 10
hmax = 450
h = pyroms_toolbox.change(h, '<', hmin, hmin)
# constrain maximum depth
h = pyroms_toolbox.change(h, '>', hmax, hmax)
# fix depth of land points
h[water == 0] = hmin

# ------------------------------------------------------------------------
# smooth bathymetry
rx0_max = 0.35
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, water)
print 'Max Roughness value is: ', RoughMat.max()
h = bathy_smoother.bathy_smoothing.smoothing_Positive_rx0(water, h, rx0_max)
RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(h, water)
print 'Max Roughness value is: ', RoughMat.max()

# ------------------------------------------------------------------------
# constrain water depth of coastal cells. This is critical for p-source
# style runoff.
hmin0 = 3
idx = []
idy = []
sign = []
rdir = []
river = []
# read in coast cells
river_num = 0
fin = open('river_cells.txt', 'r')
for line in fin:
    llist = line.rstrip('\n').split(',')
    if int(llist[0]) == -1:
        river_num += 1
        irdir = int(llist[1])
        isign = int(llist[2])
    elif int(llist[0]) == -10:
        break
    else:
        river.append(river_num)
        idx.append(int(llist[0]))
        idy.append(int(llist[1]))
        rdir.append(irdir)
        sign.append(isign)

idx = np.array(idx)
idy = np.array(idy)
rdir = np.array(rdir)
sign = np.array(sign)
river = np.array(river)

h[idx, idy] = hmin0

# smooth fake estuaries
def local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.3):

    h1 = h[xmin:xmax, ymin:ymax]
    msk1 = water[xmin:xmax, ymin:ymax]

    # smooth bathymetry
    hs = bathy_smoother.bathy_smoothing.smoothing_Negative_rx0(msk1, h1, rx0_max)
    h[xmin:xmax, ymin:ymax] = hs

for i in range(1, river_num+1):
    msk = river == i
    iidx = idx[msk]
    iidy = idy[msk]
    irdir = rdir[msk][0]
    isign = sign[msk][0]
    if (isign == 1) & (irdir == 0):
        xmin = iidx[0]; xmax = iidx[-1]+1
        ymin = iidy[0]; ymax = iidy[-1]+5
    elif (isign == 1) & (irdir == 1):
        xmin = iidx[0]; xmax = iidx[-1]+5
        ymin = iidy[0]; ymax = iidy[-1]+1
    elif (isign == -1) & (irdir == 0):
        xmin = iidx[0]; xmax = iidx[-1]+1
        ymin = iidy[0]-4; ymax = iidy[-1]+1
    elif (isign == -1) & (irdir == 1):
        xmin = iidx[0]-4; xmax = iidx[-1]+1
        ymin = iidy[0]; ymax = iidy[-1]+1

    local_smooth(h, water, xmin, xmax, ymin, ymax, rx0_max=0.2)

# smooth again
h = bathy_smoother.bathy_smoothing.smoothing_Negative_rx0(water, h, rx0_max)

# ------------------------------------------------------------------------
# # insure that depth is always deeper than hmin (again)
# h = pyroms_toolbox.change(h, '<', hmin, hmin0)
# # constrain maximum depth
# h = pyroms_toolbox.change(h, '>', hmax, hmax)
# # # fix depth of land points
# h[water == 0] = hmin0

# ------------------------------------------------------------------------
# redesign the vertical coordinate
theta_b = 2.0
theta_s = 8.0
Tcline = 10
N = 40
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

# ------------------------------------------------------------------------
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
# write grid file
pyroms.grid.write_ROMS_grid(grd, filename = out_dir + 'grd/' + grd_name + '_grd.nc')
