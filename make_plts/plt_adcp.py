'''
plot adcp station data
'''

import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from ocean_toolbox import noaa_adcp
import pyroms

import read_host_info
sv = read_host_info.read_host_info()
data_dir = sv['in_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'
grd = pyroms.grid.get_ROMS_grid(grd1)
ang = grd.hgrid.angle_rho.mean()

# read data
info = {'stn' : 'SEA0839',
        'file_dir': data_dir + 'NOAA_ADCP/',
        'sl': 'l',
        'Wp_hrs': -1}

crt = noaa_adcp.get_noaa_current(info)
crt()

t = crt.ctime
Uraw = crt.u + 1j*crt.v
z = crt.z
lat = crt.info['lat']
lon = crt.info['lon']

# rotate speed vector
U = Uraw*np.exp(ang*1j)

Urawbar = Uraw.mean(axis=1)
Ubar = U.mean(axis=1)

plt.figure()
# plt.plot(np.real(Ubar), z, 'r')
# plt.plot(np.imag(Ubar), z, 'b')
# plt.gca().invert_yaxis()

plt.quiver(np.zeros(z.shape), -z, np.real(Urawbar), np.imag(Urawbar))
