import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from ocean_toolbox import noaa_adcp

info = {'stn' : 'SEA0839',
        'file_dir': '/glade/p/work/chuning/data/NOAA_ADCP/',
        'sl': 'l',
        'Wp_hrs': -1}

crt = noaa_adcp.get_noaa_current(info)
crt()

t = crt.ctime
u = crt.u
v = crt.v
z = crt.z
lat = crt.info['lat']
lon = crt.info['lon']

ubar = u.mean(axis=1)
vbar = v.mean(axis=1)

plt.figure()
plt.plot(ubar, z, 'r')
plt.plot(vbar, z, 'b')
plt.gca().invert_yaxis()
