""" plot ctd climatology. """

import numpy as np
import matplotlib.pyplot as plt
from geopy.distance import vincenty

from ocean_toolbox import ctd

stn_list = [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': stn_list,
        'clim_deep_interp': 'yes',
        'filter': 'no'}

ctd_gb = ctd.ctd(info)
ctd_gb()

# ------------ calculate distance ------------------------------------
dis = np.zeros(len(stn_list))
for i in range(1, len(dis)):
    dis[i] = vincenty(
                      (ctd_gb.climatology['lat'][i-1], ctd_gb.climatology['lon'][i-1]),
                      (ctd_gb.climatology['lat'][i], ctd_gb.climatology['lon'][i])
                     ).meters
dis = np.cumsum(dis)/1000.  # [km]

depth = ctd_gb.data['fathometer_depth'][ctd_gb.climatology['station'].astype(int)]

# ------------ make plots --------------------------------------------
plt.contourf(dis, ctd_gb.climatology['z'], ctd_gb.climatology['salt'][:, 0, :])
plt.fill_between(dis, depth,  500, facecolor='lightgrey')
plt.colorbar()
