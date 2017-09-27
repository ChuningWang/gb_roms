import numpy as np
import matplotlib.pyplot as plt

from ocean_toolbox import ctd
info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd.nc',
        'sl': 's',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
       }
c = ctd.ctd(info)
c()


