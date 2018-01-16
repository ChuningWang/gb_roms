import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from datetime import datetime

from ocean_toolbox import ctd
info = {'data_dir': '/Users/CnWang/Documents/gb_roms/ctd_raw/',
        'file_dir': '/Users/CnWang/Documents/gb_roms/',
        'file_name': 'ctd.nc',
        'sl': 'l',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'clim_station': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20],
        'clim_deep_interp': 'yes',
        'filter': 'no',
       }
c = ctd.ctd(info)
c()
t0 = nc.date2num(datetime(2016, 07, 15), 'days since 1900-01-01') 
# stn_list = [21, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
stn_list = [20, 19, 18, 17, 16, 14, 13, 4, 3, 2, 1, 0]
c.get_trans(['temp', 'salt', 'fluor'], stn_list, t0)
