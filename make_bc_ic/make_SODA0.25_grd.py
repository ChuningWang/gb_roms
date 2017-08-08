import netCDF4 as nc
import numpy as np

import read_host_info
sv = read_host_info.read_host_info()
dst_dir = sv['out_dir']
soda_dir = sv['soda_dir']

fin  = nc.Dataset(soda_dir + 'soda3.3.1_5dy_ocean_or_2000_01_03.nc', 'r')
fout = nc.Dataset(soda_dir + 'grid/SODA3_0.25deg_grid.nc', 'w')

fout.createDimension('xu_ocean', len(fin.dimensions['xu_ocean']))
fout.createDimension('xt_ocean', len(fin.dimensions['xt_ocean']))
fout.createDimension('yu_ocean', len(fin.dimensions['yu_ocean']))
fout.createDimension('yt_ocean', len(fin.dimensions['yt_ocean']))
fout.createDimension('st_ocean', len(fin.dimensions['st_ocean']))
fout.createDimension('sw_ocean', len(fin.dimensions['sw_ocean']))
fout.createDimension('st_edges_ocean', len(fin.dimensions['st_edges_ocean']))
fout.createDimension('sw_edges_ocean', len(fin.dimensions['sw_edges_ocean']))

fout.createVariable('geolon_t', 'd', ('yt_ocean', 'xt_ocean'))
fout.createVariable('geolat_t', 'd', ('yt_ocean', 'xt_ocean'))
fout.createVariable('geolon_c', 'd', ('yu_ocean', 'xu_ocean'))
fout.createVariable('geolat_c', 'd', ('yu_ocean', 'xu_ocean'))
fout.createVariable('st_ocean', 'd', ('st_ocean'))
fout.createVariable('sw_ocean', 'd', ('sw_ocean'))
fout.createVariable('st_edges_ocean', 'd', ('st_edges_ocean'))
fout.createVariable('sw_edges_ocean', 'd', ('sw_edges_ocean'))
fout.createVariable('coriolis_param', 'd', ('yt_ocean', 'xt_ocean'))
fout.createVariable('kmt', 'd', ('yt_ocean', 'xt_ocean'))
fout.createVariable('ht', 'd', ('yt_ocean', 'xt_ocean'))
fout.createVariable('kmu', 'd', ('yu_ocean', 'xu_ocean'))

lont, latt = np.meshgrid(fin.variables['xt_ocean'][:], fin.variables['yt_ocean'][:])
lonu, latu = np.meshgrid(fin.variables['xu_ocean'][:], fin.variables['yu_ocean'][:])

fout.variables['geolon_t'][:] = lont
fout.variables['geolat_t'][:] = latt
fout.variables['geolon_c'][:] = lonu
fout.variables['geolat_c'][:] = latu
fout.variables['st_ocean'][:] = fin.variables['st_ocean'][:]
fout.variables['sw_ocean'][:] = fin.variables['sw_ocean'][:]
fout.variables['st_edges_ocean'][:] = fin.variables['st_edges_ocean'][:]
fout.variables['sw_edges_ocean'][:] = fin.variables['sw_edges_ocean'][:]

# fout.Dataset['coriolis_param'][:] = fin.Dataset['coriolis_param'][:]
# fout.Dataset['kmt'][:] = fin.Dataset['kmt'][:]
# fout.Dataset['ht'][:] = fin.Dataset['ht'][:]
# fout.Dataset['kmu'][:] = fin.Dataset['kmu'][:]
# 
# fout.Dataset['kmt']._FillValue = fin.Dataset['kmt']._FillValue
# fout.Dataset['ht']._FillValue = fin.Dataset['ht']._FillValue
# fout.Dataset['kmu']._FillValue = fin.Dataset['kmu']._FillValue

fin.close()
fout.close()
