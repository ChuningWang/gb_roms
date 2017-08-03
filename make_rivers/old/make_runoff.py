import subprocess
import numpy as np
import netCDF4 as nc
from datetime import datetime

import pyroms
import pyroms_toolbox

tstep = 0

# load hydrology data from D. Hill
# data_dir = '/Volumes/R1/ROMS/hydrology/GOA/'
# filelst = subprocess.check_output(['ls', data_dir]).replace('/n', ' ').split()
fh = nc.Dataset('../data/GB_discharge.nc', 'r')
t = fh.variables['t'][tstep]
lon = fh.variables['lon'][:]
lat = fh.variables['lat'][:]
d = fh.variables['discharge'][tstep, :, :]
coast_cells = fh.variables['coast_cells'][:]
fh.close()

print("Loading GB grid:")
# load GB grid object
grd = pyroms.grid.get_ROMS_grid('GB')




# # load 2-dimentional interannual discharge data 
# # from 1948-2007. See Dai and Trenberth (2002) and Dai et al. (2009)
# print 'Load discharge data'
# nc_data = netCDF.Dataset('/Volumes/R2/Data/Dai-Trenberth_runoff/runoff.daitren.iaf.101020.nc', 'r')
# data = nc_data.variables['runoff'][:]
# # time with leap year
# #time = np.zeros(36)
# time = np.zeros(data.shape[0])
# t0 = 17531 #12/31/1947
# #t0 = 37255 #12/31/2001
# tidx=0
# #for year in range(2002,2005):
# for year in range(1980,2008):
#     if year%4 == 0:
#       daysinmonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#     else :
#       daysinmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#     for month in range(12):
#       if tidx == 0:
#         time[tidx] = t0 + daysinmonth[month]/2.
#       else:
#         if month == 0:
#           time[tidx] = time[tidx-1] + daysinmonth[month]/2. + 31/2.
#         else :
#           time[tidx] = time[tidx-1] + daysinmonth[month]/2. + daysinmonth[month-1]/2.
#       tidx = tidx + 1
# 
# # define some variables
# wts_file = 'remap_weights_daitren_to_NWA_conservative_nomask.nc'
# nt = data.shape[0]
# Mp, Lp = grd.hgrid.mask_rho.shape
# spval = -1e30
# runoff_raw = np.zeros((Mp,Lp))
# runoff = np.zeros((Mp,Lp))
# rspread = 6
# 
# print("Creating runoff file:")
# # create runoff file
# #runoff_file = 'runoff_NWA_daitren_inter_annual_2002-2004.nc'
# runoff_file = 'runoff_NWA_daitren_inter_annual_1948-2007.nc'
# nc = netCDF.Dataset(runoff_file, 'w', format='NETCDF3_CLASSIC')
# nc.Description = 'D. Hill Glacier Bay hydrology map, since 1981'
# nc.Author = 'make_runoff.py'
# nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
# nc.title = 'D. Hill hydrology map'
# 
# print("Creating variables:")
# # creat dimensions and variables
# nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
# nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
# nc.createDimension('runoff_time', None)
# 
# nc.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
# nc.variables['lon_rho'].long_name = 'longitude of RHO-points'
# nc.variables['lon_rho'].units = 'degree_east'
# nc.variables['lon_rho'].field = 'lon_rho, scalar'
# nc.variables['lon_rho'][:] = grd.hgrid.lon_rho
# 
# nc.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
# nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
# nc.variables['lat_rho'].units = 'degree_north'
# nc.variables['lat_rho'].field = 'lat_rho, scalar'
# nc.variables['lat_rho'][:] = grd.hgrid.lat_rho
# 
# nc.createVariable('runoff_time', 'f8', ('runoff_time'))
# nc.variables['runoff_time'].long_name = 'time'
# nc.variables['runoff_time'].units = 'days since 1900-01-01 00:00:00'
# nc.variables['runoff_time'][:] = time
# 
# nc.createVariable('Runoff_raw', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
# nc.variables['Runoff_raw'].long_name = 'Dai_Trenberth River Runoff raw'
# nc.variables['Runoff_raw'].missing_value = str(spval)
# nc.variables['Runoff_raw'].units = 'kg/s/m^2'
# 
# nc.createVariable('Runoff', 'f8', ('runoff_time', 'eta_rho', 'xi_rho'))
# nc.variables['Runoff'].long_name = 'Dai_Trenberth River Runoff'
# nc.variables['Runoff'].missing_value = str(spval)
# nc.variables['Runoff'].units = 'kg/s/m^2'
# 
# 
# # get littoral (here 4 cells width)
# width = 4
# idx = []
# idy = []
# maskl = grd.hgrid.mask_rho.copy()
# for w in range(width):
#     lit = pyroms_toolbox.get_littoral(maskl)
#     idx.extend(lit[0])
#     idy.extend(lit[1])
#     maskl[lit] = 0
# 
# littoral_idx = (np.array(idx), np.array(idy))
# maskl = np.zeros(grd.hgrid.mask_rho.shape)
# maskl[littoral_idx] = 1
# 
# mask_idx = np.where(grd.hgrid.mask_rho == 0)
# 
# nct=0
# print "\n nct = ", nct
# #for t in range(672,672+36):  # year 2002 to 2004
# for t in range(nt):
#     print 'Remapping runoff for time %f' %time[nct]
#     # conservative horizontal interpolation using scrip
#     runoff_raw = pyroms.remapping.remap(data[t,:,:], wts_file, \
#                                            spval=spval)
#     idx = np.where(runoff_raw != 0)
#     runoff = pyroms_toolbox.move_runoff(runoff_raw, \
#                   np.array(idx).T + 1, np.array(littoral_idx).T + 1, maskl, \
#                   grd.hgrid.x_rho, grd.hgrid.y_rho, grd.hgrid.dx, grd.hgrid.dy)
# 
#     # spread the runoff within the littoral band
#     runoff_spread = np.zeros((Mp,Lp))
#     idx = np.where(runoff != 0)
#     for p in range(np.size(idx,1)):
#         j = range(max(0,idx[0][p]-rspread), min(Mp-1,idx[0][p]+rspread+1))
#         i = range(max(0,idx[1][p]-rspread), min(Lp-1,idx[1][p]+rspread+1))
#         ji = np.meshgrid(j,i)
#         sidx = np.where(maskl[ji] == 1)
#         nbpt = np.size(sidx) / 2
#         rpt = runoff[idx[0][p],idx[1][p]] * grd.hgrid.dx[idx[0][p],idx[1][p]] * grd.hgrid.dy[idx[0][p],idx[1][p]]
#         rpt = rpt / nbpt
#         for pa in range(nbpt):
#             pai = sidx[0][pa] + ji[1].min()
#             paj = sidx[1][pa] + ji[0].min()
#             runoff_spread[paj, pai] = runoff_spread[paj, pai] + \
#                            rpt / (grd.hgrid.dx[paj, pai] * grd.hgrid.dy[paj, pai])
# 
# 
#     # spval
#     runoff_spread[mask_idx] = spval
# 
#     # write data in destination file
#     nc.variables['Runoff'][nct] = runoff_spread
#     nc.variables['Runoff_raw'][nct] = runoff_raw
# 
#     nct = nct + 1
#     print "\n nct = ", nct
# 
# # close netcdf file
# nc.close()
