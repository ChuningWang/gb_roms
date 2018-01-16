'''
this is an example to generate a ROMS tide file for TPXO8-atlas tidal product. Till this moment I couldn't found
a good python script to do so, so this script is created for this purpose. TPXO8-atlas is a high resolution global
model of ocean tides. it is available here:

http://volkov.oce.orst.edu/tides/tpxo8_atlas.html

download the 8 major constitutes with 1/30 deg resotuion. Other constitutes are also supported, but haven't been
tested. Only netCDF format is supported in this version, so make sure you download the right files.

The script 'make_remap_weights_file.py' is used to generate weight file between two grids, in this example the
Glacier Bay and TPXO8 grid files. TPXO8 grid file is kept in the same directory as data files:

/Volumes/R1/scratch/chuning/data/tpxo8nc/grid_tpxo8atlas_30_v1.nc

Run this script before 'make_tide.py' to generate the weight files (remap_*.nc).

To use the remap function a new set of grid is created (CGrid_TPXO8). I followed the structure of GLORY CGrid 
from the PYROMS packaage itself, so that it is consistent with PYROMS.

---------------------------------------------------------------------------------------------------------------------
tidal parameters are converted to elliptical parameters using ap2ep, origionally written in MATLAB and later
translated to python package PyFVCOM. See

CGrid_TPXO8/tidal_ellipse.py for more details.

-Chuning Wang, Dec 18 2016

'''

import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
import sys

import pyroms
import pyroms_toolbox

import CGrid_TPXO8

from datetime import datetime, timedelta
from tide_astro import t_vuf

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']

if len(sys.argv)>1:
    grd1 = sys.argv[-1]
else:
    grd1 = 'GB_lr'

# read ROMS grid
dst_grd = pyroms.grid.get_ROMS_grid(grd1)
lat = dst_grd.hgrid.lat_rho
lon = dst_grd.hgrid.lon_rho
msk = dst_grd.hgrid.mask_rho
eta, xi = lat.shape

# read TPXO8 files
pth_tpxo = in_dir + "tpxo8nc/"
src_grd = CGrid_TPXO8.get_nc_CGrid_TPXO8(pth_tpxo+"grid_tpxo8atlas_30_v1.nc")
src_grd_lr = CGrid_TPXO8.get_nc_CGrid_TPXO8(pth_tpxo+"grid_tpxo8_atlas6.nc", name='TPXO8atlas6', xrange=(870, 910), yrange=(1320, 1380))
missing_value = src_grd.missing_value

# -------------------------------------------------------------------------
# tidal constituents
consts1 =['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2'] 
# consts2 = ['MF']
consts2 = []
consts = consts1+consts2
# consts =['M2'] 
consts_num = len(consts)

# define tidal constituents names and periods
tide_name = np.array([ list('Q1  '), list('O1  '), list('P1  '), list('K1  '),
                       list('N2  '), list('M2  '), list('S2  '), list('K2  ')])
#                        list('MF  ')])
tide_period = np.array([26.8683567047119, 25.8193397521973, 24.0658893585205, 23.9344692230225, 
                        12.6583499908447, 12.420599937439, 12, 11.9672346115112])
#                         13.66079*24])

nodal_corr = 1

# -------------------------------------------------------------------------
# read weights file
wts_file_t = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_t_to_rho.nc' 
wts_file_u = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_u_to_rho.nc' 
wts_file_v = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_v_to_rho.nc' 

wts_file_t_lr = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_t_to_rho.nc' 
wts_file_u_lr = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_u_to_rho.nc' 
wts_file_v_lr = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_v_to_rho.nc' 

# grid sub-sample
xrange = src_grd.xrange
yrange = src_grd.yrange

xrange_lr = src_grd_lr.xrange
yrange_lr = src_grd_lr.yrange

# set output file
out_file = dst_grd.name + '_tides_otps.nc'

# -------------------------------------------------------------------------
# preparation for nodal correction
if nodal_corr == 1:
    d0 = datetime(1899, 12, 31, 12, 00, 00)
    # jd = [datetime(2000, 07, 01, 00, 00, 00), datetime(2001, 07, 01, 00, 00, 00)]
    jd = [datetime(2008, 01, 01, 00, 00, 00)]
    d = np.array([(jd[i]-d0).total_seconds()/timedelta(days=1).total_seconds() for i in range(len(jd))])
    # ltype = 'nodal'
    fh = nc.Dataset('tide_consts.nc', 'r')
    name = fh.variables['const_name'][:]
    fh.close()
    ju = 6
    lat = 58

    # idx = np.zeros(consts_num)*np.NaN
    idx = list()
    for i in range(consts_num):
        idx.append(np.where(np.all(name == tide_name[i], axis=1))[0][0])

    V, U, F = t_vuf(d, idx, lat)
    V = V*360  # convert phase to degree
    U = U*360  # convert phase to degree

# -------------------------------------------------------------------------
# initiate variables
hamp = np.zeros((consts_num, eta, xi))*np.nan
hpha = np.zeros((consts_num, eta, xi))*np.nan
cmax = np.zeros((consts_num, eta, xi))*np.nan
cmin = np.zeros((consts_num, eta, xi))*np.nan
cang = np.zeros((consts_num, eta, xi))*np.nan
cpha = np.zeros((consts_num, eta, xi))*np.nan

k = -1
for cst in consts:
    k = k+1

    if cst in consts1: 
        fid = nc.Dataset(pth_tpxo+'hf.'+cst.lower()+'_tpxo8_atlas_30c_v1.nc', 'r')
        hRe = fid.variables['hRe'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        hIm = fid.variables['hIm'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif cst in consts2:
        fid = nc.Dataset(pth_tpxo+'hf.'+cst.lower()+'_tpxo8_atlas_6.nc', 'r')
        hRe = fid.variables['hRe'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        hIm = fid.variables['hIm'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
    fid.close()

    if cst in consts1:
        fid = nc.Dataset(pth_tpxo+'uv.'+cst.lower()+'_tpxo8_atlas_30c_v1.nc', 'r')
        uRe = fid.variables['uRe'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        uIm = fid.variables['uIm'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        vRe = fid.variables['vRe'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        vIm = fid.variables['vIm'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif cst in consts2:
        fid = nc.Dataset(pth_tpxo+'uv.'+cst.lower()+'_tpxo8_atlas_6.nc', 'r')
        uRe = fid.variables['uRe'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        uIm = fid.variables['uIm'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        vRe = fid.variables['vRe'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        vIm = fid.variables['vIm'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
    fid.close()

    hRe = hRe*1.0
    hIm = hIm*1.0
    uRe = uRe*1.0
    uIm = uIm*1.0
    vRe = vRe*1.0
    vIm = vIm*1.0

    # unit convertion
    hRe = hRe/1000.  # mm to m
    hIm = hIm/1000.  # mm to m
    uRe = uRe/10000.  # cm2s-1 to m2s-1
    uIm = uIm/10000.  # cm2s-1 to m2s-1
    vRe = vRe/10000.  # cm2s-1 to m2s-1
    vIm = vIm/10000.  # cm2s-1 to m2s-1

    # convert vertical transport (uRe, vRe) to velocity
    if cst in consts1:
        zu = src_grd.z_u
        zv = src_grd.z_v
    elif cst in consts2:
        zu = src_grd_lr.z_u
        zv = src_grd_lr.z_v

    uRe = uRe/zu  # ms-1
    uIm = uIm/zu  # ms-1
    vRe = vRe/zv  # ms-1
    vIm = vIm/zv  # ms-1

    if cst in consts1:
        # mask invalid data
        hRe = np.ma.masked_where(~src_grd.mask_t, hRe)
        hIm = np.ma.masked_where(~src_grd.mask_t, hIm)
        uRe = np.ma.masked_where(~src_grd.mask_u, uRe)
        uIm = np.ma.masked_where(~src_grd.mask_u, uIm)
        vRe = np.ma.masked_where(~src_grd.mask_v, vRe)
        vIm = np.ma.masked_where(~src_grd.mask_v, vIm)

        # remap h, u, v
        hRe = pyroms.remapping.remap(hRe, wts_file_t, spval=missing_value)
        hIm = pyroms.remapping.remap(hIm, wts_file_t, spval=missing_value)
        uRe = pyroms.remapping.remap(uRe, wts_file_u, spval=missing_value)
        uIm = pyroms.remapping.remap(uIm, wts_file_u, spval=missing_value)
        vRe = pyroms.remapping.remap(vRe, wts_file_v, spval=missing_value)
        vIm = pyroms.remapping.remap(vIm, wts_file_v, spval=missing_value)

    elif cst in consts2:
        # mask invalid data
        hRe = np.ma.masked_where(~src_grd_lr.mask_t, hRe)
        hIm = np.ma.masked_where(~src_grd_lr.mask_t, hIm)
        uRe = np.ma.masked_where(~src_grd_lr.mask_u, uRe)
        uIm = np.ma.masked_where(~src_grd_lr.mask_u, uIm)
        vRe = np.ma.masked_where(~src_grd_lr.mask_v, vRe)
        vIm = np.ma.masked_where(~src_grd_lr.mask_v, vIm)

        # remap h, u, v
        hRe = pyroms.remapping.remap(hRe, wts_file_t_lr, spval=missing_value)
        hIm = pyroms.remapping.remap(hIm, wts_file_t_lr, spval=missing_value)
        uRe = pyroms.remapping.remap(uRe, wts_file_u_lr, spval=missing_value)
        uIm = pyroms.remapping.remap(uIm, wts_file_u_lr, spval=missing_value)
        vRe = pyroms.remapping.remap(vRe, wts_file_v_lr, spval=missing_value)
        vIm = pyroms.remapping.remap(vIm, wts_file_v_lr, spval=missing_value)

    hamp[k, :, :] = (hRe**2+hIm**2)**0.5  # mm
    hpha[k, :, :] = np.arctan(-hIm/hRe)/np.pi*180.  # deg
    uamp = (uRe**2+uIm**2)**0.5  # mm
    upha = np.arctan(-uIm/uRe)/np.pi*180.  # deg
    vamp = (vRe**2+vIm**2)**0.5  # mm
    vpha = np.arctan(-vIm/vRe)/np.pi*180.  # deg

    # -------------------------------------------------------------------------
    # nodal correction
    if nodal_corr==1:
        hamp[k, :, :] = hamp[k, :, :]*F[k]
        hpha[k, :, :] = hpha[k, :, :]-U[k]
        uamp = uamp*F[k]
        upha = upha-U[k]
        vamp = vamp*F[k]
        vpha = vpha-U[k]

    # -------------------------------------------------------------------------
    # convert ap to ep
    cmax[k, :, :], ecc, cang[k, :, :], cpha[k, :, :], w = CGrid_TPXO8.tidal_ellipse.ap2ep(uamp, upha, vamp, vpha)
    cmin[k, :, :] = cmax[k, :, :]*ecc

# -------------------------------------------------------------------------
# mask land grid points
hamp[:, dst_grd.hgrid.mask_rho==0] = missing_value
hpha[:, dst_grd.hgrid.mask_rho==0] = missing_value
cmax[:, dst_grd.hgrid.mask_rho==0] = missing_value
cmin[:, dst_grd.hgrid.mask_rho==0] = missing_value
cang[:, dst_grd.hgrid.mask_rho==0] = missing_value
cpha[:, dst_grd.hgrid.mask_rho==0] = missing_value

# -------------------------------------------------------------------------
savedata = 1
if savedata == 1:
    # write tidal information to nc file
    # -------------------------------------------------------------------------
    # create nc file
    fh = nc.Dataset(out_dir + 'frc/' + out_file, 'w')
    fh.createDimension('namelen', 4)
    fh.createDimension('tide_period', consts_num)
    fh.createDimension('eta_rho', eta)
    fh.createDimension('xi_rho', xi)

    fh.history = 'Tides from TPXO8'
    import time
    fh.creation_date = time.strftime('%c')
    fh.Type = 'ROMS Tidal Forcing File'
    fh.Title = 'Forcing for' + dst_grd.name + 'domain'
    fh.grid_file = dst_grd.name
    fh.Source = 'OTPS'

    name_nc = fh.createVariable('tide_name', 'c', ('tide_period', 'namelen'))

    period_nc = fh.createVariable('tide_period', 'd', ('tide_period'), fill_value=-9999)
    period_nc.field = 'tide_period, scalar'
    period_nc.long_name = 'tidal angular period'
    period_nc.units = 'hours'

    lat_nc = fh.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'))
    lat_nc.field = 'lat_rho, scalar'
    lat_nc.long_name = 'latitude of RHO-points'
    lat_nc.units = 'degree north'

    lon_nc = fh.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'))
    lon_nc.field = 'lon_rho, scalar'
    lon_nc.long_name = 'longitude of RHO-points'
    lon_nc.units = 'degree east'

    msk_nc = fh.createVariable('mask_rho', 'd', ('eta_rho', 'xi_rho'))
    msk_nc.long_name = 'mask on RHO-points'
    msk_nc.option_0 = 'land'
    msk_nc.option_1 = 'water'

    Eamp_nc = fh.createVariable('tide_Eamp', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Eamp_nc.field = 'tide_Eamp, scalar'
    Eamp_nc.long_name = 'tidal elevation amplitude'
    Eamp_nc.units = 'meter'

    Ephase_nc = fh.createVariable('tide_Ephase', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Ephase_nc.field = 'tide_Ephase, scalar'
    Ephase_nc.long_name = 'tidal elevation phase angle'
    Ephase_nc.units = 'degrees, time of maximum elevation with respect chosen time orgin'

    Cmax_nc = fh.createVariable('tide_Cmax', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cmax_nc.field = 'tide_Cmax, scalar'
    Cmax_nc.long_name = 'maximum tidal current, ellipse semi-major axis'
    Cmax_nc.units = 'meter second-1'

    Cmin_nc = fh.createVariable('tide_Cmin', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cmin_nc.field = 'tide_Cmin, scalar'
    Cmin_nc.long_name = 'minimum tidal current, ellipse semi-minor axis'
    Cmin_nc.units = 'meter second-1'

    Cangle_nc = fh.createVariable('tide_Cangle', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cangle_nc.field = 'tide_Cangle, scalar'
    Cangle_nc.long_name = 'tidal current inclination angle'
    Cangle_nc.units = 'degrees between semi-major axis and East'

    Cphase_nc = fh.createVariable('tide_Cphase', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cphase_nc.field = 'tide_Cphase, scalar'
    Cphase_nc.long_name = 'tidal current phase angle'
    Cphase_nc.units = 'degrees, time of maximum velocity'

    # -------------------------------------------------------------------------
    # write data
    name_nc[:, :] = tide_name
    period_nc[:] = tide_period

    lat_nc[:, :] = lat
    lon_nc[:, :] = lon
    msk_nc[:, :] = msk
    Eamp_nc[:, :, :] = hamp
    Ephase_nc[:, :, :] = hpha
    Cmax_nc[:, :, :] = cmax
    Cmin_nc[:, :, :] = cmin
    Cangle_nc[:, :, :] = cang
    Cphase_nc[:, :, :] = cpha

    fh.close()
