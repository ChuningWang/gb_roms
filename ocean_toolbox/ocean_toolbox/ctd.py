'''
standard SEABIRD CTD data wrapper
'''

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.dates as mdates
import netCDF4 as nc
from datetime import datetime
# import matplotlib.cm as cm
# from matplotlib.mlab import griddata
# from mpl_toolkits.basemap import Basemap
# from datetime import datetime, timedelta

import os
import sys
from glob import glob
import fnmatch
import re
import pdb

# ------------------------------------------------------------
class ctd(object):
    ''' main class '''

    def __init__(self, info):
        self.info = info
        if 'zlev' not in self.info.keys():
            self.info['zlev'] = 500
        self.varnames = {'salt':  ['sal00'],
                         'temp':  ['tv290C', 't090C'],
                         'pre':   ['prdM', 'prSM'],
                         'rho':   ['sigma-t00'],
                         'o2':    ['sbeox0ML/L'],
                         'fluor': ['flECO-AFL', 'wetStar'],
                         'tur':   ['turbWETntu0', 'obs', 'upoly0'],
                         'par':   ['par']}

    def __call__(self):
        # self.tidy_cnv()
        self.read_ctd()
        self.save_data()

    def tidy_cnv(self):
        ''' tidy .cnv files to get rid of some formatting issues. trashed with the new data reader. '''

        self.flist = self.recursive_glob(pattern = '*.cnv')
        for fname in self.flist:
            f = open(fname, 'r')
            text = f.read()
            f.close()
            str1 = 'number of voltages sampled'
            str2 = 'nquan'
            n1 = text.find(str1) + 32
            n2 = text.find(str2) - 3
            text = text[:n1] + text[n2:]
            f = open(fname, 'w')
            f.write(text)
            f.close()
        return None

    def read_ctd(self):
        ''' read, interpolate and combine ctd casts. '''

        self.flist = self.recursive_glob(pattern = '*.cnv')
        self.data = {}
        pr_num = len(self.flist)  # number of total entry
        # create variables
        z = np.arange(self.info['zlev'])
        self.data['z'] = z
        info_list = ['station', 'lat', 'lon', 'time']
        for i in info_list:
            self.data[i] = []
        for var in self.info['var']:
            self.data[var] = []

        for fname in self.flist:
            print('parsing '+fname)
            # read CTD cast
            info, data = self.read_cnv(fname)

            # parse header info
            for i in info_list:
                if i=='time':
                    time_str = info['date'] + ' ' + info['time']
                    ctime = nc.date2num(datetime.strptime(time_str, '%Y/%m/%d %H:%M'), 'days since 1900-01-01 00:00:00')
                    self.data[i].append(ctime)
                else:
                    self.data[i].append(info[i])

            # check data flag
            msk = data[:, -1] == 0
            data = data[msk, :]

            # find cast profile depth info
            idx = info['depSM_index']
            pr_depth = data[:, idx]

            # iterate to process data
            for var in self.info['var']:
                varlist = self.varnames[var]
                for varname in varlist:
                    if (varname + '_index') in info.keys():
                        idx = info[varname + '_index']
                pr = self.interp(pr_depth, data[:, idx])
                self.data[var].append(pr)

        for i in info_list:
            self.data[i] = np.array(self.data[i])
        for var in self.info['var']:
            self.data[var] = np.array(self.data[var]).T

    def read_cnv(self, fname):
        ''' read single CTD cast cnv file. '''

        # read raw text
        f = open(fname, 'r')
        text = f.readlines()
        f.close()

        # split header and data
        text = [x.strip() for x in text]
        text_split = text.index('*END*')
        header = text[:text_split]
        cast_data = text[text_split+1:]

        # parse header
        cast_info = dict()
        for line in header:
            if 'Station:' in line:
                cast_info['station'] = float(line.split(':')[-1])
            if 'Latitude:' in line:
                cast_info['lat'] = float(line.split(':')[-1])
            if 'Longitude:' in line:
                cast_info['lon'] = float(line.split(':')[-1])
            if 'Date GMT:' in line:
                cast_info['date'] = line.split(':')[-1].replace(' ', '')
            if 'Time GMT:' in line:
                cast_info['time'] = line.split(':')[-2] + ':' + line.split(':')[-1].replace(' ', '')
            if 'Fathometer Depth:' in line:
                cast_info['fathometer_depth'] = float(line.split(':')[-1])
            if 'Cast Target Depth:' in line:
                cast_info['cast_target_depth'] = float(line.split(':')[-1])
            # variable names
            for i in range(15):
                name = 'name ' + str(i) + ' ='
                if name in line:
                    var_info = line.split(' = ')[1].split(':')
                    cast_info[var_info[0] + '_index'] = i
                    cast_info[var_info[0] + '_unit'] = var_info[1]
                    # cast_info['var' + str(i)] = var_info[0]
                    # cast_info['unit' + str(i)] = var_info[1]

        # parse data
        cast_data = [line.split() for line in cast_data]
        cast_data = np.array(cast_data).astype(float)
        return cast_info, cast_data

    def interp(self, pr_z, pr):
        ''' interpolate ctd profiles to standard depth. '''

        z = self.data['z']
        pr_new = np.NaN*np.zeros(self.info['zlev'])
        msk = (z >= pr_z.min()) & (z <= pr_z.max())
        pr_new[msk] = interp1d(pr_z, pr)(z[msk])
        return pr_new

    def recursive_glob(self, pattern = '*'):
        ''' Search recursively for files matching a specified pattern. '''

        matches = []
        for root, dirnames, filenames in os.walk(self.info['data_dir']):
            for filename in fnmatch.filter(filenames, pattern):
                matches.append(os.path.join(root, filename))

        return matches

    def save_data(self):
        ''' save data to netCDF file for easy access later. '''

        # open new nc file
        fh = nc.Dataset(self.info['file_dir'] + self.info['file_name'], 'w')
        fh.title = 'Glacier Bay data collection'
        fh.createDimension('z', self.info['zlev'])
        fh.createDimension('cast')

        # write data
        fh.createVariable['z', 'd', ('z')]
        fh.variables['z'][:] = self.data['z']
        info_list = ['station', 'lat', 'lon', 'time']
        for i in info_list:
            fh.createVariable(i, 'd', ('cast'))
            fh.variables[i][:] = self.data[i]

        for var in self.info['var']:
            fh.createVariable(var, 'd', ('z', 'cast'))
            fh.variables[var][:] = self.data[var]

        fh.close()

    # def rd_ctd(ctd_dir):
    #     '''
    #     Read CTD data collected in Glacier Bay National Park.
    #     Usage
    #        Input
    #            ctd_dir -- path of CTD NetCDF file
    #        Output
    #             ctd    -- A python dict variable containing
    #                       - d => depth
    #                       - p => pressure
    #                       - t => temperature
    #                       - s => salinity
    #                       - o2 => oxygen concentration
    #                       - f => chlorophyll fluorescence
    #                       - tb => water turbidity
    #                       - par => insolation

    #     The NetCDF file is generated with the MATLAB script rdctd_gb.m and get_ctd.m.
    #     There is a module that can load CTD data into python workspace, however, it is
    #     not as good as the MATLAB script I generated for earlier work.
    #            Chuning Wang, 2015/05/27
    #     '''

    #     import netCDF4 as nc

    #     # read CTD data
    #     print 'Reading CTD NetCDF file...'
    #     fh = nc.Dataset(ctd_dir,mode='r')

    #     mt = fh.variables['mtime'][:]
    #     mt_vec = fh.variables['mtime_vector'][:]
    #     stn = fh.variables['station'][:]
    #     lat = fh.variables['latitude'][:]
    #     lon = fh.variables['longitude'][:]
    #     d = fh.variables['depth'][:]
    #     p = fh.variables['pressure'][:]
    #     t = fh.variables['temperature'][:]
    #     s = fh.variables['salinity'][:]
    #     rho = fh.variables['density'][:]
    #     o2 = fh.variables['oxygen'][:]
    #     f = fh.variables['fluorescence'][:]
    #     tb = fh.variables['turbidity'][:]
    #     par = fh.variables['par'][:]

    #     # mask out NaNs
    #     print 'Masking invalid data...'
    #     d = np.ma.masked_where(np.isnan(d),d)
    #     p = np.ma.masked_where(np.isnan(p),p)
    #     t = np.ma.masked_where(np.isnan(t),t)
    #     s = np.ma.masked_where(np.isnan(s),s)
    #     rho = np.ma.masked_where(np.isnan(rho),rho)
    #     o2 = np.ma.masked_where(np.isnan(o2),o2)
    #     f = np.ma.masked_where(np.isnan(f),f)
    #     tb = np.ma.masked_where(np.isnan(tb),tb)
    #     par = np.ma.masked_where(np.isnan(par),par)

    #     # convert matlab time to python datetime
    #     print 'Converting MATLAB time to Python datetime...'
    #     pyt = np.array([datetime(int(mt_vec[0,i]),
    #         int(mt_vec[1,i]),
    #         int(mt_vec[2,i]),
    #         int(mt_vec[3,i]),
    #         int(mt_vec[4,i]),
    #         int(mt_vec[5,i])
    #         ) for i in range(mt.size)])

    #     # averge Lat & Lon for each station
    #     print 'Getting station coordinates...'
    #     lat_stn = np.empty(24)
    #     lon_stn = np.empty(24)
    #     lat_stn[:] = np.nan
    #     lon_stn[:] = np.nan

    #     for i in range(int(max(stn))):
    #         lat_stn[i] = np.mean(lat[stn==i])
    #         lon_stn[i] = np.mean(lon[stn==i])

    #     ctd = { 'mt':      mt,
    #             'mt_vec':  mt_vec,
    #             'pyt':     pyt,
    #             'stn':     stn,
    #             'lat':     lat,
    #             'lon':     lon,
    #             'lat_stn': lat_stn,
    #             'lon_stn': lon_stn,
    #             'd':       d,
    #             'p':       p,
    #             't':       t,
    #             's':       s,
    #             'rho':     rho,
    #             'o2':      o2,
    #             'f':       f,
    #             'tb':      tb,
    #             'par':     par
    #             }
    #     return ctd

    # # ------------------------------------------------------------
    # def plthov_ctd(t, z, data, rho=-1, depth = 200, fig = -1, ax = -1):
    #     # Plot hovmuller diagram of Glacier Bay CTD data.
    #     # Input
    #     #   t - python datetime
    #     #   z - depth (1D)
    #     #   data - variable to plot
    #     #   depth (optional) - maximum depth of y axis
    #     #   fig (optional) - figure handle
    #     #   ax (optional) - axes handle
    #     # Output
    #     #   fig - figure handle
    #     #   ax - axes handle
    #     #      Chuning Wang 2016/05/27

    #     import matplotlib.pyplot as plt

    #     if fig==-1 and ax==-1:
    #         fig = plt.figure()
    #         ax  = fig.add_subplot(111)
    #     elif fig!=-1 and ax==-1:
    #         ax = fig.add_subplot(111)
    #         print 'Using input figure handle...'
    #     elif fig!=-1 and ax!=-1:
    #         print 'Using input figure and axes handle...'
    #     else:
    #         print 'Please specify fig when ax is specified!!!'
    #         fig = plt.gcf()

    # #     plt.pcolormesh(t,z,data,cmap=plt.cm.Blues)
    #     plt.pcolormesh(mdates.date2num(t),z,data)
    #     plt.xlabel('Time')
    #     plt.ylabel('Depth [m]')
    # #     plt.yscale('log')
    #     plt.ylim(0, depth)
    #     plt.xlim(min(t),max(t))
    #     plt.clim(np.min(data),np.max(data))
    #     plt.gca().invert_yaxis()
    #     plt.colorbar()
    #     # reset xticks
    #     years = mdates.YearLocator()
    #     ax = plt.gca()
    #     ax.xaxis.set_major_locator(years)
    #     fig.autofmt_xdate()

    #     if np.any(rho!=-1):
    #         # contour density
    #         clevs = np.arange(1020, 1030, 1)
    #         rhoc = plt.contour(mdates.date2num(t), z, rho, clevs, colors='k', linewidths=.2)
    #         clevs = clevs[::4]
    #         rhoc = plt.contour(mdates.date2num(t), z, rho, clevs, colors='k', linewidths=.4)
    #         plt.clabel(rhoc, fontsize=3)

    #     # plt.show(block=False)
    #     return fig, ax

    # # ------------------------------------------------------------
    # def cal_ctd_climatology(ctd):
    #     # Calculate climatology profiles of each station with Glacier Bay CTD data.
    #     # QC
    #     msk = (ctd['d'] > 180) & (ctd['s'] < 30)

    #     z = ctd['d'][:, 0]
    #     zz = z.size
    #     # Initiate
    #     p   = np.empty([12,24,zz]); p[:]   = np.NaN
    #     t   = np.empty([12,24,zz]); t[:]   = np.NaN
    #     s   = np.empty([12,24,zz]); s[:]   = np.NaN
    #     o2  = np.empty([12,24,zz]); o2[:]  = np.NaN
    #     f   = np.empty([12,24,zz]); f[:]   = np.NaN
    #     tb  = np.empty([12,24,zz]); tb[:]  = np.NaN
    #     par = np.empty([12,24,zz]); par[:] = np.NaN

    #     mm = ctd['mt_vec'][1,:]
    #     # Loop through stations & months
    #     for i in range(12):
    #         for j in range(24):
    #             msk = (mm == i) & (ctd['stn'] == j)
    #             p[i, j, :]   = np.nanmean(ctd['p'][:, msk],   axis=1)
    #             t[i, j, :]   = np.nanmean(ctd['t'][:, msk],   axis=1)
    #             s[i, j, :]   = np.nanmean(ctd['s'][:, msk],   axis=1)
    #             o2[i, j, :]  = np.nanmean(ctd['o2'][:, msk],  axis=1)
    #             f[i, j, :]   = np.nanmean(ctd['f'][:, msk],   axis=1)
    #             tb[i, j, :]  = np.nanmean(ctd['tb'][:, msk],  axis=1)
    #             par[i, j, :] = np.nanmean(ctd['par'][:, msk], axis=1)

    #     ctd_climatology = { 'd':   z,
    #                         'p':   p,
    #                         't':   t,
    #                         's':   s,
    #                         'o2':  o2,
    #                         'f':   f,
    #                         'tb':  tb,
    #                         'par': par,
    #                         'lat': ctd['lat_stn'],
    #                         'lon': ctd['lon_stn']
    #                         }
    #     return ctd_climatology

    # # ------------------------------------------------------------
    # def pltscatter_ctd(t, z, data, rho=-1, clim = 'auto', depth = 200, fig=-1, ax=-1):
    #     # Similar to plthov_ctd

    #     import matplotlib.pyplot as plt

    #     if fig==-1 and ax==-1:
    #         fig = plt.figure()
    #         ax  = fig.add_subplot(111)
    #     elif fig!=-1 and ax==-1:
    #         ax = fig.add_subplot(111)
    #         print 'Using input figure handle...'
    #     elif fig!=-1 and ax!=-1:
    #         print 'Using input figure and axes handle...'
    #     else:
    #         print 'Please specify fig when ax is specified!!!'
    #         fig = plt.gcf()

    #     ax.set_axis_bgcolor((.85, .85, .85))
    #     if clim == 'auto':
    #         cmax = np.max(data)
    #         cmin = np.min(data)
    #     else:
    #         cmax = clim[1]
    #         cmin = clim[0]
    #     data[data>cmax] = cmax
    #     data[data<cmin] = cmin
    #     steps = 20
    #     dc = (cmax-cmin)/steps

    #     mt = mdates.date2num(t)
    #     # for i in range(steps):
    #     #     ii = np.where((data>=cmin+i*dc) & (data<=cmin+(i+1)*dc))
    #     #     plt.plot(t[ii[1]],z[ii[0]],'s',c=cmap(i*255/steps),mec='none',ms=3)

    #     #     # for j in range(np.shape(ii)[1]):
    #     #     #     pp = patches.Rectangle((mdates.date2num(t[ii[1][j]])-3,z[ii[0][j]]),
    #     #     #                            7,
    #     #     #                            1,
    #     #     #                            facecolor=cmap(i*255/steps),
    #     #     #                            edgecolor='none'
    #     #     #          )
    #     #     #     ax.add_patch(pp)
    #     tt, zz = np.meshgrid(mt,z)
    #     tt = tt.flatten()
    #     zz = zz.flatten()
    #     dd = data.flatten()
    #     plt.scatter(tt, zz, s=8,
    #             c=dd, vmin=cmin, vmax=cmax, cmap=cm.get_cmap('RdBu_r'),
    #             marker='s', edgecolor='none')
    #     plt.ylim(0, depth)
    #     plt.xlim(min(t),max(t))
    #     plt.gca().invert_yaxis() 
    #     plt.clim(cmin,cmax)
    #     plt.colorbar()
    #     # reset xticks
    #     years = mdates.YearLocator()
    #     ax.xaxis.set_major_locator(years)
    #     fig.autofmt_xdate()

    #     if np.any(rho!=-1):
    #         # contour density
    #         clevs = np.arange(1020, 1030, 1)
    #         rhoc = plt.contour(mdates.date2num(t), z, rho, clevs, colors='w', linewidths=.4)
    #         clevs = clevs[::4]
    #         rhoc = plt.contour(mdates.date2num(t), z, rho, clevs, colors='w', linewidths=.8)
    #         plt.clabel(rhoc, fontsize=5)

    #     plt.show(block=False)

    # # ------------------------------------------------------------
    # def plttrans_ctd(lat_stn, lon_stn, z, data, clim = 'auto', depth = 200):

    #     import matplotlib.pyplot as plt
    #     from geopy.distance import vincenty

    #     if clim == 'auto':
    #         cmax = np.max(data)
    #         cmin = np.min(data)
    #     else:
    #         cmax = clim[1]
    #         cmin = clim[0]

    #     stn_tr1 = range(13)+range(21,22)
    #     dis1 = np.zeros(lat_stn[stn_tr1].size)
    #     for i in range(1,lat_stn[stn_tr1].size):
    #         dis1[i] = vincenty(
    #                 (lat_stn[stn_tr1[i-1]], lon_stn[stn_tr1[i-1]]),
    #                 (lat_stn[stn_tr1[i]],   lon_stn[stn_tr1[i]])
    #                 ).meters
    #         dis1 = np.cumsum(dis1)
    #     dis1 = dis1/1000 # [km]

    #     stn_tr2 = range(5)+range(13,21)
    #     dis2 = np.zeros(lat_stn[stn_tr2].size)
    #     for i in range(1,lat_stn[stn_tr2].size):
    #         dis2[i] = vincenty(
    #                 (lat_stn[stn_tr2[i-1]], lon_stn[stn_tr2[i-1]]),
    #                 (lat_stn[stn_tr2[i]],   lon_stn[stn_tr2[i]])
    #                 ).meters
    #         dis2 = np.cumsum(dis2)
    #     dis2 = dis2/1000 # [km]

    #     data1 = data[:, stn_tr1]
    #     data2 = data[:, stn_tr2]

    #     if clim == 'auto':
    #         cmax = np.max(data)
    #         cmin = np.min(data)
    #     else:
    #         cmax = clim[1]
    #         cmin = clim[0]

    #     fig, (ax1, ax2) = plt.subplots(2)
    #     # fig = plt.figure()
    #     # ax0 = fig.add_subplot(211)

    #     plt.sca(ax1)
    #     plt.pcolormesh(dis1,z,data1)
    #     # plt.xlabel('Distance [m]')
    #     plt.ylabel('Depth [m]')
    #     # plt.yscale('log')
    #     plt.ylim(0, depth)
    #     plt.xlim(0, 120)
    #     plt.clim(cmin,cmax)
    #     plt.gca().invert_yaxis()
    #     plt.gca().xaxis.tick_top()
    #     plt.xticks(dis1,["%02d"%i for i in stn_tr1],rotation=90)
    #     plt.grid('on')

    #     for i in range(len(dis1)-1):
    #         pr = (data1[:, i]-np.min(data1[:, i]))/(np.max(data1[:, i])-np.min(data1[:, i]))*(dis1[i+1]-dis1[i])
    #         plt.plot(pr+dis1[i], z, 'k', lw=2)

    #     # ax1 = fig.add_subplot(212)

    #     plt.sca(ax2)
    #     pcl = plt.pcolormesh(dis2,z,data2)
    #     plt.xlabel('Station Along Channel')
    #     plt.ylabel('Depth [m]')
    #     # plt.yscale('log')
    #     plt.ylim(0, depth)
    #     plt.xlim(0, 120)
    #     plt.clim(cmin,cmax)
    #     plt.gca().invert_yaxis()
    #     plt.xticks(dis2,["%02d"%i for i in stn_tr2],rotation=90)
    #     plt.grid('on')

    #     for i in range(len(dis2)-1):
    #         pr = (data2[:, i]-np.min(data2[:, i]))/(np.max(data2[:, i])-np.min(data2[:, i]))*(dis2[i+1]-dis2[i])
    #         plt.plot(pr+dis2[i], z, 'k', lw=2)

    #     # Plot ruler
    #     plt.plot(np.array([105, 115]),depth*0.73*np.ones(2),'k',lw=2)
    #     plt.text(110,depth*0.7,'10 km',ha='center')

    #     # Fig Settings
    #     fig.subplots_adjust(hspace = 0)
    #     # Add Colorbar
    #     cbar_cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
    #     fig.colorbar(pcl, cax = cbar_cax)
    #     plt.show(block=False)

    #     return dis1, dis2, stn_tr1, stn_tr2, ax1, ax2
    # # ------------------------------------------------------------
    # def get_cruise(ctd):
    #     # Find indices for each cruise
    #     dt = np.diff(np.floor(mdates.date2num(ctd['pyt'])))
    #     dt = np.hstack([0, dt])
    #     k1 = np.squeeze(np.where(dt>=7))
    #     k1 = np.hstack([0, k1])
    #     k2 = np.hstack([k1[1:]-1, np.size(dt)+1])
    #     k3 = np.squeeze(np.where((k2-k1)>15))
    #     return k1, k2, k3
    # # ------------------------------------------------------------
    # def get_topo(ctd):
    #     # Find bottom topography
    #     topo = np.zeros(np.int(np.max(ctd['stn'])))
    #     dd = ctd['d']
    #     dd.mask = ctd['t'].mask
    #     for i in range(np.int(np.max(ctd['stn']))): 
    #         topo[i] = np.max(dd[:,ctd['stn']==i])
    #     del dd
    #     return topo
    # # ------------------------------------------------------------
    # def pltmapview_ctd(lat_stn, lon_stn, data,
    #         d_avg=15, clevs='auto', cborient='horizontal', method='scatter',
    #         fig=-1, ax=-1):

    #     import matplotlib.pyplot as plt

    #     # average over the top d_avg meters
    #     if np.shape(data.shape)[0]==2:
    #         data = np.nanmean(data[:d_avg,:],axis=0)
    #     msk = ~np.isnan(data)

    #     k = np.any(msk)
    #     if k:
    #         lon_stn2 = lon_stn[msk]
    #         lat_stn2 = lat_stn[msk]
    #         data = data[msk]

    #         if clevs == 'auto':
    #             clevs = np.linspace(np.min(data),np.max(data),11)

    #         if fig==-1 and ax==-1:
    #             fig = plt.figure()
    #             ax  = fig.add_subplot(111)
    #         elif fig!=-1 and ax==-1:
    #             ax = fig.add_subplot(111)
    #             print 'Using input figure handle...'
    #         elif fig!=-1 and ax!=-1:
    #             print 'Using input figure and axes handle...'
    #         else:
    #             print 'Please specify fig when ax is specified!!!'
    #             fig = plt.gcf()

    #         m = Basemap(llcrnrlon=-137.5,llcrnrlat=58,urcrnrlon=-135.5,urcrnrlat=59.25,
    #                 projection='stere',lat_0=58,lon_0=-137.5,
    #                 resolution='h'
    #                 )

    #         m.drawcoastlines(linewidth=.3)
    #         # m.fillcontinents(color='lightgrey')

    #         m.drawmeridians(np.arange(-137.5,-135.5,0.5), linewidth=.3, labels=[0,0,0,1], fontsize=6)
    #         m.drawparallels(np.arange(58,59.25,0.25), linewidth=.3, labels=[1,0,0,0], fontsize=6)

    #         if method=='contourf':
    #             numcols, numrows = 100, 100
    #             xi = np.linspace(lon_stn2.min(), lon_stn2.max(), numcols)
    #             yi = np.linspace(lat_stn2.min(), lat_stn2.max(), numrows)
    #             xi, yi = np.meshgrid(xi, yi)

    #             zi = griddata(lon_stn2, lat_stn2, data, xi, yi)
    #             xi, yi = m(xi,yi)

    #             ci = m.contourf(xi,yi,zi,clevs)

    #             # draw stations
    #             x, y = m(lon_stn[msk],lat_stn[msk])
    #             m.plot(x,y,'ok',markersize=3)
    #             x, y = m(lon_stn[~msk],lat_stn[~msk])
    #             m.plot(x,y,'o',color='grey',markersize=3)

    #         elif method=='scatter':
    #             x, y = m(lon_stn[msk], lat_stn[msk])
    #             ci = m.scatter(x, y, s=30,
    #                     c=data, vmin=clevs[0], vmax=clevs[-1],
    #                     marker='o', edgecolor='none')
    #             x, y = m(lon_stn[~msk],lat_stn[~msk])
    #             m.scatter(x, y, s=30, marker='o', color='k', edgecolor='none')

    #         # add colorbar
    #         cbar = fig.colorbar(ci,orientation=cborient)

    #         cbar.ax.tick_params(labelsize=8) 
    #         # cbar.ax.update_ticks()

    #         plt.show(block=False)

    #         return k
