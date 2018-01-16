""" NOAA sea level prediction data parser. """

import datetime as dt

import numpy as np
import netCDF4 as nc

import urllib2
import xmltodict

class get_noaa_slpred(object):
    """ noaa sea level prediction data parser. """

    def __init__(self, info):
        self.info = info
        self.get_metadata()

        if 'prod' not in self.info.keys():
            self.info['prod'] = 'predictions'
        if 'units' not in self.info.keys():
            self.info['units'] = 'metric'
        if 'datum' not in self.info.keys():
            self.info['datum'] = 'STND'
        if 'tz' not in self.info.keys():
            self.info['tz'] = 'lst'
        if 'fmt' not in self.info.keys():
            self.info['fmt'] = 'csv'
        if 't0' not in self.info.keys():
            self.info['t0'] = '1900-01-01 00:00:00'
        if 'sl' not in self.info.keys():
            self.info['sl'] = 'l'

        self.info['file_name'] = self.info['stn'] + '_' + \
            str(self.info['year']) + '.nc'

        return None

    def __call__(self):
        self.make_url()
        if self.info['sl'] == 's':
            self.wget_slpred()
            self.save_data()
        else:
            self.load_data()

        return None

    def get_metadata(self):
        """ get metadata information. """

        self.info['stn_url'] = 'http://tidesandcurrents.noaa.gov/mdapi/v0.6/webapi/stations/' + \
                               self.info['stn'] + '.xml'

        page = urllib2.urlopen(self.info['stn_url'])
        s_info = xmltodict.parse(page.read())
        s_info = s_info['Stations']['Station']
        self.info['stn_name'] = s_info['name']
        self.info['lat'] = float(s_info['lat'])
        self.info['lon'] = float(s_info['lng'])
        self.info['bdate'] = str(self.info['year']) + '0101'
        self.info['edate'] = str(self.info['year'] + 1) + '0101'

        return None

    def make_url(self):
        """ make downloading url.  """

        self.url = []
        for mm in range(12):
            bdate = "%04d%02d%02d" % (self.info['year'], mm+1, 1)
            if mm == 11:
                edate = "%04d%02d%02d" % (self.info['year']+1, 1, 1)
            else:
                edate = "%04d%02d%02d" % (self.info['year'], mm+2, 1)

            self.url.append('https://tidesandcurrents.noaa.gov/api/datagetter?' +
                            'begin_date=' + bdate +
                            '&end_date=' + edate +
                            '&station=' + self.info['stn'] +
                            '&product=' + self.info['prod'] +
                            '&units=' + self.info['units'] +
                            '&datum=' + self.info['datum'] +
                            '&time_zone=' + self.info['tz'] +
                            '&application=web_services' +
                            '&format=' + self.info['fmt'])

        return None

    def wget_slpred(self):
        """ get sea level prediction data with python wget. """

        for kk, url in enumerate(self.url):
            print 'Acquiring data from ' + url + '...'
            page = urllib2.urlopen(url)
            web = page.read().split('\n')
            web.pop(0)
            web.pop(-1)
            ct = len(web)
            for i in range(ct):
                web[i] = web[i].split(',')

            time = np.array([dt.datetime.strptime(web[i][0], '%Y-%m-%d %H:%M') for i in range(ct)])
            ctime = nc.date2num(time, 'days since '+self.info['t0'])
            slpred = np.array([float(web[i][1]) for i in range(ct)])

            if kk == 0:
                self.ctime = ctime
                self.slpred = slpred
            else:
                self.ctime = np.concatenate((self.ctime, ctime))
                self.slpred = np.concatenate((self.slpred, slpred))

        self.dt_hrs = self._get_dt(self.ctime)

        return None

    def _get_dt(self, ctime):
        """ Calculate data time interval. """

        ctime_hrs = np.round(ctime*24.*10.)/10.
        val, idx = np.unique(np.diff(ctime_hrs), return_inverse=True)
        dt_hrs = val[np.argmax(np.bincount(idx))]

        return dt_hrs

    def save_data(self):
        """ save data to NetCDF file. """

        fh = nc.Dataset(self.info['file_dir'] + self.info['file_name'], 'w')

        fh.title = 'NOAA tidal prediction collection'
        fh.station = self.info['stn']
        fh.station_name = self.info['stn_name']
        fh.lat = self.info['lat']
        fh.lon = self.info['lon']
        fh.tz = self.info['tz']
        fh.station_url = self.info['stn_url']
        fh.bdate = self.info['bdate']
        fh.edate = self.info['edate']
        fh.product = self.info['prod']
        fh.units = self.info['units']
        fh.t0 = self.info['t0']

        fh.createDimension('t')

        fh.createVariable('time', 'd', ('t'))
        fh.createVariable('slpred', 'd', ('t'))

        fh.variables['time'][:] = self.ctime
        fh.variables['slpred'][:] = self.slpred
        fh.close()

        return None

    def load_data(self):
        """ load data from NetCDF file. """

        fh = nc.Dataset(self.info['file_dir'] + self.info['file_name'], 'r')
        self.ctime = fh.variables['time'][:]
        self.slpred = fh.variables['slpred'][:]
        fh.close()

        self.dt_hrs = self._get_dt(self.ctime)

        return None
