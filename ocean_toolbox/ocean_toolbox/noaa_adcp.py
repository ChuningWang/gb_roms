import datetime as dt
import numpy as np
import urllib2
import xmltodict
from scipy.interpolate import interp1d
from scipy.signal import buttord, butter, filtfilt
import netCDF4 as nc

class get_noaa_current(object):
    def __init__(self, info):
        self.info = info
        self.get_metadata()
        if 'bdate' not in self.info.keys():
            bdate = self.info['deployed'][:4] + \
                    self.info['deployed'][5:7] + \
                    self.info['deployed'][8:10]
            self.info['bdate'] = bdate
        if 'edate' not in self.info.keys():
            edate = self.info['retrieved'][:4] + \
                    self.info['retrieved'][5:7] + \
                    self.info['retrieved'][8:10]
            edate = dt.datetime.strptime(edate, '%Y%m%d') + dt.timedelta(1)
            edate = edate.strftime('%Y%m%d')
            self.info['edate'] = edate
        if 'prod' not in self.info.keys():
            self.info['prod'] = 'currents'
        if 'units' not in self.info.keys():
            self.info['units'] = 'metric'
        if 'tz' not in self.info.keys():
            self.info['tz'] = 'gmt'
        if 'fmt' not in self.info.keys():
            self.info['fmt'] = 'csv'
        if 't0' not in self.info.keys():
            self.info['t0'] = '1900-01-01 00:00:00'
        if 'sl' not in self.info.keys():
            self.info['sl'] = 'l'
        if 'Wp_hrs' not in self.info.keys():
            self.info['Wp_hrs'] = -1
        self.info['file_name'] = self.info['stn'] + '_' + self.info['bdate'] + \
                                                    '_' + self.info['edate'] + '.nc'

    def __call__(self):
        print 'Formating download urls'
        self.time_splitter()
        self.make_url()
        if self.info['sl'] == 's':
            print 'Acquiring data from NOAA server'
            self.get_current()
        if self.info['sl'] == 'l':
            print 'Loading data locally'
            self.load_data()
        self.compute_uv()
        if self.info['Wp_hrs'] > 0:
            self.filter()
        if self.info['sl'] == 's':
            self.save_data()

        return None

    def get_metadata(self):
        self.info['stn_url'] = 'http://tidesandcurrents.noaa.gov/mdapi/v0.6/webapi/stations/' + \
                               self.info['stn'] + '.xml'

        page = urllib2.urlopen(self.info['stn_url'])
        s_info = xmltodict.parse(page.read())
        s_info = s_info['Stations']['Station']
        self.info['stn_name'] = s_info['name']
        self.info['lat'] = float(s_info['lat'])
        self.info['lon'] = float(s_info['lng'])
        self.info['deployed'] = s_info['deployed']
        self.info['retrieved'] = s_info['retrieved']
        self.info['tz_offset'] = float(s_info['timezone_offset'])
        self.info['dply_url'] = s_info['deployments']['@self']
        self.info['bins_url'] = s_info['bins']['@self']

        page = urllib2.urlopen(self.info['dply_url'])
        d_info = xmltodict.parse(page.read())
        d_info = d_info['Deployments']
        self.info['depth'] = float(d_info['depth'])

        page = urllib2.urlopen(self.info['bins_url'])
        b_info = xmltodict.parse(page.read())
        b_info = b_info['Bins']['Bin']
        self.info['zlevs'] = len(b_info)
        z = np.array([float(b_info[i]['depth']) for i in range(self.info['zlevs'])])
        self.z = z

    def time_splitter(self):
        tseg_length = 30.
        tb = dt.datetime.strptime(self.info['bdate'], '%Y%m%d')
        te = dt.datetime.strptime(self.info['edate'], '%Y%m%d')
        tdelta = te-tb
        tslice = int(np.ceil(tdelta.days/tseg_length))
        self.bdate_list = []
        self.edate_list = []
        for i in range(tslice):
            tbi = dt.datetime.strftime(tb+i*dt.timedelta(30), '%Y%m%d')
            if i == tslice-1:
                tei = self.info['edate']
            else:
                tei = dt.datetime.strftime(tb+(i+1)*dt.timedelta(tseg_length)-dt.timedelta(1), \
                                           '%Y%m%d')

            self.bdate_list.append(tbi)
            self.edate_list.append(tei)

        return None

    def make_url(self):
        self.url = [['https://tidesandcurrents.noaa.gov/api/datagetter?' + \
                     'begin_date=' + self.bdate_list[tslice] + \
                     '&end_date=' + self.edate_list[tslice] + \
                     '&station=' + self.info['stn'] + \
                     '&product=' + self.info['prod'] + \
                     '&bin=' + str(bin) + \
                     '&units=' + self.info['units'] + \
                     '&time_zone=' + self.info['tz'] + \
                     '&application=web_services' + \
                     '&format=' + self.info['fmt'] for tslice in range(len(self.bdate_list))] \
                                                   for bin in range(1, self.info['zlevs']+1)]

        return None

    def wget_current(self, idx, tslice):
        page = urllib2.urlopen(self.url[idx][tslice])
        web = page.read().split('\n')
        web.pop(0)
        web.pop(-1)
        ct = len(web)
        for i in range(ct):
            web[i] = web[i].split(',')

        time = np.array([dt.datetime.strptime(web[i][0], '%Y-%m-%d %H:%M') for i in range(ct)])
        speed = np.array([float(web[i][1]) for i in range(ct)])
        dir = np.array([float(web[i][2]) for i in range(ct)])

        # unit conversion
        speed = speed/100  # m s-1

        return time, speed, dir

    def get_dt(self, ctime):
        '''
        Calculate data time interval.
        '''
        ctime_hrs = np.round(ctime*24.*10.)/10.
        val, idx = np.unique(np.diff(ctime_hrs), return_inverse=True)
        dt_hrs = val[np.argmax(np.bincount(idx))]
        return dt_hrs

    def get_current(self):
        for idx in range(len(self.url)):
            print 'Acquiring data from ' + self.url[idx][0] + '...'
            for tslice in range(len(self.url[idx])):
                time0, speed0, dir0 = self.wget_current(idx, tslice)
                ctime0 = nc.date2num(time0, 'days since '+self.info['t0'])
                if tslice == 0:
                    ctime = ctime0.copy()
                    speed = speed0.copy()
                    dir = dir0.copy()
                else:
                    ctime = np.concatenate((ctime, ctime0))
                    speed = np.concatenate((speed, speed0))
                    dir = np.concatenate((dir, dir0))

            if idx == 0:
                # interpolate onto a new time grid
                # remove the first and last few records to make interpolation easier
                self.dt_hrs = self.get_dt(ctime)
                self.ctime = np.arange(ctime[5], ctime[-5], self.dt_hrs/24.)

            speed = interp1d(ctime, speed)(self.ctime)
            dir = interp1d(ctime, dir)(self.ctime)

            if idx == 0:
                self.speed = speed
                self.dir = dir
            else:
                self.speed = np.vstack((self.speed, speed))
                self.dir = np.vstack((self.dir, dir))

        return None

    def compute_uv(self):
        self.u = self.speed*np.sin((self.dir-0.)*(np.pi/180))
        self.v = self.speed*np.cos((self.dir-0.)*(np.pi/180))

        self.uraw = self.u.copy()
        self.vraw = self.v.copy()
        self.speedraw = self.speed.copy()
        self.dirraw = self.dir.copy()

        return None

    def filter(self):

        wp = int(self.info['Wp_hrs']/self.dt_hrs)
        b = np.ones(wp)/float(wp)
        a = 1

        self.u = filtfilt(b, a, self.u)
        self.v = filtfilt(b, a, self.v)
        U = self.u + 1j*self.v
        self.speed = np.abs(U)
        self.dir = 90-np.angle(U)*180/np.pi

        return None

    def save_data(self):
        fh = nc.Dataset(self.info['file_dir'] + self.info['file_name'], 'w')

        fh.title = 'NOAA ADCP data collection'
        fh.station = self.info['stn']
        fh.station_name = self.info['stn_name']
        fh.lat = self.info['lat']
        fh.lon = self.info['lon']
        fh.deploy_time = self.info['deployed']
        fh.retrieve_time = self.info['retrieved']
        fh.tz = self.info['tz']
        fh.tz_offset = self.info['tz_offset']
        fh.station_url = self.info['stn_url']
        fh.dply_url = self.info['dply_url']
        fh.bins_url = self.info['bins_url']
        fh.bdate = self.info['bdate']
        fh.edate = self.info['edate']
        fh.depth = self.info['depth']
        fh.product = self.info['prod']
        fh.units = self.info['units']
        fh.t0 = self.info['t0']
        fh.bins = self.info['zlevs']

        fh.createDimension('z', self.info['zlevs'])
        fh.createDimension('t')

        fh.createVariable('z', 'd', ('z'))
        fh.createVariable('time', 'd', ('t'))
        fh.createVariable('speedraw', 'd', ('z', 't'))
        fh.createVariable('dirraw', 'd', ('z', 't'))
        fh.createVariable('uraw', 'd', ('z', 't'))
        fh.createVariable('vraw', 'd', ('z', 't'))
        fh.createVariable('speed', 'd', ('z', 't'))
        fh.createVariable('dir', 'd', ('z', 't'))
        fh.createVariable('u', 'd', ('z', 't'))
        fh.createVariable('v', 'd', ('z', 't'))

        fh.variables['z'][:] = self.z
        fh.variables['time'][:] = self.ctime
        fh.variables['speedraw'][:] = self.speedraw
        fh.variables['dirraw'][:] = self.dirraw
        fh.variables['uraw'][:] = self.uraw
        fh.variables['vraw'][:] = self.vraw
        fh.variables['speed'][:] = self.speed
        fh.variables['dir'][:] = self.dir
        fh.variables['u'][:] = self.u
        fh.variables['v'][:] = self.v
        fh.close()

        return None

    def load_data(self):
        fh = nc.Dataset(self.info['file_dir'] + self.info['file_name'], 'r')
        self.z = fh.variables['z'][:]
        self.ctime = fh.variables['time'][:]
        self.speed = fh.variables['speedraw'][:]
        self.dir = fh.variables['dirraw'][:]
        fh.close()

        self.dt_hrs = self.get_dt(self.ctime)

        return None
