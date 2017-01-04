import sys
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

# import gb_toolbox

vlevel = 3

fh = nc.Dataset('../../data/ctd.nc', 'r')
mt = fh.variables['mtime'][:]
stn = fh.variables['station'][:]
t = fh.variables['temperature'][:]
s = fh.variables['salinity'][:]
fh.close()

msk = stn==12
mt = mt[msk]
t = t[:, msk]
s = s[:, msk]

pyt = [datetime.fromordinal(int(i)) + timedelta(days=i%1) - timedelta(days=366) for i in mt]
yday = np.array([i.timetuple().tm_yday for i in pyt])

# lowess smoothing
yday = np.concatenate((yday-366, yday, yday+366))
t = np.concatenate((t, t, t), axis=1)
s = np.concatenate((s, s, s), axis=1)

yd = range(366)
tt = np.zeros((vlevel, 366))
ss = np.zeros((vlevel, 366))

for i in range(vlevel):
    ll = lowess(t[i, :], yday, frac=0.05)
    tt[i, :] = interp1d(ll[:, 0], ll[:, 1])(yd)
    ll = lowess(s[i, :], yday, frac=0.05)
    ss[i, :] = interp1d(ll[:, 0], ll[:, 1])(yd)

# lowess smoothing again
for i in range(vlevel):
    ll = lowess(np.concatenate((t[i, :], t[i, :], t[i, :])), np.concatenate((yday-366, yday, yday+366)), frac=0.1)
    tt[i, :] = interp1d(ll[:, 0], ll[:, 1])(yd)
    ll = lowess(np.concatenate((s[i, :], s[i, :], s[i, :])), np.concatenate((yday-366, yday, yday+366)), frac=0.1)
    ss[i, :] = interp1d(ll[:, 0], ll[:, 1])(yd)

ttime = yd
temp = tt.mean(axis=0)
salt = np.zeros(366)

mkplt = 1
if mkplt==1:
    import matplotlib.pyplot as plt
    # plt.plot(yday, t[0:5, :].T, '.')
    # plt.plot(ll[:, 0], ll[:, 1], '.k')
    plt.plot(yd, tt.T, '.')
    plt.show(block=False)

savefile = 1
if savefile==1:

    outfile = sys.argv[1]

    # create file with all the objects
    out = nc.Dataset(outfile, 'a', format='NETCDF3_64BIT')

    out.createDimension('river_tracer_time', len(ttime))

    times = out.createVariable('river_tracer_time', 'f8', ('river_tracer_time'))
    times.units = 'day'
    times.cycle_length = 365.25
    times.long_name = 'river tracer time'

    temp = out.createVariable('river_temp', 'f8', ('river_tracer_time'))
    temp.long_name = 'river runoff potential temperature'
    temp.units = 'Celsius'
    temp.time = 'river_tracer_time'

    salt = out.createVariable('river_salt', 'f8', ('river_tracer_time'))
    salt.long_name = 'river runoff salinity'
    salt.time = 'river_tracer_time'

    out.variables['river_tracer_time'][:] = ttime
    out.variables['river_temp'][:] = temp
    out.variables['river_salt'][:] = salt

    out.close()
