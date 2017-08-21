import csv, glob
import numpy as np
import netCDF4 as nc
from datetime import datetime

def read_station(stn_file):
    fh = open(stn_file, 'rb')
    reader = csv.reader(fh, delimiter=',')
    rn = 0

    data = []

    for row in reader:
        # Save header row.
        if rn == 0:
            header = row
        else:
            data.append(row)
        rn += 1
                  
    fh.close()

    data = np.array(data)
    stn = data[0, 0]
    name = data[0, 1]
    lat = float(data[0, 2])
    lon = float(data[0, 3])
    lev = data[0, 4]
    date = data[:, 5]
    ctime = np.array([nc.date2num(datetime.strptime(date[i], '%Y-%m-%d'), 'days since 1900-01-01 00:00:00') for i in range(len(date))])
    awnd = data[:, 6]
    awnd = awnd.astype(np.float)
    prcp = data[:, 7]
    prcp = prcp.astype(np.float)
    snow = data[:, 8]
    snow = snow.astype(np.float)
    # tair = data[:, 9]
    # tair = tair.astype(np.float)
    tmax = data[:, 10]
    tmax = tmax.astype(np.float)
    tmin = data[:, 11]
    tmin = tmin.astype(np.float)
    data = {'stn': stn, 'name': name, 'lat': lat, 'lon': lon, 'lev': lev, 
            'ctime': ctime, 'awnd': awnd, 'prcp': prcp, 'snow': snow, 'tmax': tmax, 'tmin': tmin}
    return data
