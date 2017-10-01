import numpy as np
from wodpy import wod
from scipy.interpolate import interp1d
import netCDF4 as nc
from matplotlib import path
import matplotlib.pyplot as plt

c0 = [(-134.90, 57.95),
      (-135.50, 57.95),
      (-135.50, 58.50),
      (-134.90, 58.50)]

# c0 = [(-135.6, 58.10),
#       (-136.8, 58.25),
#       (-136.1, 58.40),
#       (-135.7, 58.25),
#       (-135.6, 58.45),
#       (-135.0, 58.35),
#       (-135.0, 58.05)]

# c0 = [(-137.50, 57.75),
#       (-137.50, 59.25),
#       (-134.50, 59.25),
#       (-134.50, 57.25)]

p = path.Path(c0)

time = []
yearday = []
lat = []
lon = []
salt = []
temp = []
zlev = 500
z = np.arange(zlev)

fid = open('/Users/CnWang/Documents/gb_roms/CTDS7513')

n = 0  # counter
t = True
while t:
    pr = wod.WodProfile(fid)
    lon0, lat0 = pr.longitude(), pr.latitude()
    if p.contains_point((lon0, lat0)):
        zpr = pr.z()
        if len(zpr)>1:
            zmin, zmax = zpr[0], zpr[-1]
            zmsk = (z>=zmin) & (z<=zmax)
            spr = np.zeros(zlev)*np.NaN
            spr[zmsk] = interp1d(pr.z(), pr.s())(z[zmsk])
            tpr = np.zeros(zlev)*np.NaN
            tpr[zmsk] = interp1d(pr.z(), pr.t())(z[zmsk])
            salt.append(spr)
            temp.append(tpr)
            time.append(nc.date2num(pr.datetime(), 'days since 1900-01-01'))
            yearday.append(pr.datetime().timetuple().tm_yday)
            lat.append(lat0)
            lon.append(lon0)
    t = not pr.is_last_profile_in_file(fid)
    n += 1

salt = np.array(salt)
temp = np.array(temp)
lon = np.array(lon)
lat = np.array(lat)
time = np.array(time)
yearday = np.array(yearday)

t = np.array([-46, 46, 138, 229, 321, 412])
season_msk = np.zeros((len(time), 6))
season_msk[:, 0] = yearday>275
season_msk[:, 1] = yearday<=92
season_msk[:, 2] = (yearday>92) & (yearday<=183)
season_msk[:, 3] = (yearday>183) & (yearday<=275)
season_msk[:, 4] = yearday>275
season_msk[:, 5] = yearday<=92

salt_season = np.zeros((6, zlev))
temp_season = np.zeros((6, zlev))

for i in range(6):
    salt_season[i, :] = np.nanmean(salt[season_msk[:, i]==1, :], axis=0)
    temp_season[i, :] = np.nanmean(temp[season_msk[:, i]==1, :], axis=0)

# process the profile
salt_season[:, z > 250] = np.NaN
temp_season[:, z > 250] = np.NaN
