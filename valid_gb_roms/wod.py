import numpy as np
from wodpy import wod
from matplotlib import path
import matplotlib.pyplot as plt

c0 = [(-136.4, 58.05),
      (-136.6, 58.30),
      (-136.9, 58.25),
      (-136.8, 57.95)]

p = path.Path(c0)

lat = []
lon = []

fid = open('/Users/CnWang/Documents/gb_roms/CTDS7513')

n = 0  # counter
t = True
while t:
    pr = wod.WodProfile(fid)
    lon0, lat0 = pr.longitude(), pr.latitude()
    if p.contains_point((lon0, lat0)):
        lat.append(lat0)
        lon.append(lon0)
        plt.plot(pr.s(), -pr.z())
        # print pr.z()
    t = not pr.is_last_profile_in_file(fid)
    n += 1
