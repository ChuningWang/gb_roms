import numpy as np
import pandas as pd
import xarray as xr
from glob import glob
from scipy.signal import buttord, butter, filtfilt
import ttide

import read_host_info
sv = read_host_info.read_host_info()
in_dir = sv['in_dir']
out_dir = sv['out_dir']
model_dir = sv['model_dir']

model = 'tmpdir_GB-CIRC/outputs/2008/'

flist = glob(model_dir+model+'*.nc')
flist = flist[0:2]

data = []

for f in flist:
    fh = xr.open_dataset(f, decode_times=False)

# dt = 24*float(fh['time'][1]-fh['time'][0])
# dt = round(dt*10)/10
# 
# # barotropic
# fh['ubar'] = fh['uraw'].mean(dim='z')
# fh['vbar'] = fh['vraw'].mean(dim='z')
# tfit_v = ttide.t_tide(np.array(fh['ubar']+1j*fh['vbar']), dt)
# fh.close()
