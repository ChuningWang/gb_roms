import numpy as np
from gb_toolbox import gb_current

stn_list = ['SEA0845', 'SEA0846', 'SEA0847', 'SEA0848', 'SEA0849', 'SEA0850']
bdate_list = ['20080810']
edate_list = ['20080910']

# stn_list = ['SEA1008', 'SEA1009', 'SEA1010']
# bdate_list = ['20100525', '20100625', '20100725']
# edate_list = ['20100625', '20100725', '20100815']

for stn in stn_list:
    for i in range(len(bdate_list)):
        bdate = bdate_list[i]
        edate = edate_list[i]
        info = {'stn' : stn,
                'bdate' : bdate,
                'edate' : edate,
                'filename': '/glade/p/work/chuning/data/NOAA_ADCP/'+stn+'_'+bdate+'_'+edate+'.nc',
                'sl': 's',
                'Wp_hrs': 2}

        crt = gb_current.get_noaa_current(info)
        crt()
