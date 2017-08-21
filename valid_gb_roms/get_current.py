import numpy as np
from gb_toolbox import gb_current

# stn_list = ['SEA1002', 'SEA1003', 'SEA1004' ,'SEA1008', 'SEA1009', 'SEA1010']
# bdate = '20100725'
# edate = '20100815'

stn_list = ['SEA0845', 'SEA0846', 'SEA0847', 'SEA0848', 'SEA0849', 'SEA0850']
bdate = '20080810'
edate = '20080910'

for stn in stn_list:
    info = {'stn' : stn,
            'bdate' : bdate,
            'edate' : edate,
            'filename': '/glade/p/work/chuning/data/NOAA_ADCP/'+stn+'_'+bdate+'_'+edate+'.nc',
            'sl': 's',
            'Wp_hrs': 2}

    crt = gb_current.get_noaa_current(info)
    crt()
