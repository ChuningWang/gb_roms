from ocean_toolbox import noaa_adcp

stn_list = ['SEA0845', 'SEA0846', 'SEA0847', 'SEA0848', 'SEA0849', 'SEA0850']

for stn in stn_list:
    info = {'stn' : stn,
            'file_dir': '/glade/p/work/chuning/data/NOAA_ADCP/',
            'sl': 's',
            'Wp_hrs': 2}

    crt = noaa_adcp.get_noaa_current(info)
    crt()
