from ocean_toolbox import noaa_adcp

stn_list = ['SEA0839']

for stn in stn_list:
    info = {'stn' : stn,
            'file_dir': '/glade/p/work/chuning/data/NOAA_ADCP/',
            'sl': 's',
            'Wp_hrs': 2}

    crt = noaa_adcp.get_noaa_current(info)
    crt()
