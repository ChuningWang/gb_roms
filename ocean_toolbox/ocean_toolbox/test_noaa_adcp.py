from ocean_toolbox import noaa_adcp

stn = 'SEA1009'
info = {'stn' : stn,
        'file_dir': '~/Documents/gb_roms/NOAA_ADCP/',
        'sl': 's',
        'Wp_hrs': 2}

crt = noaa_adcp.get_noaa_current(info)
crt()
