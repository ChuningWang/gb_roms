from ocean_toolbox import noaa_slpred

my_year = 2008
# my_month = 01
# stn_list = ['9452294']

# for stn in stn_list:

info = {'stn': '9452294',
        'year': my_year,
        # 'month' : my_month,
        'file_dir': '/Users/CnWang/Documents/gb_roms/NOAA_slpred/',
        'sl': 'l'}

slp = noaa_slpred.get_noaa_slpred(info)
slp()
