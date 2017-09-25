import numpy as np

from ocean_toolbox import ctd
info = {'data_dir': './'}
info = {'data_dir': '/Users/CnWang/Downloads/ctd/',
        'var': ['salt', 'temp', 'o2', 'rho', 'pre', 'fluor', 'tur', 'par'],
        'file_dir': './',
        'file_name': 'ctd.nc'}
c = ctd.ctd(info)
c()

# # read raw text
# f = open('/Users/CnWang/Downloads/ctd/cnv_2016/1610_5_0254_02_c_f_a_m_l_d_B.cnv', 'r')
# text = f.readlines()
# f.close()
# 
# # split header and data
# text = [x.strip() for x in text]
# text_split = text.index('*END*')
# header = text[:text_split]
# data = text[text_split+1:]
# 
# # process header
# cast_info = dict()
# for line in header:
#     if 'Station:' in line:
#         cast_info['station'] = float(line.split(':')[-1])
#     if 'Latitude:' in line:
#         cast_info['lat'] = float(line.split(':')[-1])
#     if 'Longitude:' in line:
#         cast_info['lon'] = float(line.split(':')[-1])
#     if 'Date GMT:' in line:
#         cast_info['date'] = line.split(':')[-1]
#     if 'Time GMT:' in line:
#         cast_info['time'] = line.split(':')[-1]
#     if 'Fathometer Depth:' in line:
#         cast_info['fathometer_depth'] = float(line.split(':')[-1])
#     if 'Cast Target Depth:' in line:
#         cast_info['cast_target_depth'] = float(line.split(':')[-1])
#     # variable names
#     for i in range(15):
#         name = 'name ' + str(i) + ' ='
#         if name in line:
#             var_info = line.split(' = ')[1].split(':')
#             cast_info['var' + str(i)] = var_info[0]
#             cast_info['unit' + str(i)] = var_info[1]
# 
# # process data
# data = [line.split() for line in data]
# data = np.array(data).astype(float)


# f = open('./9908_2_0047_04_c_f_a_m_l_d_B.cnv', 'r')
# f = open('./9908_2_0047_04_c_f_a_m_l_d_B.cnv', 'r')
# f = open('./0103_1_0065_09_c_f_a_m_l_d_B.cnv', 'r')
# text = f.read()
# f.close()
# str1 = 'number of voltages sampled'
# str2 = 'nquan'
# n1 = text.find(str1) + 32
# n2 = text.find(str2) - 3
# aa = text[n1:n2]
# text = text[:n1] + text[n2:]
# f = open('./aaa.cnv', 'w')
# f.write(text)
# f.close()
