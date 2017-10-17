import csv
import socket

def read_host_info():
    hostname = socket.gethostname()
    print('hostname: ' +hostname)
    hostname = hostname[:3]

    if hostname in ['ysl', 'che', 'yel', 'gey', 'pro', 'cal']:
        file_in = '/glade/u/home/chuning/git/gb_roms/Modules/yellowstone.info'
    elif hostname in ['aln']:
        file_in = '/Users/chuning/git/gb_roms/Modules/alnilam.info'
    elif hostnam ['Cod', 'nwk']:
        file_in = '/Users/CnWang/git/gb_roms/Modules/cod.info'

    f = open(file_in, 'rb')
    reader = csv.reader(f)
    sv = dict(reader)
    f.close()
    return sv
