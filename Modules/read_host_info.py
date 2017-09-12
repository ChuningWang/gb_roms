import csv
import socket

def read_host_info():
    hostname = socket.gethostname()
    print('hostname: ' +hostname)
    hostname = hostname[:3]

    if hostname=='ysl' or hostname=='yel' or hostname=='gey' or hostname=='pro' or hostname=='cal':
        file_in = '/glade/u/home/chuning/git/gb_roms/Modules/yellowstone.info'
    elif hostname=='aln':
        file_in = '/Users/chuning/git/gb_roms/Modules/alnilam.info'
    elif hostname=='Cod' or hostname=='nwk':
        file_in = '/Users/CnWang/git/gb_roms/Modules/cod.info'

    f = open(file_in, 'rb')
    reader = csv.reader(f)
    sv = dict(reader)
    f.close()
    return sv
