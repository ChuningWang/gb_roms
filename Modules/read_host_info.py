import csv
import socket

def read_host_info():
    hostname = socket.gethostname()
    print('hostname: ' +hostname)
    hostname = hostname[:5]

    if hostname=='yslog' or hostname=='yello' or hostname=='geyse' or hostname=='prong' or hostname=='calde':
        file_in = '/glade/u/home/chuning/git/gb_roms/Modules/yellowstone.info'
    elif hostname=='alnil':
        file_in = '/Users/chuning/git/gb_roms/Modules/alnilam.info'

    f = open(file_in, 'rb')
    reader = csv.reader(f)
    sv = dict(reader)
    f.close()
    return sv
