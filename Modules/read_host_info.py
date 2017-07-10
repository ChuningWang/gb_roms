import csv
import socket

def read_host_info():
    hostname = socket.gethostname()[:5]
    print('hostname: ' +hostname)

    if hostname=='yslog' or hostname=='yello' or hostname=='geyse' or hostname=='calde':
        file_in = '/glade/u/home/chuning/git/gb_roms/Modules/yellowstone.info'
    elif hostname=='alnil':
        file_in = '/Users/chuning/git/gb_roms/Modules/alnilam.info'

    f = open(file_in, 'rb')
    reader = csv.reader(f)
    sv = dict(reader)
    f.close()
    return sv
