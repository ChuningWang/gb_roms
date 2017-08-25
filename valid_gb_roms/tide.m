clear; clc; close all

in_file = '/glade/p/work/chuning/data/NOAA_ADCP/SEA0847_20080810_20080910.nc';

z = ncread(in_file, 'z');
t = ncread(in_file, 'time');
u = ncread(in_file, 'uraw');
v = ncread(in_file, 'vraw');

U = u+1i*v
U = tide_filter(U, t, 1);
Utr = mean(U, 2);
Ucl = U-repmat(Utr, [1, length(z)]);

% barotropic
[tidestruc, pout] = t_tide(Utr, 'interval',0.1);

% baroclinic
% for i=1:length(z)
%     [tidestruc, pout] = t_tide(Ucl, 'interval',0.1);
% end
