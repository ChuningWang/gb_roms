clear; clc; close all

in_dir = '/glade/scratch/chuning/tmpdir_GB-CIRC/outputs/2008/';
out_dir = '/glade/p/work/chuning/gb_roms/tides/';
flist = dir([in_dir '*sta*.nc']);

save([out_dir 'Tide_model_stn.mat'])

tlist = {'Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2', 'MF'};

btr = struct;
bcl = struct;

for i=1:length(tlist)
    tname = tlist{i};
    eval(['btr.', tname, '=nan(4, length(flist));'])
end

lat = zeros(size(flist));
lon = zeros(size(flist));

t = [];
ubar = [];
vbar = [];
u = [];
v = [];

for i=1:length(flist)

    in_file = [in_dir flist(i).name];
    if i==1
        model_name = flist(i).name(1:7);
        h = ncread(in_file, 'h');
        lat = ncread(in_file, 'lat_rho');
        lon = ncread(in_file, 'lon_rho');
        ang = ncread(in_file, 'angle');
    end

    t = [t; ncread(in_file, 'ocean_time')];
    ubar = [ubar ncread(in_file, 'ubar')];
    vbar = [vbar ncread(in_file, 'vbar')];
    u = cat(3, u, ncread(in_file, 'u'));
    v = cat(3, v, ncread(in_file, 'u'));
end

[zz, ss, tt] = size(u);

t = t/(24*60*60);
Ubar = ubar+1i*vbar;
Ubar = Ubar.*exp(1i*repmat(ang, [1, tt]));
U = u+1i*v;
U = U.*exp(1i*repmat(ang', [zz, 1, tt]));

ts = zeros(ss, 8);
ts = zeros(ss, tt);

stn_name = {'SEA0845', 'SEA0846', 'SEA0847', 'SEA0848', 'SEA0849', 'SEA0850', 'SEA1002', 'SEA1003', 'SEA1004'};

for j=1:length(tlist)
    tname = tlist{j};
    eval(['btr.', tname, '=nan(4, length(stn_name));'])
end

for i=1:ss
    [ts, pout] = t_tide(Ubar(i, :), 'interval', 0.1);
    eval(['btr.Utide_' stn_name{i}, '= pout;'])
    for j=1:length(tlist)
        tname = tlist{j};
        idx = strmatch(tname, ts.name);
        if ~isempty(idx)
            eval(['btr.', tname, '(1, i) = ts.tidecon(idx, 1);'])
            eval(['btr.', tname, '(2, i) = ts.tidecon(idx, 3);'])
            eval(['btr.', tname, '(3, i) = ts.tidecon(idx, 5);'])
            eval(['btr.', tname, '(4, i) = ts.tidecon(idx, 7);'])
        end
    end
end

save([out_dir, 'Tide_model_stn.mat'], 'btr', 'lat', 'lon', '-append')

% U = tide_filter(U, t, 1);
%     Ubar = mean(U, 2);
% 
%     for j=1:length(tlist)
%         tname = tlist{j};
%         eval(['bcl.', tname, '_', stn_name, '=nan(4, length(z));'])
%     end
% 
%     for j=1:length(tlist)
%         tname = tlist{j};
%         idx = strmatch(tname, ts.name);
%         if ~isempty(idx)
%             eval(['btr.', tname, '(1, i) = ts.tidecon(idx, 1);'])
%             eval(['btr.', tname, '(2, i) = ts.tidecon(idx, 3);'])
%             eval(['btr.', tname, '(3, i) = ts.tidecon(idx, 5);'])
%             eval(['btr.', tname, '(4, i) = ts.tidecon(idx, 7);'])
%         end

%     % baroclinic
%     Ucl = U-repmat(Ubar, [1, length(z)]);
%     Uit = zeros(size(Ucl));
%     Urs = zeros(size(Ucl));
%     fmaj_cl = zeros(size(z));
%     fmin_cl = zeros(size(z));
%     finc_cl = zeros(size(z));
%     fpha_cl = zeros(size(z));
% 
%     for j=1:length(z)
%         [ts, pout] = t_tide(Ucl(:, j), 'interval',0.1);
%         for k=1:length(tlist)
%             tname = tlist{k};
%             idx = strmatch(tname, ts.name);
%             if ~isempty(idx)
%                 eval(['bcl.', tname, '_', stn_name, '(1, j) = ts.tidecon(idx, 1);'])
%                 eval(['bcl.', tname, '_', stn_name, '(2, j) = ts.tidecon(idx, 3);'])
%                 eval(['bcl.', tname, '_', stn_name, '(3, j) = ts.tidecon(idx, 5);'])
%                 eval(['bcl.', tname, '_', stn_name, '(4, j) = ts.tidecon(idx, 7);'])
%             end
%         end
%         Uit(:, j) = pout;
%         Urs(:, j) = Ucl(:, j)-pout;
%     end
% 
%     eval(['bcl.', stn_name, '_Ucl = Ucl;'])
%     eval(['bcl.', stn_name, '_Uit = Uit;'])
%     eval(['bcl.', stn_name, '_Urs = Urs;'])
%     eval(['bcl.', stn_name, '_z = z;'])
% 
% end
% 
% save([out_dir, 'Tide_stn_lr.mat'], 'btr', 'bcl', 'lat', 'lon', '-append')
