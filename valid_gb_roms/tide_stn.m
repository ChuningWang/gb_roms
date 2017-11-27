clear; clc; close all

in_dir = '/glade/p/work/chuning/data/NOAA_ADCP/';
out_dir = '/glade/p/work/chuning/gb_roms/tides/';
flist = dir([in_dir '*.nc']);

save([out_dir 'Tide_stn.mat'])

tlist = {'Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2', 'MF'};

btr = struct;
bcl = struct;

for i=1:length(tlist)
    tname = tlist{i};
    eval(['btr.', tname, '=nan(4, length(flist));'])
end

lat = zeros(size(flist));
lon = zeros(size(flist));


for i=1:length(flist)
    in_file = [in_dir flist(i).name];
    stn_name = flist(i).name(1:7);
    z = ncread(in_file, 'z');
    t = ncread(in_file, 'time');
    u = ncread(in_file, 'uraw');
    v = ncread(in_file, 'vraw');
    lat(i) = ncreadatt(in_file, '/', 'lat');
    lon(i) = ncreadatt(in_file, '/', 'lon');

    U = u+1i*v;
    U = tide_filter(U, t, 3);
    Ubar = mean(U, 2);

    [ts, pout] = t_tide(Ubar, 'interval',0.1);
    eval(['btr.time_' stn_name, '= t;'])
    eval(['btr.Utide_' stn_name, '= pout;'])
    for j=1:length(tlist)
        tname = tlist{j};
        eval(['bcl.', tname, '_', stn_name, '=nan(4, length(z));'])
    end

    % save barotropic tide velocity to ncfiles
    nccreate([out_dir stn_name '.nc'], 'time', ...
             'dimensions', {'time', length(t)});
    nccreate([out_dir stn_name '.nc'], 'lat')
    nccreate([out_dir stn_name '.nc'], 'lon')
    nccreate([out_dir stn_name '.nc'], 'utide', ...
             'dimensions', {'time', length(t)});
    nccreate([out_dir stn_name '.nc'], 'vtide', ...
             'dimensions', {'time', length(t)});
    nccreate([out_dir stn_name '.nc'], 'ures', ...
             'dimensions', {'time', length(t)});
    nccreate([out_dir stn_name '.nc'], 'vres', ...
             'dimensions', {'time', length(t)});
    ncwrite([out_dir stn_name '.nc'], 'time', t)
    ncwrite([out_dir stn_name '.nc'], 'lat', lat(i))
    ncwrite([out_dir stn_name '.nc'], 'lon', lon(i))
    ncwrite([out_dir stn_name '.nc'], 'utide', real(pout))
    ncwrite([out_dir stn_name '.nc'], 'vtide', imag(pout))
    ncwrite([out_dir stn_name '.nc'], 'ures', real(Ubar - pout))
    ncwrite([out_dir stn_name '.nc'], 'vres', imag(Ubar - pout))

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

    % baroclinic
    Ucl = U-repmat(Ubar, [1, length(z)]);
    Uit = zeros(size(Ucl));
    Urs = zeros(size(Ucl));
    fmaj_cl = zeros(size(z));
    fmin_cl = zeros(size(z));
    finc_cl = zeros(size(z));
    fpha_cl = zeros(size(z));

    for j=1:length(z)
        [ts, pout] = t_tide(Ucl(:, j), 'interval',0.1);
        for k=1:length(tlist)
            tname = tlist{k};
            idx = strmatch(tname, ts.name);
            if ~isempty(idx)
                eval(['bcl.', tname, '_', stn_name, '(1, j) = ts.tidecon(idx, 1);'])
                eval(['bcl.', tname, '_', stn_name, '(2, j) = ts.tidecon(idx, 3);'])
                eval(['bcl.', tname, '_', stn_name, '(3, j) = ts.tidecon(idx, 5);'])
                eval(['bcl.', tname, '_', stn_name, '(4, j) = ts.tidecon(idx, 7);'])
            end
        end
        Uit(:, j) = pout;
        Urs(:, j) = Ucl(:, j)-pout;
    end

    eval(['bcl.time_' stn_name, '= t;'])
    eval(['bcl.', stn_name, '_Ucl = Ucl;'])
    eval(['bcl.', stn_name, '_Uit = Uit;'])
    eval(['bcl.', stn_name, '_Urs = Urs;'])
    eval(['bcl.', stn_name, '_z = z;'])

end

save([out_dir, 'Tide_stn.mat'], 'btr', 'bcl', 'lat', 'lon', '-append')
