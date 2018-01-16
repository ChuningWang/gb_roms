clear; close all

% first, concatenate ubar, vbar from netCDF files
in_dir = '/glade/scratch/chuning/tmpdir_GB-ref/outputs/2008/';
grd_file = '/glade/p/work/chuning/gb_roms/grd/GlacierBay_lr_grd.nc';
out_file = '/glade/p/work/chuning/gb_roms/tides/Tide_model_lr.nc';
flist = dir([in_dir '*his*nc']);
flist = flist(2:end);
% flist = flist(end-29:end);

% stride
dd = 1;
ts = 12;

interval = 24/ts;

% read all dimensions
nci = ncinfo([in_dir flist(1).name]);
nci = nci.Dimensions;
for i=1:length(nci)
	eval([nci(i).Name '=' num2str(nci(i).Length) ';']);
end

% read lon lat
lat = ncread(grd_file, 'lat_rho');
lon = ncread(grd_file, 'lon_rho');
msk = ncread(grd_file, 'mask_rho');
h = ncread(grd_file, 'h');
ang = ncread(grd_file, 'angle');

lat = lat(1:dd:end, 1:dd:end);
lon = lon(1:dd:end, 1:dd:end);
msk = msk(1:dd:end, 1:dd:end);
h = h(1:dd:end, 1:dd:end);
ang = ang(1:dd:end, 1:dd:end);

[xi_s, eta_s] = size(h);

% total number of entries
tt = length(flist);

% create variables
t = zeros(tt*ts, 1);
Ubar = zeros(xi_s, eta_s, tt*ts);

for i=1:tt
	disp(['Loading time step ' num2str(i)])
    t((i-1)*ts+1:(i-1)*ts+ts) = ncread([in_dir flist(i).name], 'ocean_time');
	ubar = ncread([in_dir flist(i).name], 'ubar');
	vbar = ncread([in_dir flist(i).name], 'vbar');
    ubar(isnan(ubar)) = 0;
    vbar(isnan(vbar)) = 0;
    ubar = 0.5*(ubar(1:end-1, :, :) + ubar(2:end, :, :));
    vbar = 0.5*(vbar(:, 1:end-1, :) + vbar(:, 2:end, :));
    Ub = ubar(:, 2:end-1, :) + 1i*vbar(2:end-1, :, :);
    Ubar(2:end-1, 2:end-1, (i-1)*ts+1:(i-1)*ts+ts) = Ub;
end

Ubar = Ubar(1:dd:end, 1:dd:end, :);

% rotate velocity vector
Ubar = Ubar.*exp(1j*repmat(ang, [1, 1, tt*ts]));
Utide = nan(size(Ubar));

% tidal analysis
tlist = {'Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2', 'MF'};
btr = struct;
for i=1:length(tlist)
	tname = tlist{i};
	eval(['btr.', tname, '=nan(4, xi_s, eta_s);'])
end

for i=1:xi_s
	for j=1:eta_s
		if any(Ubar(i, j, :) ~= 0)
			[tts, pout] = t_tide(Ubar(i, j, :), 'interval', interval, ...
                                 'error', 'wboot', 'output', 'none');
			Utide(i, j, :) = pout;
			% for k=1:length(tlist)
			% 	tname = tlist{k};
			% 	idx = strmatch(tname, tts.name);
            %     if ~isempty(idx)
            %         eval(['btr.', tname, '(1, i, j) = tts.tidecon(idx, 1);'])
            %         eval(['btr.', tname, '(2, i, j) = tts.tidecon(idx, 3);'])
            %         eval(['btr.', tname, '(3, i, j) = tts.tidecon(idx, 5);'])
            %         eval(['btr.', tname, '(4, i, j) = tts.tidecon(idx, 7);'])
            %     end
			% end
		end
	end
end

% save(out_file, 'btr', 'Utide', 'lat', 'lon', 'h', 'ang', 'msk')
nccreate(out_file, 'time', 'dimension', ...
         {'t', length(t)})
nccreate(out_file, 'utide', 'dimension', ...
         {'y', size(Utide, 1), 'x', size(Utide, 2), 't', size(Utide, 3)})
nccreate(out_file, 'vtide', 'dimension', ...
         {'y', size(Utide, 1), 'x', size(Utide, 2), 't', size(Utide, 3)})
nccreate(out_file, 'ures', 'dimension', ...
         {'y', size(Utide, 1), 'x', size(Utide, 2), 't', size(Utide, 3)})
nccreate(out_file, 'vres', 'dimension', ...
         {'y', size(Utide, 1), 'x', size(Utide, 2), 't', size(Utide, 3)})
ncwrite(out_file, 'time', t)
ncwrite(out_file, 'utide', real(Utide))
ncwrite(out_file, 'vtide', imag(Utide))
ncwrite(out_file, 'ures', real(Ubar - Utide))
ncwrite(out_file, 'vres', imag(Ubar - Utide))
