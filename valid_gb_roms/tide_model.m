clear; clc; close all

% first, concatenate ubar, vbar from netCDF files
in_dir = '/glade/scratch/chuning/tmpdir_GB-ref/outputs/2008/';
grd_file = '/glade/p/work/chuning/gb_roms/grd/GlacierBay_lr_grd.nc';
out_file = '/glade/p/work/chuning/gb_roms/tides/Tide_model_lr.mat';
flist = dir([in_dir '*his*nc']);
% flist = flist(end-29:end);

save(out_file)

% stride
dd = 1;
ts = 24;

% read all dimensions
nci = ncinfo([in_dir flist(1).name]);
nci = nci.Dimensions;
for i=1:length(nci)
	eval([nci(i).Name '=' num2str(nci(i).Length) ';']);
end

% read lon lat
lat = ncread(grd_file, 'lat_psi');
lon = ncread(grd_file, 'lon_psi');
msk = ncread(grd_file, 'mask_psi');
h = ncread(grd_file, 'h');
ang = ncread(grd_file, 'angle');
h = 0.25*(h(1:end-1, 1:end-1)+h(1:end-1, 2:end)+h(2:end, 1:end-1)+h(2:end, 2:end));
ang = 0.25*(ang(1:end-1, 1:end-1)+ang(1:end-1, 2:end)+ang(2:end, 1:end-1)+ang(2:end, 2:end));

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
Ubar = zeros(xi_psi, eta_psi, tt*ts);
Ubar = Ubar(1:dd:end, 1:dd:end, :);

Utide = nan(xi_psi, eta_psi, tt*ts);
Utide = Utide(1:dd:end, 1:dd:end, :);

for i=1:tt
	disp(['Loading time step ' num2str(i)])
    t((i-1)*ts+1:(i-1)*ts+ts) = ncread([in_dir flist(i).name], 'ocean_time');
	ubar = ncread([in_dir flist(i).name], 'ubar');
	vbar = ncread([in_dir flist(i).name], 'vbar');
	ubar = 0.5*(ubar(:, 1:end-1, :)+ubar(:, 2:end, :));
	vbar = 0.5*(vbar(1:end-1, :, :)+vbar(2:end, :, :));
	Ub = ubar(1:dd:end, 1:dd:end, :) + 1i*vbar(1:dd:end, 1:dd:end, :);
	Ub = Ub.*exp(1j*repmat(ang, [1, 1, ts]));
	Ubar(:, :, (i-1)*ts+1:(i-1)*ts+ts) = Ub;
end

Ubar = permute(Ubar, [3, 1, 2]);
Utide = permute(Utide, [3, 1, 2]);

% tidal analysis
tlist = {'Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2', 'MF'};
btr = struct;
for i=1:length(tlist)
	tname = tlist{i};
	eval(['btr.', tname, '=nan(4, xi_s, eta_s);'])
end

for i=1:xi_s
	for j=1:eta_s
		if all(~isnan(Ubar(:, i, j)))
			[ts, pout] = t_tide(Ubar(:, i, j), 'interval',1);
			Utide(:, i, j) = pout;
			for k=1:length(tlist)
				tname = tlist{k};
				idx = strmatch(tname, ts.name);
                if ~isempty(idx)
                    eval(['btr.', tname, '(1, i, j) = ts.tidecon(idx, 1);'])
                    eval(['btr.', tname, '(2, i, j) = ts.tidecon(idx, 3);'])
                    eval(['btr.', tname, '(3, i, j) = ts.tidecon(idx, 5);'])
                    eval(['btr.', tname, '(4, i, j) = ts.tidecon(idx, 7);'])
                end
			end
		end
	end
end

save(out_file, 'btr', 'Utide', 'lat', 'lon', 'h', 'ang', 'msk', '-append')
