clc; clear; close all

mod = load('/glade/p/work/chuning/gb_roms/tides/Tide_model_hr.mat');
stn = load('/glade/p/work/chuning/gb_roms/tides/Tide_stn.mat');
tag = 'hr'

dd = 14;

tlist = {'Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2', 'MF'};

for i=1:length(tlist)
    tname = tlist{i};
    disp(tname)
    eval(['fmaj_m = squeeze(mod.btr.' tname '(1, floor(dd/2):dd:end, floor(dd/2):dd:end));'])
    eval(['fmin_m = squeeze(mod.btr.' tname '(2, floor(dd/2):dd:end, floor(dd/2):dd:end));'])
    eval(['finc_m = squeeze(mod.btr.' tname '(3, floor(dd/2):dd:end, floor(dd/2):dd:end));'])
    eval(['fpha_m = squeeze(mod.btr.' tname '(4, floor(dd/2):dd:end, floor(dd/2):dd:end));'])

    eval(['fmaj_s = squeeze(stn.btr.' tname '(1, :));'])
    eval(['fmin_s = squeeze(stn.btr.' tname '(2, :));'])
    eval(['finc_s = squeeze(stn.btr.' tname '(3, :));'])
    eval(['fpha_s = squeeze(stn.btr.' tname '(4, :));'])

    figure('visible','off')
    m_proj('Lambert Conformal Conic', 'lat', [58 59.2], 'lon', [-137.5 -135])
    m_gshhs_h('color',[.5, .5, .5]);
    m_grid('tickdir','in');
    if strcmp(tname, 'MF')
        m_ellipse(mod.lon(floor(dd/2):dd:end, floor(dd/2):dd:end),mod.lat(floor(dd/2):dd:end, floor(dd/2):dd:end),...
                  fmaj_m, fmin_m, finc_m, fpha_m,...
                  500,'line','linewidth',0.5,'linestyle','-')
        m_ellipse(stn.lon,stn.lat,...
                  fmaj_s, fmin_s, finc_s, fpha_s,...
                  500,'line','linewidth',1.0,'color','k')
    else
        m_ellipse(mod.lon(floor(dd/2):dd:end, floor(dd/2):dd:end),mod.lat(floor(dd/2):dd:end, floor(dd/2):dd:end),...
                  fmaj_m, fmin_m, finc_m, fpha_m,...
                  50,'line','linewidth',0.1,'linestyle','-')
        m_ellipse(stn.lon,stn.lat,...
                  fmaj_s, fmin_s, finc_s, fpha_s,...
                  50,'line','linewidth',0.1,'color','k')
    end
    print(gcf, ['/glade/p/work/chuning/gb_roms/figs/tides_cmp/' tname '_' tag], '-dpng', '-r600')
    close
end
