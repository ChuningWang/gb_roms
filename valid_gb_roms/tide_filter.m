function df = tide_filter(d, t, Wp_hrs)
    dt = (t(2)-t(1))*24.;  % hours
    samplefreq = 24./dt;  % rad per day
    stopfreq = 24./Wp_hrs;
    Wp = stopfreq*2./samplefreq;
    Ws = 2.*Wp;

    [n, Wn] = buttord(Wp, Ws, 3, 60);
    [b, a] = butter(n, Wn);
    df = filtfilt(b, a, d);
end
