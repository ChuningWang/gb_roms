# Read buoy data downloaded from NDBC webpage using wget_buoy.py

# --------------------------------------------------------
# Station Info
# BLTA2 (58.455 N 135.889 W)
# Site elevation: 5 m above mean sea level
# Air temp height: 5.5 m above site elevation
# Anemometer height: 9.1 m above site elevation
# Barometer elevation: 4.8 m above mean sea level

# ELFA2 (58.193 N 136.347 W)
# Sea temp depth: 1.1 m below MLLW

# CSPA2 (58.199 N 136.640 W)
# Site elevation: 25 m above mean sea level
# Air temp height: 1.5 m above site elevation
# Anemometer height: 10 m above site elevation
# Barometer elevation: 26.5 m above mean sea level

# 46083 (58.300 N 137.997 W)
# Site elevation: sea level
# Air temp height: 4 m above site elevation
# Anemometer height: 5 m above site elevation
# Barometer elevation: sea level
# Sea temp depth: 1 m below water line
# Water depth: 136 m
# Watch circle radius: 215 yards

# --------------------------------------------------------
# Notes
# WDIR    Wind direction (the direction the wind is coming from in degrees clockwise from true N) during the same period used for WSPD. See Wind Averaging Methods
# WSPD    Wind speed (m/s) averaged over an eight-minute period for buoys and a two-minute period for land stations. Reported Hourly. See Wind Averaging Methods.
# GST Peak 5 or 8 second gust speed (m/s) measured during the eight-minute or two-minute period. The 5 or 8 second period can be determined by payload, See the Sensor Reporting, Sampling, and Accuracy section.
# WVHT    Significant wave height (meters) is calculated as the average of the highest one-third of all of the wave heights during the 20-minute sampling period. See the Wave Measurements section.
# DPD Dominant wave period (seconds) is the period with the maximum wave energy. See the Wave Measurements section.
# APD Average wave period (seconds) of all waves during the 20-minute period. See the Wave Measurements section.
# MWD The direction from which the waves at the dominant period (DPD) are coming. The units are degrees from true North, increasing clockwise, with North as 0 (zero) degrees and East as 90 degrees. See the Wave Measurements section.
# PRES    Sea level pressure (hPa). For C-MAN sites and Great Lakes buoys, the recorded pressure is reduced to sea level using the method described in NWS Technical Procedures Bulletin 291 (11/14/80). ( labeled BAR in Historical files)
# ATMP    Air temperature (Celsius). For sensor heights on buoys, see Hull Descriptions. For sensor heights at C-MAN stations, see C-MAN Sensor Locations
# WTMP    Sea surface temperature (Celsius). For buoys the depth is referenced to the hull's waterline. For fixed platforms it varies with tide, but is referenced to, or near Mean Lower Low Water (MLLW).
# DEWP    Dewpoint temperature taken at the same height as the air temperature measurement.
# VIS Station visibility (nautical miles). Note that buoy stations are limited to reports from 0 to 1.6 nmi.
# PTDY    Pressure Tendency is the direction (plus or minus) and the amount of pressure change (hPa)for a three hour period ending at the time of observation. (not in Historical files)
# TIDE    The water level in feet above or below Mean Lower Low Water (MLLW).

# ----------------------------------------------------------------------
def wget_buoy(buoy, pth, yr_range, wind=-1):
    import urllib2
    print 'Downloading data to '+pth
    for yr in range(yr_range[0],yr_range[1]+1):
        if wind!=1:
            url = ('http://www.ndbc.noaa.gov/view_text_file.php?filename='+
                   buoy+'h'+"%04d"%yr+'.txt.gz&dir=data/historical/stdmet/')
            try:
                response = urllib2.urlopen(url)
                webContent = response.read()

                print 'Downloading year '+"%04d"%yr+'...'
                f = open(pth+'b'+buoy+'_'+"%04d"%yr+'.txt', 'w')
                f.write(webContent)
                f.close
            except urllib2.HTTPError, e:
                print 'No data in year '+"%04d"%yr+'!!!'
        else:
            url = ('http://www.ndbc.noaa.gov/view_text_file.php?filename='+
                   buoy+'c'+"%04d"%yr+'.txt.gz&dir=data/historical/cwind/')
            try:
                response = urllib2.urlopen(url)
                webContent = response.read()

                print 'Downloading year '+"%04d"%yr+'...'
                f = open(pth+'w'+buoy+'_'+"%04d"%yr+'.txt', 'w')
                f.write(webContent)
                f.close
            except urllib2.HTTPError, e:
                print 'No data in year '+"%04d"%yr+'!!!'


# ----------------------------------------------------------------------
def rd_buoy(pth,stn,mthd = 'bin_avg1day'):

    import csv, glob
    import numpy as np
    from datetime import datetime, timedelta
    import pdb
    import matplotlib.dates as mdates
    
    ft2m = 0.3048 # Feet to Meter

    # stn = 'blta2'
    print 'Getting file names...'
    flist = glob.glob(pth+'b'+stn.lower()+'*.txt')

    # Initiate
    print 'Initiating...'
    YY = []; MM = []; DD = []; hh = []
    wdir = []; wspd = []; gst = []
    wvht = []; dpd = []; apd = []; mwd = []
    pres = []; atmp = []; wtmp = []; dewp = []
    vis = []; tide = [] # tidal height in feet

    for ff in flist:
        print 'Loading '+ff+'...'
        csvData = []
        csvHeader = []
        csvFileObj = open(ff)
        readerObj = csv.reader(csvFileObj, delimiter=' ', skipinitialspace=True)
        for row in readerObj:
            # pdb.set_trace()
            if row[0][0][0] != '2':
                csvHeader.append(row)    # read header
                continue
            csvData.append(row)
        csvFileObj.close()

        csvData = np.asarray(csvData)
        csvData = csvData.astype(np.float)
        csvHeader = csvHeader[0]

        for i in range(np.size(csvHeader)):
            if csvHeader[i]=='#YY' or csvHeader[i]=='YYYY':
                YY = np.concatenate([YY, csvData[:,i]])
            elif csvHeader[i]=='MM':
                MM = np.concatenate([MM, csvData[:,i]])
            elif csvHeader[i]=='DD':
                DD = np.concatenate([DD, csvData[:,i]])
            elif csvHeader[i]=='hh':
                hh = np.concatenate([hh, csvData[:,i]])
            elif csvHeader[i]=='WD' or csvHeader[i]=='WDIR' or csvHeader[i]=='DIR':
                wdir = np.concatenate([wdir, csvData[:,i]])
            elif csvHeader[i]=='WSPD' or csvHeader[i]=='SPD':
                wspd = np.concatenate([wspd, csvData[:,i]])
            elif csvHeader[i]=='GST' or csvHeader[i]=='GSP':
                gst = np.concatenate([gst, csvData[:,i]])
            elif csvHeader[i]=='WVHT':
                wvht = np.concatenate([wvht, csvData[:,i]])
            elif csvHeader[i]=='DPD' or csvHeader[i]=='DOMPD':
                dpd = np.concatenate([dpd, csvData[:,i]])
            elif csvHeader[i]=='APD' or csvHeader[i]=='AVP':
                apd = np.concatenate([apd, csvData[:,i]])
            elif csvHeader[i]=='MWD':
                mwd = np.concatenate([mwd, csvData[:,i]])
            elif csvHeader[i]=='PRES' or csvHeader[i]=='BAR':
                pres = np.concatenate([pres, csvData[:,i]])
            elif csvHeader[i]=='ATMP':
                atmp = np.concatenate([atmp, csvData[:,i]])
            elif csvHeader[i]=='WTMP':
                wtmp = np.concatenate([wtmp, csvData[:,i]])
            elif csvHeader[i]=='DEWP':
                dewp = np.concatenate([dewp, csvData[:,i]])
            elif csvHeader[i]=='VIS':
                vis = np.concatenate([vis, csvData[:,i]])
            elif csvHeader[i]=='TIDE':
                tide = np.concatenate([tide, csvData[:,i]])

    # pdb.set_trace()
    # Change tidal heights from feet to meter
    tide = tide*ft2m

    # Convert YYMMDDhh to python datetime
    print 'Converting time vector to Python datetime...'
    pyt = np.array([datetime(int(YY[i]),
                             int(MM[i]),
                             int(DD[i]),
                             int(hh[i]),
                             0, 0
                            ) for i in range(YY.size)])

    # Mask out invalid data
    print 'Setting invalid data to NaN...'
    wdir[wdir==999]  = np.nan
    wspd[wspd==99]   = np.nan
    gst[gst==99]     = np.nan
    wvht[wvht==99]   = np.nan
    dpd[dpd==99]     = np.nan
    apd[apd==99]     = np.nan
    mwd[mwd==999]    = np.nan
    pres[pres==9999] = np.nan
    atmp[atmp==999]  = np.nan
    wtmp[wtmp==999]  = np.nan
    dewp[dewp==999]  = np.nan
    vis[vis==99]     = np.nan
    tide[tide==99]   = np.nan
    
    if mthd == 'fill':
        # Fill up data gaps
        print 'Filling up data gaps...'
        # Insert some NaN into data gaps
        dt = timedelta(hours=1) # sampling interval
        difft = np.diff(pyt)
        pp = np.squeeze(np.where(difft!=dt))
        kk = 0
        for i in pp:
            nt = int(difft[i].total_seconds()/3600)-1
            pyt_in = np.array([pyt[i+kk]+(j+1)*dt for j in range(nt)]) 
            # pdb.set_trace()
            pyt  = np.concatenate((pyt[:i+kk+1], pyt_in, pyt[i+kk+1:]))
            wdir = np.concatenate([wdir[:i+kk+1], np.ones(nt)*np.nan, wdir[i+kk+1:]])
            wspd = np.concatenate([wspd[:i+kk+1], np.ones(nt)*np.nan, wspd[i+kk+1:]])
            gst  = np.concatenate([gst[:i+kk+1],  np.ones(nt)*np.nan,  gst[i+kk+1:]])
            wvht = np.concatenate([wvht[:i+kk+1], np.ones(nt)*np.nan, wvht[i+kk+1:]])
            dpd  = np.concatenate([dpd[:i+kk+1],  np.ones(nt)*np.nan,  dpd[i+kk+1:]])
            apd  = np.concatenate([apd[:i+kk+1],  np.ones(nt)*np.nan,  apd[i+kk+1:]])
            mwd  = np.concatenate([mwd[:i+kk+1],  np.ones(nt)*np.nan,  mwd[i+kk+1:]])
            pres = np.concatenate([pres[:i+kk+1], np.ones(nt)*np.nan, pres[i+kk+1:]])
            atmp = np.concatenate([atmp[:i+kk+1], np.ones(nt)*np.nan, atmp[i+kk+1:]])
            wtmp = np.concatenate([wtmp[:i+kk+1], np.ones(nt)*np.nan, wtmp[i+kk+1:]])
            dewp = np.concatenate([dewp[:i+kk+1], np.ones(nt)*np.nan, dewp[i+kk+1:]])
            vis  = np.concatenate([vis[:i+kk+1],  np.ones(nt)*np.nan,  vis[i+kk+1:]])
            tide = np.concatenate([tide[:i+kk+1], np.ones(nt)*np.nan, tide[i+kk+1:]])
            kk = kk+nt
    elif mthd == 'bin_avg1day':
        print 'Bin-averaging data with 1 day interval'
        dt = timedelta(days=1) # bin average interval
        nt = (pyt[-1].date()-pyt[0].date()).days+1
        pyt2 = np.array([pyt[0].date()+i*dt for i in range(nt)])
        # Initiate
        wdir2 = np.empty(nt); wspd2 = np.empty(nt); gst2  = np.empty(nt)
        wvht2 = np.empty(nt); dpd2  = np.empty(nt); apd2  = np.empty(nt)
        mwd2  = np.empty(nt); pres2 = np.empty(nt); atmp2 = np.empty(nt)
        wtmp2 = np.empty(nt); dewp2 = np.empty(nt); vis2  = np.empty(nt)
        tide2 = np.empty(nt)
        wdir2[:] = np.nan; wspd2[:] = np.nan; gst2[:]  = np.nan
        wvht2[:] = np.nan; dpd2[:]  = np.nan; apd2[:]  = np.nan
        mwd2[:]  = np.nan; pres2[:] = np.nan; atmp2[:] = np.nan
        wtmp2[:] = np.nan; dewp2[:] = np.nan; vis2[:]  = np.nan
        tide2[:] = np.nan;
        # pdb.set_trace()
        mt = mdates.date2num(pyt)
        mt2 = mdates.date2num(pyt2)
        for i in range(np.size(pyt2[:-1])):
            # pdb.set_trace()
            msk = (mt>=mt2[i]) & (mt<mt2[i+1])
            wdir2[i] = np.nanmean(wdir[msk])
            wspd2[i] = np.nanmean(wspd[msk])
            gst2[i]  = np.nanmean(gst[msk])
            wvht2[i] = np.nanmean(wvht[msk])
            dpd2[i]  = np.nanmean(dpd[msk])
            apd2[i]  = np.nanmean(apd[msk])
            mwd2[i]  = np.nanmean(mwd[msk])
            pres2[i] = np.nanmean(pres[msk])
            atmp2[i] = np.nanmean(atmp[msk])
            wtmp2[i] = np.nanmean(wtmp[msk])
            dewp2[i] = np.nanmean(dewp[msk])
            vis2[i]  = np.nanmean(vis[msk])
            tide2[i] = np.nanmean(tide[msk])

        pyt  = pyt2
        wdir = wdir2
        wspd = wspd2
        gst  = gst2
        wvht = wvht2
        dpd  = dpd2
        apd  = apd2
        mwd  = mwd2
        pres = pres2
        atmp = atmp2
        wtmp = wtmp2
        dewp = dewp2
        vis  = vis2
        tide = tide2
        
    # pdb.set_trace()
    # Mask out invalid data
    print 'Masking invalid data...'
    wdir = np.ma.masked_invalid(wdir)
    wspd = np.ma.masked_invalid(wspd)
    gst  = np.ma.masked_invalid(gst)
    wvht = np.ma.masked_invalid(wvht)
    dpd  = np.ma.masked_invalid(dpd)
    apd  = np.ma.masked_invalid(apd)
    mwd  = np.ma.masked_invalid(mwd)
    pres = np.ma.masked_invalid(pres)
    atmp = np.ma.masked_invalid(atmp)
    wtmp = np.ma.masked_invalid(wtmp)
    dewp = np.ma.masked_invalid(dewp)
    vis  = np.ma.masked_invalid(vis)
    tide = np.ma.masked_invalid(tide)
    
    buoy = {'stn':  stn,
            'pyt':  pyt,
            'wdir': wdir,
            'wspd': wspd,
            'gst':  gst,
            'wvht': wvht,
            'dpd':  dpd,
            'apd':  apd,
            'mwd':  mwd,
            'pres': pres,
            'atmp': atmp,
            'wtmp': wtmp,
            'dewp': dewp,
            'vis':  vis,
            'tide': tide
           }

    return buoy

# ----------------------------------------------------------------------
def smooth2a(matrixIn,Nr,Nc=-1):
    import numpy as np
    from scipy.sparse import spdiags
    import pdb

    if Nc==-1:
        Nc = Nr
    
    if hasattr(matrixIn,'mask'):
        matrixIn[matrixIn.mask] = None

    # matrixIn = np.matrix(matrixIn)
    matrixIn = np.matrix(matrixIn).T

    # pdb.set_trace()
    row = matrixIn.shape[0]
    col = matrixIn.shape[1]
    eL = spdiags(np.ones([2*Nr+1,row]),np.array(range(-Nr,Nr+1)),row,row).todense()
    eR = spdiags(np.ones([2*Nc+1,col]),np.array(range(-Nc,Nc+1)),col,col).todense()
    
    A = np.isnan(matrixIn)
    matrixIn[A]=0

    nrmlize = eL*(~A)*eR
    # pdb.set_trace()
    matrixOut = eL*matrixIn*eR
    matrixOut = np.divide(matrixOut,nrmlize)
    matrixOut = np.squeeze(np.array(matrixOut))
    matrixOut = matrixOut.T
    return matrixOut

