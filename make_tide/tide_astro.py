import numpy as np
from scipy import io
from scipy import sparse

def t_astron(d):
    
    D = d/10000.

    # Compute astronomical constants at time d1
    args = np.vstack((np.ones(len(d)), d, D*D, D**3))

    # These are the coefficients of the formulas in the Explan. Suppl.
    sc  = np.array([ 270.434164,13.1763965268,-0.0000850, 0.000000039]);
    hc  = np.array([ 279.696678, 0.9856473354, 0.00002267,0.000000000]);
    pc  = np.array([ 334.329556, 0.1114040803,-0.0007739,-0.00000026]);
    npc = np.array([-259.183275, 0.0529539222,-0.0001557,-0.000000050]);

    # first coeff was 281.220833 in Foreman but Expl. Suppl. has 44.
    ppc = np.array([ 281.220844, 0.0000470684, 0.0000339, 0.000000070]);

    # Compute the parameters; we only need the factional part of the cycle.
    astro = np.fmod(np.dot(np.vstack((sc, hc, pc, npc, ppc)), args)/360.0, 1);

    # Compute lunar time tau, based on fractional part of solar day.
    # We add the hour angle to the longitude of the sun and subtract the
    # longitude of the moon.
    tau = np.fmod(d+0.5, 1) + astro[1, :] - astro[0, :];
    astro = np.vstack((tau, astro));

    # Compute rates of change.
    dargs = np.vstack((np.zeros(len(d)), np.ones(len(d)), 2.0e-4*D, 3.0e-4*D*D))

    ader = np.dot(np.vstack((sc, hc, pc, npc, ppc)), dargs)/360.0;
    dtau = 1.0 + ader[1, :] - ader[0, :]
    ader = np.vstack((dtau, ader));

    return astro, ader

def t_vuf(ltype, d, ju, lat, consts_file = 't_tide/t_constituents.mat'):

    consts = io.loadmat(consts_file)
    doodson = consts['const']['doodson'][0, 0].squeeze()
    semi = consts['const']['semi'][0, 0].squeeze()
    isat = consts['const']['isat'][0, 0].squeeze()
    ishallow = consts['const']['ishallow'][0, 0].squeeze()
    nshallow = consts['const']['nshallow'][0, 0].squeeze()
    rr = consts['sat']['amprat'][0, 0].squeeze()
    ilatfac = consts['sat']['ilatfac'][0, 0].squeeze()
    deldood = consts['sat']['deldood'][0, 0].squeeze()
    phcorr = consts['sat']['phcorr'][0, 0].squeeze()
    iconst = consts['sat']['iconst'][0, 0].squeeze()
    shallow_iname = consts['shallow']['iname'][0, 0].squeeze()
    shallow_coef = consts['shallow']['coef'][0, 0].squeeze()

    # Calculate astronomical arguments at mid-point of data time series.
    astro, ader = t_astron(d)

    v = np.fmod(np.dot(doodson, astro).squeeze()+semi, 1)

    if (abs(lat) < 5):
        lat = np.sign(lat)*5

    # Satellite amplitude ratio adjustment for latitude. 
    slat = np.sin(np.pi*lat/180)

    j = np.where(ilatfac==1)  # latitude correction for diurnal constituents
    rr[j] = rr[j]*0.36309*(1.0-5.0*slat*slat)/slat

    j = np.where(ilatfac==2)  # latitude correction for semi-diurnal constituents
    rr[j] = rr[j]*2.59808*slat

    # Calculate nodal amplitude and phase corrections.
    uu = np.fmod(np.dot(deldood, astro[3:6]).squeeze() + phcorr, 1)

    # Sum up all of the satellite factors for all satellites.

    nsat = len(iconst)
    nfreq = len(isat)

    fsum = (1 + np.sum(
        sparse.csr_matrix(
            ((rr*np.exp(1j*2*np.pi*uu)).squeeze(), (np.arange(nsat), iconst.squeeze()-1)), 
            shape=(nsat, nfreq)).todense(), 
        axis = 0))

    fsum = np.asarray(fsum).squeeze()

    f = abs(fsum)
    u = np.angle(fsum)/(2*np.pi)

    # Compute amplitude and phase corrections for shallow water constituents. 
    for k in np.where(np.isfinite(ishallow.squeeze()))[0].tolist():
        ik = (ishallow[k] + np.arange(nshallow[k]) - 1).tolist()
        f[k] = np.prod(f[shallow_iname[ik]-1]**abs(shallow_coef[ik]))
        u[k] = np.sum(u[shallow_iname[ik]-1]*shallow_coef[ik])
        v[k] = np.sum(v[shallow_iname[ik]-1]*shallow_coef[ik])

    f=f[ju]
    u=u[ju]
    v=v[ju]
    return v, u, f

# from datetime import datetime, timedelta
# d0 = datetime(1899, 12, 31, 12, 00, 00)
# # jd = [datetime(2000, 07, 01, 00, 00, 00), datetime(2001, 07, 01, 00, 00, 00)]
# jd = [datetime(2000, 07, 01, 00, 00, 00)]
# d = np.array([(jd[i]-d0).total_seconds()/timedelta(days=1).total_seconds() for i in range(len(jd))])
# 
# ltype = 'nodal'
# ju = 6
# lat = 58
# 
# v, u, f = t_vuf(ltype, d, ju, lat)
