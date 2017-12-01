import numpy as np
from scipy import io
from scipy import sparse
import netCDF4 as nc

def t_astron(d):

    D = d/10000.

    # Compute astronomical constants at time d1
    args = np.vstack((np.ones(len(d)), d, D*D, D**3))

    # These are the coefficients of the formulas in the Explan. Suppl.
    sc  = np.array([ 270.434164, 13.1763965268, -0.0000850,  0.000000039])
    hc  = np.array([ 279.696678,  0.9856473354, 0.00002267,  0.000000000])
    pc  = np.array([ 334.329556,  0.1114040803, -0.0007739,  -0.00000026])
    npc = np.array([-259.183275,  0.0529539222, -0.0001557, -0.000000050])

    # first coeff was 281.220833 in Foreman but Expl. Suppl. has 44.
    ppc = np.array([ 281.220844, 0.0000470684, 0.0000339, 0.000000070])

    # Compute the parameters; we only need the factional part of the cycle.
    astro = np.fmod(np.dot(np.vstack((sc, hc, pc, npc, ppc)), args)/360.0, 1)

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

def t_getconsts(outfile, consts_file = 't_tide/t_constituents.mat'):

    consts = 146

    c = io.loadmat(consts_file)

    fh = nc.Dataset(outfile, 'w')
    fh.createDimension('namelen', 4)
    fh.createDimension('doodson', 6)
    fh.createDimension('deldood', 3)
    fh.createDimension('consts', consts)
    fh.createDimension('sat', 162)
    fh.createDimension('shallow', 251)

    fh.Title = 'Constants from t_tide'
    fh.Source = 't_tide/t_consts.mat'

    const_name = fh.createVariable('const_name', 'c', ('consts', 'namelen'))
    const_freq = fh.createVariable('const_freq', 'd', ('consts'))
    const_kmpr = fh.createVariable('const_kmpr', 'c', ('consts', 'namelen'))
    const_ikmpr = fh.createVariable('const_ikmpr', 'd', ('consts'))
    const_df = fh.createVariable('const_df', 'd', ('consts'))
    const_doodson = fh.createVariable('const_doodson', 'd', ('consts', 'doodson'))
    const_semi = fh.createVariable('const_semi', 'd', ('consts'))
    const_isat = fh.createVariable('const_isat', 'd', ('consts'))
    const_nsat = fh.createVariable('const_nsat', 'd', ('consts'))
    const_ishallow = fh.createVariable('const_ishallow', 'd', ('consts'))
    const_nshallow = fh.createVariable('const_nshallow', 'd', ('consts'))
    const_doodsonamp = fh.createVariable('const_doodsonamp', 'd', ('consts'))
    const_doodsonspecies = fh.createVariable('const_doodsonspecies', 'd', ('consts'))

    sat_deldood = fh.createVariable('sat_deldood', 'd', ('sat', 'deldood'))
    sat_phcorr = fh.createVariable('sat_phcorr', 'd', ('sat'))
    sat_amprat = fh.createVariable('sat_amprat', 'd', ('sat'))
    sat_ilatfac = fh.createVariable('sat_ilatfac', 'd', ('sat'))
    sat_iconst = fh.createVariable('sat_iconst', 'd', ('sat'))

    shallow_iconst = fh.createVariable('shallow_iconst', 'd', ('shallow'))
    shallow_coef = fh.createVariable('shallow_coef', 'd', ('shallow'))
    shallow_iname = fh.createVariable('shallow_iname', 'd', ('shallow'))

    for i in range(consts):

        const_name[i, :] = list(c['const']['name'][0, 0][i])
        const_kmpr[i, :] = list(c['const']['kmpr'][0, 0][i])

    const_freq[:] = c['const']['freq'][0, 0].squeeze()
    const_ikmpr[:] = c['const']['ikmpr'][0, 0].squeeze()
    const_df[:] = c['const']['df'][0, 0].squeeze()
    const_doodson[:, :] = c['const']['doodson'][0, 0].squeeze()
    const_semi[:] = c['const']['semi'][0, 0].squeeze()
    const_isat[:] = c['const']['isat'][0, 0].squeeze()
    const_nsat[:] = c['const']['nsat'][0, 0].squeeze()
    const_ishallow[:] = c['const']['ishallow'][0, 0].squeeze()
    const_nshallow[:] = c['const']['nshallow'][0, 0].squeeze()
    const_doodsonamp[:] = c['const']['doodsonamp'][0, 0].squeeze()
    const_doodsonspecies[:] = c['const']['doodsonspecies'][0, 0].squeeze()

    sat_deldood[:, :] = c['sat']['deldood'][0, 0].squeeze()
    sat_phcorr[:] = c['sat']['phcorr'][0, 0].squeeze()
    sat_amprat[:] = c['sat']['amprat'][0, 0].squeeze()
    sat_ilatfac[:] = c['sat']['ilatfac'][0, 0].squeeze()
    sat_iconst[:] = c['sat']['iconst'][0, 0].squeeze()

    shallow_iconst[:] = c['shallow']['iconst'][0, 0].squeeze()
    shallow_coef[:] = c['shallow']['coef'][0, 0].squeeze()
    shallow_iname[:] = c['shallow']['iname'][0, 0].squeeze()

    fh.close()

def t_vuf(d, ju, lat, consts_file = 'tide_consts.nc', ltype = 'nodal'):

    fh = nc.Dataset('tide_consts.nc', 'r')
    doodson = fh.variables['const_doodson'][:]
    semi = fh.variables['const_semi'][:]
    isat = fh.variables['const_isat'][:]
    ishallow = fh.variables['const_ishallow'][:]
    nshallow = fh.variables['const_nshallow'][:]

    rr = fh.variables['sat_amprat'][:]
    ilatfac = fh.variables['sat_ilatfac'][:]
    deldood = fh.variables['sat_deldood'][:]
    phcorr = fh.variables['sat_phcorr'][:]
    iconst = fh.variables['sat_iconst'][:]

    shallow_iname = fh.variables['shallow_iname'][:]
    shallow_coef = fh.variables['shallow_coef'][:]
    fh.close()

    # Calculate astronomical arguments at mid-point of data time series.
    astro, ader = t_astron(d)

    v = np.fmod(np.dot(doodson, astro).squeeze()+semi, 1)

    if (abs(lat) < 5):
        lat = np.sign(lat)*5

    # Satellite amplitude ratio adjustment for latitude. 
    slat = np.sin(np.pi*lat/180)

    j = np.where(ilatfac == 1)  # latitude correction for diurnal constituents
    rr[j] = rr[j]*0.36309*(1.0-5.0*slat*slat)/slat

    j = np.where(ilatfac == 2)  # latitude correction for semi-diurnal constituents
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

        f[k] = np.prod(f[shallow_iname[ik].astype(int)-1]**abs(shallow_coef[ik]))
        u[k] = np.sum(u[shallow_iname[ik].astype(int)-1]*shallow_coef[ik])
        v[k] = np.sum(v[shallow_iname[ik].astype(int)-1]*shallow_coef[ik])

    f = f[ju]
    u = u[ju]
    v = v[ju]
    return v, u, f

t_getconsts('tide_consts.nc')

