import numpy as np
import pyroms
from CGrid_TPXO8 import CGrid_TPXO8


def get_nc_CGrid_Hill(grdfile, name='TPXO8', \
                       xrange=(4350, 4550), yrange=(6600, 6900), missing_value=-9999):
    """
    grd = get_nc_CGrid_Hill(grdfile)

    Load Cgrid object for Hill's hydrology model from netCDF file
    """

    nc = pyroms.io.Dataset(grdfile)

    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]

    nc.close()

    # # land mask
    # h_msk = zh!=0

    # longitude from -180 to 180
    lonh[lonh>180] = lonh[lonh>180]-360

    lathh, lonhh = np.meshgrid(lath, lonh)

    # generate tpxo8 grid
    # xrange = [4400, 4600]
    # yrange = [6600, 6900]
    return CGrid_TPXO8(lonhh, lathh, msk, missing_value, 'Hill', xrange, yrange)

