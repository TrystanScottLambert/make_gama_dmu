"""
Module of transformational functions to convert and manage units of properties to better
match the older version of the GAMA DMU
"""

import numpy as np
import astropy.units as u
import pandas

from config import astropy_cosmo


def convert_angular_to_physical_sep(
    angular_separations_deg: np.ndarray[float], redshift: float
) -> np.ndarray[float]:
    """
    Converts the angular separation in DEG to the physical separation in Kpc/h
    """
    val = (
        angular_separations_deg
        * u.deg
        / astropy_cosmo.arcsec_per_kpc_comoving(redshift)
    )
    return val.to(u.kpc).value

def make_ra_positive(ra_column: pandas.core.series.Series) -> np.ndarray[float]:
    "Adds 360 degress to any negative values"
    arr = np.array(ra_column)
    negative = np.where(arr < 0)
    arr[negative] += 360
    return arr
