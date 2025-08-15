"""
Sanity check module. Comparing physically derived properties to G3Cv10 for conistancy.
"""

import numpy as np
import pylab as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd



def read_fits(fits_name: str) -> pd.DataFrame:
    """
    Reads in the fits file and returns a data frame
    """
    data = Table.read(fits_name, format="fits")
    return data.to_pandas()


def bin_data(array_1: np.ndarray, array_2: np.ndarray, title: str) -> None:
    """
    Automatically bins the data to the same resolution and normalizes.
    """
    n_bin = 100
    bins = np.linspace(np.min(array_2), np.max(array_2), n_bin)
    plt.hist(array_1, bins=bins, density=True, histtype="step")
    plt.hist(array_2, bins=bins, density=True, histtype="step")
    plt.yscale('log')
    plt.title(title)
    plt.show()


def check_central_accuracy(
    ra_1: pd.Series, dec_1: pd.Series, ra_2: pd.Series, dec_2: pd.Series, title: str
) -> None:
    "Does a simple cross match and works out the angular offset which is then displayed."
    ra_1, dec_1 = np.array(ra_1), np.array(dec_1)
    ra_2, dec_2 = np.array(ra_2), np.array(dec_2)

    c1 = SkyCoord(ra=ra_1 * u.deg, dec=dec_1 * u.deg)
    c2 = SkyCoord(ra=ra_2 * u.deg, dec=dec_2 * u.deg)
    _, d2d, _ = c1.match_to_catalog_sky(c2)
    plt.hist(d2d.to(u.arcsec), bins=100)
    plt.xlabel("Offset arcseconds")
    plt.ylabel("counts")
    plt.title(title)
    plt.show()


gama_10 = read_fits("gama_catalogs/gama_v10_group_cats/G3CFoFGroupv10.fits")
gama_10 = gama_10[gama_10["IterCenRA"] > 100]

gama_11 = read_fits("G3CFoFGroup.fits")
gama_11 = gama_11[gama_11["GAMARegion"] != b"g23"]


things_to_check = [
    "IterCenRA",
    "IterCenDec",
    "CenRA",
    "CenDec",
    "BCGRA",
    "BCGDec",
    "Rad50",
    "Rad100",
    "Rad1Sig",
    "VelDisp",
    # "VelDispErr",
    "MassProxy",
    "TotFluxProxy",
    "TotRmag",
    "MassA",
    "LumB",
    "MassAfunc",
    "LumBfunc",
    "Nfof",
    "Zfof",
]

centers = ["IterCen", "Cen", "BCG"]

# for center in centers:
#     ra = f'{center}RA'
#     dec = f'{center}Dec'
#     check_central_accuracy(gama_10[ra], gama_10[dec], gama_11[ra], gama_11[dec], center)


for thing in things_to_check:
    array_10 = np.array(gama_10[thing])
    array_11 = np.array(gama_11[thing])
    bin_data(array_10, array_11, thing)
